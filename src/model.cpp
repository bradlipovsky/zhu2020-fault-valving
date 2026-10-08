// Independent implementation from the equations in Zhu et al. (2020).
// No source or state from the authors' simulator or the earlier attempt is used.
#include <fftw3.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Vec = std::vector<double>;
using State = std::array<Vec,4>; // slip deficit [m], psi, normalized k*, excess p [Pa]
constexpr double YEAR=365.25*86400., PI=3.14159265358979323846;
constexpr double gamma_ars=1.-0.7071067811865475244;
constexpr double delta_ars=1.-1./(2.*gamma_ars);

struct Parameters {
    int n=32768, threads=4;
    double Ly=500000, Lz=500000, mu=32.4e9, radiation=4.68e6;
    double Vp=1e-9, V0=1e-6, f0=.6, dc=.002;
    double rho=1000, gravity=9.8, viscosity=1e-4, storage=1e-11;
    double normal_gradient=22000, stress_scale=30e6;
    double kmin=1e-19, kmax=1e-15, L=1, T=1e8, influx=3e-9;
    double a_surface=.01, a_15km=.03, a_deep_slope=3e-6;
    double ab_shallow=-.01, ab_corner=13500, vw_bottom=17000;
    double perturbation=1e-4, perturb_depth=12000, perturb_width=1000;
    double initial_locking_depth=0, initial_transition_width=1500, initial_velocity_fraction=1e-4;
    double years=200, tolerance=1e-3, dt_max=1e6;
    double output_years=.025, output_seconds=.5, output_depth=30000;
    double output_spacing=30, seismic_threshold=1e-3;
    bool fixed_pressure=false;
};

Parameters read_config(const std::string& path) {
    Parameters p;
    std::map<std::string,double*> d={
        {"Ly",&p.Ly},{"Lz",&p.Lz},{"mu",&p.mu},{"radiation",&p.radiation},
        {"Vp",&p.Vp},{"V0",&p.V0},{"f0",&p.f0},{"dc",&p.dc},
        {"rho",&p.rho},{"gravity",&p.gravity},{"viscosity",&p.viscosity},
        {"storage",&p.storage},{"normal_gradient",&p.normal_gradient},
        {"stress_scale",&p.stress_scale},{"kmin",&p.kmin},{"kmax",&p.kmax},
        {"L",&p.L},{"T",&p.T},{"influx",&p.influx},
        {"a_surface",&p.a_surface},{"a_15km",&p.a_15km},{"a_deep_slope",&p.a_deep_slope},
        {"ab_shallow",&p.ab_shallow},{"ab_corner",&p.ab_corner},{"vw_bottom",&p.vw_bottom},
        {"perturbation",&p.perturbation},{"perturb_depth",&p.perturb_depth},
        {"perturb_width",&p.perturb_width},{"years",&p.years},
        {"initial_locking_depth",&p.initial_locking_depth},
        {"initial_transition_width",&p.initial_transition_width},
        {"initial_velocity_fraction",&p.initial_velocity_fraction},
        {"tolerance",&p.tolerance},{"dt_max",&p.dt_max},
        {"output_years",&p.output_years},{"output_seconds",&p.output_seconds},
        {"output_depth",&p.output_depth},{"output_spacing",&p.output_spacing},
        {"seismic_threshold",&p.seismic_threshold}};
    std::ifstream f(path); if(!f) throw std::runtime_error("Cannot open config: "+path);
    std::string line;
    while(std::getline(f,line)) {
        line=line.substr(0,line.find('#'));
        auto e=line.find('='); if(e==std::string::npos) continue;
        std::string key=line.substr(0,e);
        key.erase(std::remove_if(key.begin(),key.end(),[](unsigned char c){return std::isspace(c);}),key.end());
        double v=std::stod(line.substr(e+1));
        if(key=="n") p.n=int(v);
        else if(key=="threads") p.threads=int(v);
        else if(key=="fixed_pressure") p.fixed_pressure=(v!=0);
        else if(d.count(key)) *d[key]=v;
        else throw std::runtime_error("Unknown config key: "+key);
    }
    if(p.n<8 || p.Lz<=0 || p.Ly<=0 || p.T<=0 || p.dc<=0 || p.tolerance<=0)
        throw std::runtime_error("Invalid grid, time, or length scale");
    return p;
}

State zeros(int n) { State x; for(auto& v:x) v.resize(n,0.); return x; }
double harmonic(double x,double y) { return 2*x*y/(x+y); }
double asinh_exp(double x) { return x>30 ? x+std::log(2.) : std::asinh(std::exp(x)); }
double friction(double v,double psi,double a,const Parameters& p) {
    if(v==0) return 0;
    return std::copysign(a*asinh_exp(std::log(std::abs(v)/(2*p.V0))+psi/a),v);
}
double velocity(double tau,double psi,double a,double normal,const Parameters& p) {
    if(tau==0) return 0;
    double target=std::abs(tau), hi=std::log(target/p.radiation);
    double ratio=target/(normal*a);
    double log_sinh=ratio>30?ratio-std::log(2.):std::log(std::sinh(ratio));
    double lo=-740, w=std::min(hi,std::log(2*p.V0)+log_sinh-psi/a);
    w=std::max(lo,w);
    for(int j=0;j<60;j++) {
        double v=std::exp(w), arg=w-std::log(2*p.V0)+psi/a;
        double f=normal*a*asinh_exp(arg)+p.radiation*v-target;
        if(std::abs(f)<1e-8+target*2e-13) return std::copysign(v,tau);
        if(f>0) hi=w; else lo=w;
        double slope=normal*a*(arg>30?1.:1./std::sqrt(1.+std::exp(-2*arg)))+p.radiation*v;
        double candidate=w-f/slope;
        w=(std::isfinite(candidate)&&candidate>lo&&candidate<hi)?candidate:.5*(lo+hi);
    }
    return std::copysign(std::exp(w),tau);
}

struct Model {
    Parameters p;
    int n;
    double h;
    Vec z,a,b,prestress,buffer,spectral,stiffness,vel,tau;
    Vec permeability,face,lower,diagonal,upper,rhs,iterate,next;
    fftw_plan forward,backward;
    Model(Parameters parameters):p(parameters),n(p.n),h(p.Lz/n),
        z(n),a(n),b(n),prestress(n),buffer(n),spectral(n),stiffness(n),vel(n),tau(n),
        permeability(n),face(n+1),lower(n),diagonal(n),upper(n),rhs(n),iterate(n),next(n) {
        omp_set_num_threads(p.threads);
        for(int i=0;i<n;i++) {
            z[i]=(i+.5)*h;
            a[i]=z[i]<=15000?p.a_surface+(p.a_15km-p.a_surface)*z[i]/15000:
                p.a_15km+p.a_deep_slope*(z[i]-15000);
            double ab=z[i]<=p.ab_corner?p.ab_shallow:
                p.ab_shallow+(-p.ab_shallow)*(z[i]-p.ab_corner)/(p.vw_bottom-p.ab_corner);
            b[i]=a[i]-ab;
            if(a[i]<=0 || b[i]<=0) throw std::runtime_error("Nonpositive friction parameter");
            double wave=PI*i/p.Lz;
            stiffness[i]=i==0?p.mu/(2*p.Ly):.5*p.mu*wave/std::tanh(wave*p.Ly);
        }
        forward=fftw_plan_r2r_1d(n,buffer.data(),spectral.data(),FFTW_REDFT10,FFTW_ESTIMATE);
        backward=fftw_plan_r2r_1d(n,spectral.data(),buffer.data(),FFTW_REDFT01,FFTW_ESTIMATE);
        if(!forward||!backward) throw std::runtime_error("FFTW plan failed");
    }
    ~Model(){fftw_destroy_plan(forward);fftw_destroy_plan(backward);}
    double kstar(double u) const { return p.kmin+(p.kmax-p.kmin)*u; }
    double normal(int i,double excess) const {return (p.normal_gradient-p.rho*p.gravity)*z[i]-excess;}
    double k_value(int i,double u,double excess) const {
        return p.kmin+(p.kmax-p.kmin)*u*std::exp(-normal(i,excess)/p.stress_scale);
    }
    void faces(const Vec& u,const Vec& excess) {
        #pragma omp parallel for schedule(static) if(n>4096)
        for(int i=0;i<n;i++) permeability[i]=k_value(i,u[i],excess[i]);
        face[0]=harmonic(kstar(u[0]),permeability[0]);
        for(int i=1;i<n;i++) face[i]=harmonic(permeability[i-1],permeability[i]);
        face[n]=permeability[n-1];
    }
    State initial() {
        State x=zeros(n);
        double u=(p.Vp/p.L)/(p.Vp/p.L+1/p.T);
        for(int i=0;i<n;i++) {
            x[2][i]=u;
            double low=(i==0?0:x[3][i-1]), high=(p.normal_gradient-p.rho*p.gravity)*z[i];
            for(int k=0;k<70;k++) {
                double mid=.5*(low+high), current=k_value(i,u,mid);
                double left=i==0?kstar(u):k_value(i-1,u,x[3][i-1]);
                double flux=harmonic(left,current)*(mid-(i==0?0:x[3][i-1]))/(p.viscosity*h*(i==0?.5:1.));
                if(flux>p.influx) high=mid; else low=mid;
            }
            x[3][i]=.5*(low+high);
            x[1][i]=p.f0-b[i]*std::log(p.Vp/p.V0);
            double initial_v=p.Vp;
            if(p.initial_locking_depth>0)
                initial_v*=p.initial_velocity_fraction+(1-p.initial_velocity_fraction)*.5*
                    (1+std::tanh((z[i]-p.initial_locking_depth)/p.initial_transition_width));
            prestress[i]=normal(i,x[3][i])*friction(initial_v,x[1][i],a[i],p)+p.radiation*initial_v;
        }
        for(int i=0;i<n;i++) x[1][i]-=p.perturbation*std::exp(-.5*std::pow((z[i]-p.perturb_depth)/p.perturb_width,2));
        return x;
    }
    void elastic(const Vec& slip_deficit) {
        std::copy(slip_deficit.begin(),slip_deficit.end(),buffer.begin());
        fftw_execute(forward);
        for(int i=0;i<n;i++) spectral[i]*=stiffness[i];
        fftw_execute(backward);
        for(int i=0;i<n;i++) tau[i]=prestress[i]-buffer[i]/(2*n);
    }
    bool mechanical(const State& x,State& f) {
        elastic(x[0]);
        int bad=0;
        #pragma omp parallel for reduction(|:bad) schedule(static) if(n>4096)
        for(int i=0;i<n;i++) {
            double N=normal(i,x[3][i]);
            if(N<=0 || !std::isfinite(N) || x[2][i]<0 || x[2][i]>1 || !std::isfinite(x[1][i])) {bad=1;continue;}
            double v=velocity(tau[i],x[1][i],a[i],N,p), speed=std::abs(v);
            vel[i]=v;
            f[0][i]=v-p.Vp;
            f[1][i]=b[i]/p.dc*(p.V0*std::exp((p.f0-x[1][i])/b[i])-speed);
            f[2][i]=speed/p.L*(1-x[2][i])-x[2][i]/p.T;
            if(!std::isfinite(f[1][i])) bad=1;
        }
        return !bad;
    }
    bool pressure(const Vec& base,const Vec& u,double dt,Vec& result) {
        if(p.fixed_pressure || dt==0) {result=base;return true;}
        iterate=base;
        double factor=dt/(p.storage*p.viscosity*h*h);
        for(int k=0;k<50;k++) {
            faces(u,iterate);
            for(int i=0;i<n;i++) {
                double l=factor*face[i]*(i==0?2:1), r=i==n-1?0:factor*face[i+1];
                lower[i]=i==0?0:-l; upper[i]=-r; diagonal[i]=1+l+r;
                rhs[i]=base[i];
            }
            rhs[n-1]+=dt*p.influx/(p.storage*h);
            for(int i=1;i<n;i++) {
                double m=lower[i]/diagonal[i-1];
                diagonal[i]-=m*upper[i-1]; rhs[i]-=m*rhs[i-1];
            }
            next[n-1]=rhs[n-1]/diagonal[n-1];
            for(int i=n-2;i>=0;i--) next[i]=(rhs[i]-upper[i]*next[i+1])/diagonal[i];
            double change=0;
            for(int i=0;i<n;i++) change=std::max(change,std::abs(next[i]-iterate[i]));
            if(!std::isfinite(change)) return false;
            if(change<1e-3) {result=next;return true;}
            // Fixed-point damping changes only the nonlinear solve, not the equation.
            double relaxation=k<8?1.:.5;
            for(int i=0;i<n;i++) iterate[i]+=relaxation*(next[i]-iterate[i]);
        }
        return false;
    }
    bool step(const State& x,double dt,State& out,State& stage,State& f0,State& f1) {
        if(!mechanical(x,f0)) return false;
        for(int j=0;j<3;j++) for(int i=0;i<n;i++) stage[j][i]=x[j][i]+gamma_ars*dt*f0[j][i];
        if(!pressure(x[3],stage[2],gamma_ars*dt,stage[3])) return false;
        if(!mechanical(stage,f1)) return false;
        for(int j=0;j<3;j++) for(int i=0;i<n;i++) out[j][i]=x[j][i]+dt*(delta_ars*f0[j][i]+(1-delta_ars)*f1[j][i]);
        Vec base(n);
        for(int i=0;i<n;i++) base[i]=x[3][i]+(1-gamma_ars)/gamma_ars*(stage[3][i]-x[3][i]);
        return pressure(base,out[2],gamma_ars*dt,out[3]);
    }
    double error(const State& coarse,const State& fine) const {
        double result=0;
        for(int i=0;i<n;i++) {
            result=std::max(result,std::abs(fine[0][i]-coarse[0][i])/(3*p.tolerance*p.dc));
            result=std::max(result,std::abs(fine[1][i]-coarse[1][i])/(3*p.tolerance*a[i]));
            result=std::max(result,std::abs(fine[2][i]-coarse[2][i])/(3*p.tolerance*(.01+std::abs(fine[2][i]))));
            result=std::max(result,std::abs(fine[3][i]-coarse[3][i])/(3*p.tolerance*p.stress_scale));
        }
        return result;
    }
    Vec flux(const State& x) {
        faces(x[2],x[3]); Vec q(n+1);
        q[0]=face[0]*x[3][0]*2/(p.viscosity*h);
        for(int i=1;i<n;i++) q[i]=face[i]*(x[3][i]-x[3][i-1])/(p.viscosity*h);
        q[n]=p.influx;return q;
    }
};

// Kennedy and Carpenter, NASA/TM-2001-211038, Appendix C.
// Coefficients checked against the SUNDIALS 7.7 written Butcher-table documentation.
struct ARK4 {
    Model& m;
    std::array<State,6> rate;
    State stage;
    Vec base;
    static constexpr double ae[6][6]={
        {0,0,0,0,0,0}, {.5,0,0,0,0,0},
        {13861./62500,6889./62500,0,0,0,0},
        {-116923316275./2393684061468,-2731218467317./15368042101831,9408046702089./11113171139209,0,0,0},
        {-451086348788./2902428689909,-2682348792572./7519795681897,12662868775082./11960479115383,3355817975965./11060851509271,0,0},
        {647845179188./3216320057751,73281519250./8382639484533,552539513391./3454668386233,3354512671639./8306763924573,4040./17871,0}};
    static constexpr double ai[6][6]={
        {0,0,0,0,0,0}, {.25,.25,0,0,0,0},
        {8611./62500,-1743./31250,.25,0,0,0},
        {5012029./34652500,-654441./2922500,174375./388108,.25,0,0},
        {15267082809./155376265600,-71443401./120774400,730878875./902184768,2285395./8070912,.25,0},
        {82889./524892,0,15625./83664,69875./102672,-2260./8211,.25}};
    static constexpr double high[6]={82889./524892,0,15625./83664,69875./102672,-2260./8211,.25};
    static constexpr double low[6]={4586570599./29645900160,0,178811875./945068544,814220225./1159782912,-3700637./11593932,61727./225920};
    ARK4(Model& model):m(model),stage(zeros(m.n)),base(m.n) {for(auto& r:rate)r=zeros(m.n);}
    bool step(const State& x,double dt,State& out,double& error) {
        for(int s=0;s<6;s++) {
            for(int j=0;j<4;j++)for(int i=0;i<m.n;i++) {
                double change=0;
                for(int k=0;k<s;k++)change+=(j==3?ai[s][k]:ae[s][k])*rate[k][j][i];
                stage[j][i]=x[j][i]+dt*change;
            }
            if(s>0) {
                base=stage[3];
                if(!m.pressure(base,stage[2],ai[s][s]*dt,stage[3]))return false;
            }
            if(!m.mechanical(stage,rate[s]))return false;
            if(!m.p.fixed_pressure) {
                Vec q=m.flux(stage);
                for(int i=0;i<m.n;i++)rate[s][3][i]=(q[i+1]-q[i])/(m.p.storage*m.h);
            }
        }
        error=0;
        for(int j=0;j<4;j++)for(int i=0;i<m.n;i++) {
            double increment=0,estimate=0;
            for(int s=0;s<6;s++) {
                increment+=high[s]*rate[s][j][i];
                estimate+=(high[s]-low[s])*rate[s][j][i];
            }
            out[j][i]=x[j][i]+dt*increment;
            double scale=j==0?m.p.dc:j==1?m.a[i]:j==2?.01+std::abs(out[j][i]):m.p.stress_scale;
            error=std::max(error,std::abs(dt*estimate)/(m.p.tolerance*scale));
            if(!std::isfinite(out[j][i]))return false;
            if(j==2 && (out[j][i]<0 || out[j][i]>1))return false;
            if(j==3 && m.normal(i,out[j][i])<=0)return false;
        }
        return true;
    }
};

template<typename T> void put(std::ostream& f,const T& x){f.write(reinterpret_cast<const char*>(&x),sizeof(T));}
template<typename T> void get(std::istream& f,T& x){f.read(reinterpret_cast<char*>(&x),sizeof(T));}

struct Output {
    Model& m;
    std::string dir;
    std::ofstream fields,history;
    std::vector<int> indices;
    double last_time=-1e100,last_snapshot_speed=0;
    int64_t snapshots=0;
    Output(Model& model,std::string directory,bool resume):m(model),dir(directory) {
        std::filesystem::create_directories(dir);
        int stride=std::max(1,int(std::round(m.p.output_spacing/m.h)));
        for(int i=0;i<m.n&&m.z[i]<=m.p.output_depth;i+=stride) indices.push_back(i);
        fields.open(dir+"/fields.bin",std::ios::binary|(resume?std::ios::app:std::ios::trunc));
        history.open(dir+"/history.csv",resume?std::ios::app:std::ios::trunc);
        if(!fields||!history) throw std::runtime_error("Cannot open output files");
        if(!resume) {
            fields.write("ZHUIND01",8); uint64_t nz=indices.size();put(fields,nz);
            for(int i:indices) put(fields,m.z[i]);
            history<<"time_s,dt_s,step,vmax_m_s,z_vmax_m,min_effective_pa,surface_flux_m_s,error";
            for(int depth:{5000,10000,15000,20000})
                for(auto s:{"velocity","slip","effective","permeability","flux"})history<<","<<s<<"_"<<depth;
            history<<"\n";
        }
        history<<std::setprecision(17);
    }
    void record(const State& x,double t,double dt,int64_t step,double err,bool force=false) {
        State f=zeros(m.n); if(!m.mechanical(x,f)) throw std::runtime_error("Invalid accepted state");
        Vec q=m.flux(x);
        auto it=std::max_element(m.vel.begin(),m.vel.end(),[](double x,double y){return std::abs(x)<std::abs(y);});
        int imax=int(it-m.vel.begin()); double vmax=std::abs(*it), minN=1e100;
        for(int i=0;i<m.n;i++)minN=std::min(minN,m.normal(i,x[3][i]));
        history<<t<<","<<dt<<","<<step<<","<<vmax<<","<<m.z[imax]<<","<<minN<<","<<q[0]<<","<<err;
        for(int depth:{5000,10000,15000,20000}) {
            int i=std::min(m.n-1,int(depth/m.h));
            history<<","<<m.vel[i]<<","<<x[0][i]+m.p.Vp*t<<","<<m.normal(i,x[3][i])<<","<<m.permeability[i]<<","<<.5*(q[i]+q[i+1]);
        }
        history<<"\n";
        double spacing=vmax>m.p.seismic_threshold?m.p.output_seconds:m.p.output_years*YEAR;
        bool transition=(vmax>m.p.seismic_threshold)!=(last_snapshot_speed>m.p.seismic_threshold);
        bool rapid_change=std::abs(std::log10(std::max(vmax,1e-30)/std::max(last_snapshot_speed,1e-30)))>.2;
        if(force || t-last_time>=spacing || transition || rapid_change) {
            put(fields,t);put(fields,step);
            for(int k=0;k<6;k++)for(int i:indices) {
                double value=k==0?x[0][i]+m.p.Vp*t:k==1?m.vel[i]:k==2?m.normal(i,x[3][i]):
                    k==3?m.permeability[i]:k==4?m.kstar(x[2][i]):.5*(q[i]+q[i+1]);
                float v=float(value);put(fields,v);
            }
            last_time=t;last_snapshot_speed=vmax;snapshots++;
        }
    }
    void checkpoint(const State& x,double t,double dt,int64_t step,int64_t rejected) {
        fields.flush();history.flush();
        std::ofstream f(dir+"/checkpoint.tmp",std::ios::binary);
        put(f,t);put(f,dt);put(f,step);put(f,rejected);
        for(auto& v:x) f.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));
        f.close();std::filesystem::rename(dir+"/checkpoint.tmp",dir+"/checkpoint.bin");
    }
};

void run(const std::string& config,const std::string& directory,bool resume) {
    Parameters p=read_config(config); Model m(p); State x=m.initial();
    State fine=zeros(p.n);ARK4 solver(m);
    double t=0,dt=10; int64_t steps=0,rejected=0;
    std::filesystem::create_directories(directory);
    if(resume) {
        std::ifstream f(directory+"/checkpoint.bin",std::ios::binary);
        get(f,t);get(f,dt);get(f,steps);get(f,rejected);
        for(auto& v:x) f.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(double));
        if(!f) throw std::runtime_error("Invalid checkpoint");
    } else std::filesystem::copy_file(config,directory+"/config.cfg",std::filesystem::copy_options::overwrite_existing);
    Output out(m,directory,resume);
    if(!resume)out.record(x,t,0,steps,0,true);
    auto start=std::chrono::steady_clock::now(),last=start;
    double end=p.years*YEAR;
    while(t<end) {
        dt=std::min({dt,p.dt_max,end-t});
        double error=1e20;
        bool okay=solver.step(x,dt,fine,error);
        if(!okay)error=1e20;
        if(!std::isfinite(error))error=1e20;
        if(error<=1) {
            x.swap(fine);t+=dt;steps++;
            out.record(x,t,dt,steps,error,t==end);
        } else rejected++;
        double multiplier=error==0?2.:std::clamp(.9*std::pow(error,-.25),.15,2.);
        dt*=multiplier;
        if(dt<1e-10 || t+dt==t)throw std::runtime_error("Step below 1e-10 s; model or nonlinear solver failed at year "+std::to_string(t/YEAR));
        auto now=std::chrono::steady_clock::now();
        if(std::chrono::duration<double>(now-last).count()>20) {
            double elapsed=std::chrono::duration<double>(now-start).count();
            std::cout<<std::setprecision(8)<<"year="<<t/YEAR<<" steps="<<steps<<" rejected="<<rejected<<" dt="<<dt<<" elapsed_s="<<elapsed<<std::endl;
            out.checkpoint(x,t,dt,steps,rejected);last=now;
        }
    }
    out.checkpoint(x,t,dt,steps,rejected);
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream meta(directory+"/completed.json");
    meta<<std::setprecision(15)<<"{\n  \"completed\": true, \"time_s\": "<<t<<", \"years\": "<<t/YEAR
        <<", \"accepted_steps\": "<<steps<<", \"rejected_steps\": "<<rejected<<", \"elapsed_s\": "<<elapsed
        <<", \"n\": "<<p.n<<", \"dz_m\": "<<m.h<<", \"tolerance\": "<<p.tolerance<<"\n}\n";
    std::cout<<"COMPLETED years="<<t/YEAR<<" steps="<<steps<<" elapsed_s="<<elapsed<<std::endl;
}

void require(bool condition,const std::string& message) {if(!condition)throw std::runtime_error("TEST FAILED: "+message);}
void test() {
    Parameters p;p.n=256;p.threads=1;p.perturbation=0;Model m(p);State x=m.initial();
    Vec q=m.flux(x);double qerr=0;
    for(double v:q)qerr=std::max(qerr,std::abs(v/p.influx-1));
    require(qerr<1e-7,"steady Darcy flux");
    State f=zeros(p.n);require(m.mechanical(x,f),"initial state admissible");
    double verr=0;for(double v:m.vel)verr=std::max(verr,std::abs(v/p.Vp-1));
    require(verr<1e-8,"steady sliding friction solve");
    for(double N:{1e4,1e6,50e6})for(double v:{-1.,-1e-12,1e-15,1e-9,1e-3,1.,10.}) {
        double tau=N*friction(v,.72,.02,p)+p.radiation*v;
        require(std::abs(velocity(tau,.72,.02,N,p)/v-1)<1e-7,"friction inversion");
    }
    Vec slip(p.n);for(int i=0;i<p.n;i++)slip[i]=std::cos(7*PI*m.z[i]/p.Lz);
    m.elastic(slip);double elastic_error=0;
    for(int i=0;i<p.n;i++)elastic_error=std::max(elastic_error,std::abs((m.prestress[i]-m.tau[i])/m.stiffness[7]-slip[i]));
    require(elastic_error<1e-10,"finite-strip elastic eigenmode");
    State out=zeros(p.n),stage=zeros(p.n),f0=zeros(p.n),f1=zeros(p.n);
    require(m.step(x,1e5,out,stage,f0,f1),"stationary IMEX step");
    require(m.error(x,out)<1e-5,"steady state preservation");
    ARK4 solver(m);double estimate;
    require(solver.step(x,1e5,out,estimate),"stationary production ARK4 step");
    require(m.error(x,out)<1e-5,"production ARK4 steady state preservation");
    std::cout<<std::setprecision(12)<<"steady_flux_relative_error="<<qerr<<"\nsteady_velocity_relative_error="<<verr<<"\nelastic_mode_error="<<elastic_error<<"\n";
    // Constant-k diffusion: Dirichlet at the surface, no flux at the base.
    // The exact mode decays as exp(-D [pi/(2H)]^2 t).
    double previous=0;
    for(int n:{32,64,128}) {
        Parameters d=p;d.n=n;d.Lz=1000;d.Ly=1000;d.influx=0;d.kmin=1e-15;d.kmax=1e-15;
        d.normal_gradient=22000;Model md(d);Vec u(n,0),e(n),first(n),second(n),base(n);
        for(int i=0;i<n;i++)e[i]=1e6*std::sin(PI*md.z[i]/(2*d.Lz));
        double D=d.kmin/(d.viscosity*d.storage),tend=1e5,dt=2000.*32/n;
        double balance_error=0;
        for(double t=0;t<tend-.5*dt;t+=dt) {
            require(md.pressure(e,u,gamma_ars*dt,first),"diffusion stage 1");
            for(int i=0;i<n;i++)base[i]=e[i]+(1-gamma_ars)/gamma_ars*(first[i]-e[i]);
            require(md.pressure(base,u,gamma_ars*dt,second),"diffusion stage 2");
            double change=0;for(int i=0;i<n;i++)change+=(second[i]-e[i])*md.h*d.storage;
            double q1=d.kmin*first[0]*2/(d.viscosity*md.h),q2=d.kmin*second[0]*2/(d.viscosity*md.h);
            balance_error=std::max(balance_error,std::abs(change+dt*((1-gamma_ars)*q1+gamma_ars*q2)));
            e.swap(second);
        }
        double err=0;for(int i=0;i<n;i++)err=std::max(err,std::abs(e[i]/1e6-std::sin(PI*md.z[i]/(2*d.Lz))*std::exp(-D*std::pow(PI/(2*d.Lz),2)*tend)));
        std::cout<<"diffusion_n="<<n<<" relative_error="<<err<<" storage_balance_m="<<balance_error<<"\n";
        require(balance_error<1e-12,"discrete fluid mass conservation");
        if(previous>0)require(previous/err>3.5,"second-order diffusion convergence");
        previous=err;
    }
    // Constant velocity limit of k* and state: direct equations, independent exact solutions.
    double speed=1e-8,u0=.05,psi0=.75,duration=1e6;
    double rate=speed/p.L+1/p.T,ueq=(speed/p.L)/rate;
    double exact_u=ueq+(u0-ueq)*std::exp(-rate*duration);
    double theta0=p.dc/p.V0*std::exp((psi0-p.f0)/.02);
    double exact_theta=p.dc/speed+(theta0-p.dc/speed)*std::exp(-speed*duration/p.dc);
    double exact_psi=p.f0+.02*std::log(exact_theta*p.V0/p.dc);
    previous=0;
    for(int steps:{100,200,400}) {
        double u=u0,psi=psi0,dt=duration/steps;
        auto ratepsi=[&](double s){return .02/p.dc*(p.V0*std::exp((p.f0-s)/.02)-speed);};
        for(int j=0;j<steps;j++) {
            double fu=speed/p.L*(1-u)-u/p.T,fs=ratepsi(psi);
            double u1=u+gamma_ars*dt*fu,s1=psi+gamma_ars*dt*fs;
            u+=dt*(delta_ars*fu+(1-delta_ars)*(speed/p.L*(1-u1)-u1/p.T));
            psi+=dt*(delta_ars*fs+(1-delta_ars)*ratepsi(s1));
        }
        double err=std::max(std::abs(u-exact_u),std::abs(psi-exact_psi));
        std::cout<<"constitutive_steps="<<steps<<" error="<<err<<"\n";
        if(previous>0)require(previous/err>3.5,"second-order state and permeability convergence");
        previous=err;
    }
    std::cout<<"ALL EQUATION AND LIMIT TESTS PASSED\n";
}

void laws(const std::string& config,const std::string& directory) {
    Parameters p=read_config(config);Model m(p);State x=m.initial();
    std::filesystem::create_directories(directory);
    std::ofstream stress(directory+"/stress_law.csv"),heal(directory+"/healing_law.csv"),
        slip(directory+"/slip_law.csv"),profile(directory+"/steady_profile.csv");
    for(auto f:{&stress,&heal,&slip,&profile})*f<<std::setprecision(15);
    stress<<"effective_pa,kstar_1e15,kstar_1e16,kstar_1e17,kstar_1e18\n";
    heal<<"time_s,T_1e8,T_1e9,T_1e10\n";slip<<"slip_m,kstar_m2\n";
    profile<<"depth_m,permeability_m2,pressure_pa,effective_pa,normal_pa,hydrostatic_pa,lithostatic_pa,a,b,kstar_m2\n";
    for(int j=0;j<=400;j++) {
        double N=j*1e6;stress<<N;
        for(double ks:{1e-15,1e-16,1e-17,1e-18})stress<<","<<p.kmin+(ks-p.kmin)*std::exp(-N/p.stress_scale);
        stress<<"\n";
        double t=j*YEAR;heal<<t;
        for(double T:{1e8,1e9,1e10})heal<<","<<p.kmin+(p.kmax-p.kmin)*std::exp(-t/T);
        heal<<"\n";
        double displacement=j*.01;slip<<displacement<<","<<p.kmax-(p.kmax-p.kmin)*std::exp(-displacement/p.L)<<"\n";
    }
    for(int i=0;i<m.n && m.z[i]<=30000;i++) {
        double z=m.z[i];
        profile<<z<<","<<m.k_value(i,x[2][i],x[3][i])<<","<<x[3][i]+p.rho*p.gravity*z<<","<<m.normal(i,x[3][i])
            <<","<<p.normal_gradient*z<<","<<p.rho*p.gravity*z<<","<<2700*p.gravity*z<<","<<m.a[i]<<","<<m.b[i]<<","<<m.kstar(x[2][i])<<"\n";
    }
}

int main(int argc,char** argv) {
    try {
        if(argc==2&&std::string(argv[1])=="--test") {test();return 0;}
        if(argc==4&&std::string(argv[1])=="--laws") {laws(argv[2],argv[3]);return 0;}
        if(argc<3) {std::cerr<<"Usage: valving CONFIG OUTPUT_DIRECTORY [--resume]\n";return 2;}
        run(argv[1],argv[2],argc>3&&std::string(argv[3])=="--resume");return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
