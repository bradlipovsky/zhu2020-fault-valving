// Zhu et al. (2020), equations (1)--(9), all quantities in SI units.
// Exact spectral elastic Dirichlet-to-Neumann map; conservative hydraulic FV;
// implicit midpoint + step doubling for adaptive second-order coupled evolution.
#include <fftw3.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using Vec = std::vector<double>;
constexpr double pi=3.14159265358979323846, year=365.0*86400;
constexpr double mu=32.4e9, damping=4.68e6, vp=1e-9, v0=1e-6;
constexpr double f0=.6, dc=.002, storage=1e-11, viscosity=1e-4;
constexpr double rhog=9800, normal_gradient=22050, stress_scale=30e6;
constexpr double kfloor=1e-19, L=1;

double interp(double x, const Vec& xp, const Vec& fp) {
    if(x<=xp.front()) return fp.front();
    if(x>=xp.back()) return fp.back();
    size_t i=std::upper_bound(xp.begin(),xp.end(),x)-xp.begin()-1;
    return fp[i]+(fp[i+1]-fp[i])*(x-xp[i])/(xp[i+1]-xp[i]);
}
// log(asinh(exp(x))) is not needed: evaluate asinh(exp(x)) without overflow.
double asinh_exp(double x) { return x>30 ? x+std::log(2.0) : std::asinh(std::exp(x)); }
double friction(double v,double psi,double a) {
    if(v==0) return 0;
    return std::copysign(a*asinh_exp(std::log(std::abs(v)/(2*v0))+psi/a),v);
}
double velocity(double tau,double psi,double a,double effective) {
    if(tau==0) return 0;
    double target=std::abs(tau), lo=-745, hi=std::log(target/damping);
    double wstatic=target/(effective*a);
    double logsinh=wstatic>30 ? wstatic-std::log(2.0) : std::log(std::sinh(wstatic));
    double x=std::min(hi,std::log(2*v0)-psi/a+logsinh);
    x=std::max(lo,x);
    for(int it=0;it<60;++it) {
        double v=std::exp(x), w=x-std::log(2*v0)+psi/a;
        double F=damping*v+effective*a*asinh_exp(w)-target;
        if(std::abs(F)<1e-8+target*2e-13) return std::copysign(v,tau);
        if(F>0) hi=x; else lo=x;
        double slope=damping*v+effective*a/(w>30 ? 1.0 : std::sqrt(1+std::exp(-2*w)));
        double nx=x-F/slope;
        x=(std::isfinite(nx) && nx>lo && nx<hi) ? nx : (lo+hi)/2;
    }
    throw std::runtime_error("friction solve did not converge");
}

struct State { Vec d,psi,k,h; explicit State(size_t n):d(n),psi(n),k(n),h(n){} };
// k in State is k-star / 1e-15; h = p-rho*g*z; d = slip - Vp*t.
struct Model {
    int n; double H,W,dz,T,q0; bool coupled;
    Vec z,a,b,sigma,klo,khi,ss_h,tau0,buf,eigen,kwork,lower,diag,upper,rhs,old;
    std::vector<State> rk;
    fftw_plan transform;
    Model(int intervals,double healing,double height,double width,bool feedback):
      n(intervals+1),H(height),W(width),dz(H/intervals),T(healing),
      q0(healing==1e7?3.3764e-10:3e-9),coupled(feedback),
      z(n),a(n),b(n),sigma(n),klo(n),khi(n),ss_h(n),tau0(n),buf(n),eigen(n),
      kwork(n),lower(n),diag(n),upper(n),rhs(n),old(n) {
        // Original input files taper both slip-dependent bounds to a constant
        // from 30 to 60 km; pressure-dependent floor remains 1e-19 m^2.
        double kdeep=T==1e7?1e-17:(T==1e8?9.1e-17:(T==1e9?5e-16:9.091e-16));
        for(int i=0;i<n;++i) {
            z[i]=i*dz;
            a[i]=interp(z[i],{0,14900,27600,60000},{.0105,.03,.07,.173});
            b[i]=interp(z[i],{0,13600,14900,27600,60000},{.02,.0378,.0356,.0375,.0424});
            sigma[i]=normal_gradient*z[i];
            klo[i]=interp(z[i],{0,30000,60000,H},{kfloor,kfloor,kdeep,kdeep})/1e-15;
            khi[i]=interp(z[i],{0,30000,60000,H},{1e-15,1e-15,kdeep,kdeep})/1e-15;
            double w=pi*i/H;
            eigen[i]=i==0 ? -mu/(2*W) : -mu*w/(2*std::tanh(w*W));
        }
        transform=fftw_plan_r2r_1d(n,buf.data(),buf.data(),FFTW_REDFT00,FFTW_ESTIMATE);
        if(!transform) throw std::runtime_error("FFTW plan failed");
        for(int j=0;j<7;++j)rk.emplace_back(n);
    }
    ~Model(){fftw_destroy_plan(transform);}
    void elastic(const Vec& d,Vec& tau) {
        std::copy(d.begin(),d.end(),buf.begin()); fftw_execute(transform);
        for(int i=0;i<n;++i) buf[i]*=eigen[i]/(2*(n-1));
        fftw_execute(transform);tau=buf;
    }
    double perm(int i,double h,double ks) const {
        return kfloor+(ks*1e-15-kfloor)*std::exp(-(sigma[i]-rhog*z[i]-h)/stress_scale);
    }
    void permeability(const State& s,Vec& out) const {
        for(int i=0;i<n;++i) out[i]=perm(i,s.h[i],s.k[i]);
    }
    void initialize(State& s,double perturb) {
        for(int i=0;i<n;++i) {
            s.k[i]=(vp/L*khi[i]+klo[i]/T)/(vp/L+1/T);
            s.psi[i]=f0-b[i]*std::log(vp/v0);
        }
        // Steady discrete Darcy flux: each face satisfies q=q0 exactly.
        for(int i=1;i<n;++i) {
            double left=perm(i-1,s.h[i-1],s.k[i-1]);
            double low=s.h[i-1],high=low+q0*viscosity*dz/kfloor;
            for(int j=0;j<90;++j) {
                double mid=(low+high)/2;
                double k=.5*(left+perm(i,mid,s.k[i]));
                if(k*(mid-s.h[i-1])/dz/viscosity>q0) high=mid; else low=mid;
            }
            s.h[i]=(low+high)/2;
        }
        ss_h=s.h;
        for(int i=0;i<n;++i) {
            double eff=std::max(1e6,sigma[i]-rhog*z[i]-s.h[i]);
            tau0[i]=eff*friction(vp,s.psi[i],a[i])+damping*vp;
            // Specified small state perturbation starts the unstable solution.
            s.psi[i]-=perturb*std::exp(-std::pow((z[i]-14000)/1000,2));
        }
    }
    void rates(const State& s,State& r,Vec& v) {
        elastic(s.d,v);
        for(int i=0;i<n;++i) {
            double eff=std::max(1e6,sigma[i]-rhog*z[i]-s.h[i]);
            v[i]=velocity(tau0[i]+v[i],s.psi[i],a[i],eff);
            // The paper assumes a fixed slip direction, V>=0. Abort if violated.
            if(v[i]<0) throw std::runtime_error("negative slip velocity outside model assumptions");
            r.d[i]=v[i]-vp;
            r.psi[i]=b[i]/dc*(v0*std::exp((f0-s.psi[i])/b[i])-v[i]);
            r.k[i]=v[i]/L*(khi[i]-s.k[i])-(s.k[i]-klo[i])/T;
        }
    }
    void load_profile(State& s,const std::string& path) {
        std::ifstream in(path);if(!in)throw std::runtime_error("cannot read initial profile");
        std::string line;std::getline(in,line);
        Vec zz,ps,ks,pp,tt;
        while(std::getline(in,line)) {
            std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
            double z1,p1,k1,h1,t1;if(!(row>>z1>>p1>>k1>>h1>>t1))throw std::runtime_error("invalid profile row");
            zz.push_back(z1);ps.push_back(p1);ks.push_back(k1);pp.push_back(h1);tt.push_back(t1);
        }
        for(int i=0;i<n;++i) {
            s.d[i]=0;s.psi[i]=interp(z[i],zz,ps);s.k[i]=interp(z[i],zz,ks)/1e-15;
            s.h[i]=interp(z[i],zz,pp)-rhog*z[i];tau0[i]=interp(z[i],zz,tt);
        }
        ss_h=s.h;
    }
    bool all_rates(const State& s,State& r,Vec& v,double dt) {
        rates(s,r,v);
        std::fill(r.h.begin(),r.h.end(),0);
        if(!coupled)return true;
        double largest=0;
        for(int i=0;i<n;++i) {
            kwork[i]=perm(i,s.h[i],s.k[i]);
            if(kwork[i]<=0 || !std::isfinite(kwork[i]))return false;
            largest=std::max(largest,kwork[i]);
        }
        if(dt*largest/(storage*viscosity*dz*dz)>.4)return false;
        double previous=.5*(kwork[0]+kwork[1])*(s.h[1]-s.h[0])/(viscosity*dz);
        for(int i=1;i<n-1;++i) {
            double next=.5*(kwork[i]+kwork[i+1])*(s.h[i+1]-s.h[i])/(viscosity*dz);
            r.h[i]=(next-previous)/(storage*dz);previous=next;
        }
        r.h[n-1]=2*(q0-previous)/(storage*dz);
        return true;
    }
    // Dormand--Prince 5(4) for subsecond steps: all four equations are explicit.
    // The hydraulic CFL condition is checked at every stage. This avoids the
    // cost of repeated implicit solves during effectively undrained ruptures.
    bool dopri(const State& s,double dt,State& high,State& low,State& stage,Vec& v) {
        static const double A[7][7]={
            {0},{1./5},{3./40,9./40},{44./45,-56./15,32./9},
            {19372./6561,-25360./2187,64448./6561,-212./729},
            {9017./3168,-355./33,46732./5247,49./176,-5103./18656},
            {35./384,0,500./1113,125./192,-2187./6784,11./84}};
        static const double B[7]={5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};
        for(int j=0;j<7;++j) {
            stage=s;
            for(int l=0;l<j;++l)if(A[j][l]!=0)for(int i=0;i<n;++i) {
                stage.d[i]+=dt*A[j][l]*rk[l].d[i];stage.psi[i]+=dt*A[j][l]*rk[l].psi[i];
                stage.k[i]+=dt*A[j][l]*rk[l].k[i];stage.h[i]+=dt*A[j][l]*rk[l].h[i];
            }
            for(int i=0;i<n;++i)if(stage.k[i]<=0)return false;
            if(!all_rates(stage,rk[j],v,dt))return false;
        }
        high=stage;low=s;
        for(int j=0;j<7;++j)if(B[j]!=0)for(int i=0;i<n;++i) {
            low.d[i]+=dt*B[j]*rk[j].d[i];low.psi[i]+=dt*B[j]*rk[j].psi[i];
            low.k[i]+=dt*B[j]*rk[j].k[i];low.h[i]+=dt*B[j]*rk[j].h[i];
        }
        return true;
    }
    // Nonlinear backward-Euler hydraulic stage with conservative face fluxes.
    bool pressure(const Vec& h0,const Vec& ks,double dt,Vec& h) {
        if(!coupled){h=ss_h;return true;}
        h=h0;
        for(int it=0;it<80;++it) {
            old=h;
            for(int i=0;i<n;++i) {
                kwork[i]=perm(i,h[i],ks[i]);
                if(!std::isfinite(kwork[i]) || kwork[i]<=0) return false;
            }
            diag[0]=1; upper[0]=0;rhs[0]=0;
            for(int i=1;i<n;++i) {
                double c=dt/(storage*viscosity*dz*dz);
                double kl=.5*(kwork[i-1]+kwork[i]), kr=i<n-1?.5*(kwork[i]+kwork[i+1]):0;
                if(i==n-1) c*=2; // half control volume at flux boundary
                lower[i]=-c*kl; upper[i]=-c*kr;diag[i]=1+c*(kl+kr);
                rhs[i]=h0[i]+(i==n-1 ? 2*dt*q0/(storage*dz) : 0);
            }
            for(int i=1;i<n;++i) {
                double m=lower[i]/diag[i-1];diag[i]-=m*upper[i-1];rhs[i]-=m*rhs[i-1];
            }
            h[n-1]=rhs[n-1]/diag[n-1];
            double err=0;
            for(int i=n-2;i>=0;--i) h[i]=(rhs[i]-upper[i]*h[i+1])/diag[i];
            for(int i=0;i<n;++i) err=std::max(err,std::abs(h[i]-old[i]));
            if(err<.01) return true; // Pa, far below time integration tolerance
        }
        return false;
    }
    bool midpoint(const State& s,double dt,State& out,State& mid,State& r,Vec& v) {
        rates(s,r,v);
        for(int i=0;i<n;++i) {
            mid.d[i]=s.d[i]+dt/2*r.d[i]; mid.psi[i]=s.psi[i]+dt/2*r.psi[i];
            mid.k[i]=s.k[i]+dt/2*r.k[i];
            if(mid.k[i]<=0 || !std::isfinite(mid.psi[i])) return false;
        }
        if(!pressure(s.h,mid.k,dt/2,mid.h)) return false;
        rates(mid,r,v);
        for(int i=0;i<n;++i) {
            out.d[i]=s.d[i]+dt*r.d[i];out.psi[i]=s.psi[i]+dt*r.psi[i];
            out.k[i]=s.k[i]+dt*r.k[i];out.h[i]=coupled?2*mid.h[i]-s.h[i]:ss_h[i];
            if(out.k[i]<=0 || !std::isfinite(out.psi[i])) return false;
        }
        return true;
    }
    // Fully implicit midpoint removes the shallow elastic relaxation CFL limit.
    // Local aging-law states are eliminated analytically from the Newton matrix;
    // a preconditioned CG solve applies the elastic Jacobian through FFTs.
    bool implicit_midpoint(const State& s,double dt,State& out,State& mid,State& r,Vec& v) {
        const double h=dt/2;
        rates(s,r,v);
        Vec current=v, tau(n), D(n), residual(n), delta(n), cg_r(n), cg_z(n), direction(n), action(n);
        for(int newton=0;newton<14;++newton) {
            for(int i=0;i<n;++i) {
                mid.d[i]=s.d[i]+h*(current[i]-vp);
                double ps=s.psi[i];bool state_converged=false;
                for(int j=0;j<30;++j) {
                    double vv=v0*std::exp((f0-ps)/b[i]);
                    double step=(ps-s.psi[i]-h*b[i]/dc*(vv-current[i]))/(1+h*vv/dc);
                    ps-=step;
                    if(!std::isfinite(ps))return false;
                    if(std::abs(step)<1e-13){state_converged=true;break;}
                }
                if(!state_converged)return false;
                mid.psi[i]=ps;
                mid.k[i]=(s.k[i]+h*(current[i]*khi[i]/L+klo[i]/T))/(1+h*(current[i]/L+1/T));
            }
            if(!pressure(s.h,mid.k,h,mid.h))return false;
            elastic(mid.d,tau);
            double norm=0;
            for(int i=0;i<n;++i) {
                double eff=std::max(1e6,sigma[i]-rhog*z[i]-mid.h[i]);
                residual[i]=tau0[i]+tau[i]-damping*current[i]-eff*friction(current[i],mid.psi[i],a[i]);
                norm=std::max(norm,std::abs(residual[i]));
                double w=std::log(current[i]/(2*v0))+mid.psi[i]/a[i];
                double factor=w>30?1:1/std::sqrt(1+std::exp(-2*w));
                double psiv=-h*b[i]/dc/(1+h*v0/dc*std::exp((f0-mid.psi[i])/b[i]));
                D[i]=damping+eff*factor*(a[i]/current[i]+psiv);
                if(D[i]<=0 || !std::isfinite(D[i]))return false;
            }
            if(norm<1e-4) {
                for(int i=0;i<n;++i) {
                    out.d[i]=2*mid.d[i]-s.d[i];out.psi[i]=2*mid.psi[i]-s.psi[i];
                    out.k[i]=2*mid.k[i]-s.k[i];out.h[i]=coupled?2*mid.h[i]-s.h[i]:ss_h[i];
                    if(out.k[i]<=0)return false;
                }
                v=current;return true;
            }
            double rz=0;
            for(int i=0;i<n;++i) {
                delta[i]=0;cg_r[i]=residual[i];cg_z[i]=cg_r[i]/(D[i]+h*mu/dz);
                direction[i]=cg_z[i];rz+=(i==0 || i==n-1?.5:1)*cg_r[i]*cg_z[i];
            }
            bool solved=false;
            for(int cg=0;cg<100;++cg) {
                elastic(direction,action);double pap=0;
                for(int i=0;i<n;++i){action[i]=D[i]*direction[i]-h*action[i];pap+=(i==0 || i==n-1?.5:1)*direction[i]*action[i];}
                if(pap<=0 || !std::isfinite(pap))return false;
                double alpha=rz/pap,maxr=0;
                for(int i=0;i<n;++i){delta[i]+=alpha*direction[i];cg_r[i]-=alpha*action[i];maxr=std::max(maxr,std::abs(cg_r[i]));}
                if(maxr<std::max(1e-6,norm*1e-7)){solved=true;break;}
                double rznew=0;
                for(int i=0;i<n;++i){cg_z[i]=cg_r[i]/(D[i]+h*mu/dz);rznew+=(i==0 || i==n-1?.5:1)*cg_r[i]*cg_z[i];}
                double beta=rznew/rz;
                for(int i=0;i<n;++i)direction[i]=cg_z[i]+beta*direction[i];
                rz=rznew;
            }
            if(!solved)return false;
            double fraction=1;
            for(int i=0;i<n;++i)if(delta[i]<0)fraction=std::min(fraction,-.8*current[i]/delta[i]);
            for(int i=0;i<n;++i)current[i]+=fraction*delta[i];
        }
        return false;
    }
};

void require(bool pass,const std::string& msg) { if(!pass)throw std::runtime_error("TEST FAILED: "+msg); }
void tests() {
    Model m(1024,1e8,500e3,500e3,true);State s(m.n),r(m.n);m.initialize(s,0);
    Vec d(m.n),tau(m.n),v(m.n);double err=0;
    for(int i=0;i<m.n;++i)d[i]=std::cos(7*pi*m.z[i]/m.H);
    m.elastic(d,tau);
    for(int i=0;i<m.n;++i)err=std::max(err,std::abs(tau[i]-m.eigen[7]*d[i])/std::abs(m.eigen[7]));
    require(err<1e-11,"elastic cosine mode");std::cout<<"elastic_mode_relative_error "<<err<<'\n';
    std::fill(d.begin(),d.end(),2.0);m.elastic(d,tau);
    require(std::abs(tau[12]+mu/m.W)<1e-7,"elastic uniform loading");
    err=0;
    for(double vv: {1e-20,1e-9,1e-5,1.0,10.0}) {
        double tt=damping*vv+50e6*friction(vv,.6,.015);
        err=std::max(err,std::abs(velocity(tt,.6,.015,50e6)/vv-1));
    }
    require(err<1e-9,"friction inversion");std::cout<<"friction_relative_error "<<err<<'\n';
    Vec h(m.n);require(m.pressure(s.h,s.k,1e7,h),"steady hydraulic solver");err=0;
    for(int i=0;i<m.n;++i)err=std::max(err,std::abs(h[i]-s.h[i]));
    require(err<.01,"steady flux preservation");std::cout<<"steady_pressure_error_Pa "<<err<<'\n';
    // Constant k manufactured eigenmode, Dirichlet top and zero-flux bottom.
    m.q0=0; std::fill(s.k.begin(),s.k.end(),kfloor/1e-15);
    double w=pi/(2*m.H),dt=1e8,D=kfloor/(storage*viscosity);
    for(int i=0;i<m.n;++i)s.h[i]=1e6*std::sin(w*m.z[i]);
    require(m.pressure(s.h,s.k,dt,h),"diffusion solver");
    double eig=4*D/(m.dz*m.dz)*std::pow(std::sin(w*m.dz/2),2);err=0;
    for(int i=0;i<m.n;++i)err=std::max(err,std::abs(h[i]-s.h[i]/(1+dt*eig)));
    require(err<1e-6,"manufactured backward Euler eigenmode");
    double mass=0;for(int i=1;i<m.n;++i)mass+=(h[i]-s.h[i])*m.dz*(i==m.n-1?.5:1)*storage;
    double surface=kfloor/viscosity*(h[1]-h[0])/m.dz;
    require(std::abs(mass+dt*surface)<1e-10,"discrete fluid conservation");
    std::cout<<"diffusion_error_Pa "<<err<<"\nmass_balance_error_m "<<std::abs(mass+dt*surface)<<'\n';
    // Uniform spring-slider limit: the spectral zero mode gives stiffness mu/2W.
    // Compare implicit midpoint integrations to a refined numerical reference.
    Model q(64,1e8,500e3,500e3,true);State init(q.n);q.initialize(init,0);q.q0=0;
    for(int i=0;i<q.n;++i) {
        q.sigma[i]=50e6+rhog*q.z[i];q.a[i]=.015;q.b[i]=.02;
        init.h[i]=0;init.k[i]=kfloor/1e-15;q.klo[i]=q.khi[i]=init.k[i];
        init.psi[i]=f0-q.b[i]*std::log(vp/v0)-.001;
        q.tau0[i]=50e6*friction(vp,init.psi[i]+.001,q.a[i])+damping*vp;
    }
    auto evolve=[&](int count) {
        State x=init,y(q.n),mid(q.n),rr(q.n);Vec vv(q.n);
        for(int j=0;j<count;++j) {
            require(q.implicit_midpoint(x,1e5/count,y,mid,rr,vv),"implicit midpoint convergence");
            std::swap(x,y);
        }
        return x.d[0];
    };
    double exact=evolve(64),e1=std::abs(evolve(1)-exact),e2=std::abs(evolve(2)-exact);
    require(e1/e2>3.7 && e1/e2<4.3,"midpoint second-order time convergence");
    std::cout<<"midpoint_error_ratio "<<e1/e2<<"\nPASS\n";
    auto evolve_rk=[&](int count) {
        State x=init,y(q.n),low(q.n),stage(q.n);Vec vv(q.n);
        for(int j=0;j<count;++j) {
            require(q.dopri(x,1e6/count,y,low,stage,vv),"Dormand-Prince stages");
            std::swap(x,y);
        }
        return x.d[0];
    };
    exact=evolve_rk(64);e1=std::abs(evolve_rk(1)-exact);e2=std::abs(evolve_rk(2)-exact);
    require(e1/e2>20 && e1/e2<50,"fifth-order time convergence");
    std::cout<<"dopri_error_ratio "<<e1/e2<<"\nPASS\n";
}

int main(int argc,char**argv) try {
    int intervals=16384,stride=10;double T=1e8,years=120,H=578952.681034,rtol=1e-4,maxdt=3e5,perturb=1e-4;
    bool coupled=true,append=false;double start_year=0;
    std::string output="results/coupled",profile="",state_file="";
    for(int i=1;i<argc;++i) {
        std::string a=argv[i];
        if(a=="--test"){tests();return 0;}
        if(a=="--fixed"){coupled=false;continue;}
        if(a=="--append"){append=true;continue;}
        if(i+1>=argc)throw std::runtime_error("missing value for "+a);
        std::string v=argv[++i];
        if(a=="--n")intervals=std::stoi(v);else if(a=="--T")T=std::stod(v);
        else if(a=="--years")years=std::stod(v);else if(a=="--height")H=std::stod(v);
        else if(a=="--rtol")rtol=std::stod(v);else if(a=="--max-dt")maxdt=std::stod(v);
        else if(a=="--perturb")perturb=std::stod(v);else if(a=="--output")output=v;
        else if(a=="--profile")profile=v;
        else if(a=="--state")state_file=v;else if(a=="--start-year")start_year=std::stod(v);
        else if(a=="--stride")stride=std::stoi(v);else throw std::runtime_error("unknown argument "+a);
    }
    if(intervals<16 || rtol<=0 || years<=0 || stride<1 || H<=60000 || maxdt<=0)
        throw std::runtime_error("invalid configuration");
    if(T!=1e7 && T!=1e8 && T!=1e9 && T!=1e10)
        throw std::runtime_error("T must be 1e7, 1e8, 1e9, or 1e10 s for the supplied depth profiles");
    if(start_year<0 || start_year>=years || (start_year>0 && state_file.empty()) || (append && state_file.empty()))
        throw std::runtime_error("invalid continuation arguments");
    std::filesystem::create_directories(output);
    Model m(intervals,T,H,500e3,coupled);State s(m.n),full(m.n),half(m.n),fine(m.n),mid(m.n),r(m.n);
    Vec v(m.n),k(m.n);m.initialize(s,perturb);
    if(!profile.empty())m.load_profile(s,profile);
    if(!state_file.empty()) {
        std::ifstream state(state_file,std::ios::binary);
        for(auto p:{&s.d,&s.psi,&s.k,&s.h})state.read(reinterpret_cast<char*>(p->data()),m.n*8);
        if(!state || state.peek()!=EOF)throw std::runtime_error("checkpoint size does not match grid");
    }
    // Archive top 30 km on the computational nodes. Float64 time and fields.
    int nz=std::upper_bound(m.z.begin(),m.z.end(),30000)-m.z.begin();
    if(append) {
        std::ifstream previous(output+"/fields.bin",std::ios::binary|std::ios::ate);
        auto bytes=previous.tellg();std::streamoff width=8*(1+7*nz);
        if(bytes<width || bytes%width!=0)throw std::runtime_error("invalid existing output for append");
        previous.seekg(-width,std::ios::end);double last_time;previous.read(reinterpret_cast<char*>(&last_time),8);
        if(std::abs(last_time-start_year*year)>1e-5)throw std::runtime_error("continuation time does not match existing output");
    }
    std::ofstream meta(output+"/metadata.json");
    meta<<std::setprecision(17)<<"{\n\"n\": "<<m.n<<", \"nz\": "<<nz<<", \"dz\": "<<m.dz
        <<", \"height_m\": "<<H<<", \"width_m\": "<<m.W<<", \"T_s\": "<<T<<", \"q0_m_s\": "<<m.q0
        <<", \"rtol\": "<<rtol<<", \"max_dt_s\": "<<maxdt<<", \"years\": "<<years
        <<", \"coupled\": "<<(coupled?"true":"false")<<", \"initial_state_perturbation\": "<<perturb
        <<", \"continuation_start_year\": "<<start_year
        <<", \"initial_profile\": \""<<profile<<"\", \"fields\": [\"slip_m\",\"velocity_m_s\",\"effective_stress_Pa\",\"permeability_m2\",\"flux_m_s\",\"psi\",\"kstar_m2\"]}\n";meta.close();
    auto mode=std::ios::out|(append?std::ios::app:std::ios::trunc);
    std::ofstream bin(output+"/fields.bin",mode|std::ios::binary),history(output+"/history.csv",mode);
    if(!append)history<<"time_s,dt_s,max_velocity_m_s,accepted,rejected,error\n";
    history<<std::setprecision(15);
    auto write=[&](double t) {
        m.rates(s,r,v);m.permeability(s,k);
        bin.write(reinterpret_cast<char*>(&t),8);
        for(int f=0;f<7;++f)for(int j=0;j<nz;++j) {
            double val=0;
            if(f==0)val=s.d[j]+vp*t;
            if(f==1)val=v[j];
            if(f==2)val=std::max(1e6,m.sigma[j]-rhog*m.z[j]-s.h[j]);
            if(f==3)val=k[j];
            if(f==4)val=.5*(k[j]+k[j+1])*(s.h[j+1]-s.h[j])/m.dz/viscosity;
            if(f==5)val=s.psi[j];
            if(f==6)val=s.k[j]*1e-15;
            bin.write(reinterpret_cast<char*>(&val),8);
        }
    };
    double t=start_year*year,dt=100,stop=years*year,last_report=-1,last_save=t;long accepted=0,rejected=0;
    auto start=std::chrono::steady_clock::now();if(!append)write(t);
    while(t<stop) {
        dt=std::min({dt,maxdt,stop-t});
        bool ok=false,use_rk=dt<1;double divisor=use_rk?1:3,exponent=use_rk?-.2:-1./3;
        try {
            if(use_rk) {
                ok=m.dopri(s,dt,fine,full,mid,v);
            } else if(dt>3000) {
                ok=m.implicit_midpoint(s,dt,full,mid,r,v) && m.implicit_midpoint(s,dt/2,half,mid,r,v)
                     && m.implicit_midpoint(half,dt/2,fine,mid,r,v);
            } else {
                ok=m.midpoint(s,dt,full,mid,r,v) && m.midpoint(s,dt/2,half,mid,r,v)
                     && m.midpoint(half,dt/2,fine,mid,r,v);
            }
        } catch(const std::exception&) {ok=false;}
        double error=0;
        if(ok)for(int i=0;i<m.n;++i) {
            if(!std::isfinite(fine.d[i]+fine.psi[i]+fine.k[i]+fine.h[i]) ||
               !std::isfinite(full.d[i]+full.psi[i]+full.k[i]+full.h[i])) {
                error=INFINITY;break;
            }
            error=std::max(error,std::abs(fine.d[i]-full.d[i])/(divisor*rtol*dc));
            error=std::max(error,std::abs(fine.psi[i]-full.psi[i])/(divisor*rtol*m.b[i]));
            error=std::max(error,std::abs(fine.k[i]-full.k[i])/(divisor*rtol*(1e-3+std::abs(fine.k[i]))));
            error=std::max(error,std::abs(fine.h[i]-full.h[i])/(divisor*rtol*stress_scale));
        }
        if(!ok || !std::isfinite(error) || error>1) {
            ++rejected;dt*=ok?std::max(.1,.8*std::pow(error,exponent)):.25;
            if(dt<1e-10)throw std::runtime_error("time step underflow at t="+std::to_string(t));
            continue;
        }
        std::swap(s,fine);t+=dt;++accepted;
        bool save=accepted%stride==0 || t-last_save>=.02*year || t>=stop;
        if(save){write(t);last_save=t;}
        if(accepted%100==0 || t>=stop)history<<t<<','<<dt<<','<<*std::max_element(v.begin(),v.end())<<','<<accepted<<','<<rejected<<','<<error<<'\n';
        double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(elapsed-last_report>=30 || t>=stop) {
            std::cout<<std::setprecision(7)<<"t_year="<<t/year<<" dt="<<dt<<" Vmax="<<*std::max_element(v.begin(),v.end())
                     <<" accepted="<<accepted<<" rejected="<<rejected<<" wall_s="<<elapsed<<std::endl;
            last_report=elapsed;bin.flush();history.flush();
        }
        dt*=std::clamp(.9*std::pow(std::max(error,1e-8),exponent),.5,2.0);
    }
    std::ofstream checkpoint(output+"/final_state.bin",std::ios::binary);
    for(auto p:{&s.d,&s.psi,&s.k,&s.h})checkpoint.write(reinterpret_cast<const char*>(p->data()),m.n*8);
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
