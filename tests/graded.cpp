// Checks specific to Galerkin elasticity and conservative flow on a graded mesh.
#define main model_program_entry
#include "../src/model.cpp"
#undef main

int main() {
    try {
        std::cout<<std::setprecision(12);
        for(double surface_ratio:{1.,8.}) {
        Parameters p;p.n=32768;p.threads=1;p.fine_depth=30000;p.perturbation=0;p.surface_ratio=surface_ratio;
        Model m(p);State initial=m.initial();
        std::fill(m.prestress.begin(),m.prestress.end(),0.);
        Vec one(m.n,1),x(m.n),y(m.n);
        m.elastic(one);double constant=0;
        for(double tau:m.tau)constant=std::max(constant,std::abs(tau/m.stiffness[0]+1));
        require(constant<1e-10,"graded constant elastic mode");
        for(int i=0;i<m.n;i++) {
            x[i]=std::cos(7*PI*m.z[i]/p.Lz)+.2*std::sin(53*PI*m.z[i]/p.Lz);
            y[i]=std::cos(13*PI*m.z[i]/p.Lz)+.3*std::sin(71*PI*m.z[i]/p.Lz);
        }
        m.elastic(x);Vec kx=m.tau;m.elastic(y);Vec ky=m.tau;
        double xy=0,yx=0,energy=0,absolute=0;
        for(int i=0;i<m.n;i++) {
            xy+=m.volume[i]*x[i]*ky[i];yx+=m.volume[i]*y[i]*kx[i];
            energy-=m.volume[i]*x[i]*kx[i];absolute+=m.volume[i]*std::abs(x[i]*ky[i]);
        }
        double reciprocity=std::abs(xy-yx)/absolute;
        require(reciprocity<1e-12 && energy>0,"graded reciprocity and positive elastic energy");
        Parameters pu=p;pu.fine_depth=0;pu.surface_ratio=1;Model uniform(pu);Vec localized(m.n),ulocal(uniform.n);
        auto bump=[](double z){double u=(z-12000)/8000;return std::abs(u)<1?std::exp(1-1/(1-u*u)):0.;};
        for(int i=0;i<m.n;i++)localized[i]=bump(m.z[i]);
        for(int i=0;i<uniform.n;i++)ulocal[i]=bump(uniform.z[i]);
        m.elastic(localized);uniform.elastic(ulocal);double elastic=0,scale=0;
        for(int i=0;i<m.n && m.z[i]<25000;i++) {
            if(m.z[i]<3000)continue;
            elastic=std::max(elastic,std::abs(m.tau[i]-uniform.tau[m.nodes[i]]));
            scale=std::max(scale,std::abs(uniform.tau[m.nodes[i]]));
        }
        require(elastic/scale<1e-7,"upper-fault localized elastic loading agrees with uniform mesh");
        State perturbed=initial;Vec pressure(m.n);
        for(int i=0;i<m.n;i++)perturbed[3][i]+=1e5*std::sin(PI*m.z[i]/p.Lz);
        double dt=2e5;require(m.pressure(perturbed[3],perturbed[2],dt,pressure),"graded nonlinear pressure solve");
        State final=perturbed;final[3]=pressure;Vec q=m.flux(final);
        double stored=0;
        for(int i=0;i<m.n;i++)stored+=p.storage*m.volume[i]*(pressure[i]-perturbed[3][i]);
        double budget=std::abs(stored-dt*(p.influx-q[0]));
        require(budget<1e-10,"graded nonlinear fluid mass balance");
        std::cout<<"surface_ratio="<<surface_ratio<<" graded_nodes="<<m.n<<" fourier_nodes="<<p.n
            <<"\nconstant_mode_error="<<constant<<" reciprocity_error="<<reciprocity
            <<"\nlocalized_elastic_error="<<elastic/scale<<" storage_residual_m="<<budget<<"\n";
        double previous=0;
        for(int n:{2048,4096,8192,16384}) {
            Parameters s=p;s.Lz=30000;s.n=n;s.fine_depth=3000;s.stretch_scale=5000;
            s.maximum_spacing=16*s.Lz/n;s.kmin=0;
            Model ms(s);State state=ms.initial();double error=0;
            double G=s.normal_gradient-s.rho*s.gravity;
            double ratio=s.viscosity*s.influx/(G*s.kmax*state[2][0]);
            for(int i=0;i<ms.n;i++) {
                double exact=-s.stress_scale*std::log(ratio+(1-ratio)*std::exp(-G*ms.z[i]/s.stress_scale));
                error=std::max(error,std::abs(ms.normal(i,state[3][i])-exact)/s.stress_scale);
            }
            std::cout<<"graded_steady_fourier_n="<<n<<" dynamic_nodes="<<ms.n<<" error="<<error<<"\n";
            if(previous>0)require(previous/error>3.5,"graded steady-flow spatial convergence");
            previous=error;
        }
        }
        double previous=0;
        for(int n:{256,512,1024}) {
            Parameters d;d.n=n;d.Lz=1000;d.Ly=1000;d.threads=1;d.influx=0;d.kmin=d.kmax=1e-15;
            d.fine_depth=300;d.stretch_scale=200;d.maximum_spacing=16*d.Lz/n;
            d.surface_ratio=8;d.surface_scale=100;
            Model m(d);Vec u(m.n,0),e(m.n),first(m.n),second(m.n),base(m.n);
            for(int i=0;i<m.n;i++)e[i]=1e6*std::sin(PI*m.z[i]/(2*d.Lz));
            double D=d.kmin/(d.viscosity*d.storage),end=1e5,dt=1000.*256/n;
            for(double t=0;t<end-.5*dt;t+=dt) {
                require(m.pressure(e,u,gamma_ars*dt,first),"graded diffusion stage 1");
                for(int i=0;i<m.n;i++)base[i]=e[i]+(1-gamma_ars)/gamma_ars*(first[i]-e[i]);
                require(m.pressure(base,u,gamma_ars*dt,second),"graded diffusion stage 2");e.swap(second);
            }
            double error=0;
            for(int i=0;i<m.n;i++)error=std::max(error,std::abs(e[i]/1e6-
                std::sin(PI*m.z[i]/(2*d.Lz))*std::exp(-D*std::pow(PI/(2*d.Lz),2)*end)));
            std::cout<<"graded_diffusion_fourier_n="<<n<<" dynamic_nodes="<<m.n<<" error="<<error<<"\n";
            if(previous>0)require(previous/error>3.5,"graded diffusion convergence");
            previous=error;
        }
        std::cout<<"GRADED MESH TESTS PASSED\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
