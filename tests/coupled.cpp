// Additional verification invokes the actual model operators, not copies of them.
#define main model_program_entry
#include "../src/model.cpp"
#undef main

int main() {
    try {
        Parameters p;p.n=512;p.Lz=30000;p.Ly=500000;p.threads=1;
        p.initial_locking_depth=17000;p.perturbation=.001;
        Model m(p);State initial=m.initial();
        for(int i=0;i<p.n;i++) {
            initial[3][i]+=5e5*std::sin(PI*m.z[i]/p.Lz);
            initial[2][i]*=1+.1*std::sin(PI*m.z[i]/p.Lz);
        }
        Vec solved(p.n);double dt=1e6;
        require(m.pressure(initial[3],initial[2],dt,solved),"nonlinear hydraulic solve");
        State next=initial;next[3]=solved;Vec q=m.flux(next);
        double change=0;
        for(int i=0;i<p.n;i++)change+=(solved[i]-initial[3][i])*p.storage*m.h;
        double residual=std::abs(change-dt*(p.influx-q[0]));
        require(residual<1e-9,"nonlinear Darcy storage budget");
        std::cout<<std::setprecision(14)<<"nonlinear_storage_residual_m="<<residual<<"\n";
        auto integrate=[&](double timestep){
            State x=initial,out=zeros(p.n),stage=zeros(p.n),f0=zeros(p.n),f1=zeros(p.n);
            for(double t=0;t<1e6-.5*timestep;t+=timestep) {
                require(m.step(x,timestep,out,stage,f0,f1),"coupled time integration");x.swap(out);
            }
            return x;
        };
        State exact=integrate(2500);double prior=0;
        for(double timestep:{20000.,10000.,5000.}) {
            State result=integrate(timestep);double error=0;
            for(int i=0;i<p.n;i++) {
                error=std::max(error,std::abs(result[0][i]-exact[0][i])/p.dc);
                error=std::max(error,std::abs(result[1][i]-exact[1][i])/m.a[i]);
                error=std::max(error,std::abs(result[2][i]-exact[2][i]));
                error=std::max(error,std::abs(result[3][i]-exact[3][i])/p.stress_scale);
            }
            std::cout<<"coupled_dt_s="<<timestep<<" scaled_error="<<error<<"\n";
            if(prior>0)require(prior/error>3.4,"second-order coupled convergence");
            prior=error;
        }
        std::cout<<"NONLINEAR COUPLED TESTS PASSED\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
