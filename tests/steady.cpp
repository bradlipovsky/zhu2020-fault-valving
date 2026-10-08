// Compare the production steady-flow solver with a separately integrated ODE.
#define main model_program_entry
#include "../src/model.cpp"
#undef main

int main() {
    try {
        std::cout<<std::setprecision(12);
        for(double floor:{0.,1e-19}) {
            double previous=0;
            for(int n:{64,128,256,512}) {
                Parameters p;p.n=n;p.Lz=30000;p.threads=1;p.kmin=floor;p.perturbation=0;
                Model m(p);State x=m.initial();
                double G=p.normal_gradient-p.rho*p.gravity;
                double C=(p.kmax-p.kmin)*x[2][0],keq=p.viscosity*p.influx/G,D=keq-p.kmin;
                double Ninf=p.stress_scale*std::log(C/D),error=0;
                for(int i=0;i<n;i++) {
                    double exact;
                    if(floor==0) {
                        double Xeq=keq/C;
                        exact=-p.stress_scale*std::log(Xeq+(1-Xeq)*std::exp(-G*m.z[i]/p.stress_scale));
                    } else {
                        // Separating dN/dz=G-eta*q/[kmin+C*exp(-N/S)] gives
                        // z=N/G-S*keq/(G*D)*log[(C-D*exp(N/S))/(C-D)].
                        double lo=0,hi=Ninf;
                        for(int it=0;it<70;it++) {
                            double mid=.5*(lo+hi);
                            double depth=mid/G-p.stress_scale*keq/(G*D)*
                                std::log1p(-D*std::expm1(mid/p.stress_scale)/(C-D));
                            if(depth<m.z[i])lo=mid;else hi=mid;
                        }
                        exact=.5*(lo+hi);
                    }
                    error=std::max(error,std::abs(m.normal(i,x[3][i])-exact)/p.stress_scale);
                }
                std::cout<<"steady_kmin="<<floor<<" n="<<n<<" relative_pressure_error="<<error<<"\n";
                if(previous>0)require(previous/error>3.7,"second-order nonlinear steady-pressure convergence");
                previous=error;
            }
        }
        std::cout<<"ANALYTICAL NONLINEAR STEADY-FLOW TESTS PASSED\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
