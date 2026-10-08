// Compare both independent integrators from a checkpoint produced by our model.
#define main model_program_entry
#include "../src/model.cpp"
#undef main


int main(int argc,char** argv) {
    try {
        std::cout<<std::setprecision(14);
        require(argc==5,"CONFIG CHECKPOINT METHOD DURATION_SECONDS");
        Parameters p=read_config(argv[1]);Model m(p);State initial=m.initial();
        double t0,dt0;int64_t previous_steps,previous_rejects;
        std::ifstream f(argv[2],std::ios::binary);get(f,t0);get(f,dt0);get(f,previous_steps);get(f,previous_rejects);
        for(auto& v:initial)f.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(double));
        require(bool(f),"candidate checkpoint read");
        int method=std::stoi(argv[3]);double end=std::stod(argv[4]);
        State x=initial,out=zeros(m.n),coarse=zeros(m.n),half=zeros(m.n),stage=zeros(m.n),f0=zeros(m.n),f1=zeros(m.n);
        ARK4 solver(m);double t=0,dt=dt0;int steps=0,rejects=0;
        auto begin=std::chrono::steady_clock::now();
        while(t<end) {
            dt=std::min(dt,end-t);double error=1e20;
            bool okay=method==4?solver.step(x,dt,out,error):
                m.step(x,dt,coarse,stage,f0,f1)&&m.step(x,.5*dt,half,stage,f0,f1)&&m.step(half,.5*dt,out,stage,f0,f1);
            if(okay&&method==2)error=m.error(coarse,out);
            if(!okay||!std::isfinite(error))error=1e20;
            if(error<=1){x.swap(out);t+=dt;steps++;}else rejects++;
            dt*=error==0?2:std::clamp(.9*std::pow(error,method==4?-.25:-1./3),.15,2.);
            require(dt>1e-10,"candidate timestep");
        }
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        std::cout<<"method="<<method<<" t0="<<t0<<" duration="<<t<<" steps="<<steps<<" rejects="<<rejects<<" wall_s="<<seconds<<std::endl;
        std::ofstream output(std::string(argv[2])+".method"+argv[3],std::ios::binary);
        for(auto& v:x)output.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
