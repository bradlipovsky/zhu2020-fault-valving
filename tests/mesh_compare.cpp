// Compare independently generated checkpoints on a common physical depth range.
#define main model_program_entry
#include "../src/model.cpp"
#undef main

State read_state(const std::string& path,Model& model,double& time) {
    State x=model.initial();double dt;int64_t step,rejected;
    std::ifstream input(path,std::ios::binary);get(input,time);get(input,dt);get(input,step);get(input,rejected);
    for(auto& field:x)input.read(reinterpret_cast<char*>(field.data()),field.size()*sizeof(double));
    require(bool(input),"read comparison checkpoint");return x;
}
int main(int argc,char** argv) {
    try {
        if(argc!=5 && argc!=7)throw std::runtime_error("Usage: mesh_compare CONFIG_A CHECKPOINT_A CONFIG_B CHECKPOINT_B [MIN_DEPTH MAX_DEPTH]");
        Model a(read_config(argv[1])),b(read_config(argv[3]));double ta,tb;
        State x=read_state(argv[2],a,ta),y=read_state(argv[4],b,tb);
        require(ta==tb,"comparison endpoints must match exactly");
        double low=argc==7?std::stod(argv[5]):3000,high=argc==7?std::stod(argv[6]):25000;
        State rate=zeros(a.n);require(a.mechanical(x,rate),"comparison state A");
        rate=zeros(b.n);require(b.mechanical(y,rate),"comparison state B");
        std::array<double,4> error{};double velocity_error=0;
        for(int i=0;i<a.n;i++) {
            if(a.z[i]<low || a.z[i]>high)continue;
            int j=std::clamp(int(std::upper_bound(b.z.begin(),b.z.end(),a.z[i])-b.z.begin())-1,0,b.n-2);
            double w=std::clamp((a.z[i]-b.z[j])/(b.z[j+1]-b.z[j]),0.,1.);
            for(int field=0;field<4;field++)error[field]=std::max(error[field],
                std::abs(x[field][i]-((1-w)*y[field][j]+w*y[field][j+1])));
            velocity_error=std::max(velocity_error,std::abs(a.vel[i]-((1-w)*b.vel[j]+w*b.vel[j+1])));
        }
        std::cout<<std::setprecision(15)<<"{\"years\":"<<ta/YEAR<<",\"minimum_depth_m\":"<<low
            <<",\"maximum_depth_m\":"<<high<<",\"slip_difference_m\":"<<error[0]
            <<",\"state_difference\":"<<error[1]<<",\"normalized_kstar_difference\":"<<error[2]
            <<",\"pressure_difference_pa\":"<<error[3]<<",\"velocity_difference_m_per_s\":"<<velocity_error<<"}\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
