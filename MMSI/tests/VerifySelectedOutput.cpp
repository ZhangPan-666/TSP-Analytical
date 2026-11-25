#include "_TSPMMSI.h"
#include <cstdio>

TSP_MMSI::Model model(bool magnetic) {
    TSP_MMSI::Model m;m.Fields=magnetic?3:1;
    m.vertices={{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,1},{0,0,-1}};
    m.faces={{4,0,1},{4,1,2},{4,2,3},{4,3,0},{5,1,0},{5,2,1},{5,3,2},{5,0,3}};
    m.data.resize(6*(2+m.Fields));
    for(int i=0;i<6;++i)m.SetState(i,{1,2,{1,.3,-.2}});
    return m;
}
int main() {
    try {
        const double latitude[]{31,-22},longitude[]{42,71},radius[]{3,1.5};
        for(bool magnetic:{false,true}){
            TSP_MMSI::Options options;options.AbsoluteTolerance=options.RelativeTolerance=1e-16;
            options.DerivativeOrder=magnetic?2:1;
            const int fields=magnetic?10:20;const std::uint32_t mask=magnetic?14u:8u;
            std::array<std::array<double,2>,20>reference{},selected{};
            std::array<double*,20>all{},out{};
            for(int k=0;k<fields;++k){all[k]=reference[k].data();selected[k].fill(-123);
                if(mask>>k&1u)out[k]=selected[k].data();}
            TSP_MMSI::work_detail::Run(model(magnetic),2,latitude,longitude,radius,all,options,0);
            TSP_MMSI::work_detail::Run(model(magnetic),2,latitude,longitude,radius,out,options,0,mask);
            for(int k=0;k<fields;++k)for(int i=0;i<2;++i)
                if(selected[k][i]!=(mask>>k&1u?reference[k][i]:-123))
                    throw std::runtime_error("selected output mismatch or unrequested output overwritten");
            bool rejected=false;
            try{TSP_MMSI::work_detail::Run(model(magnetic),2,latitude,longitude,radius,out,options,0,0);}
            catch(const std::invalid_argument&){rejected=true;}
            if(!rejected)throw std::runtime_error("empty selection accepted for nonempty observations");
            std::array<double*,20>empty{};
            TSP_MMSI::work_detail::Run(model(magnetic),0,nullptr,nullptr,nullptr,empty,options,0,mask);
        }
        std::puts("PASS selected output, legacy default, null unused pointers, empty input, invalid selection");return 0;
    }catch(const std::exception&e){std::fprintf(stderr,"FAIL %s\n",e.what());return 1;}
}
