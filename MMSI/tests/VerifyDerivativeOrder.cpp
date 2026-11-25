// Build with -I..; absolute/relative targets remain 1e-16.
#ifndef TSP_ORDER_HEADER
#define TSP_ORDER_HEADER "_TSPMMSI.h"
#endif
#include TSP_ORDER_HEADER
#include <cstdio>
#include <thread>

namespace {
void require(bool good,const char*message) {if(!good)throw std::runtime_error(message);}
TSP_MMSI::Model source(int kind) {
    TSP_MMSI::Model m;m.Fields=kind==2?3:1;
    m.vertices={{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,1},{0,0,-1}};
    m.faces={{4,0,1},{4,1,2},{4,2,3},{4,3,0},{5,1,0},{5,2,1},{5,3,2},{5,0,3}};
    m.data.resize(m.vertices.size()*(2+m.Fields));
    for(std::size_t i=0;i<m.vertices.size();++i){const auto&v=m.vertices[i];
        const double top=kind==0?2:2+.02*v[0]+.01*v[2];
        m.SetState(i,{1,top,{1,.3,-.2}});}
    return m;
}
int components(int order) {return order==3?20:order==2?10:order==1?4:1;}
void analytic() {
    const TSP::Prism prism{{{{1,0,0},{0,1,0},{0,0,1}}},1,2,1};
    const std::array<TSP::Vec,5>points{{{.29,-.27,3.1},{.1,.2,.3},{.8,.9,.9},{.3,.4,45},{0,0,3}}};
    const auto axes=TSP::LocalAxes(-23,-71);
    for(std::size_t p=0;p<points.size();++p)for(int order=0;order<=3;++order){
        const tsp_math::JetScope scope((std::uint32_t(1)<<components(order))-1);
        TSP::Options options;options.absolute_tolerance=options.relative_tolerance=1e-16;
        const auto cart=TSP::GravityCartesian(prism,points[p],options);
        const auto rotated=TSP::RotateGravity(cart,axes);
        const auto local=TSP::GravityLocalCartesian(prism,points[p],options);
        const auto magnetic=TSP::MagneticCartesian(prism,points[p],{.3,-.2,.7},options);
        std::printf("A %zu %d",p,order);
        for(int k=0;k<components(order);++k)std::printf(" %.17g %.17g %.17g",cart[k],rotated[k],local[k]);
        const int n=order==3?10:order==2?4:order==1?1:0;
        for(int k=0;k<n;++k)std::printf(" %.17g",magnetic[k]);
        std::puts("");
    }
}
void shell() {
    const auto axes=TSP_MMSI::detail::LocalAxes(31,42);
    for(int kind=0;kind<3;++kind)for(int order=0;order<=3;++order){
        TSP_MMSI::Options options;options.AbsoluteTolerance=options.RelativeTolerance=1e-16;options.DerivativeOrder=order;
        TSP_MMSI::FastShellIntegrator engine(source(kind),options,1);
        for(int p=0;p<3;++p){
            const double radius=p==0?3:p==1?1.5:.4;
            const auto point=TSP_MMSI::detail::Scale(axes[2],radius);
            TSP_MMSI::FastDiagnostics d;
            const auto field=engine.Fields(point,&d);
            require(tsp_math::active_jet_mask==tsp_math::JET_FULL_MASK,"field mask leaked");
            require(d.QuadratureConverged==(d.ErrorRatio<=1),"false convergence");
            std::printf("M %d %d %d %llu %zu %zu %zu %d %d %d %.17g",kind,order,p,
                (unsigned long long)d.KernelNodes,d.AnalyticalTriangles,d.AnalyticalEdges,d.Refinements,
                int(d.QuadratureConverged),int(d.LimitReached),int(d.NearApproximation),d.ErrorRatio);
            for(int f=0;f<(kind==2?3:1);++f){std::array<double,20>part{};
                for(int k=0;k<components(order);++k)part[k]=field[f*20+k];
                const auto local=TSP_MMSI::detail::Rotate(part,axes,components(order));
                for(int k=0;k<components(order);++k)std::printf(" %.17g %.17g",part[k],local[k]);
                for(int k=components(order);k<20;++k)require(field[f*20+k]==0,"unrequested field computed");
            }
            std::puts("");std::fflush(stdout);
        }
    }
}
void restoration() {
    TSP_MMSI::Options o;o.AbsoluteTolerance=o.RelativeTolerance=1e-16;
    TSP_MMSI::FastShellIntegrator engine(source(0),o,1);
    const auto point=TSP_MMSI::detail::Scale(TSP_MMSI::detail::Direction(31,42),3);
    const auto full=engine.Fields(point);
    {
        const tsp_math::JetScope outer(15);
        const auto actual=engine.Fields(point);
#ifndef TSP_ORDER_BASELINE
        require(actual==full,"request order inherited a caller's lower-order mask");
#else
        (void)actual;(void)full;
#endif
        require(tsp_math::active_jet_mask==15,"nested mask not restored");
        bool rejected=false;
        try {engine.Fields(TSP_MMSI::detail::Scale(TSP_MMSI::detail::Direction(31,42),2));}
        catch(const std::domain_error&){rejected=true;}
        require(rejected&&tsp_math::active_jet_mask==15,"exception mask not restored");
    }
    std::exception_ptr error;
    std::thread worker([&]{try{
        require(tsp_math::active_jet_mask==tsp_math::JET_FULL_MASK,"worker initial mask");
        const tsp_math::JetScope scope(1023);
        engine.Fields(point);
        require(tsp_math::active_jet_mask==1023,"worker mask not restored");
    }catch(...){error=std::current_exception();}});
    worker.join();if(error)std::rethrow_exception(error);
    require(tsp_math::active_jet_mask==tsp_math::JET_FULL_MASK,"main thread mask changed");
}
}
int main() {
    try {analytic();shell();restoration();std::puts("PASS derivative-order regression");return 0;}
    catch(const std::exception&e){std::fprintf(stderr,"FAIL: %s\n",e.what());return 1;}
}
