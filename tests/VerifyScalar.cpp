#include "_TSP.h"
#include <chrono>
#include <cstdio>
#include <iomanip>
#include <random>
using tsp_math::J;
using tsp_math::Real;
TSP::detail::EllipticValues evaluate(const J&u,const J&m,const J&n,bool scalar) {
    if(scalar)return TSP::detail::EvaluateElliptic(u,m,n);
    using namespace TSP;
    auto elliptic=[&](const J&theta,bool complete) {
        return detail::CarlsonElliptic(complete?J(tsp_math::PI/2):theta,m,n,complete);
    };
    const int periods=static_cast<int>(std::floor((u.value()+PI/2)/PI));
    J v=u-periods*tsp_math::PI;
    if(std::abs(std::abs(v.value())-PI/2)<1e-10) {
        const int sign=v.value()>0?1:-1;
        auto out=elliptic(v,true);J d=v-sign*tsp_math::PI/2;
        const J root=tsp_math::sqrt(1-m),den=1-n;
        out.F=sign*out.F+d/root-m*tsp_math::pow(d,3)/(6*tsp_math::pow(root,3));
        out.E=sign*out.E+d*root+m*tsp_math::pow(d,3)/(6*root);
        out.P=sign*out.P+d/(den*root)-tsp_math::pow(d,3)*(n/den+m/(2*(1-m)))/(3*den*root);
        if(periods) {auto c=elliptic(v,true);out.F+=2*periods*c.F;out.E+=2*periods*c.E;out.P+=2*periods*c.P;}
        return out;
    }
    auto out=elliptic(v,false);
    if(periods) {auto c=elliptic(v,true);out.F+=2*periods*c.F;out.E+=2*periods*c.E;out.P+=2*periods*c.P;}
    return out;
}
int main() {
    try {
        std::setvbuf(stdout,nullptr,_IONBF,0);
        std::mt19937_64 rng(20261028);std::uniform_real_distribution<double>uniform(0,1);
        std::ofstream output("seed_samples.csv");output<<std::setprecision(32);
        output<<"case,amplitude,m,n,field,coefficient,reference,scalar,scaled_error\n";
        double worst=0;int count=0;
        auto sample=[&](double t,double x,double y) {
            const double m=1-std::exp(x),n=1-std::exp(y);
            const auto u=J::variable(t,0),mj=J::variable(m,1),nj=J::variable(n,2);
            const auto reference=evaluate(u,mj,nj,false),actual=evaluate(u,mj,nj,true);
            const J a[3]={reference.F,reference.E,reference.P},b[3]={actual.F,actual.E,actual.P};
            for(int f=0;f<3;++f)for(int k=0;k<20;++k) {
                const double error=std::abs(static_cast<double>(a[f].a[k]-b[f].a[k]))/(1+std::abs(static_cast<double>(a[f].a[k])));
                worst=std::max(worst,error);
                output<<count<<','<<t<<','<<m<<','<<n<<','<<f<<','<<k<<','<<a[f].a[k].value<<','<<b[f].a[k].value<<','<<error<<'\n';
                if(!std::isfinite(error)||error>1e-16) {
                    std::fprintf(stderr,"FAIL seed case=%d t=%.17g X=%g Y=%g f=%d k=%d error=%g\n",count,t,x,y,f,k,error);
                    throw std::runtime_error("scalar seed derivative mismatch");
                }
            }
            for(auto mask:{15u,1023u}) {
                tsp_math::JetScope scope(mask);
                const auto partial=evaluate(u,mj,nj,true);
                const J fields[3]={partial.F,partial.E,partial.P};
                for(int f=0;f<3;++f)for(int k=0;k<20;++k)if(tsp_math::jet_active(k))
                    if(std::abs(static_cast<double>(fields[f].a[k]-b[f].a[k]))>1e-16*(1+std::abs(static_cast<double>(b[f].a[k]))))
                        throw std::runtime_error("masked derivative recurrence mismatch");
            }
            ++count;
        };
        for(int i=0;i<120;++i)sample((.0001+.9998*uniform(rng))*TSP::PI/2,.011+21.8*uniform(rng),-21.8+43.6*uniform(rng));
        for(double t:{0.,1e-6,.01,.4,TSP::PI/2-1e-6,TSP::PI/2,3.4,-.7})
            for(auto xy:std::array<std::pair<double,double>,5>{{{.02,.03},{6.,8.},{20.,-21.},{21.,21.},{10.,10.}}})sample(t,xy.first,xy.second);
        for(double t:{0.,.4,TSP::PI/2,3.4})for(double x:{.005,13.8675})
            for(double delta:{-1e-9,0.,1e-9})sample(t,x,std::log(.001)+delta);
        std::printf("PASS seeds: %d complete/principal/period cases, 60 coefficients each, max scaled difference %.6g\n",count,worst);
        if(!output)throw std::runtime_error("cannot write scalar seed checks");
        std::puts("ALL SCALAR TESTS PASSED");return 0;
    }catch(const std::exception&e){std::fprintf(stderr,"FAIL %s\n",e.what());return 1;}
}
