#include "_TSP.h"
#include <iostream>
#include <iomanip>
using tsp_math::J;
using tsp_math::Real;

int main() {
    try {
        Real worst=0;int count=0;
        for(double radius:{.3,.7,1.,1.4,2.3})for(double height:{-.999999,-.6,.2,.99,.999999})
            for(double angle:{-2.8,-.2,.01,.6,2.8}) {
                const Real a(.5),h0(height),k0=tsp_native_math::sqrt(1-h0*h0);
                const J h(h0),k(k0),R(radius),phi=J::variable(angle,0),c=h*tsp_math::cos(phi);
                auto H=[&](const Real&s) {
                    const J D=tsp_math::sqrt((R-s)*(R-s)+2*R*s*(1-c));
                    return tsp_math::pow(D,3)/3+R*c*((s-R*c)*D+R*R*(1-c*c)*tsp_math::log(s-R*c+D))/2;
                };
                auto pole=[&](const Real&s) {
                    const J z=s-R;return (z.a[0]>=0?1:-1)*(tsp_math::pow(z,3)/3+R*z*z/2);
                };
                const J jp=pole(Real(1))-pole(a),jm=R*(1-a*a)/2+(1-a*a*a)/3;
                J primitive=TSP::detail::RegularEdge(a,phi,h,k,R);
                const int signs=(R.a[0]<=1?1:-1)-(R.a[0]<=a?1:-1);
                if(signs) {
                    const J half=phi/2;
                    const J A=tsp_math::atan2((1+h)*tsp_math::sin(half),k*tsp_math::cos(half));
                    primitive-=signs*tsp_math::pow(R,3)*k*k*A/2;
                }
                const J kernel=k*((H(Real(1))-H(a))/(1-c*c)-jp/(2*(1-c))-jm/(2*(1+c)));
                for(int degree=0;degree<3;++degree) {
                    const int i=degree==0?0:degree==1?1:4,j=degree==0?1:degree==1?4:10;
                    const Real error=tsp_native_math::abs((degree+1)*primitive.a[j]-kernel.a[i])/(1+tsp_native_math::abs(kernel.a[i]));
                    if(error>worst)worst=error;
                    if(!error.isfinite()||error>Real("1e-16")) {
                        std::cerr<<"identity failed R="<<radius<<" h="<<height<<" phi="<<angle<<" degree="<<degree<<" error="<<error.value<<'\n';
                        return 1;
                    }
                    ++count;
                }
            }
        std::cout<<std::setprecision(18)<<"PASS regular primitive: "<<count<<" coefficients; maximum normalized residual "<<worst.value<<'\n';
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
