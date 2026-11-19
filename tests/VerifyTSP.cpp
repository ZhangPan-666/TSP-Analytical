#include "_TSP.h"
#include <chrono>
#include <cstdio>
#include <iomanip>
#include <random>
#include <sstream>
using namespace TSP;
namespace {
struct Rule {std::vector<double>x,w;explicit Rule(int n) {
    x.resize(n);w.resize(n);
    for(int i=0;i<(n+1)/2;++i) {
        double z=std::cos(PI*(i+.75)/(n+.5)),der=0;
        for(int k=0;k<40;++k) {double p0=1,p1=z;for(int j=2;j<=n;++j){double t=((2*j-1)*z*p1-(j-1)*p0)/j;p0=p1;p1=t;}
            der=n*(z*p1-p0)/(z*z-1);double next=z-p1/der;if(std::abs(next-z)<3e-16){z=next;break;}z=next;}
        x[i]=(1-z)/2;x[n-1-i]=(1+z)/2;w[i]=w[n-1-i]=1/((1-z*z)*der*der);
    }
}};
void compensated(double value,double&sum,double&error) {const double t=sum+value;error+=std::abs(sum)>=std::abs(value)?(sum-t)+value:(value-t)+sum;sum=t;}
void check(bool ok,const char*message) {if(!ok)throw std::runtime_error(message);}
Vec polar(double t,double l) {return {std::sin(t)*std::cos(l),std::sin(t)*std::sin(l),std::cos(t)};}
Prism generic() {return {{{polar(.50,0),polar(.63,.85),polar(1.05,.39)}},1,2,1};}
// Independent 3-D Gauss rule using affine planar triangle projection. It does
// not call Geometry, jets, elliptic integrals or the surface evaluator.
Gravity volume(const Prism&prism,Vec p,int order) {
    const auto a=detail::unit(prism.vertices[0]),b=detail::unit(prism.vertices[1]),c=detail::unit(prism.vertices[2]);
    const double det=std::abs(detail::dot(a,detail::cross(b,c)));
    Rule rule(order);Gravity out{},corr{};
    for(int i=0;i<order;++i)for(int j=0;j<order;++j) {
        const double t=rule.x[i],v=rule.x[j];Vec z{};
        for(int k=0;k<3;++k)z[k]=(1-t)*(1-v)*a[k]+t*b[k]+(1-t)*v*c[k];
        const double ln=detail::norm(z);const Vec u=detail::scale(z,1/ln);
        const double angular=det/(ln*ln*ln)*(1-t)*rule.w[i]*rule.w[j];
        for(int h=0;h<order;++h) {
            const double s=prism.r1+(prism.r2-prism.r1)*rule.x[h];
            const double weight=angular*rule.w[h]*(prism.r2-prism.r1)*s*s*prism.density;
            Vec q{};for(int k=0;k<3;++k)q[k]=s*u[k]-p[k];
            const double d=detail::norm(q),iv=1/d,iv3=iv*iv*iv,iv5=iv3*iv*iv,iv7=iv5*iv*iv;
            Gravity f{};f[0]=iv;for(int k=0;k<3;++k)f[k+1]=q[k]*iv3;
            for(int k=4;k<20;++k) {
                const auto e=tsp_math::ex[k];int ax[3],n=0;
                for(int r=0;r<e.x;++r)ax[n++]=0;for(int r=0;r<e.y;++r)ax[n++]=1;for(int r=0;r<e.z;++r)ax[n++]=2;
                const int x=ax[0],y=ax[1];
                if(n==2)f[k]=3*q[x]*q[y]*iv5-(x==y)*iv3;
                else {int z=ax[2];f[k]=15*q[x]*q[y]*q[z]*iv7-3*((x==y)*q[z]+(x==z)*q[y]+(y==z)*q[x])*iv5;}
            }
            for(int k=0;k<20;++k)compensated(weight*f[k],out[k],corr[k]);
        }
    }
    for(int k=0;k<20;++k)out[k]+=corr[k];return out;
}
Magnetic dipoles(const Prism&prism,Vec p,Vec m,int order) {
    const Vec a=detail::unit(prism.vertices[0]),b=detail::unit(prism.vertices[1]),c=detail::unit(prism.vertices[2]);
    const double det=std::abs(detail::dot(a,detail::cross(b,c)));Rule rule(order);Magnetic out{},err{};
    constexpr int axes[6][2]={{0,0},{0,1},{1,1},{2,0},{2,1},{2,2}};
    for(int i=0;i<order;++i)for(int j=0;j<order;++j) {
        const double t=rule.x[i],v=rule.x[j];Vec z{};for(int k=0;k<3;++k)z[k]=(1-t)*(1-v)*a[k]+t*b[k]+(1-t)*v*c[k];
        const double length=detail::norm(z);Vec u=detail::scale(z,1/length);
        for(int h=0;h<order;++h) {
            const double s=prism.r1+(prism.r2-prism.r1)*rule.x[h];
            const double w=det/std::pow(length,3)*(1-t)*rule.w[i]*rule.w[j]*rule.w[h]*(prism.r2-prism.r1)*s*s;
            Vec q{};for(int k=0;k<3;++k)q[k]=s*u[k]-p[k];const double d=detail::norm(q),r3=std::pow(d,-3),r5=std::pow(d,-5),r7=std::pow(d,-7),mq=detail::dot(m,q);
            Magnetic f{};f[0]=mq*r3;for(int k=0;k<3;++k)f[k+1]=3*q[k]*mq*r5-m[k]*r3;
            for(int k=0;k<6;++k){int x=axes[k][0],y=axes[k][1];f[k+4]=15*q[x]*q[y]*mq*r7-3*(m[x]*q[y]+m[y]*q[x]+(x==y)*mq)*r5;}
            for(int k=0;k<10;++k)compensated(w*f[k],out[k],err[k]);
        }
    }
    for(int k=0;k<10;++k)out[k]+=err[k];return out;
}
template<std::size_t N> double error(const std::array<double,N>&a,const std::array<double,N>&b) {double d=0;for(int k=0;k<N;++k)d=std::max(d,std::abs(a[k]-b[k]));return d;}
std::ofstream csv;
void exterior(const char*name,const Prism&p,Vec x,Method method=Method::Automatic,double tol=2e-8) {
    Options opt;opt.method=method;Diagnostics diag;
    const auto start=std::chrono::steady_clock::now();const auto a=GravityCartesian(p,x,opt,&diag);
    const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    const auto b=volume(p,x,28),c=volume(p,x,40);
    double d=error(a,c),conv=error(b,c);
    std::printf("%-24s method=%s max_error=%.3e reference_convergence=%.3e ms=%.2f\n",name,diag.method==Method::Elliptic?"elliptic":"surface",d,conv,ms);
    check(conv<tol/4,"independent volume reference did not converge");check(d<tol,"gravity component verification failed");
    for(int k=0;k<20;++k)csv<<name<<','<<k<<','<<x[0]<<','<<x[1]<<','<<x[2]<<','<<a[k]<<','<<c[k]<<','<<a[k]-c[k]<<','<<conv<<','<<ms<<','<<(diag.method==Method::Elliptic?"elliptic":"surface")<<"\n";
}
void ell_functions() {
    using namespace tsp_math;
    // m=nu=0 is exactly F=E=Pi=u, including amplitudes beyond one period.
    for(double u:{-4.,-TSP::PI/2,-.3,0.,.4,TSP::PI/2,4.}) {
        auto f=detail::EvaluateElliptic(J::variable(u,0),J(0),J(0));
        for(auto v:{f.F,f.E,f.P}){check(std::abs(v.value()-u)<1e-14,"elliptic amplitude continuation");check(std::abs(static_cast<double>(v.a[1])-1)<1e-14,"elliptic amplitude derivative");}
    }
    // Parameter coalescence m=nu and amplitude=0 must retain their jets.
    auto m=J::variable(-.4,1),n=J::variable(-.4,2),u=J::variable(0,0);
    const auto f=detail::EvaluateElliptic(u,m,n);check(std::abs(static_cast<double>(f.P.a[1])-1)<1e-15,"zero amplitude jet");
    std::puts("elliptic branch and coalescent parameter identities: PASS");
}
void public_api() {
    auto p=generic();double f1=90-.5*180/PI,f2=90-.63*180/PI,f3=90-1.05*180/PI,l1=0,l2=.85*180/PI,l3=.39*180/PI,r1=1,r2=2,rho=2;
    double lat=61,lon=-25,r=3.4;Gravity got{};
    TSP_GravityEstimate(1,&f1,&f2,&f3,&l1,&l2,&l3,&r1,&r2,&rho,1,&lat,&lon,&r,&got[0],&got[1],&got[2],&got[3],&got[4],&got[5],&got[7],&got[6],&got[8],&got[9],&got[10],&got[11],&got[12],&got[14],&got[13],&got[16],&got[17],&got[15],&got[18],&got[19]);
    p.density=2;auto expected=RotateGravity(volume(p,work_detail::observation(lat,lon,r),36),LocalAxes(lat,lon));check(error(got,expected)<1e-8,"public NEU gravity interface");
    double mx=.3,my=-.2,mz=.7;Magnetic magnetic{};
    TSP_MagneticEstimate(1,&f1,&f2,&f3,&l1,&l2,&l3,&r1,&r2,&mx,&my,&mz,1,&lat,&lon,&r,&magnetic[0],&magnetic[1],&magnetic[2],&magnetic[3],&magnetic[4],&magnetic[5],&magnetic[6],&magnetic[7],&magnetic[8],&magnetic[9]);
    p.density=1;auto global=work_detail::source_magnetization(p,mx,my,mz);auto axes=LocalAxes(lat,lon);auto independent=dipoles(p,work_detail::observation(lat,lon,r),global,36);
    Magnetic local{};local[0]=independent[0];for(int i=0;i<3;++i)for(int k=0;k<3;++k)local[i+1]+=axes[i][k]*independent[k+1];
    int ij[6][2]={{0,0},{0,1},{1,1},{2,0},{2,1},{2,2}};double H[3][3]={{independent[4],independent[5],independent[7]},{independent[5],independent[6],independent[8]},{independent[7],independent[8],independent[9]}};
    for(int a=0;a<6;++a)for(int i=0;i<3;++i)for(int j=0;j<3;++j)local[a+4]+=axes[ij[a][0]][i]*axes[ij[a][1]][j]*H[i][j];
    check(error(magnetic,local)<1e-8,"public NEU magnetic interface");
    double vz,bx,by,bz;TSP_VzEstimate(1,&f1,&f2,&f3,&l1,&l2,&l3,&r1,&r2,&rho,1,&lat,&lon,&r,&vz);
    TSP_MagneticVectorEstimate(1,&f1,&f2,&f3,&l1,&l2,&l3,&r1,&r2,&mx,&my,&mz,1,&lat,&lon,&r,&bx,&by,&bz);
    check(std::abs(vz-got[3])<1e-14&&std::abs(bx-magnetic[1])+std::abs(by-magnetic[2])+std::abs(bz-magnetic[3])<1e-14,"selected APIs");
    std::puts("array API / NEU rotation / source magnetization frame: PASS");
}
void local_precision() {
    Prism p{{Vec{-0.29619813272602386,0.8137976813493738,0.49999999999999994},Vec{-0.428583650351056,0.74232865770136425,0.51503807491005416},Vec{-0.31157834524108396,0.69981642136369815,0.64278760968653925}},6341.2,6371.2,2670};const Vec point{-2315.4911752135808,4971.4229713279174,3242.785290465918};
    const Gravity expected{133846678.86319067,579311.72183865658,-64328.632337306604,-254019.09110253843,-191.23949486430249,-8.8000226376464319,-54755.818634398987,-369.56706949368953,6368.3264096562389,560.80656435799199,-5233243.4760315977,608766.32489537634,0.18665053237831047,-70825.533271330991,0.0039106286725060242,5304069.009302929,8241.1530464753505,0.29026362518064863,-617007.4779418516,-0.47691415755895911};const auto local=GravityLocalCartesian(p,point);
    for(int k=0;k<20;++k)check(std::abs(local[k]-expected[k])<=1e-16*(1+std::abs(expected[k])),"unrounded metre-clearance local tensor");
    const Magnetic magneticExpected{160.58570568352292,5.1828258305912867,-11.698385736611467,-9.7994574061883508,5.6555567580222155,-165.41491350368005,567.99633434624923,-736.32698378050088,1243.2295783036404,-573.6518911042715};const auto magnetic=MagneticCartesian(p,point,{.3,-.2,.7});
    for(int k=0;k<10;++k)check(magnetic[k]==magneticExpected[k],"unrounded magnetic contraction");
    for(const auto&angle:std::array<std::array<double,2>,6>{{{0,0},{90,123},{-90,78},{32,114},{0,180},{89.99999999,-179}}}) {
        const auto axes=detail::local_axes(angle[0],angle[1]);
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)check(tsp_native_math::abs(detail::rdot(axes[i],axes[j])-(i==j?1:0))<detail::Real("1e-29"),"two-word local-axis orthonormality");
    }
    for(const auto&point:std::array<Vec,3>{{{0,0,3},{0,0,-3},{0,0,0}}}) {
        const auto axes=detail::local_axes(detail::convert(point));
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)check(tsp_native_math::abs(detail::rdot(axes[i],axes[j])-(i==j?1:0))<detail::Real("1e-29"),"Cartesian pole/origin frame");
    }
    std::puts("two-word local frame / near-zero third derivative / magnetic contraction: PASS");
}
void internal_and_errors() {
    auto p=generic();Vec n=detail::unit(detail::add(detail::add(p.vertices[0],p.vertices[1]),p.vertices[2]));
    Options a;a.method=Method::Elliptic;
    auto x=detail::scale(n,1.5);auto e=GravityCartesian(p,x,a);
    check(std::abs(e[4]+e[7]+e[9]+4*PI)<1e-8,"interior Poisson trace");
    // A rotated model and observer must give rotated global tensors.
    const std::array<Vec,3>axes{{{0,1,0},{0,0,1},{1,0,0}}};auto pr=p;
    for(auto&v:pr.vertices)v={v[1],v[2],v[0]};Vec xp{x[1],x[2],x[0]};
    check(error(GravityCartesian(pr,xp),RotateGravity(e,axes))<2e-8,"rotation covariance");
    auto q=p;q.density=-2;auto z=GravityCartesian(q,x);for(int k=0;k<20;++k)check(std::abs(z[k]+2*e[k])<2e-9,"density linearity");
    bool caught=false;try{GravityCartesian(p,detail::scale(n,2));}catch(const std::domain_error&){caught=true;}check(caught,"physical boundary must throw");
    caught=false;try{auto bad=p;bad.vertices[1]=bad.vertices[0];GravityCartesian(bad,x);}catch(const std::invalid_argument&){caught=true;}check(caught,"invalid triangle must throw");
    double xx=0,yy=0,zz=0;TSP_VzEstimate(0,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,1,&xx,&yy,&zz,&xx);check(xx==0,"empty model must produce zero");
    std::puts("interior Poisson / covariance / errors: PASS");
}
}
int main() {
    std::setvbuf(stdout,nullptr,_IONBF,0);
    try {
        csv.open("verification_components.csv");csv<<std::setprecision(17)<<"case,component,x,y,z,computed,volume_reference,residual,reference_convergence,runtime_ms,method\n";
        local_precision();ell_functions();auto p=generic();
        exterior("generic exterior",p,{.29,-.27,3.1},Method::Elliptic);
        exterior("inner cavity",p,{.25,-.1,.3},Method::Elliptic);
        exterior("radial shell exterior",p,{0,0,1.3},Method::Elliptic);
        auto n=detail::unit(detail::add(detail::add(p.vertices[0],p.vertices[1]),p.vertices[2]));
        exterior("north angular cap",p,detail::scale(n,3.4),Method::Elliptic);
        exterior("south angular cap",p,detail::scale(n,-3.4),Method::Elliptic);
        exterior("R equals r1",p,{0,0,1},Method::Elliptic);
        exterior("R equals r2",p,{0,0,2},Method::Elliptic);
        exterior("center",p,{0,0,0});
        auto mid=detail::unit(detail::add(p.vertices[0],p.vertices[1]));exterior("h equals 1",p,detail::scale(mid,3));
        auto N=detail::unit(detail::cross(p.vertices[0],p.vertices[1]));exterior("h equals 0",p,detail::scale(N,3));
        auto pi_miss=detail::unit(detail::add(mid,detail::scale(N,.005)));
        exterior("table Pi outside",p,detail::scale(pi_miss,3),Method::Elliptic);
        auto fe_miss=detail::unit(detail::add(detail::scale(mid,-1),detail::scale(N,.02)));
        exterior("table FE outside",p,detail::scale(fe_miss,2),Method::Elliptic);
        Prism octant{{{{1,0,0},{0,1,0},{0,0,1}}},1,2,1};
        auto north=GravityCartesian(octant,{0,0,3}),south=GravityCartesian(octant,{0,0,-3});
        check(std::abs(north[0]-1.5270566268062245)<2e-11&&std::abs(south[0]-.9164043259858369)<2e-11,"elementary axial benchmarks");
        const Vec M{.3,-.2,.7},obs{.29,-.27,3.1};auto magnetic=MagneticCartesian(p,obs,M);auto direct=dipoles(p,obs,M,36);
        std::printf("magnetic direct dipole error: %.3e\n",error(magnetic,direct));check(error(magnetic,direct)<1e-8,"independent magnetic dipoles");
        std::mt19937_64 rng(20261022);std::uniform_real_distribution<double>u(-1,1);
        for(int i=0;i<6;++i){Vec q{u(rng),u(rng),u(rng)};q=detail::scale(detail::unit(q),3+i*.2);std::string name="random "+std::to_string(i);exterior(name.c_str(),p,q);}
        public_api();internal_and_errors();csv.close();std::puts("ALL TSP TESTS PASSED (seed 20261022)");return 0;
    }catch(const std::exception&e){std::fprintf(stderr,"VALIDATION FAILED: %s\n",e.what());return 1;}
}
