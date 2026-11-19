#define TSP_MPI
#include "_TSP.h"
int main(int argc,char**argv){return TSP::work_detail::run_mpi<3>(argc,argv);}
