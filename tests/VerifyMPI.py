"""MPI tests with 1/2/4 ranks, including more ranks than observers.
Requires the serial and MPI CMake targets. No user input/results are replaced.
"""
from pathlib import Path
import argparse
import math, os, struct, subprocess, tempfile
from VerifyWorkflow import ROOT, executable, inputs, run

def launch(launcher,exe,directory,ranks,good=True):
    env=os.environ.copy()
    env['OMP_NUM_THREADS']='1'
    p=subprocess.run([launcher,'-n',str(ranks),str(exe.resolve())],cwd=directory,env=env,capture_output=True,text=True,timeout=60)
    assert (p.returncode==0)==good, f'{exe}, ranks={ranks}: {p.returncode}\n{p.stdout}\n{p.stderr}'

def unpack(data):
    n=struct.unpack_from('<i',data)[0]
    return n,struct.unpack_from('<'+'d'*((len(data)-4)//8),data,4)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, default=ROOT / 'build', help='Directory containing the built programs')
    parser.add_argument('--mpiexec', default='mpiexec', help='MPI launcher whose runtime matches the build')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='tsp-mpi-') as tmp:
        d=Path(tmp)
        for magnetic in [False,True]:
            param='TFM.ForPar' if magnetic else 'TFG.ForPar'
            result='Result_TFM.dat' if magnetic else 'Result_TFG.dat'
            full='Work1_MagneticEstimate' if magnetic else 'Work1_GravityEstimate'
            selected='Work2_MagneticEstimateVectorMPI' if magnetic else 'Work2_GravityEstimateVzMPI'
            full_exe = executable(args.bin_dir, full + 'MPI')
            selected_exe = executable(args.bin_dir, selected)
            (d/param).write_bytes(inputs(magnetic,True,include_protected=True))
            run(executable(args.bin_dir, full),d)
            reference=(d/result).read_bytes()
            n,values=unpack(reference)
            for ranks in [1,2,4]:
                launch(args.mpiexec,full_exe,d,ranks)
                assert (d/result).read_bytes()==reference, 'MPI full output differs from serial'
                launch(args.mpiexec,selected_exe,d,ranks)
                count,partial=unpack((d/result).read_bytes())
                assert count==n and partial[:3*n]==values[:3*n]
                kept=[1,2,3] if magnetic else [3]
                for field in range(10 if magnetic else 20):
                    actual=partial[(field+3)*n:(field+4)*n]
                    expected=values[(field+3)*n:(field+4)*n]
                    if field in kept:
                        assert all(abs(a-b)<=1e-16*(1+abs(b)) for a,b in zip(actual,expected)), 'selected MPI component differs'
                    else:
                        assert all(map(math.isnan,actual)), 'uncomputed field must be NaN'
            (d/param).write_bytes(struct.pack('<ii',0,0))
            launch(args.mpiexec,full_exe,d,4)
            assert (d/result).read_bytes()==struct.pack('<i',0)
            (d/param).write_bytes(b'bad')
            launch(args.mpiexec,full_exe,d,2,False)
            print(full+'MPI: 1/2/4 ranks, serial byte equality, selected fields, empty ranks, invalid input PASS')
    print('ALL MPI TESTS PASSED')

if __name__=='__main__':main()
