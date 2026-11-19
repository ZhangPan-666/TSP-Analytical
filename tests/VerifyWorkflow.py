"""Standard-library binary workflow tests for the CMake executables.
Uses temporary directories and never replaces user parameter/result files.
Pass --write-examples to generate the two documented sample inputs in C++Codes.
"""
from pathlib import Path
import argparse
import math, struct, subprocess, tempfile, os

ROOT = Path(__file__).resolve().parents[1]

def inputs(magnetic=False, two_models=False, include_protected=False):
    f = [90-x*180/math.pi for x in [.50,.63,1.05]]
    model = [[v] for v in (*f,0,.85*180/math.pi,.39*180/math.pi,1,2)]
    model += [[.3],[-.2],[.7]] if magnetic else [[1.]]
    if two_models:
        model = [x+x for x in model]
    coords = [[-25.,0.],[61.,90.],[3.4,3.]]
    if include_protected:
        # Opposite an edge at R=R2 exercises the radial-endpoint chart.
        a=(math.sin(.5),0.,math.cos(.5))
        b=(math.sin(.63)*math.cos(.85),math.sin(.63)*math.sin(.85),math.cos(.63))
        unit=lambda v: [x/math.sqrt(sum(y*y for y in v)) for x in v]
        mid=unit([a[k]+b[k] for k in range(3)])
        normal=unit([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
        q=unit([-mid[k]+.02*normal[k] for k in range(3)])
        coords[0].append(math.degrees(math.atan2(q[1],q[0])))
        coords[1].append(math.degrees(math.asin(q[2])))
        coords[2].append(2.)
    pack = lambda v: struct.pack('<'+'d'*len(v),*v)
    return (struct.pack('<i',len(model[0]))+b''.join(map(pack,model))
            +struct.pack('<i',len(coords[0]))+b''.join(map(pack,coords)))

def executable(directory, name):
    return directory / (name + ('.exe' if os.name == 'nt' else ''))


def run(exe, directory, good=True):
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '2'
    proc = subprocess.run([str(exe.resolve())], cwd=directory, env=env, capture_output=True, text=True, timeout=60)
    if (proc.returncode==0) != good:
        raise AssertionError(f'{exe}: return={proc.returncode}\n{proc.stdout}\n{proc.stderr}')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, default=ROOT / 'build', help='Directory containing the built programs')
    parser.add_argument('--write-examples', action='store_true', help='Regenerate both example input files in the repository root')
    args = parser.parse_args()
    if args.write_examples:
        for mag in [False,True]:
            (ROOT/('TFM.ForPar' if mag else 'TFG.ForPar')).write_bytes(inputs(mag))
        print('Example TFG.ForPar and TFM.ForPar created')
        return
    with tempfile.TemporaryDirectory(prefix='tsp-workflow-') as temp:
        d=Path(temp)
        for mag in [False,True]:
            param='TFM.ForPar' if mag else 'TFG.ForPar'
            result='Result_TFM.dat' if mag else 'Result_TFG.dat'
            exe='Work1_MagneticEstimate' if mag else 'Work1_GravityEstimate'
            serial_exe = executable(args.bin_dir, exe)
            openmp_exe = executable(args.bin_dir, exe + 'OpenMP')
            (d/param).write_bytes(inputs(mag))
            run(serial_exe,d)
            serial=(d/result).read_bytes()
            n=struct.unpack_from('<i',serial)[0]
            assert n==2 and len(serial)==4+n*(13 if mag else 23)*8
            values=struct.unpack_from('<'+'d'*((len(serial)-4)//8),serial,4)
            assert all(map(math.isfinite,values))
            if openmp_exe.is_file():
                run(openmp_exe,d)
                assert (d/result).read_bytes()==serial, 'OpenMP differs from serial'
            (d/param).write_bytes(inputs(mag,include_protected=True))
            run(serial_exe,d)
            protected_result=(d/result).read_bytes()
            if openmp_exe.is_file():
                run(openmp_exe,d)
                assert (d/result).read_bytes()==protected_result, 'OpenMP protected Carlson branch differs'
            (d/param).write_bytes(inputs(mag,True))
            run(serial_exe,d)
            doubled=(d/result).read_bytes()
            x=struct.unpack_from('<'+'d'*((len(serial)-4)//8),doubled,4)
            assert all(abs(a-2*b)<1e-16*(1+abs(a)) for a,b in zip(x[3*n:],values[3*n:])), 'model superposition'
            (d/param).write_bytes(struct.pack('<i',0)+struct.pack('<i',0))
            run(serial_exe,d)
            assert (d/result).read_bytes()==struct.pack('<i',0), 'empty dataset'
            (d/param).write_bytes(b'bad')
            run(serial_exe,d,False)
            (d/param).write_bytes(inputs(mag)+b'extra')
            run(serial_exe,d,False)
            print(exe+': binary layout, superposition, empty data, malformed input PASS')
            if openmp_exe.is_file():
                print(exe+': OpenMP byte comparison PASS')
    print('ALL WORKFLOW TESTS PASSED (analytic production)')

if __name__=='__main__':
    main()
