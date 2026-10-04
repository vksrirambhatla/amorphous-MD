"""Expensive fixed-cell first-principles melt/quench sampling using ASE and QE forces."""
import numpy as np
from common import *

def main():
    from ase.io import read, write
    from ase import units
    from ase.calculators.espresso import Espresso, EspressoProfile
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
    from ase.optimize import BFGS
    p=parser(__doc__); p.add_argument('structure',help='Periodic CIF/extxyz of a sensible molecular cell')
    p.add_argument('pseudopotentials'); p.add_argument('--pseudo-dir',required=True)
    p.add_argument('--command',default='pw.x',help='e.g. mpiexec -n 8 pw.x')
    p.add_argument('--schedule',type=positive,nargs='+',required=True,help='Temperature stages K, e.g. 600 500 400 300')
    p.add_argument('--steps-per-stage',type=int,required=True)
    p.add_argument('--timestep-fs',type=positive,default=0.5)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--ecutwfc',type=positive,default=80)
    p.add_argument('--ecutrho',type=positive,default=640)
    a=p.parse_args()
    if a.steps_per_stage<1: p.error('steps-per-stage must be positive')
    atoms=read(a.structure)
    if not all(atoms.pbc) or len(atoms)<2: raise ValueError('Need a periodic multi-atom cell')
    out=output_dir(a).resolve()
    if (out/'md.traj').exists(): raise FileExistsError('Use fresh output directory')
    pp=read_json(a.pseudopotentials); ppdir=Path(a.pseudo_dir).resolve()
    for s in set(atoms.get_chemical_symbols()):
        if s not in pp or not (ppdir/pp[s]).is_file(): raise FileNotFoundError(f'Missing PP: {s}')
    atoms.calc=Espresso(profile=EspressoProfile(command=a.command,pseudo_dir=str(ppdir)),
        directory=str(out/'qe'),pseudopotentials=pp,kpts=(1,1,1),
        input_data={'control':{'tprnfor':True},'system':{'ecutwfc':a.ecutwfc,'ecutrho':a.ecutrho,
                    'input_dft':'PBE','vdw_corr':'grimme-d3','occupations':'fixed'},
                    'electrons':{'conv_thr':1e-8,'mixing_beta':0.3}})
    rng=np.random.default_rng(a.seed)
    MaxwellBoltzmannDistribution(atoms,temperature_K=a.schedule[0],rng=rng); Stationary(atoms)
    dyn=Langevin(atoms,a.timestep_fs*units.fs,temperature_K=a.schedule[0],
        friction=0.01/units.fs,rng=rng,trajectory=str(out/'md.traj'),logfile=str(out/'md.log'))
    for t in a.schedule:
        dyn.set_temperature(temperature_K=t); dyn.run(a.steps_per_stage)
        write(out/f'stage_{t:g}K.extxyz',atoms)
    opt=BFGS(atoms,logfile=str(out/'quench.log')); converged=opt.run(fmax=0.03,steps=500)
    write(out/'quenched.cif',atoms)
    dump(out/'metadata.json',dict(converged=bool(converged),final_energy_cell_eV=float(atoms.get_potential_energy()),
        temperatures_K=a.schedule,steps_per_stage=a.steps_per_stage,seed=a.seed,
        warning='Trajectory sampling only; prove melting and retained chemistry, assess RDF/diffusion and independent replicas/density/size effects. Short quench does not establish an amorphous phase or glass stability. High T may decompose organics.'))
    if not converged: raise RuntimeError('Final quench did not converge')

if __name__ == '__main__': main()
