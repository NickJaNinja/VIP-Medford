from espresso import espresso
from ase import Atoms
from ase.build import molecule
from ase.io import write
import sys

NUM_PW = sys.argv[1]

CO_molecule = molecule('CO')
CO_molecule.center(10)

CO_molecule.calc = espresso(pw=NUM_PW,
                                dw=NUM_PW*10,
                                kpts=(1,1,1),
                                xc='PBE',
                                outdir='e_CO',
                                convergence={'energy':1e-6,
                                                 'mixing':0.05,
                                                 'mixing_mode':'local-TF',
                                                 'maxsteps':1000,
                                                 'diag':'cg'})


e_CO = CO_molecule.get_potential_energy()
print("CO potential energy: ", str(e_CO))

write('CO_calc.traj', CO_molecule)