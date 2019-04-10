import sys

from espresso import espresso
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import molecule
from ase.spacegroup import crystal
from ase.visualize import view
from ase.io import write
from ase.io import read

NUM_PW = sys.argv[1]

CO_molecule = molecule('CO')
CO_molecule.center(10)

CO_molecule.calc = espresso(pw=NUM_PW,
                                dw=4500,
                                kpts=(1,1,1),
                                xc='PBE',
                                outdir='e_CO',
                                convergence={'energy':1e-6,
                                                 'mixing':0.05,
                                                 'mixing_mode':'local-TF',
                                                 'maxsteps':1000,
                                                 'diag':'cg'})


#dyn = QuasiNewton(CO_molecule, trajectory='CO_calc.traj')
#dyn.run(fmax=0.05)

e_CO = CO_molecule.get_potential_energy()
print("CO potential energy: ", e_CO)

write('CO_calc.traj', CO_molecule)