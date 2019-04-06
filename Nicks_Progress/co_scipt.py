from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import surface, fcc111, add_adsorbate
from ase.spacegroup import crystal
from ase.visualize import view
from ase.io import write
from ase.io import read

CO_molecule = molecule('CO')
CO_molecule.center(10)

CO_molecule.calc = espresso(pw=400,
                                dw=4500,
                                kpts=(1,1,1),
                                xc='PBE',
                                outdir='CO_pt',
                                convergence={'energy':1e-6,
                                                 'mixing':0.05,
                                                 'mixing_mode':'local-TF',
                                                 'maxsteps':1000,
                                                 'diag':'cg'})

write('CO_calc.traj', CO_molecule)

#dyn = QuasiNewton(CO_molecule, trajectory='CO_calc.traj')
#dyn.run(fmax=0.05)

e_CO   = CO_molecule.get_potential_energy()
print("CO potential energy: ", e_CO)
