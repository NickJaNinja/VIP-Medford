import sys

from ase import Atoms
from espresso import espresso
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import molecule, surface, fcc111, add_adsorbate
from ase.spacegroup import crystal
from ase.visualize import view
from ase.io import write
from ase.io import read

CO_path = 'CO_calc.traj'
Pt_path = 'Pt_calc.traj'

NUM_PW  = sys.argv[1]
NUM_KPT = sys.argv[2]

CO_molecule = read(CO_path)
Pt_slab = read(Pt_path)

# height of the gas above the slab (in angstroms)
height = 3

constraint = FixAtoms(mask=[a.position[2] <= 11 for a in Pt_slab])
Pt_slab.set_constraint(constraint)

add_adsorbate(Pt_slab, CO_molecule, height, position=(1, 1))

Pt_slab.calc = espresso(pw=NUM_PW,
                        dw=4500,
                        kpts=(NUM_KPT,NUM_KPT,1),
                        xc='PBE',
                        outdir='CO_pt',
                        convergence={'energy':1e-6,
                                     'mixing':0.05,
                                     'mixing_mode':'local-TF',
                                     'maxsteps':1000,
                                     'diag':'cg'})

#dyn = QuasiNewton(Pt_slab, trajectory='PtCO_complete.traj')
#dyn.run(fmax=0.05)
E_slab_pt = Pt_slab.get_potential_energy()

print('System energy after: ', E_slab_pt)