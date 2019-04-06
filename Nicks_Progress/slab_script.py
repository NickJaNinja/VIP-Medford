from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import surface, fcc111, add_adsorbate
from ase.spacegroup import crystal
from ase.visualize import view
from ase.io import write
from ase.io import read

Pt_path = 'Pt_111.traj'

#Building the Platinum crystal structure
a = 3.923 #angstrom
Pt = crystal(['Pt'],basis=[(0,0,0)],spacegroup=225,cellpar=[a,a,a,90,90,90])

# Building a (111) Platinum surface slab
Pt_111 = surface(Pt, (1, 1, 1), 2)
Pt_111.center(vacuum=10,axis=2)
Pt_111_repeat = Pt_111.repeat((1,1,1))

dyn = QuasiNewton(Pt_111_repeat, trajectory=Pt_path)
dyn.run(fmax=3)

Pt_molecule = read(Pt_path)

Pt_molecule.calc = espresso(pw=400,
                            dw=4500,
                            kpts=(4,4,4),
                            xc='PBE',
                            outdir='CO_pt',
                            convergence={'energy':1e-6,
                                         'mixing':0.05,
                                         'mixing_mode':'local-TF',
                                         'maxsteps':1000,
                                         'diag':'cg'})

#dyn = QuasiNewton(Pt_molecule, trajectory='Pt_calc.traj')
#dyn.run(fmax=0.05)

e_Pt = Pt_molecule.get_potential_energy()
print("Pt potential energy: ", e_pt)

write('Pt_calc.traj', Pt_molecule)