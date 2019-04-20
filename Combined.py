from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import surface, fcc111, add_adsorbate, molecule, surface
from espresso import espresso
from ase.spacegroup import crystal
from ase.visualize import view
from ase.io import write
from ase.io import read
from ase.spacegroup import crystal
import sys

current_pw = int(500)
current_k = 8

a=3.923 #angstrom
Pt=crystal(['Pt'],basis=[(0,0,0)],spacegroup=225,cellpar=[a,a,a,90,90,90])

h = 2.6

conv_dict = {'energy':1e-6,
            'mixing':0.05,
            'mixing_mode':'local-TF',
            'maxsteps':1000,
            'diag':'cg'}

Pt_slab=surface(Pt,(1,1,1),2)
Pt_slab.center(vacuum=10,axis=2)
CO_molecule = molecule('CO')
CO_molecule.center(10)

add_adsorbate(Pt_slab, CO_molecule, h, position=(2.5, 2.5))

constraint = FixAtoms(mask=[True, True, True, True, False, False, False, False])
Pt_slab.set_constraint(constraint)

Pt_slab.calc=espresso(pw=current_pw,
                       dw=current_pw*10,
                       kpts=(current_k,current_k,1),
                       xc='PBE',
                       outdir='test_output',
                       convergence=conv_dict)
dyn = QuasiNewton(Pt_slab, trajectory='Pt+CO.traj')
dyn.run(fmax=0.05)

combined_energy = Pt_slab.get_potential_energy()

print(combined_energy)
f = open("Comb_" +str( current_pw) + "_" +str( current_k) + ".txt", "w")
f.write(str(combined_energy) + '\n')
f.close()
