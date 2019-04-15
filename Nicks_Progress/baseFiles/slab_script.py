import sys
from espresso import espresso
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.spacegroup import crystal
from ase.build import surface, fcc111, add_adsorbate
from ase.io import write
from ase.io import read

NUM_PW  = sys.argv[1]
NUM_KPT = sys.argv[2]

#Building the Platinum crystal structure
a = 3.923 #angstrom
Pt = crystal(['Pt'],basis=[(0,0,0)],spacegroup=225,cellpar=[a,a,a,90,90,90])

# Building a (111) Platinum surface slab
Pt_111 = surface(Pt, (1, 1, 1), 2)
Pt_111.center(vacuum=10,axis=2)
Pt_111_repeat = Pt_111.repeat((1,1,1))

constraint = FixAtoms(mask=[False,True,True,False,False,True,True,False])
Pt_111_repeat.set_constraint(constraint)

Pt_111_repeat.calc = espresso(pw=NUM_PW,
                            dw=NUM_PW*10,
                            kpts=(NUM_KPT,NUM_KPT,1),
                            xc='PBE',
                            spinpol=True,
                            outdir='relaxed_pt',
                            convergence={'energy':1e-6,
                                         'mixing':0.05,
                                         'mixing_mode':'local-TF',
                                         'maxsteps':1000,
                                         'diag':'david'})

dyn = QuasiNewton(Pt_111_repeat, trajectory='modified_slab.traj')
dyn.run(fmax=3)

modified_slab = read('modified_slab.traj')

conv_dict = {'energy':1e-6,
             'mixing':0.05,
             'mixing_mode':'local-TF',
             'maxsteps':1000,
             'diag':'cg'}

modified_slab.calc = espresso(pw=NUM_PW,
                        dw=NUM_PW*10,
                        kpts=(NUM_KPT,NUM_KPT,1),
                        xc='PBE',
                        outdir='pt_output',
                        convergence=conv_dict)

Pt_energy = modified_slab.get_potential_energy()

print("Slab energy: " + str(Pt_energy))

write('Pt_calc.traj', modified_slab)