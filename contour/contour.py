from pathlib import Path
from math import ceil, floor
from random import choice
from itertools import combinations, permutations, product

import MDAnalysis as mda

salida_path = Path("/home/pbarletta/labo/22/locuaz/daux/cacho/0-B_KASRSAGGALEGARGQ")

u = mda.Universe(str(salida_path / "pre_AmberPDBFixer_init_nonoverlapped.pdb"))

for a in u.select_atoms("name NA or name Na or name Cl or name CL"):
    print(a)

# Remove ions to later replace them
ions = u.select_atoms("name NA or name Na or name Cl or name CL")
v = u.atoms - ions

#
H = get_matrix(v.dimensions) * 10
inv_H = np.linalg.inv(H)
centro = np.sum(H / 2, axis=1)

waters = u.select_atoms("resname WAT or resname SOL and type O")
# Get transformed coordinates of the water oxygens
s_positions = (waters.positions - centro) @ inv_H
assert np.all(s_positions < .5) and np.all(s_positions > -.5), "Bad box transformation"

# Insert `nions` ions
nions = 17
splits = np.linspace(-.4, .4, nsplits)
corners = list(product(splits, splits, splits))
assert len(corners) > nions
used_corners = set()

idx_wats = []
center = tuple([0, 0, 0])
remaining_ions = nions
while remaining_ions != 0:
    # Get random corner
    corner = choice(corners)
    if corner in used_corners or np.allclose(corner, center):
        continue
    used_corners.add(corner)
    remaining_ions -= 1
    # Get the idx (np.argmin) of the oxygen that is closest to the corner coordinates.
    idx_wats.append(np.argmin(np.sum((s_positions - corner)**2, axis = 1)))

# These are the contouring water oxygens to replace with ions
wat_to_replace = { waters[idx] for idx in idx_wats }

# TODO: remove the selected waters, like:
# for wat_atmgroup in ...
#   u.atoms. -= wat_atmgroup

# Finally, create the universe with the Ions
uni_ions = mda.Universe.empty(
    nions, n_residues=nions, atom_resindex=[0]*nions, residue_segindex=[0]*nions, trajectory=True)

uni_ions.add_TopologyAttr("name", ["CL"] * nions)
uni_ions.add_TopologyAttr("type", ["CL"] * nions)
uni_ions.add_TopologyAttr("resname", ["CL"] * nions)
uni_ions.add_TopologyAttr("resid", list(range(nions)))

uni_ions.atoms.positions = np.array([ wat.position for wat in wat_to_replace ])

# output universe
mda.Merge(u.atoms, uni_ions.atoms).atoms.write("a.pdb")
