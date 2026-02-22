import pymol
from pymol import cmd
import numpy as np

pymol.pymol_argv = ['pymol', '-c', '-q']
pymol.finish_launching()
cmd.load("../data/ligands/df_aligned_sweep.sdf", "mol")

cmd.orient("not tail")
view = cmd.get_view()
print("View matrix:", view)

# Let's project tail center of mass
tail_com = cmd.centerofmass("tail")
anchor_com = cmd.centerofmass("not tail")

print("Model Tail COM:", tail_com)
print("Model Anchor COM:", anchor_com)

# To get screen coordinates, we apply the 3x3 rotation matrix (first 9 elements of view)
# PyMOL view array:
# 0-8: 3x3 rotation matrix (row-major)
# 9-11: translation
# 12-14: camera center
rot = np.array(view[0:9]).reshape(3,3)
tail_scr = np.dot(rot, tail_com)
anchor_scr = np.dot(rot, anchor_com)
print("Screen Tail COM:", tail_scr)
print("Screen Anchor COM:", anchor_scr)
print("Delta:", tail_scr - anchor_scr)

cmd.quit()
