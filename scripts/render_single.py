import sys
import os
import pymol
from pymol import cmd

def main():
    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()

    data_dir = os.path.join(os.path.dirname(__file__), "..", "data", "ligands")
    sdf_path = os.path.join(data_dir, "docosanyl_ferulate_conf1.sdf")
    out_png = os.path.join(data_dir, "df_single_molecule.png")
    
    cmd.load(sdf_path, "df_single")
    cmd.hide("everything")
    cmd.show("sticks", "df_single")
    cmd.set("valence", 1)
    
    # Render single molecule full bright
    cmd.color("green", "elem C")
    cmd.color("red", "elem O")
    cmd.set("stick_transparency", 0.0, "df_single")
    
    # Define ester for orientation
    cmd.select("ester", "id 22+23+24+25+26")
    
    # Orient strictly on the ester so the C=O bond is forced into the XY plane (parallel to viewing plane)
    cmd.orient("ester")
    
    # Orient algorithm might flip it arbitrary ways around Z, 
    # let's just make sure it fills the screen
    cmd.zoom("df_single", buffer=2.0)
    cmd.bg_color("white")
    
    cmd.png(out_png, width=1920, height=1080, dpi=300, ray=1)
    cmd.quit()

if __name__ == "__main__":
    main()
