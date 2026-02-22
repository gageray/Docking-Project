import sys
import os
import pymol
from pymol import cmd

def main():
    # Start PyMOL without GUI and in quiet mode
    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()

    data_dir = os.path.join(os.path.dirname(__file__), "..", "data", "ligands")
    sdf_path = os.path.join(data_dir, "df_aligned_sweep.sdf")
    out_png = os.path.join(data_dir, "df_sweep_visualization.png")
    
    print(f"Loading {sdf_path}...")
    cmd.load(sdf_path, "df_aligned_sweep")
    
    # Apply rendering settings
    cmd.set("all_states", 1)
    cmd.hide("everything")
    cmd.show("sticks", "df_aligned_sweep")
    cmd.color("atomic", "df_aligned_sweep")
    cmd.select("tail", "not id 29+30+31+32+33+34+35+36+28+27") # selection of everything but the rigid anchor
    cmd.color("gray50", "tail")
    cmd.set("stick_transparency", 0.7, "tail")
    
    # Orient the camera on the anchor. 
    # Because it is a flat ring, orient aligns its plane with the XY camera plane (face-on).
    cmd.orient("not tail")
    
    # Ensure the tail points towards positive X (to the right)
    view = cmd.get_view()
    rot = view[0:9] # 3x3 rotation matrix for Model -> Camera
    
    tail_com = cmd.centerofmass("tail")
    anchor_com = cmd.centerofmass("not tail")
    
    # Calculate screen X coordinate for tail and anchor vectors
    tail_x = rot[0]*tail_com[0] + rot[1]*tail_com[1] + rot[2]*tail_com[2]
    anchor_x = rot[0]*anchor_com[0] + rot[1]*anchor_com[1] + rot[2]*anchor_com[2]
    
    if tail_x < anchor_x:
        cmd.turn("z", 180)  # Flip it strictly in the 2D plane so the ring stays face-on, but tail points right
        
    cmd.zoom("df_aligned_sweep", buffer=2.5)
    cmd.bg_color("white")
    
    # Render and save the image
    print(f"Rendering image to {out_png}...")
    cmd.png(out_png, width=1920, height=1080, dpi=300, ray=1)
    print("Done rendering.")
    
    cmd.quit()

if __name__ == "__main__":
    main()
