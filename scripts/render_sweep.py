import sys
import os
import pymol
from pymol import cmd

def main():
    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()

    data_dir = os.path.join(os.path.dirname(__file__), "..", "data", "ligands")
    sdf_path = os.path.join(data_dir, "df_aligned_sweep.sdf")
    out_png = os.path.join(data_dir, "df_sweep_visualization.png")
    
    cmd.load(sdf_path, "df_aligned_sweep")
    cmd.set("all_states", 1)
    cmd.hide("everything")
    cmd.show("sticks", "df_aligned_sweep")
    cmd.set("valence", 1)
    
    # Define anatomical parts
    cmd.select("head", "id 27+28+29+30+31+32+33+34+35+36+22+23+24+25+26")
    cmd.select("tail", "id 1-21")
    
    # 1. Head (Ring + Ester): vivid green/red, completely solid
    cmd.color("green", "head and elem C")
    cmd.color("red", "head and elem O")
    cmd.set("stick_transparency", 0.0, "head")
    
    # 2. Sweeping Alkane Tail: palegreen, semi-transparent
    cmd.color("palegreen", "tail and elem C")
    cmd.set("stick_transparency", 0.6, "tail")
    
    # Center camera on the entire rigid head
    cmd.orient("head")
    
    view = cmd.get_view()
    rot = view[0:9]
    tail_com = cmd.centerofmass("tail")
    anchor_com = cmd.centerofmass("head")
    tail_x = rot[0]*tail_com[0] + rot[1]*tail_com[1] + rot[2]*tail_com[2]
    anchor_x = rot[0]*anchor_com[0] + rot[1]*anchor_com[1] + rot[2]*anchor_com[2]
    
    if tail_x < anchor_x:
        cmd.turn("z", 180)
        
    cmd.zoom("df_aligned_sweep", buffer=2.5)
    cmd.bg_color("white")
    
    cmd.png(out_png, width=1920, height=1080, dpi=300, ray=1)
    cmd.quit()

if __name__ == "__main__":
    main()
