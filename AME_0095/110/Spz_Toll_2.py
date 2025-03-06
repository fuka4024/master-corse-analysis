from pymol import cmd
import csv
import math

# タンパク質のロードと選択
cmd.load("C:/Users/monom/Downloads/4bv4.cif", "Spz_Toll")
cmd.select("Spz", "chain L or chain M")
cmd.select("Toll", "chain R")

def highlight_binding_sites(selection1, selection2, cutoff=4.0):
                    cmd.select("binding_sites", f"({selection1}) within {cutoff} of ({selection2})") 
highlight_binding_sites("Toll", "Spz", cutoff=4.0)
# PyMOL画像保存
cmd.color("blue", "Spz")
cmd.color("green", "Toll")
cmd.show("spheres", "binding_sites")  # スフィアで表示
cmd.color("orange", "binding_sites")
cmd.show("cartoon", "Spz or Toll")
cmd.hide("labels")
cmd.ray()
cmd.png("C:/Users/monom/Documents/spz_Toll/DmSpz_Toll.png")