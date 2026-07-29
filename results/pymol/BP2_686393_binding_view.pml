reinitialize
load data/processed/FtsZ_receptor_noH.pdb, receptor
load work/docking/BP2/BP2_686393_vina_out.pdbqt, ligand
hide everything
show cartoon, receptor
color slate, receptor
show sticks, ligand
color yelloworange, ligand
select contact_residues, (chain B and resi 185) or (chain B and resi 182) or (chain B and resi 304) or (chain B and resi 169) or (chain B and resi 186) or (chain B and resi 188) or (chain B and resi 227) or (chain B and resi 189) or (chain B and resi 225)
show sticks, contact_residues
color marine, contact_residues
show surface, byres contact_residues around 4 of ligand
set transparency, 0.55
distance contacts, ligand, contact_residues, 4.0
set dash_color, orange
set label_size, 18
set ray_opaque_background, off
bg_color white
orient ligand or contact_residues
zoom ligand or contact_residues, 4
png results/figures/BP2_686393_pymol_binding_view.png, width=1800, height=1400, dpi=250, ray=1
# Contact residue labels: GLU185B+SER182B+ARG304B+MET169B+VAL186B+LEU188B+SER227B+ASN189B+ILE225B
