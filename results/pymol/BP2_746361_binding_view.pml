reinitialize
load data/processed/FtsZ_receptor_noH.pdb, receptor
load work/docking/BP2/BP2_746361_vina_out.pdbqt, ligand
hide everything
show cartoon, receptor
color slate, receptor
show sticks, ligand
color yelloworange, ligand
select contact_residues, (chain B and resi 306) or (chain B and resi 305) or (chain B and resi 260) or (chain B and resi 189) or (chain B and resi 304) or (chain B and resi 196) or (chain B and resi 262) or (chain B and resi 225) or (chain B and resi 185)
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
png results/figures/BP2_746361_pymol_binding_view.png, width=1800, height=1400, dpi=250, ray=1
# Contact residue labels: THR306B+VAL305B+SER260B+ASN189B+ARG304B+ASP196B+ALA262B+ILE225B+GLU185B
