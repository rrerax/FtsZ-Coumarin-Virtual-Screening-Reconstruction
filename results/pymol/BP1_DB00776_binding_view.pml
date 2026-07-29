reinitialize
load data/processed/FtsZ_receptor_noH.pdb, receptor
load work/docking/BP1/BP1_DB00776_vina_out.pdbqt, ligand
hide everything
show cartoon, receptor
color slate, receptor
show sticks, ligand
color yelloworange, ligand
select contact_residues, (chain B and resi 47) or (chain B and resi 55) or (chain B and resi 57) or (chain B and resi 50) or (chain B and resi 56) or (chain B and resi 49) or (chain B and resi 39) or (chain B and resi 41) or (chain B and resi 48)
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
png results/figures/BP1_DB00776_pymol_binding_view.png, width=1800, height=1400, dpi=250, ray=1
# Contact residue labels: LEU47B+LYS55B+ASP57B+SER50B+LEU56B+MET49B+ALA39B+ASN41B+LEU48B
