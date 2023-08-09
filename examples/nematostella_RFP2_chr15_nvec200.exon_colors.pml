hide everything
bg white
show cartoon, (chain A)
set_color excolor_1, [0.65,0.81,0.89]
select exon_1, (resi 1-1) & (chain A)
color excolor_1, exon_1
set_color excolor_2, [0.12,0.47,0.71]
select exon_2, (resi 2-34) & (chain A)
color excolor_2, exon_2
set_color excolor_3, [0.7,0.87,0.54]
select exon_3, (resi 36-108) & (chain A)
color excolor_3, exon_3
set_color excolor_4, [0.2,0.63,0.17]
select exon_4, (resi 110-178) & (chain A)
color excolor_4, exon_4
set_color excolor_5, [0.98,0.6,0.6]
select exon_5, (resi 180-238) & (chain A)
color excolor_5, exon_5
set_color outphase, [0.3,0.3,0.3]
select out_of_phase, (resi 35,109,179) & (chain A)
color outphase, out_of_phase
