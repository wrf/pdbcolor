hide everything
bg white
show cartoon, (chain A or chain B)
set_color excolor_1, [0.65,0.81,0.89]
select exon_1, (resi 1-29) & (chain A or chain B)
color excolor_1, exon_1
set_color excolor_2, [0.12,0.47,0.71]
select exon_2, (resi 31-60) & (chain A or chain B)
color excolor_2, exon_2
set_color excolor_3, [0.7,0.87,0.54]
select exon_3, (resi 62-93) & (chain A or chain B)
color excolor_3, exon_3
set_color excolor_4, [0.2,0.63,0.17]
select exon_4, (resi 94-130) & (chain A or chain B)
color excolor_4, exon_4
set_color excolor_5, [0.98,0.6,0.6]
select exon_5, (resi 132-181) & (chain A or chain B)
color excolor_5, exon_5
set_color excolor_6, [0.89,0.1,0.11]
select exon_6, (resi 183-227) & (chain A or chain B)
color excolor_6, exon_6
set_color excolor_7, [0.99,0.75,0.44]
select exon_7, (resi 228-275) & (chain A or chain B)
color excolor_7, exon_7
set_color excolor_8, [1.0,0.5,0.0]
select exon_8, (resi 276-328) & (chain A or chain B)
color excolor_8, exon_8
set_color excolor_9, [0.79,0.7,0.84]
select exon_9, (resi 330-384) & (chain A or chain B)
color excolor_9, exon_9
set_color excolor_10, [0.42,0.24,0.6]
select exon_10, (resi 385-414) & (chain A or chain B)
color excolor_10, exon_10
set_color excolor_11, [0.9,0.65,0.5]
select exon_11, (resi 416-453) & (chain A or chain B)
color excolor_11, exon_11
set_color outphase, [0.3,0.3,0.3]
select out_of_phase, (resi 30,61,131,182,329,415) & (chain A or chain B)
color outphase, out_of_phase
