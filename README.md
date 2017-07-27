# pdbcolor
Python code for PyMOL to color a PDB structure based on percent identity from a multiple sequence alignment. By default, gray indicates less than 50% identity, following reverse rainbow order (so blue to red) to show increasing identity, with magenta showing 100% identity (excluding gaps or missing data). 

![percent_identity_color_scheme.png](https://github.com/wrf/pdbcolor/blob/master/percent_identity_color_scheme.png)

Colors are defined in the `color_by_identity.py` script, where triplets are RGB values of 0 to 1 (so 1,1,1 is white). Thresholds of identity are defined in the `pdb_site_identity.py` script, where the scores are calculated and printed for each of the ATOM records in the PDB file.

In a normal PDB file, the temperatureFactor or beta-factor is the second to last term, here in the first atom it is 0.82.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00  0.82           N  `

The script identifies the residue number (here alanine is the 11th residue of the protein sequence, but the first in the model), and changes the temperatureFactor for each atom of that residue. Here it is 5.00, indicating reasonably high identity (>80%) at this position, and this alanine will be colored green.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00  5.00           N  `

These values range from 1.00 to 9.00, though can be changed to represent any arbitrary value for however many colors are needed. Any other score ranging from 1 to 9 could also be used, say for dN/dS ratios, hydrophobicity, etc.

## Instructions ##
1) Generate a multiple sequence alignment, which must include the sequence of the target PDB structure. For instance, the protein [Symplectin](https://bitbucket.org/wrf/squid-transcriptomes/src) is the target, and residues will be colored based on the percentage of residues in other sequences that are identical to this sequence. Here they are aligned with the program `mafft`.

`mafft-linsi symplectin_w_outgroups.fasta > symplectin_w_outgroups.aln`

2) Using the alignment and the PDB file (here generated using the [SWISS-MODEL](https://www.swissmodel.expasy.org/) server), reassign the temperatureFactors based on conservation scores using `pdb_site_identity.py`. The target sequence is indicated using the `-s` option, and must exactly match the name given in the alignment. The script assumes the residues in the PDB file are numbered based on this sequence, meaning even if the PDB file starts with residue 20 (say the first 19 were disordered or cleaved), the sequence still must start with residue 1.

`pdb_site_identity.py -p symplectin_model.pdb -s Symplectin -a symplectin_w_outgroups.aln > symplectin_model_w_scores.pdb`

3) Open the PDB file in PyMOL:

`pymol symplectin_model_w_scores.pdb`

4) In the PyMOL console, run the command:

`run color_by_identity.py`

![symplectin_domains_by_conservation.png](https://github.com/wrf/pdbcolor/blob/master/symplectin_domains_by_conservation.png)

In the left panel, certain residues in the binding pocket are strongly conserved (shown in red, >95%), such as the catalytic triad of E-K-C (though C is green, meaning only >80% identity, as this is sometimes serine). In the right panel, the other domain is not well conserved outside of disulfide-forming cysteines.

## References ##
The colorization script was modified from the `consurf_new.py` script from the [ConSurf Server](http://consurf.tau.ac.il/2016/)
