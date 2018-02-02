# pdbcolor
Python code for PyMOL to color a PDB structure based on various parameters obtained from a multiple sequence alignment. The scripts could be modified to accept essentially any parameter that can be obtained for each residue, say for dN/dS ratios, hydrophobicity, etc. and could be changed to represent any arbitrary value series for however many colors are needed.

For all cases, the scripts work by changing a value for each ATOM record in a PDB file. In a normal PDB file, the temperatureFactor or beta-factor is the second to last term, here in the first atom it is 0.82.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00  0.82           N  `

Each script identifies the residue number (here alanine is the 11th residue of the protein sequence, but the first in the model), and changes the temperatureFactor for each atom of that residue. For instance, in the `pdb_site_identity.py` script, the value would be recoded to 5.00, indicating reasonably high identity (>80%) at this position, and this alanine will be colored green.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00  5.00           N  `

At the moment, scripts can only recode a single protein per PDB file, meaning hetero-multimers are not co-colored. This may be developed in the future.

## percent identity ##
For a target sequence, recolor by percent identity from a multiple sequence alignment. By default, gray indicates less than 50% identity, following reverse rainbow order (so blue to red) to show increasing identity, with magenta showing 100% identity (excluding gaps or missing data). 

![percent_identity_color_scheme.png](https://github.com/wrf/pdbcolor/blob/master/percent_identity_color_scheme.png)

Colors are defined in the `color_by_identity.py` script, where triplets are RGB values of 0 to 1 (so 1,1,1 is white). Thresholds of identity are defined in the `pdb_site_identity.py` script, where the scores are calculated and printed for each of the ATOM records in the PDB file.

These values range from 1.00 to 9.00, for a gradient of percentage points. The colors are meant to specifically highlight strongly conserved sites, which is why site identity below 50% is colored gray.

### Instructions ###
1) Generate a multiple sequence alignment, which must include the sequence of the target PDB structure. For instance, the protein [Symplectin](https://bitbucket.org/wrf/squid-transcriptomes/src) is the target, and residues will be colored based on the percentage of residues in other sequences that are identical to this sequence. Here they are aligned with the program `mafft`.

`mafft-linsi examples/symplectins_w_outgroups.fasta > examples/symplectins_w_outgroups.aln`

2) Using the alignment and the PDB file (here generated using the [SWISS-MODEL](https://www.swissmodel.expasy.org/) server), reassign the temperatureFactors based on conservation scores using `pdb_site_identity.py`. The target sequence is indicated using the `-s` option, and must exactly match the name given in the alignment. The script assumes the residues in the PDB file are numbered based on this sequence, meaning even if the PDB file starts with residue 20 (say the first 19 were disordered or cleaved), the sequence still must start with residue 1.

`./pdb_site_identity.py -p examples/symplectin_swissmodel_01.pdb -s Symplectin -a examples/symplectins_w_outgroups.aln > examples/symplectin_swissmodel_01_w_id.pdb`

3) Open the PDB file in PyMOL:

`pymol examples/symplectin_swissmodel_01_w_id.pdb`

4) In the PyMOL console, run the command:

`run color_by_identity.py`

![symplectin_domains_by_conservation.png](https://github.com/wrf/pdbcolor/blob/master/symplectin_domains_by_conservation.png)

In the left panel, certain residues in the binding pocket are strongly conserved (shown in red, >95%), such as the catalytic triad of E-K-C (though C is green, meaning only >80% identity, as this is sometimes serine). In the right panel, the other domain is not well conserved outside of disulfide-forming cysteines.

## RAxML site-wise likelihood ##
For an alignment and a series of phylogenetic trees with fixed topologies, [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) can produce a table of site-wise log-likelihoods for each tree topology (using the `-f G` option). The difference between two topologies can be computed for each site (the "dlnL"). While most sites are marginally different (values between -0.5 and 0.5, colored rust and mud), a small number of sites have a strong effect on the likelihood of the final tree (colored in bright pink and teal).

![likelihood_color_scheme_v1.png](https://github.com/wrf/pdbcolor/blob/master/likelihood_color_scheme_v1.png)

## heteropecilly ##
See the [heteropecilly github repo](https://github.com/wrf/heteropecilly)

![heteropecilly_color_scheme.png](https://github.com/wrf/pdbcolor/blob/master/heteropecilly_color_scheme.png)

## References ##
The colorization script was modified from the `consurf_new.py` script from the [ConSurf Server](http://consurf.tau.ac.il/2016/)
