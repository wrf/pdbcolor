# pdbcolor
Python code to color a PDB structure based on percent identity from a multiple sequence alignment

## Instructions ##
1) Generate a multiple sequence alignment, which must include the sequence of the target PDB structure. For instance, the protein [Symplectin](https://bitbucket.org/wrf/squid-transcriptomes/src) is the target, and residues will be colored based on the percentage identity of this sequence.

`mafft-linsi symplectin_w_outgroups.fasta > symplectin_w_outgroups.aln`

2) Using the alignment and the PDB file, assign the temperatureFactors to conservation scores using `pdb_site_identity.py`:

`pdb_site_identity.py -p symplectin_model.pdb -s Symplectin -a symplectin_w_outgroups.aln > symplectin_model_w_scores.pdb`

3) Open the PDB file in PyMOL:

`pymol symplectin_model_w_scores.pdb`

4) In the PyMOL console, run the command:

`run color_by_identity.py`

By default, gray indicates less than 50% identity, following reverse rainbow order (so blue to red) to show increasing identity, with magenta showing 100% identity (excluding gaps or missing data). Colors are defined in the `color_by_identity.py` script, where triplets are RGB values of 0 to 1 (so 1,1,1 is white). Thresholds of identity are actually determined in the `pdb_site_identity.py` script, where the scores are then printed in the ATOM records in the PDB file.

![symplectin_domains_by_conservation.png](https://github.com/wrf/pdbcolor/blob/master/symplectin_domains_by_conservation.png)

In the left panel, certain residues in the binding pocket are strongly conserved (shown in red, >95%), such as the catalytic triad of E-K-C (though C is green, meaning only >80% identity, as this is sometimes serine). Comparatively, the other domain is not well conserved outside of disulfide-forming cysteines.

## References ##
The colorization script was modified from the `consurf_new.py` script from the [ConSurf Server](http://consurf.tau.ac.il/2016/)
