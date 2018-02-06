# pdbcolor
Python code for [PyMOL](https://pymol.org/2/) to color a [PDB structure](http://www.rcsb.org/) based on various parameters obtained from a multiple sequence alignment. The scripts could be modified to accept essentially any parameter that can be obtained for each residue, say for dN/dS ratios, hydrophobicity, etc. and could be changed to represent any arbitrary value series for however many colors are needed.

For all cases, the scripts work by changing a value for [each ATOM record in a PDB file](http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences-and-the-pdb-format). In a normal PDB file, the temperatureFactor or beta-factor is the second to last term, here in the first atom it is 0.82.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00  0.82           N  `

Each script identifies the residue number (here alanine is the 11th residue of the protein sequence, but the first in the model), and changes the temperatureFactor for each atom of that residue. For instance, in the `pdb_site_identity.py` script, the value would be recoded to 5.00, indicating reasonably high identity (>80%) at this position, and this alanine will be colored green.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00  5.00           N  `

**At the moment, scripts can only recode a single protein per PDB file, meaning hetero-multimers are not co-colored.** If a PDB file contains multiple other proteins, metals, ligands, DNA, etc., all of these will be colored the "null" color. This may be developed in the future. Running the `color_by_` scripts within PyMOL recolors all atoms, though the color of any of these can be manually changed afterwards, perhaps to highlight ligands or specific domains.

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
For an alignment and a series of phylogenetic trees with fixed topologies, [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) can produce a table of site-wise log-likelihoods for each tree topology (using the `-f G` option). The difference between two topologies can be computed for each site (the "dlnL"), showing which sites contribute strongly to one topology or the other. These dlnL values are recoded to a hex string where a dlnL of 0 becomes exactly 8.0 to include the range of values but also permit each value to be coded as a single character (for display as a fasta sequence). 

While most sites are marginally different - values between -0.5 and 0.5 - colored mud (7) and rust (8), a small number of sites have a strong effect on the likelihood of the final tree, colored in blue/teal (0-6) for sites that favor the second tree and bright pink (9-f) for sites favoring the first tree. The scripts to make use of this assume that proteins and site-wise likelihood values derive from a protein supermatrix. Because proteins in supermatrices are often trimmed compared to the complete sequence, any missing values are colored gray as gaps. Sites that were constant in the alignment need to be recoded as such, and are then colored green in the resulting protein.

![likelihood_color_scheme_v1.png](https://github.com/wrf/pdbcolor/blob/master/likelihood_color_scheme_v1.png)

Detailed instructions are currently in the [heteropecilly github repo](https://github.com/wrf/heteropecilly)

The workflow is meant to begin from a supermatrix, meaning a lot of other data need to be generated and organized and assumes that most proteins in the supermatrix are trimmed or are incomplete. Briefly, the steps required are:

1) Starting from a protein supermatrix ( see [supermatrix repo for scripts to manipulate matrices](https://github.com/wrf/supermatrix) ), run `RAxML` to get sitewise likelihoods of the entire dataset, then convert that output to tabular using the `sitewise_ll_to_columns.py` script.
2) Because constant sites may still receive lnL values from RAxML (which is biologically meaningless) these should be recoded to a value indicating they are constant using the `sitewise_recode_constant_sites.py` [script](https://github.com/wrf/heteropecilly/blob/master/sitewise_recode_constant_sites.py). Using the matrix and the tabular likelihoods, recode constant sites and generate a new tabular output file that can be used in step 5.
3) Split the supermatrix into taxa [with this script](https://github.com/wrf/supermatrix/blob/master/split_supermatrix_to_taxa.py), where a fasta file is generated for each taxon and each protein gets a unique name based on the partition.
4) Using `blastp`, align the proteins against a protein set from a reference taxon that has crystal structures of the proteins of interest (probably human) and where the proteins in the reference set are complete, meaning untrimmed.
5) With the blast results, generate an alignment of each protein with the reference and a line of the recoded sitewise likelihoods (with `blast_to_align_pairs.py`). Files will be automatically named and contain exactly three fasta entries, the trimmed protein from the supermatrix (which will have gaps), the full protein (the protein of interest that presumably has a crystal structure), and the sitewise likelihoods that should match the trimmed protein. This alignment file is the input for `pdb_log_likelihood.py`.
6) Run `pdb_log_likelihood.py` to recode the ATOM records in the PDB file of the protein of interest.
7) View in PyMOL, and run the `color_by_likelihood.py` script in the PyMOL console.

![3d4j_w_lnl.png](https://github.com/wrf/pdbcolor/blob/master/3d4j_w_lnl.png)

For example, above is the structure [3d4j](https://www.rcsb.org/structure/3d4j) of human diphosphomevalonate decarboxylase ([MVD1_HUMAN](http://www.uniprot.org/uniprot/P53602)). The protein is a homodimer (see [Voynova et al 2008 Human mevalonate diphosphate decarboxylase: Characterization, investigation of the mevalonate diphosphate binding site, and crystal structure](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2709241/)), and two pori-sis-favoring residues (teal) can be seen involved in the dimerization interface of the two C-terminal domains. Rather unusually, these two residues (S275 and W282) are F and R, respectively, in most opisthokonts and many sponges, suggesting that sponges retained the ancestral state for these residues, which is why pori-sis is the favored topology. Potentially, the R guanidinium can interact with the pi-electron cloud of F, and may still support dimerization. The constant region (green, with sulfate bound) is the active site of the enzyme.

## heteropecilly ##
The steps are almost identical to the sitewise likelihood above, only the site-wise heteropecilly calculations are used instead. Steps 3 and 4 are the same. See the [heteropecilly github repo](https://github.com/wrf/heteropecilly) for details. 

![heteropecilly_color_scheme.png](https://github.com/wrf/pdbcolor/blob/master/heteropecilly_color_scheme.png)

## References ##
The colorization script was modified from the `consurf_new.py` script from the [ConSurf Server](http://consurf.tau.ac.il/2016/)
