# pdbcolor
Python code for [PyMOL](https://pymol.org/2/) to color a [PDB structure](http://www.rcsb.org/) based on various parameters obtained from a multiple sequence alignment. The scripts could be modified to accept essentially any parameter that can be obtained for each residue, say for dN/dS ratios, hydrophobicity, etc. and could be changed to represent any arbitrary value series for however many colors are needed. Currently, these scripts only support `.pdb` format, and **NOT** `.cif` format. I will probably implement this in the future.

Existing schemes include:
* [generic data](https://github.com/wrf/pdbcolor#generic-data) for each residue read from a tabular file or csv
* [percent identity](https://github.com/wrf/pdbcolor#percent-identity) calculated directly from the alignment
* [gene structure](https://github.com/wrf/pdbcolor#gene-structure), based on CDS features from a [GFF file](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

More specialized scripts include:
* [RAxML sitewise likelihoods](https://github.com/wrf/pdbcolor#raxml-site-wise-likelihood), based on relative substitution probabilities between multiple fixed phylogenetic trees
* [phylobayes sitewise likelihoods](https://github.com/wrf/pdbcolor#phylobayes-site-wise-likelihood), as above for RAxML, but using the program [phylobayes](https://github.com/bayesiancook/pbmpi)
* [heteropecilly calculations](https://github.com/wrf/pdbcolor#heteropecilly), based on calculations of lineage-specific amino acid substitutions from [Simion et al 2017](https://github.com/psimion/SuppData_Metazoa_2017)
* [statistical coupling analysis from Halabi 2009](https://github.com/wrf/pdbcolor/tree/master/sca), based on the binary approximation of co-evolving residues

## General usage on PDB files ##
### Making a PyMOL script ###
There are two practical ways of automatically coloring a structure in PyMOL. One is to process the structure and sequence and generate a script that can be run in PyMOL, which would include the color definitions and selection commands. These [commands](https://pymolwiki.org/index.php/Category:Commands) are exactly the same as those that would be [typed into the PyMOL console](https://pymol.org/pymol-command-ref.html). For example:

```
show cartoon
set_color colordefault, [0.75,0.75,0.58]
color colordefault, all
```

This is the standard option for [pdb_color_generic.py](https://github.com/wrf/pdbcolor/blob/master/pdb_color_generic.py) and [gff_cds_to_pymol_script.py](https://github.com/wrf/pdbcolor/blob/master/gff_cds_to_pymol_script.py). This type of output also can be generated for [pdb_site_identity.py](https://github.com/wrf/pdbcolor/blob/master/pdb_site_identity.py). This will recolor the residues when run in PyMOL, and might be easier to share with collaborators if the PDB files are large. The script is then run in PyMOL by typing:

`@/path/to/your/folder/name_of_your_script`

This strategy makes the most sense for cases using standard PDB files, from the databases. That way, another user can use a command like `fetch 5dv9` to get the PDF file and just run the PyMOL script or even the commands alone.

### Recoding the PDB file ###
The other strategy will color residues by changing a value for [each ATOM record in a PDB file](http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences-and-the-pdb-format), and then running a script of the [PyMOL API](https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha). This strategy may make more sense for modeled proteins, that one cannot use `fetch`. In a normal PDB file, the [temperatureFactor or beta-factor](http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/dealing-with-coordinates) is the second to last term, here in the first atom it is 0.82.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00  0.82           N  `

Each script identifies the residue number (here alanine is the 11th residue of the protein sequence, but the first in the model), and changes the temperatureFactor for each atom of that residue. For instance, in the `pdb_site_identity.py` example below, the value would be recoded to 81.88, indicating reasonably high identity (>80%) at this position, and this alanine will be colored green.

`ATOM      1  N   ALA A  11       1.483 183.030  20.022  1.00 81.88           N  `

Running the scripts within PyMOL recolors all atoms, though the color of any of these can be manually changed afterwards, perhaps to highlight ligands or specific domains.

**Note: the temperature factors are a measure of the confidence of the position of each atom. Higher temperature would mean less confident (typical for surface residues or side-chains). Thus, if a residue is highly conserved (such as those highlighted by the scripts below), it is advisable to also check the temperature factors prior to recoding, such as using the PyMOL command** `spectrum b`.

## generic data ##
Generic numerical data about each residue can be extracted from a tabular or csv file using the `pdb_color_generic.py` script. Provided that the PDB file only has a single chain (such as from structural simulations), very little information is then needed. The script will extract data from a desired column, generate an appropriate range, and print the PyMOL commands as a script to standard output. Eight color schemes are currently implemented (change with option `-l`), 4 varying from gray to color (for sequential intensity) and 4 are diverging to emphasize very high or low scores, for display on either white or black backgrounds. I normally view proteins on black background, but journals often require white background for figures.

**Note: the varied nature of input data makes it difficult to program a generic scheme, and can cause any number of unanticipated errors. Please open an issue if you encounter problems with your dataset.**

![generic_color_schemes_v1.png](https://github.com/wrf/pdbcolor/blob/master/svg/generic_color_schemes_v1.png)

For example, here I use the dN/dS data of [histamine receptor H1 3RZE](https://www.rcsb.org/structure/3rze), from the reanalysis by [Sydykova 2018](https://f1000research.com/articles/6-1845/v2). Their csv file contains many parameters, here only four are extracted to match up with the PDB file [provided in their supplemental data](https://github.com/clauswilke/proteinER). These are dNdS (non-synonymous/synonymous substitutions), lrt (dN/dS likelihood ratio test statistic), rsa (relative solvent accessibility), and wcn (weighted contact number for the side chain).

```
pdb_color_generic.py -c 4 -d , -p 3rze.pdb -i 3rze.map.rates_features.csv -l blue -g dnds --exclude-common-group > 3rze.color_by_dnds.pml 
pdb_color_generic.py -c 5 -d , -p 3rze.pdb  -i 3rze.map.rates_features.csv -l div2b -g lrt > 3rze.color_by_lrt.pml 
pdb_color_generic.py -c 12 -d , -p 3rze.pdb -i 3rze.map.rates_features.csv -l green -g rsa > 3rze.color_by_rsa.pml
pdb_color_generic.py -c 13 -d , -p 3rze.pdb -i 3rze.map.rates_features.csv -l div1b -g wcn > 3rze.color_by_wcn.pml 
```

![3rze_w_wcn_dnds.png](https://github.com/wrf/pdbcolor/blob/master/examples/3rze_w_wcn_dnds.png)

### Usage considerations ###
PyMOL appears to have trouble running commands with very long lists of residues. For cases where most of the values are 0 (or otherwise uninteresting) such as dN/dS or sitewise identity, only a few residues need to be highlighted. The lowest group can be omitted using the flag `-x` (`--exclude-common-group`).

Because scores are always processed from lowest to highest, for cases where all values are negative, the color schemes can be reversed with the flag `-r` (`--reverse-colors`). In this case, `-x` would again omit the bin of values closest to 0, hence the highest values. For diverging datasets and color schemes, the same option removes the middle bin (4, out of 0 to 8) by default. The common bin to omit can be changed with the option `-O` (`--exclusion-override`).

## percent identity ##
With the script `pdb_site_identity.py`, for a target sequence, recolor by percent identity from a multiple sequence alignment. The same sequential color schemes as generic (red, yellow, blue, and green) are used for the script output, gray indicates less than 50% identity. For the output script, red gradient is the default, but this can be changed with the option `--base-color`, such as `--base-color blue`.

#### Options are: ####
   * `-a` : alignment, required one or more multiple sequence alignments. Several files can be given if trying to recode a multimer/ complex.
   * `-s` : sequence, the sequence names or FASTA headers for the sequence/s in the alignment that match the structure. Normally this would be the protein seq ID that matches the DBREF field in the PDB file, often a Uniprot ID like `DUOX1_HUMAN`. More than one can be given, corresponding to multiple alignment files.
   * `-p` : PDB file to read or recode. Currently, `.pdb` format is required, though other options can bypass this.
   * `-w` : write script, turns on this option and will write the PyMOL script to the specified filename. If combined with `--force-recode` this can often be used to recolor a `.cif` file, such as for very large protein complexes.
   * `--base-color` : change the base color gradient, to red, yellow, green, blue, such as `--base-color green`. In all cases, the 100% identity group has a much brighter color.
   * `--force-recode` : if the DBREF field is missing or does not match the sequence from `-s` (say if sequence was renamed, using NCBI accessions, etc), then force the recoding anyway assuming the complete protein in the alignment and PDB file.

### Instructions ###
If recoding the PDB file, the default colors are following reverse rainbow order (so blue to red) to show increasing identity, with magenta showing 100% identity (excluding gaps or missing data).

![percent_identity_color_scheme.png](https://github.com/wrf/pdbcolor/blob/master/svg/percent_identity_color_scheme.png)

Colors are defined in the `color_by_identity.py` script, where triplets are RGB values of 0 to 1 (so 1,1,1 is white). Thresholds of identity are defined in the `pdb_site_identity.py` script, where the scores are calculated and printed for each of the ATOM records in the PDB file.

These values range from 0 to 100%, but are grouped into 9 bins (it is difficult to see the color differences beyond this many). The colors are meant to specifically highlight strongly conserved sites, which is why site identity below 50% is colored gray.

1) Generate a multiple sequence alignment, which must include the sequence of the target PDB structure. For instance, the protein [Symplectin](http://www.uniprot.org/uniprot/C6KYS2) is the target, and residues will be colored based on the percentage of residues in other sequences that are identical to this sequence. Here they are aligned with the program `mafft` using an example [dataset](https://bitbucket.org/wrf/squid-transcriptomes/src) from [Francis et al 2013 Symplectin evolved from multiple duplications in bioluminescent squid](https://peerj.com/articles/3633/).

    `mafft-linsi examples/symplectins_w_outgroups.fasta > examples/symplectins_w_outgroups.aln`

2) Using the alignment and the PDB file (here generated using the [SWISS-MODEL](https://www.swissmodel.expasy.org/) server), reassign the temperatureFactors based on identity percentage using `pdb_site_identity.py`. The target sequence is indicated using the `-s` option, and must exactly match the name given in the alignment. The script assumes the residues in the PDB file are numbered based on this sequence, meaning even if the PDB file starts with residue 20 (say the first 19 were disordered or cleaved), the sequence still must start with residue 1.

    `./pdb_site_identity.py -p examples/symplectin_swissmodel_01.pdb -s Symplectin -a examples/symplectins_w_outgroups.aln > examples/symplectin_swissmodel_01_w_id.pdb`

3) Open the PDB file in PyMOL:

    `pymol examples/symplectin_swissmodel_01_w_id.pdb`

4) In the PyMOL console, run the command:

    `run color_by_identity.py`

![symplectin_domains_by_conservation.png](https://github.com/wrf/pdbcolor/blob/master/examples/symplectin_domains_by_conservation.png)

In the left panel, certain residues in the binding pocket are strongly conserved (shown in red, >95%), such as the catalytic triad of E-K-C (though C is green, meaning only >80% identity, as this is sometimes serine). In the right panel, the other domain is not well conserved outside of disulfide-forming cysteines.

### Instructions for multiple proteins ###
For cases of heteromultimers, multiple alignments can be used. An alignment *may be given for each protein* in the PDB file, though if a PDB file contains multiple other proteins without alignments, or heteroatoms like metals, ligands, DNA, etc., all of these will be colored the "null" color.

The instructions above are used, but instead multiple alignments and sequence names are given for the `-a` and `-s` options in step 2. Here an [alignment](https://bitbucket.org/molpalmuc/sponge-oxygen) is used from [Mills et al 2018](https://doi.org/10.7554/eLife.31176) to plot site identity information onto the [structure of the heterodimer](http://www.rcsb.org/structure/4ZPR) of two mouse helix-loop-helix transcription factors, [HIF1A](http://www.uniprot.org/uniprot/Q61221) and [ARNT](http://www.uniprot.org/uniprot/P53762). The same alignment may be given multiple times as long as the sequences are found, otherwise several different alignments could be used, say for protein subfamilies or individual clades.

`pdb_site_identity.py -a all_hif1.aln all_hif1.aln -s HIF1A_MOUSE ARNT_MOUSE -p 4zpr.pdb > 4zpr_w_id.pdb`

In the PyMOL console, run `color_by_identity.py` script. Then, up to 4 different chains can receive unique color schemes of grayscales to red, green, blue, or yellow, with the `bychain=True` option.

`run color_by_identity.py`

`color_by_identity(bychain=True)`

Here the two proteins ARNT and HIF1a are colored red and green, respectively. Two conserved arginines (R102 and R30) are found interacting with the DNA helix, and two conserved leucines (L112 and L40) point to each other, possibly needed for the dimerization. There are more conserved residues in the HLH domain for ARNT than for HIF (8 vs 1), suggesting that despite the varied binding partners for ARNT/CYCLE/BMAL, including HIF, SIM, NPAS1/3, CLOCK, and AHR, it is likely that several of the residues that determine the binding motif are conserved across the ARNT-interacting region, and that the varied roles of these transcription factors are determined more by interactions with HIF and related proteins.

![4zpr_w_id.png](https://github.com/wrf/pdbcolor/blob/master/examples/4zpr_w_id.png)

### Making a script instead of recoding the PDB ###
Instead of directly recoding the PDB file, a PyMOL script can be generated to recolor the residues using the `-w` option, where `-w` is the file name of the output script.

`../pdb_site_identity.py -p 4zpr.pdb -s HIF1A_MOUSE ARNT_MOUSE -a hif_sim_npas.aln arnt_only.aln -w 4zpr_hif_arnt_color_by_id.pml`

This script would be run in the PyMOL console using:

`@4zpr_hif_arnt_color_by_id.pml`

### Conservation ###
The script also allows for an alternate measure of conservation, which can also be calculated using the option `--ct-conservation` based on the equations used in [Halabi, Rivoire et al 2009](http://dx.doi.org/10.1016/j.cell.2009.07.038). There, conservation is calculated as the sitewise frequency of amino acid *a* at position *i* over the global frequency of that amino acid. This is calculated as follows:

```
conservation = f(i,a) * ln( f(i,a)/q(a) ) + ( 1 - f(i,a) ) * ln ( (1-f(i,a)) / (1-q(a)) )
f(i,a) is frequency of amino acid 'a' at position 'i'
q(a) is background frequency of amino acid 'a' in all proteins in the alignment
```

Thus, when `f(i,a)` is 1, the first term reduces to `ln( 1 / q(a) )` and the second term reduces to zero, meaning the conservation score of a constant site is directly related to its overall frequency in the alignment. In effect, this means that constant sites of common amino acids are penalized, which belies their evolutionary importance over time. To me, this output is less intuitive for this reason than directly visualising the identity, so I never use this option. In any case, it would be run as follows:

`pdb_site_identity.py -p examples/symplectin_swissmodel_01.pdb -s Symplectin -a examples/symplectins_w_outgroups.aln --ct-conservation > examples/symplectin_swissmodel_01_w_cons.pdb`

## Gene structure ##
For multidomain proteins, it may be useful to view the protein structure with information from the gene structure, i.e. introns and exons. This may display whether a domain is composed of a single exon and how it may join with other parts of the protein. The script colors residues corresponding to exons, where 20 colors are currently available, and the colors repeat after 20 exons.

**Note: this assumes that the annotation in the GFF file is complete and accurate, as the first codon in the CDS is assumed to be residue 1. PyMOL selections are generated for each exon, even if those residues are not present in the final structure, meaning it is possible to have a selection that includes no residues.**

Selections of exons are automatically numbered in PyMOL (appearing as `exon_1`), corresponding to the first CDS feature. This might not match the actual first exon in the gene, due to alternative splicing or untranslated exons.

The GFF file is specified with the option `-g`. It may be most convenient to first pull all features of a particular gene out using `grep`. For the example below, `grep` would pull out all features of `PPYR_00001`.

`grep PPYR_00001 PPYR_OGS1.0.gff3 > ppyr_00001.gff`

The default chain is `A`, and can be changed with `-c`. For homomultimers, specify the chains as `-c AB`. For heteromultimers, one would need to run the script twice with two different genes (`-i`).

For example, using the data from the [recently published genome](https://github.com/photocyte/PPYR_OGS) of the firefly *Photinus pyralis*, the exons for gene `PPYR_00001`, which is the luciferase Luc1.

1) `gff_cds_to_pymol_script.py -g ppyr_00001.gff -i PPYR_00001-PA > color_ppyr_00001_5dv9.pml`

The `.pml` script is the standard output, and can be run in the PyMOL console using:

2) `@color_ppyr_00001_5dv9.pml`

![ppyr_00001_5dv9_w_label.png](https://github.com/wrf/pdbcolor/blob/master/examples/ppyr_00001_5dv9_w_label.png)

In the example of [firefly luciferase 5DV9](https://www.rcsb.org/structure/5dv9), the entire protein is divided into two domains, the large AMP-binding domain (blue, green and pink, exons 1-5), and the C-terminal domain (red and orange, exons 6-7). Dark grey residues indicate that the residue is split across two exons (also shown by the grey diamonds on the legend), which is evident from the phase in the GFF file. This also shows that while most exons of the first domain are out of phase, exon 5 ends with a multiple of three bases, then cleanly joins the next domain.

Note that the `grep` step is not strictly necessary. Because the search is specified by the GFF ID, which should be unique in a GFF, the protein can be found still even in a large GFF file. Here, using the [Alphafold-predicted structure of human nidogen1](https://alphafold.com/entry/P14543):

`grep NID1 ~/genomes/human/GCF_000001405.38_GRCh38.p12_genomic.gff > examples/human_NID1.gff`

`gff_cds_to_pymol_script.py -i cds10864 -g examples/human_NID1.gff > examples/human_NID1.exon_colors.pml`

For instance, pulling directly from the human genome annotation:

`/gff_cds_to_pymol_script.py -i cds10864 -g ~/genomes/human/GCF_000001405.38_GRCh38.p12_genomic.gff > examples/human_NID1.exon_colors.pml`

![human_nidogen1_colored_by_exons.png](https://github.com/wrf/pdbcolor/blob/master/examples/human_nidogen1_colored_by_exons.png)

## RAxML site-wise likelihood ##
For an alignment and a series of phylogenetic trees with fixed topologies, [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) can produce a table of site-wise log-likelihoods for each tree topology (using the `-f G` option). The difference between the top two topologies (out of 3) can be computed for each site (the "dlnL"), showing which sites contribute strongly to one topology or the other. These dlnL values are recoded to a string of numbers (0-8) where the dlnL refers to the max minus the median and favors whichever topology had the maximum likelihood. Likewise, constants are coded as `x`. This also permits each value to be coded as a single character (for display as a fasta sequence).

While most sites are marginally different - values between 0 and 0.5 - colored rust (0), mud (3) and lead (6), a small number of sites have a strong effect on the likelihood of the final tree, colored in pink/red (1,2) for sites favoring tree 1, teal/green (4,5) for tree 2, and cobalt/blue (7,8) for tree 3. The scripts to make use of this assume that proteins and site-wise likelihood values derive from a protein supermatrix. Because proteins in supermatrices are often trimmed compared to the complete sequence, any missing values are colored gray (-1) as gaps. Sites that were constant in the alignment need to be recoded as such, and are then colored orange in the resulting protein.

![likelihood_color_scheme_w_tree_v2.png](https://github.com/wrf/pdbcolor/blob/master/svg/likelihood_color_scheme_w_tree_v2.png)

The workflow is meant to begin from a supermatrix, meaning a lot of other data need to be generated and organized and assumes that most proteins in the supermatrix are trimmed or are incomplete. The entire workflow can be done on a normal laptop in a couple hours. The steps required are:

1) Starting from a protein supermatrix ( see [supermatrix repo for scripts to manipulate matrices](https://github.com/wrf/supermatrix) ), run `RAxML` to get sitewise likelihoods of the entire dataset, then convert that output to tabular using the `sitewise_ll_to_columns.py` script. For example, using the [Simion 2017 dataset](https://github.com/psimion/SuppData_Metazoa_2017):

    `raxmlHPC-PTHREADS-SSE3-8.2.11 -f G -s simion2017_97sp_401632pos_1719genes.phy -m PROTGAMMALG -z tree_97sp_CAT.rooted_combined.tre -n simion2017_fullset -T 6`

    `./sitewise_ll_to_columns.py RAxML_perSiteLLs.simion2017_fullset > RAxML_perSiteLLs.simion2017_fullset.tab`

2) Because constant sites may still receive lnL values from RAxML (which is biologically meaningless) these should be recoded to a value indicating they are constant using the `sitewise_recode_constant_sites.py` [script](https://github.com/wrf/pdbcolor/blob/master/sitewise_scripts/sitewise_recode_constant_sites.py). Using the matrix and the tabular likelihoods, recode constant sites and generate a new tabular output file that can be used in step 5.

    `sitewise_recode_constant_sites.py -a simion2017_97sp_401632pos_1719genes.fasta -l RAxML_perSiteLLs.simion2017_fullset.tab > RAxML_perSiteLLs.simion2017_fullset_const_recoded.tab`

3) Split the supermatrix into taxa [with this script](https://github.com/wrf/supermatrix/blob/master/split_supermatrix_to_taxa.py), where a fasta file is generated for each taxon and each protein gets a unique name based on the partition. Any gaps will need to be removed from the split alignments.

    `split_supermatrix_to_taxa.py -a simion2017_97sp_401632pos_1719genes.fasta.gz -p simion2017_partitions.txt -d simion2017_taxa`

4) Using `blastp`, align the proteins against a protein set from a reference taxon that has crystal structures of the proteins of interest (probably human) and where the proteins in the reference set are complete, meaning untrimmed.

    `blastp -query simion2017_taxa/Homo_sapiens.fasta.nogaps -db human_uniprot.fasta -outfmt 6 -max_target_seqs 1 > simion2017_taxa/hsapiens_vs_uniprot_blastp.tab`

5) With the blast results, generate an alignment of each protein with the reference and a line of the recoded sitewise likelihoods (with `blast_to_align_pairs.py`). Files will be automatically named and contain exactly three fasta entries, the trimmed protein from the supermatrix (which will have gaps), the full protein (the protein of interest that presumably has a crystal structure), and the sitewise likelihoods that should match the trimmed protein. This alignment file is the input for `pdb_log_likelihood.py`.

    `blast_to_align_pairs.py -b simion2017_taxa/hsapiens_vs_uniprot_blastp.tab -q simion2017_taxa/Homo_sapiens.fasta.nogaps -s human_uniprot.fasta -r simion2017_taxa/Homo_sapiens.fasta -l RAxML_perSiteLLs.simion2017_fullset_const_recoded.tab`

6) Run `pdb_log_likelihood.py` to recode the ATOM records in the PDB file of the protein of interest.
7) View in PyMOL, and run the `color_by_likelihood.py` script in the PyMOL console.

For example, in the analysis by [Shen et al 2017](https://www.nature.com/articles/s41559-017-0126) using the [Whelan2015-D1-Opisthokont](https://figshare.com/articles/Error_signal_and_the_placement_of_Ctenophora_sister_to_all_other_animals/1334306) set, [STT3B](http://www.uniprot.org/uniprot/Q8TCJ2) was identified as an outlier. STT3B is the catalytic subunit of the oligosaccharyltransferase complex. It contains 22 strong sites, where 20 favor one topology (ctenophora-sister), 1 favors another (porifera-sister), and 1 favors the third. A human cryo-EM structure is now available at [6S7T](https://www.rcsb.org/structure/6S7T) (see [Ramirez 2019](http://dx.doi.org/10.1126/science.aaz3505)), but was previously [modeled](https://swissmodel.expasy.org/repository/uniprot/Q8TCJ2) after the structure of the complex in yeast: [6EZN](https://www.rcsb.org/structure/6ezn) (see [Wild 2018](http://science.sciencemag.org/content/359/6375/545)). 

`~/git/pdbcolor/pdb_log_likelihood.py -a examples/26999-27593-STT3B_HUMAN.aln -p STT3B_HUMAN_mod_from_6ezn.pdb -s STT3B_HUMAN > STT3B_HUMAN_mod_from_6ezn_w_lnl.pdb`

The script can also be run on the cyro-EM structure, including other proteins:

`~/git/pdbcolor/pdb_log_likelihood.py -a examples/26999-27593-STT3B_HUMAN.aln examples/367375-367477-DAD1_HUMAN.aln examples/160305-160495-RPN2_HUMAN.aln examples/41614-41862-OST48_HUMAN.aln -p examples/6s7t.pdb -s STT3B_HUMAN DAD1_HUMAN RPN2_HUMAN OST48_HUMAN -w > examples/6s7t.w_lnl.pdb`

![STT3b_HUMAN_mod_from_6ezn_w_lnl_labels.png](https://github.com/wrf/pdbcolor/blob/master/examples/STT3b_HUMAN_mod_from_6ezn_w_lnl_labels.png)

The complex overall is membrane bound in the ER, where a large bundle of helices sits within the membrane. Above this (in the lumen of the ER) lies the active site, containing a block of seven constant residues (W604-Q610) and no strong sites. The only porifera-sister-favoring residue is L152, which is L or another aliphatic residue in most taxa and distant outgroups, and G/S/T in ctenophores and choanoflagellates. In yeast, this residue is in proximity to an aromatic cluster on WBP1/OSTD ([Ost48 in human](http://www.uniprot.org/uniprot/P39656)), potentially involved in hydrophobic packing. Two strong residues are found on a helix adjacent to OST2 ([Dad1 in human](http://www.uniprot.org/uniprot/P61803)) in yeast, S240 and C251. C251 points towards a relatively large cavity at the "back" of the complex. S240 is also S in most species, but A in ctenophores and most outgroups. Seven strong residues are found in a poorly-modeled region that sits between OST2 and OST3 ([Tusc3](http://www.uniprot.org/uniprot/Q13454)), where it was suggested that this region could interact directly with the translocon. One helix contains three alanines on the same face of one helix, where they are leucine, glycine, and an aliphatic residue in ctenophores and most outgroups, and alanines in most other species. Four strong residues are found within the lumenal domain of STT3b, including F702 (Y in ctenophores and outgroups), G708 (A in ctenophores and outgroups), L712 (M in most ctenophores and outgroups) and A744 (V in all sponges, S/T in ctenophores).

For some cases where the `DBREF` field and the sequence name do not match, perhaps as the gene was renamed or an old structure was used, use the `--force-recode` option to recode the atoms regardless of the gene name or chains in the `DBREF` fields. For the example alignment and protein structure `3d4j.pdb` ([downloaded here](https://www.rcsb.org/structure/3d4j)), the `DBREF` uses an old name from the yeast homolog `ERG19`, referring to *ergosterol biosynthesis 19*.

`pdb_log_likelihood.py -a examples/59546-59840-MVD1_HUMAN.aln -s MVD1_HUMAN -p 3d4j.pdb --force-recode > 3d4j_w_lnl.pdb`

![3d4j_w_lnl_w_labels.png](https://github.com/wrf/pdbcolor/blob/master/examples/3d4j_w_lnl_w_labels.png)

This protein is called diphosphomevalonate decarboxylase ([MVD1_HUMAN](http://www.uniprot.org/uniprot/P53602)) in humans. The protein is a homodimer (see [Voynova et al 2008 Human mevalonate diphosphate decarboxylase: Characterization, investigation of the mevalonate diphosphate binding site, and crystal structure](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2709241/)), and two pori-sis-favoring residues (teal) can be seen involved in the dimerization interface of the two C-terminal domains. Rather unusually, these two residues (S275 and W282) are F and R, respectively, in most opisthokonts and many sponges, suggesting that sponges retained the ancestral state for these residues, which is why pori-sis is the favored topology. Potentially, the R guanidinium can interact with the pi-electron cloud of F, and may still support dimerization. The constant region (orange, with sulfate bound) is the active site of the enzyme.

For PDB files that contain multiple proteins, additional alignments can optionally be listed with `-a`. The corresponding protein ID must be given in the same order with `-s`.

### detailed instructions for scripts ###
For details regarding the use of the above scripts involved in [sitewise likelihood calculations](https://github.com/wrf/pdbcolor#raxml-site-wise-likelihood), see the source code and instructions in the [sitewise_scripts folder](https://github.com/wrf/pdbcolor/tree/master/sitewise_scripts), and were used in [Francis 2020](https://doi.org/10.7717/peerj.8865).

## phylobayes site-wise likelihood ##
Average site-wise likelihoods can be calculated from [phylobayes](https://github.com/bayesiancook/pbmpi). The procedure of plotting these onto a structure is similar to the above instructions for [RAxML](https://github.com/wrf/pdbcolor#raxml-site-wise-likelihood), with a few differences in program operation and analysis.

Again, the workflow is meant to begin from a supermatrix, though the analyses in `phylobayes` take substantially longer than for `RAxML` (100x-200x longer), potentially several days for each tree topology. The steps required are:

1) Starting from a protein supermatrix in *relaxed phylip* format, run `pbmpi` to generate a Markov chain for each fixed tree topology. For example, on a hypothetical dataset called `dataset.phy` with two or three different fixed trees (same ones used for RAxML):

    `mpirun -np 6 ~/pbmpi/data/pb_mpi -d dataset.phy -T dataset_t1.tree -cat -gtr d1_t1_chain1`

2) Run the post-analysis program `readpb-mpi` to calculate average sitewise likelihoods on some part of the chain. This operation may take substantially longer than running the chain in the first place.

    `mpirun -np 6 ~/pbmpi/data/readpb_mpi -sitelogl -x 50 5 d1_t1_chain1`

3) Use the script `sitewise_join_phylobayes_sitelogl.py` to join each of the `.sitelogl` files from each tree into one.

    `sitewise_join_phylobayes_sitelogl.py d1_t1_chain1.sitelogl d1_t2_chain1.sitelogl d1_t3_chain1.sitelogl > d1_combined_sitelogl.tab`

4) Then proceed [as above from step 2 to step 7](https://github.com/wrf/pdbcolor#raxml-site-wise-likelihood), noting one additional difference in downstream processing. In the original analysis by [Shen et al 2017](https://www.nature.com/articles/s41559-017-0126), sites with a disproportionate influence (called "strong sites") were defined as those with a dlnL of `0.5` or greater. Based on the site-homogeneous calculations in `RAxML`, this comes out to around 1-2% of sites, roughly 3 standard deviations apart from the mean (exact is `0.456`). For `phylobayes` under CAT-GTR, substantially more transitions are considered "probable", and almost 10% of sites would have dlnL values above 0.5. Thus, to produce similar numbers of strong sites in `phylobayes`, 3 standard deviations approximates to around a dlnL of `1.0` (exact is `0.959`), which must instead be used as the the "strong site" threshold. Therefore, in `blast_to_align_pairs.py`, the "strong site threshold" option `-t` must be set to `1.0`.

    `blast_to_align_pairs.py -b dataset_taxa/hsapiens_vs_uniprot_blastp.tab -q dataset_taxa/Homo_sapiens.fasta.nogaps -s human_uniprot.fasta -r dataset_taxa/Homo_sapiens.fasta -l d1_combined_sitelogl_const_recoded.tab -t 1.0`

All other steps are done as above.

## heteropecilly ##
The steps are almost identical to the [sitewise likelihood above](https://github.com/wrf/pdbcolor#raxml-site-wise-likelihood), only the site-wise heteropecilly calculations are used instead. Steps 3 and 4 are the same, but the tabular heteropecilly information must instead be given in step 5 for `blast_to_align_pairs.py`. See the [heteropecilly github repo](https://github.com/wrf/heteropecilly) for further instructions and analysis.

![heteropecilly_color_scheme.png](https://github.com/wrf/pdbcolor/blob/master/svg/heteropecilly_color_scheme.png)

`blast_to_align_pairs.py -b simion2017_taxa/hsapiens_vs_uniprot_blastp.tab -q simion2017_taxa/Homo_sapiens.fasta.nogaps -s human_uniprot.fasta -r simion2017_taxa/Homo_sapiens.fasta -p hp_by_site_w_const.tab`

`pdb_heteropecilly.py -a examples/10543-11140-MSH2_HUMAN.aln -p examples/2o8b.pdb -s MSH2_HUMAN > 2o8b_w_hp.pdb`

Heteropecilly color scheme can be visualized within Pymol, using the command in the console `run ~/git/pdbcolor/color_by_heteropecilly.py`

Below is an example from [2o8b.pdb](https://www.rcsb.org/structure/2o8b), which is the structure of the mismatch repair protein [MSH2](http://www.uniprot.org/uniprot/P43246)/[MSH6](http://www.uniprot.org/uniprot/P52701) heterodimer (see [Warren et al 2007 Structure of the human MutSalpha DNA lesion recognition complex.](https://www.ncbi.nlm.nih.gov/pubmed/?term=17531815)). Heteropecilly scores were only calculated for MSH2, so the other protein is colored in pale gray. Gaps or missing data are dark gray, constant sites are green, and the colors follow the deciles as in the charts above. In this example, the DNA helix is colored yellow, to distinguish it from the protein. Several features are evident. Large sections of the alignment had been removed by trimming, resulting in gaps when compared to the reference protein, and dark gray regions throughout the protein (long helix connecting the "clamp" domain to the ATPase domain). Constant sites form a distinct sector at the top of the image, primarily comprising the ATPase domain, probably also involved in the interface with MSH6. Many heteropecillious sites (red) appear to occur on the surface of the protein, perhaps directly interacting with the solvent or other proteins, though the extent of this was not precisely calculated. This may mean that, in general, heteropecillious sites and lineage-specific changes are a reflection of unique interactions *between* proteins.

![2o8b_w_hp_w_labels.png](https://github.com/wrf/pdbcolor/blob/master/examples/2o8b_w_hp_w_labels.png)

# References #
The first colorization script was modified from the `consurf_new.py` script from the [ConSurf Server](http://consurf.tau.ac.il/2016/), by [Ashkenazy et al 2016](https://academic.oup.com/nar/article/44/W1/W344/2499373)

### dN-dS ###
* Sydykova, DK, et al (2018) [Measuring evolutionary rates of proteins in a structural context](https://doi.org/10.12688/f1000research.12874.2) F1000Research 2018, 6:1845.

### Conservation ###
* Halabi, N., Rivoire, O. et al (2009) [Protein Sectors: Evolutionary Units of Three-Dimensional Structure](http://dx.doi.org/10.1016/j.cell.2009.07.038). *Cell* 138 (4) 774-786.
* Francis, WR. et al (2017) [Symplectin evolved from multiple duplications in bioluminescent squid](https://doi.org/10.7717/peerj.3633). *PeerJ* 5 (1) e3633.
* Mills, DB., Francis, WR. et al (2018) [The last common ancestor of animals lacked the HIF pathway and respired in low-oxygen environments](https://doi.org/10.7554/eLife.31176). *eLife* 7: e31176.

### Gene structure ###
* Fallon, TR. et al (2018) [Firefly genomes illuminate parallel origins of bioluminescence in beetles](https://doi.org/10.1101/237586) *eLife* 7:e36495.

### Phylogenetics ###
* Voynova, NE. et al (2008) [Human mevalonate diphosphate decarboxylase: Characterization, investigation of the mevalonate diphosphate binding site, and crystal structure](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2709241/). *Archives of Biochemistry and Biophysics* 480 (1) 58-67.
* Shen, X. et al (2017) [Contentious relationships in phylogenomic studies can be driven by a handful of genes](https://www.nature.com/articles/s41559-017-0126). *Nature Ecology & Evolution* 1 1-10.
* Wild, R. et al (2018) [Structure of the yeast oligosaccharyltransferase complex gives insight into eukaryotic N-glycosylation](http://science.sciencemag.org/content/359/6375/545). *Science* 550: 1-12.
* Ramirez, AS. et al (2019) [Cryo-electron microscopy structures of human oligosaccharyltransferase complexes OST-A and OST-B](http://dx.doi.org/10.1126/science.aaz3505)) *Science* 366: 1372-1375.
* Francis, WR., and DE. Canfield (2020) [Very few sites can reshape the inferred phylogenetic tree](https://doi.org/10.7717/peerj.8865). *PeerJ* 8:e8865.

### Heteropecilly ###
* Warren, JJ. et al (2007) [Structure of the human MutSalpha DNA lesion recognition complex](https://www.ncbi.nlm.nih.gov/pubmed/?term=17531815). *Molecular Cell* 26 (4) 579-592.
* Roure, B. and H. Philippe (2011) [Site-specific time heterogeneity of the substitution process and its impact on phylogenetic inference](http://www.ncbi.nlm.nih.gov/pubmed/21235782). *BMC Evolutionary Biology* 11:17.
* Simion, P. et al (2017) [A Large and Consistent Phylogenomic Dataset Supports Sponges as the Sister Group to All Other Animals](http://www.sciencedirect.com/science/article/pii/S0960982217301999). *Current Biology* 27 (7) 958-967.

