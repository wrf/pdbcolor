# details for sitewise scripts #
**Note: many scripts require Python libraries** `Bio` **and** `numpy`

## sitewise_ll_to_columns ##
Because not all sites provide the same phylogenetic information for all clades (i.e. sites that are constant except in vertebrates would not affect relationships of other taxa), it is necessary to get information about which sites strongly support certain hypotheses of relationships. This can be done in [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html), using the `-f G` option to calculate site-wise likelihoods.

`raxmlHPC-PTHREADS-SSE3-8.2.11 -f G -s simion2017_97sp_401632pos_1719genes.phy -m PROTGAMMALG -z tree_97sp_CAT.rooted_combined.tre -n simion2017_97sp_401632pos_1719genes -T 6`

The output of this can be converted to columns, where it can be more easily used by other programs:

`./sitewise_ll_to_columns.py RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes > RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab`

## sitewise_get_strong_sites ##
The vast majority of sites do not weigh in on any given node of a tree. This makes sense as not all sites should be expected to be equally informative to all levels of a tree. Fast-evolving sites may be relevant for narrower time scales, genus- or family level, while slower evolving sites would be relevant for phylum level distinctions. Thus, sites critical for resolving a particular node can be extracted using the `sitewise_get_strong_sites.py` script. By default, sites where the absolute value of the difference is greater than 0.5 are kept. For the [Simion dataset](https://github.com/psimion/SuppData_Metazoa_2017), this yields 8106 sites out of 344990 potential sites.

`sitewise_get_strong_sites.py -a simion2017_97sp_401632pos_1719genes.aln -l RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab -o simion2017_97sp_perSiteLLs_strong.aln`

Due to differences in distribution of strong sites, if `.sitelogl` values from `phylobayes` are being used, it is advisable to set the threshold `-d` to `1.0`.

## sitewise_columns_to_fasta ##
As it may also be useful to identify which sites favor which of the several topologies, a fasta string can be generated where each site is the calculated difference between the topologies. Because differences span a range of around -4 to +4, (though sometimes up to +8), this can be converted to a single character from 0-8 (9 total characters, 3 for each possible tree). This is typically imagined as favoring one of three possible trees (here called T1, T2, and T3). Values 0-2 favor T1, 3-5 favor T2, and 6-8 favor T3.

``sitewise_columns_to_fasta.py RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab > simion2017_97sp_perSiteLLs.fasta``

The fasta string can be combined with the original alignment, and sites can be extracted as above.

``cat simion2017_97sp_401632pos_1719genes.aln simion2017_97sp_perSiteLLs.fasta > simion2017_97sp_strong_w_lnl.aln``

## sitewise_recode_constant_sites ##
Constant sites should have no difference between lnL of any topology, yet they do (probably due to which species have gaps/missing data). These sites can specifically be recoded (as `const` or `x`) so that they are excluded from some future analyses.

`sitewise_recode_constant_sites.py -a simion2017_97sp_401632pos_1719genes.aln -l RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab > RAxML_perSiteLLs.simion2017_const_recoded.tab`

The `const_recoded.tab` file can be used for `sitewise_columns_to_fasta.py` or `blast_to_align_pairs.py`, where constant sites will be coded as lowercase `x`.

## blast_to_align_pairs ##
Because of trimming steps, most proteins in a supermatrix do not represent the entire, or even the majority, of the original protein, and the identity of this protein (say the name of a gene) may be unknown. For cases where human was used, the IDs of the human proteins can be extract with `blastp`, as even trimmed proteins will have the top hit to a real human protein with almost 100% identity. Thus, individual alignments of each trimmed protein can be remade with the reference protein.

Because human is being used as the reference set, all human proteins must be extracted from the [SwissProt set](http://www.uniprot.org/downloads). Because of the standard naming scheme of Uniprot proteins, the `getAinB.py` [script](https://bitbucket.org/wrf/sequences/src) extracts all proteins with the species tag *_HUMAN*, creating a new file of only human proteins.

`getAinB.py _HUMAN uniprot_sprot.fasta -s > human_uniprot.fasta`

Then, generate a file of all human proteins used in the supermatrix. This script `split_supermatrix_to_taxa.py` can be found in the [supermatrix repo](https://github.com/wrf/supermatrix).

`split_supermatrix_to_taxa.py -a simion2017_97sp_401632pos_1719genes.fasta.gz -p simion2017_partitions.txt -d simion2017_taxa`

Because `blastp` does not allow gaps in sequences, and most sequences extracted from the supermatrix will still have a gap, gaps need to be removed. The `degapper.py` [script](https://bitbucket.org/wrf/sequences/src) removes gaps and automatically renames the output file, though `Find/Replace` in any text editor would suffice as well.

`degapper.py simion2017_taxa/Homo_sapiens.fasta`

Make the blast protein database with `makeblastdb` and then run `blastp` and report the results as tabular (`-outfmt 6`).

`makeblastdb -dbtype prot -in human_uniprot.fasta`

`blastp -query simion2017_taxa/Homo_sapiens.fasta.nogaps -db human_uniprot.fasta -outfmt 6 -max_target_seqs 1 > simion2017_taxa/hsapiens_vs_uniprot_blastp.tab`

The top hit for each protein should probably be the original reference protein. If this is not the case, then that protein probably should not be used in phylogeny in the first place. Using the blast results, new alignments can be generated for each protein from the supermatrix and the reference protein. Each file is in fasta format and contains three sequences, the protein used in the supermatrix (which has probably been trimmed), the reference protein, and a row of the log-likelihood scores (coded between 0-8). Constant scores (if used) are coded as `x`.

The folder containing all of the files is automatically named (based on current time).

`blast_to_align_pairs.py -b simion2017_taxa/hsapiens_vs_uniprot_blastp.tab -q simion2017_taxa/Homo_sapiens.fasta.nogaps -s human_uniprot.fasta -r simion2017_taxa/Homo_sapiens.fasta -l RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab`

This requires [MAFFT](http://mafft.cbrc.jp/alignment/software/source.html), though potentially could be modified to run another aligner.

## alignment_pair_stats ##
After using the above scripts, summary information about each alignment can be determined and printed into a tabular format.

`./alignment_pair_stats.py -a blast_alignments_20171009-153805/ > simion2017_align_pair_stats.tab`

Each row contains 8 columns:

* partition - partition from the supermatrix
* protID - blast hit, in this case the Uniprot ID
* trimmedLength - number of non-gap characters in this partition
* trimmedPercent - percent relative to the length of the reference protein
* span - first non-gap character to last non-gap character, though may still contain internal gaps
* spanLength - length of the span in amino acids
* spanPercent - percent relative to the length of the reference protein
* refProtLength - length of reference protein

`Homo_sapiens_81364-81595  sp|P42858|HD_HUMAN  232  0.074  (129, 2455)  2326  0.740  3142`

