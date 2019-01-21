# details for SCA from Halabi et al 2009 #
**Note: scripts require Python library** `Bio`

Scripts to recreate most analyses and figures from [Halabi, Rivoire et al 2009](http://dx.doi.org/10.1016/j.cell.2009.07.038). R code was "translated" from the [Matlab scripts in the supplement](https://www.sciencedirect.com/science/article/pii/S0092867409009635?via%3Dihub#app2). The R code follows the same strategy:

* reducing each column to a binary approximation of the most common residue in each column (as was argued by the authors that this was sufficient) and calculating the divergence, generating the figure `.D_bin.pdf`.
* comparing which sequences have the most common amino acid or not, generating the figure `.site_correlations.pdf`.
* randomizes the values in each column and runs again, to get an approximation of the threshold for random associations (run 100 times), generating the figure `.eigen_v_random.pdf`.
* plotting all sites for pairs of the top 5 eigenvectors, generating the figure `.eigenvectors.pdf`.
* coloring the same divergence plot with eigenvectors 2 and 4, as done in the original work, generating the figure `.D_bin_colored.pdf`.
* printing a table of eigenvectors, as `.vec_by_site.tab`

Below is a test example using [firefly luciferase 5DV9](https://www.rcsb.org/structure/5dv9). The "blue" sector is apparent, containing active site, but the corrections of most of the sites are weak, as is the case for "green" and "red". This could be due to several factors, firstly being the size of the protein, `>500`AAs vs `<200` for those in [Halabi, Rivoire et al 2009](http://dx.doi.org/10.1016/j.cell.2009.07.038). One [study of *de novo* structure prediction from covarying residues in multiple sequence alignments](https://doi.org/10.1073/pnas.1314045110) had given a guideline that the number of sequences should be five-fold more than the length, meaning 2500 sequences in this case. As only 800 were used, this may not be sufficient to determine covariation.

![5dv9_w_sectors.png](https://github.com/wrf/pdbcolor/blob/master/sca/5dv9_w_sectors.png)

A similar approach is used by [Hopf et al 2012](https://doi.org/10.1016/j.cell.2012.04.012), whereby the find similar sectors in the heatmap. However, they interpret these sectors differently, suggesting that the correlations may more generally involve other structural configurations, like dimerization interfaces or alternate conformations (e.g. open vs. closed).

Other criticisms of the sector theory include the fact that most proteins appear to only have a single sector and the fact that choosing the correct sector from the eigenvectors was a manual trial-and-error task (see [Tesileanu et al 2015](https://doi.org/10.1371/journal.pcbi.1004091), [code here](https://bitbucket.org/ttesileanu/multicov/src/default/)). Here, in the firefly luciferase example, looking at the [eigenvectors](https://github.com/wrf/pdbcolor/blob/master/sca/luciferase_firefly_w_outgroups.aln.c100trim.eigenvectors.pdf), it is not obvious which colors or vectors are the sector(s). Eigenvector 2 seems to contain most of the strongly covarying blue sites, but the others are not as clear.

The other problem of having a single sector appears to challenge the entire theory. The original proposition was that sectors could mediate different structure-functions, i.e. folding kinetics as versus enzymatic activity. If there is only a single sector comprising both, then such functions are not independent. A study looking at conservation of sites as a function of distance in enzymes (Jack et al 2016)[https://doi.org/10.1371/journal.pbio.1002452] showed that conservation decreases with distance from the active site, suggesting that this pattern held true for enzymes where the active site was buried, but not those where it was on the surface (around a tenth of their dataset). This argues that distinct sectors are not a major feature of most proteins.

## Instructions and example ##
Overall, the script accepts a two-line FASTA format alignment, and requires a PDB file of one of the sequences.

1) If starting from a FASTA file, align with any normal aligner (here using [MAFFT](http://mafft.cbrc.jp/alignment/software/source.html)).

    `mafft luciferase_firefly_w_outgroups.fasta > luciferase_firefly_w_outgroups.aln`

2) Indels, meaning alignment columns that are primarily gaps, will make false covariance with other indels in the same proteins (or protein family), so sites that are primarily gaps should be removed. The [`trim_alignment_by_coverage.py`](https://github.com/wrf/supermatrix/blob/master/trim_alignment_by_coverage.py) script will automatically rename the file based on the option `-c`, which is the threshold for minimum coverage. In this case, 100 is used, but this could be set higher depending on the number of sequences in the alignment. The output file is generated automatically with an appended `.c100trim`.

    `~/git/supermatrix/trim_alignment_by_coverage.py -a luciferase_firefly_w_outgroups.aln -c 100`

3) The R script can accept the alignment in a 2-line FASTA format, so the alignment needs to be converted to two-line.

    `fasta2twoline.py luciferase_firefly_w_outgroups.aln.c100trim > luciferase_firefly_w_outgroups.aln.c100trim.two`

4) Run the R script `sca_on_alignment.R`, with the `.two` file as input. This will make various output files, including a tabular text file `.vec_by_site.tab`.

    `Rscript sca_on_alignment.R luciferase_firefly_w_outgroups.aln.c100trim.two`

5) Using the trimmed alignment, the PDB file, the tabular file of eigenvalues, and the name of the target sequence, run the Python script `sca_to_pymol_script.py`. This generates a PyMOL script named directly after the PDB file (with `.sectors.py` appended) that can be run in the console, similar to the ones provided as examples from the Halabi 2009 supplement.

    `sca_to_pymol_script.py -s 5DV9_A -a sca/luciferase_firefly_w_outgroups.aln.c100trim -p sca/5dv9.pdb -e sca/luciferase_firefly_w_outgroups.aln.c100trim.vec_by_site.tab`

6) Load the PDB file and run the script in the console. Do **not** use the `run` command.

    `@5dv9.pdb.sectors.py`
