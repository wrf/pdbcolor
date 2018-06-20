# pre-processing commands are
# ~/git/supermatrix/trim_alignment_by_coverage.py -a luciferase_firefly_w_outgroups.aln -c 100
# fasta2twoline.py luciferase_firefly_w_outgroups.aln.c100trim > luciferase_firefly_w_outgroups.aln.c100trim.two

# last modified 2018-06-20

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/git/pdbcolor/sca/al_S1A_1388.mat"
#inputfile = "~/git/pdbcolor/sca/luciferase_firefly_w_outgroups.aln.c100trim.two"

### read lines as strings
alignlines = readLines(inputfile)

### split fasta header into names and sequences
fastanames = c()
sequences = c()

for (i in 1:length(alignlines)) {
	if (i%%2==1) {
		fastanames = c(fastanames, substr(alignlines[i],2,255))
	} else {
		sequences = c(sequences, alignlines[i])
	}
}

### convert each string to a list
msa = lapply(sequences, function(x) strsplit(x,""))

### convert list of lists to a matrix
msa_matrix = matrix(unlist(msa), nrow=length(alignlines)/2, byrow=TRUE)
print("dimensions are")
dim(msa_matrix)
### indicate only the 20 amino acids, ignoring gaps and X
AASET = unlist(strsplit("ACDEFGHIKLMNPQRSTVWY",""))

### get base frequencies of each amino acid, also gaps and X
basecounts = table(unlist(msa))

totalletters = sum(basecounts[match(AASET,names(basecounts))])
print(paste("total letters",totalletters))
base_freq = basecounts/totalletters
print("overall AA frequencies")
base_freq
nsites = dim(msa_matrix)[2]
nsites
ntaxa = dim(msa_matrix)[1]
ntaxa

### set up blank lists to add for each site
prevalent_aa = c()
bin_value_list = c()
d_bin_list = c()

#table_by_site = apply(msa_matrix,2,table)
numgaps = colSums(matrix(as.integer(msa_matrix=="-"),nrow=ntaxa,byrow=TRUE))

### for each site, get the most common AA, and calculate the frequency (bin_freq) and divergence (D_bin)
for (i in 1:nsites) {
	site_i = msa_matrix[,i]
	aa_freq_site_i = table(site_i)
	#nogaps_table = table_by_site[[i]][match(AASET, names(unlist(table_by_site[i])))]
	most_common = which(aa_freq_site_i==max(aa_freq_site_i,na.rm=TRUE))
	prevalent_aa = c(prevalent_aa, names(most_common)[1])
	bin_prevalence = as.integer(msa_matrix[,i]==names(most_common))
	bin_freq = sum(bin_prevalence)/ntaxa
	bin_value_list = c(bin_value_list, bin_prevalence)
	D_bin = bin_freq * log(bin_freq/base_freq[[names(most_common)[1]]]) + (1-bin_freq) * log((1-bin_freq)/(1-base_freq[[names(most_common)[1]]]))
	d_bin_list = c(d_bin_list, D_bin)
}

### generate figure of divergence by site
outputfile = gsub("([\\w/]+)\\....$","\\1.D_bin.pdf",inputfile,perl=TRUE)
pdf(file=outputfile, width=8, height=4)
plot(d_bin_list,type="h",lwd=2, xlab="Position number", ylab="Binary divergence", main=inputfile)
dev.off()

### construct binary approximation matrix
bin_matrix = matrix(bin_value_list, ncol=nsites)
### calculate the binary frequencies, and overall for background
freq_bin = colSums(bin_matrix)/ntaxa
freq_bin_bg = base_freq[match(prevalent_aa,names(basecounts))]

pca_bin = prcomp(bin_matrix)
#biplot(pca_bin)

### calculate the correlation
freq_pairs_bin = t(bin_matrix) %*% bin_matrix / ntaxa
freq_multiplied = freq_bin %*% t(freq_bin)
C_bin = freq_pairs_bin - (freq_multiplied)

### determine weights as gradient of relative entropy
W = log(freq_bin * (1-freq_bin_bg) / (freq_bin_bg*(1-freq_bin)))

### correlation statistic coupling analysis matrix
C_sca = (W %*% t(W)) * abs(C_bin)

### colorize for printing as second graph of heatmap
C_image = C_sca
C_image[C_image>0.5] = 0.5
outputfile = gsub("([\\w/]+)\\....$","\\1.site_correlations.pdf",inputfile,perl=TRUE)
pdf(file=outputfile, width=8, height=8)
par(mar=c(4,4,4,1))
image(z=abs(C_image), zlim=c(0,0.5), col=rev(rainbow(51,start=0.01,end=0.66)), main=inputfile, axes=FALSE)
axis(1,at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*nsites))
axis(2,at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*nsites))
dev.off()

### determine eigenvalues of the binary correlation matrix
eigvect = eigen(C_sca)
#hist(eigvect$values, breaks=seq(-5,max(eigvect$values)+1,0.1), ylim=c(0,20), plot=TRUE )
#hist(eigvect$values)

### randomly scramble the alignment, and recalculate eigenvalues
r_eigen_list = list()
for (i in 1:100) {
	reordered_sites = c()
	for (j in 1:nsites) {
		random_order = sample(ntaxa, ntaxa, replace=FALSE)
		reordered_sites = c(reordered_sites, bin_matrix[,j][random_order])
	}
	r_matrix = matrix(reordered_sites, nrow=ntaxa, byrow=FALSE)
	#r_freq_bin = freq_bin[random_order]
	r_freq_pairs_bin = t(r_matrix) %*% r_matrix / ntaxa
	#r_freq_multiplied = r_freq_bin %*% t(r_freq_bin)
	r_C_bin = r_freq_pairs_bin - freq_multiplied
	r_C_sca = (W %*% t(W)) * r_C_bin
	r_eigvect = eigen(r_C_sca)
	r_hist = hist(r_eigvect$values, breaks=seq(-5,30,0.1), plot=FALSE)
	r_eigen_list[[i]] = r_hist$counts
#	lines(r_hist$mids, r_hist$counts)
}

eigmin = 1

### check eigenvectors vs gaps
#plot(numgaps, eigvect$vectors[,1])
#plot(numgaps, eigvect$vectors[,2])
#plot(numgaps, eigvect$vectors[,3])
#plot(numgaps, eigvect$vectors[,4])
#plot(numgaps, eigvect$vectors[,5])

### generate figure of histogram of eigenvalues vs random alignment
outputfile = gsub("([\\w/]+)\\....$","\\1.eigen_v_random.pdf",inputfile,perl=TRUE)
pdf(file=outputfile, width=8, height=10)
par(mfrow=c(3,2))
eighist = hist(eigvect$values, breaks=seq(-5,max(eigvect$values)+1,0.1), ylim=c(0,20), xlim=c(0,30) )
top5evs = rev(which(eighist$counts > 0))[1:5]
lambdatext = c(expression(paste(lambda,1,sep="")),expression(paste(lambda,2,sep="")),expression(paste(lambda,3,sep="")),expression(paste(lambda,4,sep="")),expression(paste(lambda,5,sep="")))
text(eighist$mids[top5evs], eighist$counts[top5evs]+1, lambdatext, cex=1.2)
for (i in 1:100) {
	lines(r_hist$mids, r_eigen_list[[i]], col="#1234c366")
}
legend(15,20,legend=c("Alignment", "Random"), col=c("#000000","#1234c3"), pch=15, cex=1.5)
hist(eigvect$vectors[,1], breaks=seq(-eigmin,eigmin,0.01), ylim=c(0,40) )
r_hist = hist(r_eigvect$vectors[,1], breaks=seq(-eigmin,eigmin,0.01), plot=FALSE)
lines(r_hist$mids, r_hist$counts, col="#1234c3")
abline(v=c(-0.05,0.05), col="#cd251299", lty=2, lwd=1.2)
hist(eigvect$vectors[,2], breaks=seq(-eigmin,eigmin,0.01), ylim=c(0,40) )
r_hist = hist(r_eigvect$vectors[,2], breaks=seq(-eigmin,eigmin,0.01), plot=FALSE)
lines(r_hist$mids, r_hist$counts, col="#1234c3")
abline(v=c(-0.05,0.05), col="#cd251299", lty=2, lwd=1.2)
hist(eigvect$vectors[,3], breaks=seq(-eigmin,eigmin,0.01), ylim=c(0,40) )
r_hist = hist(r_eigvect$vectors[,3], breaks=seq(-eigmin,eigmin,0.01), plot=FALSE)
lines(r_hist$mids, r_hist$counts, col="#1234c3")
abline(v=c(-0.05,0.05), col="#cd251299", lty=2, lwd=1.2)
hist(eigvect$vectors[,4], breaks=seq(-eigmin,eigmin,0.01), ylim=c(0,40) )
r_hist = hist(r_eigvect$vectors[,4], breaks=seq(-eigmin,eigmin,0.01), plot=FALSE)
lines(r_hist$mids, r_hist$counts, col="#1234c3")
abline(v=c(-0.05,0.05), col="#cd251299", lty=2, lwd=1.2)
hist(eigvect$vectors[,5], breaks=seq(-eigmin,eigmin,0.01), ylim=c(0,40) )
r_hist = hist(r_eigvect$vectors[,5], breaks=seq(-eigmin,eigmin,0.01), plot=FALSE)
lines(r_hist$mids, r_hist$counts, col="#1234c3")
abline(v=c(-0.05,0.05), col="#cd251299", lty=2, lwd=1.2)
dev.off()

### define sectors using thresholds from Halabi 2009, and then color
threshold=.05
sec_blue = which(eigvect$vectors[,2] > lapply(abs(eigvect$vectors[,4]), function (x) max(threshold, x ) ) )
sec_red = which(eigvect$vectors[,2] < lapply(-abs(eigvect$vectors[,4]), function (x) min(-threshold, x ) ) )
sec_green = which(eigvect$vectors[,4] > lapply(abs(eigvect$vectors[,2]), function (x) max(threshold, x ) ) )
#sec_orange = which(eigvect$vectors[,3] > lapply(abs(eigvect$vectors[,2]), function (x) max(threshold, x ) ) & eigvect$vectors[,4] < 0.05 )

### make default color as white, then color special sectors
colorset = rep(c("#00000022"),nsites)
colorset[sec_red] = "#d21b00c8"
colorset[sec_blue] = "#3106d8c8"
colorset[sec_green] = "#37ae23c8"
#colorset[sec_orange] = "#ea7600c8"

### generate graph of scatterplots of vectors 1 to 5
outputfile = gsub("([\\w/]+)\\....$","\\1.eigenvectors.pdf",inputfile,perl=TRUE)
pdf(file=outputfile, width=8, height=16)
par(mfrow=c(5,2), mar=c(4,4,1,1))
plot(eigvect$vectors[,2], eigvect$vectors[,3], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
text(eigvect$vectors[,2], eigvect$vectors[,3], 1:nsites, pos=4)
plot(eigvect$vectors[,2], eigvect$vectors[,4], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
text(eigvect$vectors[,2], eigvect$vectors[,4], 1:nsites, pos=4)
plot(eigvect$vectors[,2], eigvect$vectors[,5], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
text(eigvect$vectors[,2], eigvect$vectors[,5], 1:nsites, pos=4)
plot(eigvect$vectors[,2], eigvect$vectors[,6], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
text(eigvect$vectors[,2], eigvect$vectors[,6], 1:nsites, pos=4)
plot(eigvect$vectors[,3], eigvect$vectors[,4], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
plot(eigvect$vectors[,3], eigvect$vectors[,5], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
plot(eigvect$vectors[,3], eigvect$vectors[,6], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
plot(eigvect$vectors[,4], eigvect$vectors[,5], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
plot(eigvect$vectors[,4], eigvect$vectors[,6], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
plot(eigvect$vectors[,5], eigvect$vectors[,6], xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), pch=21, bg=colorset )
abline(h=c(-0.05,0.05), v=c(-0.05,0.05), col="#00000088", lty=2, lwd=1)
dev.off()



colorset_d = rep(c("#000000"),nsites)
colorset_d[sec_red] = "#d21b00c8"
colorset_d[sec_blue] = "#3106d8c8"
colorset_d[sec_green] = "#37ae23c8"

### generate colored version of divergence graph
outputfile = gsub("([\\w/]+)\\....$","\\1.D_bin_colored.pdf",inputfile,perl=TRUE)
pdf(file=outputfile, width=8, height=4)
plot(d_bin_list,type="h",lwd=2, xlab="Position number", ylab="Binary divergence", main=inputfile, col=colorset_d)
dev.off()


alternativeplots = "plot(eigvect$vectors[,1], eigvect$vectors[,2], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,1], eigvect$vectors[,3], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,1], eigvect$vectors[,4], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,1], eigvect$vectors[,5], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,2], eigvect$vectors[,3], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,2], eigvect$vectors[,4], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,2], eigvect$vectors[,5], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,3], eigvect$vectors[,4], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,3], eigvect$vectors[,5], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )
plot(eigvect$vectors[,4], eigvect$vectors[,5], xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), pch=21, bg=colorset )"

### print the results of eigenvectors 1 to 5 as a table
outtabfile = gsub("([\\w/]+)\\....$","\\1.vec_by_site.tab",inputfile,perl=TRUE)
outtable = cbind(1:nsites, prevalent_aa, eigvect$vectors[,2],eigvect$vectors[,3],eigvect$vectors[,4],eigvect$vectors[,5])
colnames(outtable) = c("Site","AA","ev2","ev3","ev4","ev5")
write.table(outtable, outtabfile, sep="\t", row.names=FALSE)
#






#