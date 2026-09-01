## Cells with 2 alphas or 2 betas
# Tabulate per sample frequency
# Adapted from 1123, uses datnt2 (no transgenic, no thymus dat), inherited from 1.Dataprep

######### Objects needed #####
# datnt2 	from 1 Dataprep
# dat	 	from 1 Dataprep
# datc   	from 1 Dataprep (ALXcells section)
# Va, Vb 	from 1 Dataprep
# annot  	annotation table
# lev1   	cell-type vector (set to cell.types from 1.Dataprep)

###############################################################################
# NOTES on conversion from S-PLUS (.SSC) to R:
#  - "_" assignment -> "<-"
#  - trailing-comma substring(x, first,) -> substring(x, first)
#  - par(fin=...) removed (not settable in R); plots written to PDF
#  - text(..., crt=90) -> text(..., srt=90)  [S-PLUS rotation arg]
#  - image() given explicit 1:n coordinates so axis text positions match
###############################################################################

######### Load dataprep outputs and annotation #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

load(file.path(script_dir, "1. Dataprep and basics 0226.RData"))
# dat, datc, datntr, datbsl, datnt2, igts, Va, Vb, Ja, Jb, abP, abN, abP.ntr, abN.ntr,
# annot, cell.types

lev1 <- cell.types

sample_cells <- function(ids, n) {
  if (length(ids) == 0) stop("Cannot sample cells: subset is empty")
  sample(ids, n, replace = length(ids) < n)
}

######### Save all plots to PDF #####
pdf(file.path(script_dir, "5. Cells with 2 betas or 2 alphas.pdf"), width = 10, height = 10)




twochain.summary <- data.frame(rep("--",2),rep("--",2),  rep("--",2), rep(0, 2), rep(0, 2), rep(0, 2), rep(0, 2), rep(0, 2), stringsAsFactors=F)
colnames(twochain.summary) <- c("IGT", "orgs", "smps", "TotalCells", "ab", "aab", "abb", "aabb")

all.bb <- datnt2[1,]		# summary table format, all double beta cells (one or two alphas) - seeded with first row of datnt2, remove later
w <- unique(datnt2$IGT); w <- w[order(as.numeric(substring(w,4)))]		#reorder IGTs up front
#for (k in "IGT10") {
for (k in w) {		#loop through IGTs
	x <- datnt2[datnt2$IGT==k,]
	smps <- unique(x$sample.name)
	orgs <- rep("-", length(smps)); for (i in 1:length(orgs)) orgs[i] <- unique(x$organ[x$sample.name==smps[i]])
	res <- matrix(0, length(smps), 5); colnames(res) <- c("TotalCells", "ab", "aab", "abb", "aabb")
for (j in (1:length(smps))) {
y <- x[x$sample.name==smps[j],]
z <- y$alpha.functionality=="productive" & y$beta.functionality=="productive"
res[j,1] <- nrow(y)
res[j,2] <- sum(z & y$alpha2.vgene=="" & y$beta2.vgene=="")	#only productive
res[j,3] <- sum(z & y$alpha2.functionality=="productive" & y$beta2.vgene=="")
res[j,4] <- sum(z & y$alpha2.vgene=="" & y$beta2.functionality=="productive")
res[j,5] <- sum(z & y$alpha2.functionality=="productive" & y$beta2.functionality=="productive")
all.bb <- rbind(all.bb, y[y$beta2.functionality=="productive",])
}	#end j
z <- data.frame(rep(k, length(smps)), orgs, smps, res,  stringsAsFactors=F); colnames(z)[1] <- "IGT"
twochain.summary <- rbind(twochain.summary, z)
}		#end k
twochain.summary <- twochain.summary[-(1:2),]
z <- rowSums(twochain.summary[,5:8])>0
twochain.summary <- twochain.summary[z,]
pct.twoalpha <- round(100*twochain.summary$aab /rowSums(twochain.summary[,5:8]), digits=1)
pct.twobeta <- round(100*twochain.summary$abb/rowSums(twochain.summary[,5:8]), digits=1)
twochain.summary <- data.frame(twochain.summary, pct.twoalpha, pct.twobeta, stringsAsFactors=F)

all.bb <- all.bb[-1,]

#plots
z <- twochain.summary$TotalCells>100
plot(twochain.summary[z,9], twochain.summary[z,10], pch=".", xlab="Pct of cells with 2 alpha chains", ylab="Pct of cells with 2 beta chains", xlim=c(0,16), ylim=c(0,16))
lines(c(0,16),c(0,16))
plot(twochain.summary[z,6], twochain.summary[z,7], pch=".", xlab="Nb of cells with 2 alpha chains", ylab="Nb of cells with 2 beta chains", xlim=c(0,200), ylim=c(0,200))
lines(c(0,200),c(0,200))
plot(twochain.summary[z,9], twochain.summary[z,10], pch=".", xlab="Pct of cells with 2 alpha chains", ylab="Pct of cells with 2 beta chains", xlim=c(0,16), ylim=c(0,16))
text(twochain.summary[z,9]+0.2, twochain.summary[z,10], paste(twochain.summary$IGT[z], twochain.summary$smps[z], sep="."), cex=0.5, adj=0)
lines(c(0,16),c(0,16))
plot(twochain.summary[z,9], twochain.summary[z,10], pch=".", xlab="Pct of cells with 2 alpha chains", ylab="Pct of cells with 2 beta chains", xlim=c(0,16), ylim=c(0,16))
text(twochain.summary[z,9]+0.2, twochain.summary[z,10], twochain.summary$orgs[z], cex=0.5, adj=0)
lines(c(0,16),c(0,16))

plot(twochain.summary[z,5], twochain.summary[z,7], pch=".", xlab="Nb of total cells with ab pair", ylab="Nb of cells with 2 beta chains", xlim=c(0,4000), ylim=c(0,600))
text(twochain.summary[z,5]+50, twochain.summary[z,7], twochain.summary$orgs[z], cex=0.7, adj=0)


## Compare read depths
z <- sample(1:nrow(datnt2), 20000)
plot(datnt2$alpha.reads[z], datnt2$alpha2.reads[z], pch=".", log="xy", xlim=c(1, 50000), ylim=c(1,50000))
lines(c(1,50000),c(1,50000))
plot(datnt2$beta.reads[z], datnt2$beta2.reads[z], pch=".", log="xy", xlim=c(1, 50000), ylim=c(1,50000))
lines(c(1,50000),c(1,50000))


## Export lists of cells for doublet detection
ab <- datnt2$alpha.functionality=="productive" & datnt2$beta.functionality=="productive"
bb <- datnt2$beta.functionality=="productive" & datnt2$beta2.functionality=="productive"
aa <- datnt2$alpha.functionality=="productive" & datnt2$alpha2.functionality=="productive"
aabb <- aa & bb
chain.cell.type <- data.frame(datnt2$IGT, datnt2$IGT.cellID, ab, aa, bb, aabb, stringsAsFactors=F)
chain.cell.type <- chain.cell.type[chain.cell.type$ab==T,]

## Frequency of second chain relative to first
vb.mat <- matrix(0, length(Vb), length(Vb))
colnames(vb.mat) <- rownames(vb.mat) <- Vb
n <- rankofsame.in2ndbeta <- rep(0, length(Vb))
w <- datnt2[bb==T,]
for ( i in 1:length(Vb)) {
	x <- w$beta.vgene==Vb[i]
	n[i] <- sum(x)
	for ( j in 1:length(Vb)) {
		y <- w$beta2.vgene[x]==Vb[j]
		vb.mat[i,j] <- round(100*sum(y)/sum(x), digits=1)
	}
}
	pct <- round(100*n/sum(n),digits=1)
for (i in 1:length(Vb)) rankofsame.in2ndbeta[i] <- (1:length(Vb))[order(vb.mat[,i])==i]
vb.mat <- data.frame(n, pct, vb.mat, stringsAsFactors=F)
vb.mat.freq <- vb.mat[,-(1:2)]/rowMeans(vb.mat[,-(1:2)])


par(mfrow=c(1,1))
plot(vb.mat$pct, colMeans(vb.mat[,-(1:2)]), pch=".", xlim=c(0,10), ylim=c(0,10), xlab="Vb frequency in primary beta.vgene", ylab="Vb frequency in beta2.vgene")
text(vb.mat$pct+0.2, colMeans(vb.mat[,-(1:2)]), Vb, cex=0.8, adj=0)
lines(c(0,10), c(0,10))
x <- as.matrix(vb.mat[,-(1:2)])
image(seq_len(nrow(x)), seq_len(ncol(x)), x, axes=FALSE)
text(0, 1:length(Vb), Vb, adj=1, cex=0.7)
text(1:length(Vb), 0, Vb, adj=1, cex=0.7, srt=90)
z <- vb.mat$pct>2
test <- t(round(log2(t(x[z,z])/colMeans(x)[z]), digits=2))
image(seq_len(nrow(x)), seq_len(ncol(x)), x/rowMeans(x), axes=FALSE)
text(0, 1:length(Vb), Vb, adj=1, cex=0.7)
text(1:length(Vb), 0, Vb, adj=1, cex=0.7, srt=90)

plot(1:length(Vb), rankofsame.in2ndbeta, pch=16, ylab="Same chain in beta2 TCR, rank", xaxt="n", xlab="")
z <- rankofsame.in2ndbeta==length(Vb)
points((1:length(Vb))[z], rankofsame.in2ndbeta[z], pch=16, col=3)
text(1:length(Vb), 0, Vb, adj=1, cex=0.7, srt=90)
title(main="Rank of the same TRBV being also used in second productive TCRb rearrangement", cex.main=1)



## same, alpha
va.mat <- matrix(0, length(Va), length(Va))
colnames(va.mat) <- rownames(va.mat) <- Va
n <- rep(0, length(Va))

w <- datnt2[aa==T,]
for ( i in 1:length(Va)) {
	x <- w$alpha.vgene==Va[i]
	n[i] <- sum(x)
	for ( j in 1:length(Va)) {
		y <- w$alpha2.vgene[x]==Va[j]
		va.mat[i,j] <- round(100*sum(y)/sum(x), digits=1)
	}
}
	pct <- round(100*n/sum(n),digits=1)
va.mat <- data.frame(n, pct, va.mat, stringsAsFactors=F)

plot(va.mat$pct, colMeans(va.mat[,-(1:2)]), pch=".", xlim=c(0,8), ylim=c(0,8), xlab="Va frequency in primary alpha.vgene", ylab="Va frequency in alpha2.vgene")
text(va.mat$pct+0.2, colMeans(va.mat[,-(1:2)]), Va, cex=0.8, adj=0)
lines(c(0,8), c(0,8))

x <- table(datnt2$alpha.vgene[aa])
z <- is.element(Va, names(x)[x>400])		#limit the matrix plot to TRAVs with primary >400)
xz <- as.matrix(va.mat[,-(1:2)])[z,z]
image(seq_len(nrow(xz)), seq_len(ncol(xz)), xz, axes=FALSE)
text(0, 1:length(Va[z]), Va[z], adj=1, cex=0.5)
text(1:length(Va[z]), 0, Va[z], adj=1, cex=0.5, srt=90)

## Cell-types in which double aa and bb live			>> No particular representation (low in gdt, of course)
#uses the standard cell types in annot$level1, in lev1 vector
w <- datnt2$alpha.functionality=="productive" & datnt2$beta.functionality=="productive"
res <- matrix(0, length(lev1), 5); colnames(res) <- c("N", "Naa", "Nbb", "Pctaa", "Pctbb")
for(i in 1:length(lev1)) {
	z <- is.element(datnt2$IGT.cellID, annot$IGT.cellID[annot$level1==lev1[i]]) & w
	res[i,1] <- sum(z)
	res[i,2] <- sum(z&aa)
	res[i,3] <- sum(z&bb)
}
res[,4] <- round(100*res[,2]/res[,1], digits=2)
res[,5] <- round(100*res[,3]/res[,1], digits=2)
doubleTCRinLineages <- data.frame(lev1, res, stringsAsFactors=F)

# frequency of double-cells in CD4 subclusters		>> exactly the same distribution in CD4 or CD8, slight difference in Tregs
z <- is.element(annot$IGT.cellID, datnt2$IGT.cellID) & annot$level1=="CD4"
x <- table(annot$level2[z])
all4 <- round(100*x/sum(x), digits=1)
z <- is.element(annot$IGT.cellID, datnt2$IGT.cellID[aa]) & annot$level1=="CD4"
x <- table(annot$level2[z])
aa4 <- round(100*x/sum(x), digits=1)
z <- is.element(annot$IGT.cellID, datnt2$IGT.cellID[bb]) & annot$level1=="CD4"
x <- table(annot$level2[z])
bb4 <- round(100*x/sum(x), digits=1)



## Look for dual-b in DP thymocytes
# could use DPs from other IGTs, but largest and best defined is IGT71/72
k <- "IGT71/72"
	x <- dat[is.element(dat$IGT, c("IGT71", "IGT72")),]
	smps <- unique(x$sample.name)
	orgs <- rep("-", length(smps)); for (i in 1:length(orgs)) orgs[i] <- unique(x$organ[x$sample.name==smps[i]])
	res <- matrix(0, length(smps), 5); colnames(res) <- c("TotalCells", "ab", "aab", "abb", "aabb")
for (j in (1:length(smps))) {
y <- x[x$sample.name==smps[j],]
w <- y$alpha.functionality=="productive" & y$beta.functionality=="productive"
res[j,1] <- nrow(y)
res[j,2] <- sum(w & y$alpha2.vgene=="" & y$beta2.vgene=="")
res[j,3] <- sum(w & y$alpha2.functionality=="productive" & y$beta2.vgene=="")
res[j,4] <- sum(w & y$alpha2.vgene=="" & y$beta2.functionality=="productive")
res[j,5] <- sum(w & y$alpha2.functionality=="productive" & y$beta2.functionality=="productive")
}	#end j
DPchain.summary <- data.frame(rep(k, length(smps)), orgs, smps, res, round(res[,-1]/rowSums(res[,-1]), digits=2), stringsAsFactors=F)
colnames(DPchain.summary)[(ncol(DPchain.summary)-3):ncol(DPchain.summary)] <- paste("Pct", colnames(res)[-1], sep=".")

#5b calculate single and double frequencies, Va and Vb for table
#uses datnt2

# For table
print(round(100*sum(is.element(datnt2$alpha.vgene, c("TRAV14-1", "TRAV14-2", "TRAV14-3", "TRAV14D-1")))/sum(datnt2$alpha.functionality=="productive"), digits=2))
print(round(100*sum(is.element(datnt2$alpha.vgene, c("TRAV4-2", "TRAV4-3", "TRAV4-3")))/sum(datnt2$alpha.functionality=="productive"), digits=2))
x <- (is.element(datnt2$alpha.vgene, c("TRAV14-1", "TRAV14-2", "TRAV14-3", "TRAV14D-1")) & is.element(datnt2$alpha2.vgene, c("TRAV4-2", "TRAV4-3", "TRAV4-3"))) | (is.element(datnt2$alpha2.vgene, c("TRAV14-1", "TRAV14-2", "TRAV14-3", "TRAV14D-1")) & is.element(datnt2$alpha.vgene, c("TRAV4-2", "TRAV4-3", "TRAV4-3")))
print(round(100*sum(x)/sum(datnt2$alpha.functionality=="productive"), digits=3))


x <- (is.element(datnt2$beta.vgene, c("TRBV13-1", "TRBV13-2", "TRBV13-3"))) | (is.element(datnt2$beta2.vgene, c("TRBV13-1", "TRBV13-2", "TRBV13-3")))
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=2))
x <- datnt2$beta.vgene=="TRBV17" | datnt2$beta2.vgene=="TRBV17"
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=2))
x <- datnt2$beta.vgene=="TRBV1" | datnt2$beta2.vgene=="TRBV1"
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=2))
x <- datnt2$beta.vgene=="TRBV31" | datnt2$beta2.vgene=="TRBV31"
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=2))
x <- datnt2$beta.vgene=="TRBV16" | datnt2$beta2.vgene=="TRBV16"
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=2))
x <- datnt2$beta.vgene=="TRBV19" | datnt2$beta2.vgene=="TRBV19"
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=2))


x <- (is.element(datnt2$beta.vgene, c("TRBV13-1", "TRBV13-2", "TRBV13-3")) & datnt2$beta2.vgene=="TRBV17") | (is.element(datnt2$beta2.vgene, c("TRBV13-1", "TRBV13-2", "TRBV13-3")) & datnt2$beta.vgene=="TRBV17")
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=3))

x <- (datnt2$beta.vgene=="TRBV1" & datnt2$beta2.vgene=="TRBV31") | (datnt2$beta2.vgene=="TRBV1" & datnt2$beta.vgene=="TRBV31")
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=3))

x <- (datnt2$beta.vgene=="TRBV16" & datnt2$beta2.vgene=="TRBV19") | (datnt2$beta2.vgene=="TRBV16" & datnt2$beta.vgene=="TRBV19")
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=3))


x <- (is.element(datnt2$beta.vgene, c("TRBV13-1", "TRBV13-2", "TRBV13-3")) & datnt2$beta2.vgene=="TRBV19") | (is.element(datnt2$beta2.vgene, c("TRBV13-1", "TRBV13-2", "TRBV13-3")) & datnt2$beta.vgene=="TRBV19")
print(round(100*sum(x)/sum(datnt2$beta.functionality=="productive"), digits=3))

## From datc, lists of ab, aab, abb and aabb cells to verify UMI counts; in datc, all alpha and betas are productive
ALXcells <- matrix("-", 16000, 4); colnames(ALXcells) <- c( "ab", "aab", "abb", "aabb")
ALXcells[,1] <- sample_cells(datc$IGT.cellID[which(datc$alpha2.functionality!="productive" & datc$beta2.functionality!="productive")], 16000)
ALXcells[,2] <- sample_cells(datc$IGT.cellID[which(datc$alpha2.functionality=="productive" & datc$beta2.functionality!="productive")], 16000)
ALXcells[,3] <- sample_cells(datc$IGT.cellID[which(datc$alpha2.functionality!="productive" & datc$beta2.functionality=="productive")], 16000)
aabb.ids <- datc$IGT.cellID[which(datc$alpha2.functionality=="productive" & datc$beta2.functionality=="productive")]
z <- length(aabb.ids)
ALXcells[seq_len(min(z, nrow(ALXcells))),4] <- aabb.ids[seq_len(min(z, nrow(ALXcells)))]

dev.off()
message("Plots saved to: ", file.path(script_dir, "5. Cells with 2 betas or 2 alphas.pdf"))
