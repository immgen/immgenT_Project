# 3 Reccurring betas, first observed in the 1A2B 22 scripts
#reuses same start
# Use datcntr for primary analysis to avoid clonal amplification and transgenic artefacts
# convert to dt, remove MAIT and NKT

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

# datcntr alias: use clean datc (see conversion notes above)
datcntr <- datc

# paired.igts: same-sample dual-lane IGT pairs (from prior analyses). If absent,
# the filter that drops paired-IGT-only recurrents is skipped.
paired.igts <- read.csv(file.path(script_dir, "paired.igts.csv"))

if (!exists("paired.igts")) {
  paired.igts <- NULL
  message("NOTE: paired.igts not found; skipping paired-IGT filter on recB.")
}

z <- (datcntr$alpha.vgene=="TRAV11" & datcntr$alpha.jgene=="TRAJ18") | (datcntr$alpha.vgene=="TRAV1" & datcntr$alpha.jgene=="TRAJ33")
dt <- datcntr[z==F,]

IGTsample <- paste(dt$IGT, dt$sample.name, sep=".")

dt$organ <- gsub("lymph node", "LN", dt$organ, fixed = TRUE)

#first narrow down on single-chain clonotypes (full, nt-level) that are duplicated across the entire dataset

allA <- paste(dt$alpha.vgene, dt$alpha.junction.nt, dt$alpha.junction, dt$alpha.jgene, sep=".")
dupA <- allA[duplicated(allA)]
allB <- paste(dt$beta.vgene, dt$beta.junction.nt, dt$beta.junction, dt$beta.jgene, sep=".")
dupB <- allB[duplicated(allB)]

length(dupA)

length(dupB)

######### Save all plots to PDF #####
pdf(file.path(script_dir, "3. Recurring betas r0725.pdf"), width = 10, height = 10)

par(mfrow = c(2, 2))

x <- table(dupA)
y <- sort(as.vector(x) + 1)
plot(seq_along(y), y, log = "y", xlim = c(0, 30000), ylim = c(1, 500),
     xlab = "rank", ylab = "count + 1", pch = 16, cex = 0.5)
title(main = "alpha nt.junction duplications")
text(0, 400,
     paste(sum(x[x > 4]), "occurrences of", sum(x > 4),
           "recurring alphas (5 or more)\namong", nrow(dt), "cells"),
     adj = 0, cex = 0.7)

x <- table(dupB)
y <- sort(as.vector(x) + 1)
plot(seq_along(y), y, log = "y", xlim = c(0, 30000), ylim = c(1, 500),
     xlab = "rank", ylab = "count + 1", pch = 16, cex = 0.5)
title(main = "beta nt.junction duplications")
text(0, 400,
     paste(sum(x[x > 4]), "occurrences of", sum(x > 4),
           "recurring betas (5 or more)\namong", nrow(dt), "cells"),
     adj = 0, cex = 0.7)

# Duplicated alpha and beta clonotypes tend to associate
sum(duplicated(allB) & duplicated(allA))
# but none correspond to full alpha+beta duplications
sum(dupAB <- duplicated(paste(allA,allB, sep=".")))


# Parse intra- vs inter-IGT duplications #here, only consider duplicated betas seen more than 5 times
x <- table(dupB)
x <- unlist(names(x)[x>4])
res <- matrix("-", length(x), 6)
rescnt <- matrix(0, length(x), 5)
resctyp <- matrix(0, length(x), length(cell.types))
colnames(resctyp) <- cell.types
colnames(res) <- c("dupB", "Vb", "Jb", "IGTs", "organ", "celltype")
colnames(rescnt) <- c("NbofCells", "NbofIGT", "Nnuc", "NbofAlphas", "MainAlphaOcc")

for (i in 1:length(x)) {
	z <- allB==x[i]
	res[i,1] <- x[i]
	res[i,2] <- unique(dt$beta.vgene[z])
	res[i,3] <- unique(dt$beta.jgene[z])
	res[i,4] <- paste(substring(unique(dt$IGT[z]),4), collapse=".")
	res[i,5] <- paste(unique(substring(dt$organ[z],1,5)), collapse="/")
	rescnt[i,1] <- sum(z)
	rescnt[i,2] <- length(unique(dt$IGT[z]))
	rescnt[i,3] <- sum(unique(nchar(dt$beta.vd.nregion[z])), unique(nchar(dt$beta.dj.nregion[z])), na.rm=T)
		w <- rev(sort(table(allA[z])))
	rescnt[i,4] <- length(w)
	rescnt[i,5] <- w[1]
		w <- as.character(annot$level1[is.element(annot$IGT.cellID, dt$IGT.cellID[z])])
	for (j in 1:length(cell.types)) resctyp[i,j] <- sum(w==cell.types[j])
		w <- table(w)
	res[i,6] <- paste(paste(names(w), w, sep=":"), collapse="/")

} #end i
recB <- data.frame(res, rescnt, resctyp, stringsAsFactors=F)
recB <- recB[rev(order(recB$NbofCells)),]

#trim the single IGT occurences
recB <- recB[recB$NbofIGT>1,]
if (!is.null(paired.igts)) {
	recB <- recB[is.element(recB$IGTs, gsub("IGT", "", paired.igts[,3]))==F,]
}

#count times each cell.type dominates
row_maxs <- function(m) apply(m, 1, max)
for (j in 1:length(cell.types)) {
	z <- colnames(recB)==cell.types[j]
	x <- sum(recB[,z]==row_maxs(recB[, cell.types]) & recB[,z]>3)
	print(paste(cell.types[j], ": ", x, sep=""))
}

#frequency in DPs
z <- is.element(dat$IGT, c("IGT72", "IGT73")) & dat$beta.functionality=="productive"
w <- paste(dat$beta.vgene[z], dat$beta.junction.nt[z], dat$beta.junction[z], dat$beta.jgene[z], sep=".")
pctThyDP <- rep(0, nrow(recB))
for (i in 1:nrow(recB)) pctThyDP[i] <- sum(w==recB$dupB[i])
pctThyDP <- round(100*pctThyDP/sum(z), digits=3)
pctAllT <- round(100*recB$NbofCells/nrow(dt), digits=3)
# S-PLUS par(fin=c(8,8)) is not settable the same way in R; reset layout instead
par(mfrow=c(1,1))
plot(pctAllT+runif(nrow(recB), -0.0003, 0.0003), pctThyDP+runif(nrow(recB), -0.0003, 0.0003), xlim=c(0, 0.02), ylim=c(0, 0.02), xlab="recB in all T data, %", ylab="recB in DP thymocytes IGT72/73, %")
title(main="Recurrent beta joins in thymus vs periphery")

recB <- data.frame(recB, pctAllT, pctThyDP, stringsAsFactors=F)

# add OLGA Pgen values
# 0226 SummaryTable already carries beta.Pgen.nt (replaces pGenAllData lookup)
Pgen.beta <- rep(0, nrow(recB))
for (i in 1:nrow(recB)) {
	z <- (1:nrow(dt))[allB==recB$dupB[i]][1]
	Pgen.beta[i] <- dt$beta.Pgen.nt[z]
}
recB <- data.frame(recB[,1:3], Pgen.beta, recB[,4:ncol(recB)], stringsAsFactors=F)


#select out CD8aa dominants
recB.8aa <- recB[recB$CD8aa>3 & recB$CD8aa==row_maxs(recB[, cell.types]),]
recB.8aa <- recB.8aa[order(recB.8aa$Jb, recB.8aa$Vb),]
recB.other <- recB[is.element(recB$dupB, recB.8aa$dupB)==F,]
recB.other <- recB.other[order(recB.other$Jb, recB.other$Vb),]
recB.4 <- recB[recB$CD4>3 & recB$CD4==row_maxs(recB[, cell.types]),]
recB.4 <- recB.4[rev(order(recB.4$CD4)),]


#distribution of Vb and Jb regions involved - do without CD8aa recBs, which have very skewed (Newbury) TCRs

# background Vb/Jb frequencies (%) in dt (not saved from earlier scripts)
Vb.freq <- sapply(Vb, function(v) 100 * mean(dt$beta.vgene == v, na.rm = TRUE))
Jb.freq <- sapply(Jb, function(j) 100 * mean(dt$beta.jgene == j, na.rm = TRUE))

recB.other <- recB[is.element(recB$dupB, recB.8aa$dupB)==F,]
x <- rep(0, length(Vb))
for (i in 1:length(Vb)) x[i] <- sum(recB.other$Vb==Vb[i])
plot(Vb.freq, 100*x/sum(x), pch=".", ylab="Vb frequency in non-CD8aa recurring beta joins")
text(Vb.freq+0.1, 100*x/sum(x), Vb, adj=0, cex=0.9)

x <- rep(0, length(Jb))
for (i in 1:length(Jb)) x[i] <- sum(recB.other$Jb==Jb[i])
plot(Jb.freq, 100*x/sum(x), pch=".", ylab="Jb frequency in non-CD8aa recurring beta joins")
text(Jb.freq+0.1, 100*x/sum(x), Jb, adj=0, cex=0.9)

#how many times does each VbJb pair occur ?
x <- rev(sort(table(paste(recB.other$Vb, recB.other$Jb, sep="."))))
plot(x, xlab="VbJb pairs", ylab="Nb of occurences among non-CD8aa recurring TCRb joins")
text((1:length(x))[1:10]+3, x[1:10]+runif(10,-0.3, 0.3), names(x)[1:10], cex=0.9, adj=0)

##closer look at the dominant VbJb combinations
#frequencies of different joins involving 13-2 and 2-7
z <- recB$Vb == "TRBV13-2" & recB$Jb == "TRBJ2-7"
recB.132.27 <- recB[z, ]
z <- dt$beta.vgene == "TRBV13-2" & dt$beta.jgene == "TRBJ2-7" & dt$beta.vgene != ""
x <- allB[z]
y <- table(x)
w <- rep(0, length(y))
for (i in 1:length(w)) w[i] <- unique(nchar(dt$beta.junction[allB == names(y)[i]]))

yv <- rev(sort(as.vector(y)))
plot(seq_along(yv), yv, xlab = "13-2 x 2-7 junctions", ylab = "Nb of occurences",
     pch = 16, cex = 0.5)
title(main = "All Vb13-2/Jb2-7 junctions")

plot(w + runif(length(y), -0.3, 0.3), y + runif(length(y), -0.2, 0.2),
     xlab = "Junction length", ylab = "Nb of occurences")
title(main = "All Vb13-2/Jb2-7 junctions, AA length vs nb of occurences")

#frequencies of different joins involving 3 and 2-7
z <- recB$Vb == "TRBV3" & recB$Jb == "TRBJ2-7"
recB.3.27 <- recB[z, ]
z <- dt$beta.vgene == "TRBV3" & dt$beta.jgene == "TRBJ2-7" & dt$beta.vgene != ""
x <- allB[z]
y <- table(x)
w <- rep(0, length(y))
for (i in 1:length(w)) w[i] <- unique(nchar(dt$beta.junction[allB == names(y)[i]]))

yv <- rev(sort(as.vector(y)))
plot(seq_along(yv), yv, xlab = "3 x 2-7 junctions", ylab = "Nb of occurences",
     pch = 16, cex = 0.5)
title(main = "All Vb3/Jb2-7 junctions")

plot(w + runif(length(y), -0.3, 0.3), y + runif(length(y), -0.2, 0.2),
     xlab = "Junction length", ylab = "Nb of occurences")
title(main = "All Vb3/Jb2-7 junctions, AA length vs nb of occurences")

dev.off()
message("Plots saved to: ", file.path(script_dir, "3. Recurring betas r0725.pdf"))
