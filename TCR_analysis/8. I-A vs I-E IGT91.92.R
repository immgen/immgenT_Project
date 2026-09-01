## To analyze the TCRs in the A, E and KO , adapted for IGT91/92 in 0125 release, moved to 0725
# inherits complete "dat" and "annot" from 1.Dataprep

######### Load dataprep outputs #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

load(file.path(script_dir, "1. Dataprep and basics 0226.RData"))
# dat, annot, Va, Vb, Ja, Jb

GSE297097.mde.incremental.level1.All.data <- read.csv(
  file.path(script_dir, "GSE297097_mde_incremental_level1_All_data.csv"),
  stringsAsFactors = FALSE, row.names = 1
)
GSE297097.mde.incremental.level1.All.data$X <- GSE297097.mde.incremental.level1.All.data$mdeTOTALVI_1
GSE297097.mde.incremental.level1.All.data$Y <- GSE297097.mde.incremental.level1.All.data$mdeTOTALVI_2
GSE297097.mde.incremental.level1.All.data$IGT.cellID <- rownames(GSE297097.mde.incremental.level1.All.data)

pdf(file.path(script_dir, "8. I-A vs I-E IGT91.92.pdf"), width = 10, height = 10)

datae <- dat[is.element(dat$IGT, c("IGT91", "IGT92")),]  # & is.element(dat$sample.name, c("A", "E", "O")),]

# Reduce to productive only
z <- datae$alpha.functionality=="productive" & datae$beta.functionality=="productive"
sum(z)
datae <- datae[z,]

datae$sample.name <- gsub("_", ".", datae$sample.name, fixed = TRUE)

# sample.type ### CAUTION: HARDCODED TO FIT 91/92 sample names
samp.type <- rep("A", nrow(datae))
samp.type[regexpr("E.", datae$sample.name)>0] <- "E"
samp.type[regexpr("EA", datae$sample.name)>0] <- "EA"
samp.type[regexpr("O.", datae$sample.name)>0] <- "KO"
samp.types <- sort(unique(samp.type))


## Are there shared clonotypes (protein level, both chains) between A and E datasets?
x <- paste(datae$alpha.vgene, datae$alpha.junction, datae$beta.vgene, datae$beta.junction, sep=".")
sum(is.element(x[samp.type=="A"], x[samp.type=="E"]))
## RESULT: zero clonotypes shared between A and E samples


## iNKT in A and E
z <- datae$alpha.vgene=="TRAV11" & datae$alpha.jgene=="TRAJ18" & samp.type=="E"
sum(z)
sum(samp.type=="E")
table(datae$beta.vgene[z])
table(datae$alpha.junction[z])


###

## Make TRAV and TRBV matrices, one column for each sample.name
Va.ae9192 <- matrix(0, length(Va), length(samp.types))
Vb.ae9192 <- matrix(0, length(Vb), length(samp.types))
rownames(Va.ae9192) <- Va; colnames(Va.ae9192) <- samp.types
rownames(Vb.ae9192) <- Vb; colnames(Vb.ae9192) <- samp.types

for (i in 1:length(Va)) {
  for (j in 1:length(samp.types)) {Va.ae9192[i,j] <- sum(datae$alpha.vgene==Va[i] & samp.type==samp.types[j])}
}

for (i in 1:length(Vb)) {
  for (j in 1:length(samp.types)) {Vb.ae9192[i,j] <- sum(datae$beta.vgene==Vb[i] & samp.type==samp.types[j])}
}
Va.ae9192.freq <- round(t(t(Va.ae9192)/colSums(Va.ae9192)*100), digits=1); Va.ae9192.freq <- Va.ae9192.freq[,order(samp.types)]
Vb.ae9192.freq <- round(t(t(Vb.ae9192)/colSums(Vb.ae9192)*100), digits=1); Vb.ae9192.freq <- Vb.ae9192.freq[,order(samp.types)]


par(mfrow=c(1,1))
x <- rowMeans(Va.ae9192[,samp.types=="A", drop=FALSE]); x <- round(x/sum(x)*100, digits=1)
y <- rowMeans(Va.ae9192[,samp.types=="E", drop=FALSE]); y <- round(y/sum(y)*100, digits=1)
plot(x,y,pch=16, xlim=c(0,8), ylim=c(0,8), xlab="Va Freq in A", ylab="Va Freq in E"); title(main="IGT91.92")
x <- rowMeans(Vb.ae9192[,samp.types=="A", drop=FALSE]); x <- round(x/sum(x)*100, digits=1)
y <- rowMeans(Vb.ae9192[,samp.types=="E", drop=FALSE]); y <- round(y/sum(y)*100, digits=1)
y[8] <- 0.6	# te separate TRBV15 and 16
plot(x,y,pch=16, xlab="Vb Freq in A", ylab="Vb Freq in E", xlim=c(0,14), ylim=c(0,14)); title(main="IGT91.92")
plot(x,y,pch=".", xlab="Vb Freq in A", ylab="Vb Freq in E", xlim=c(0,14), ylim=c(0,14)); title(main="IGT91.92")
text(x+0.2,y,rownames(Vb.ae9192), cex=0.7, adj=0)

## mva-like plots
#E vs A
x <- rowMeans(Va.ae9192[,samp.types=="A", drop=FALSE]); x <- round(x/sum(x)*100, digits=1)
y <- rowMeans(Va.ae9192[,samp.types=="E", drop=FALSE]); y <- round(y/sum(y)*100, digits=1)
plot((x+y)/2, (y+0.05)/(x+0.05), pch=".", log="y", ylim=c(0.05,20), xlab="MeanVa frequency", ylab="Va E/A ratio")
text(((x+y)/2)+0.05, (y+0.05)/(x+0.05), rownames(Va.ae9192), cex=0.7, adj=0)
lines(c(0,6), c(1,1)); title(main="Differential I-E/I-A Va usage 91.92")
#same, Vb
x <- rowMeans(Vb.ae9192[,samp.types=="A", drop=FALSE]); x <- round(x/sum(x)*100, digits=1)
y <- rowMeans(Vb.ae9192[,samp.types=="E", drop=FALSE]); y <- round(y/sum(y)*100, digits=1)
y[8] <- 0.35
plot((x+y)/2, (y+0.05)/(x+0.05), pch=".", log="y", ylim=c(0.05,20), xlab="Mean Vb frequency", ylab="Vb E/A ratio")
text(((x+y)/2)+0.05, (y+0.05)/(x+0.05), rownames(Vb.ae9192), cex=0.7, adj=0)
lines(c(0,12), c(1,1)); title(main="Differential I-E/I-A Vb usage 91.92")


### Permutation test of the Va and Vb usage differences
#Va
x <- rowMeans(Va.ae9192[,samp.types=="A", drop=FALSE]); x <- round(x/sum(x)*100, digits=1)
y <- rowMeans(Va.ae9192[,samp.types=="E", drop=FALSE]); y <- round(y/sum(y)*100, digits=1)
topVa <- rownames(Va.ae9192)[((x+y)/2)>2]

z <- samp.type=="A" | samp.type=="E"
w <- samp.type[z]
X <- datae$alpha.vgene[z]
x <- table(X[w=="A"]); x <- x/sum(x); x <- x[match(topVa, names(x), nomatch=0)]
y <- table(X[w=="E"]); y <- y/sum(y); y <- y[match(topVa, names(y), nomatch=0)]
topVa.earat9192 <- y/x

iter <- 100
res <- matrix(0, length(topVa), iter)
for (k in 1:iter){
  z <- sample(w)
  x <- table(X[z=="A"]); x <- x/sum(x); x <- x[match(topVa, names(x), nomatch=0)]
  y <- table(X[z=="E"]); y <- y/sum(y); y <- y[match(topVa, names(y), nomatch=0)]
  res[,k] <- y/x
}
w <- pnorm(topVa.earat9192, mean=apply(res, 1, mean), sd=apply(res, 1, sd))
w <- pmin(w, 1-w)
w <- w*length(topVa)		#Bonferroni correction
w <- (-log10(w))
w[w<1.4] <- 0
topVa.earat9192 <- data.frame(topVa, as.numeric(round(topVa.earat9192, digits=3)), round(w, digits=2), stringsAsFactors=F)
topVa.earat9192 <- topVa.earat9192[order(topVa.earat9192[,2]),]
colnames(topVa.earat9192) <- c("TRAV", "EtoA ratio", "-log10(padj)")

#Same, Vb
x <- rowMeans(Vb.ae9192[,samp.types=="A", drop=FALSE]); x <- round(x/sum(x)*100, digits=1)
y <- rowMeans(Vb.ae9192[,samp.types=="E", drop=FALSE]); y <- round(y/sum(y)*100, digits=1)
topVb <- rownames(Vb.ae9192)[((x+y)/2)>2]
z <- samp.type=="A" | samp.type=="E"
w <- samp.type[z]
X <- datae$beta.vgene[z]
x <- table(X[w=="A"]); x <- x/sum(x); x <- x[match(topVb, names(x), nomatch=0)]
y <- table(X[w=="E"]); y <- y/sum(y); y <- y[match(topVb, names(y), nomatch=0)]
topVb.earat9192 <- y/x
iter <- 100
res <- matrix(0, length(topVb), iter)
for (k in 1:iter){
  z <- sample(w)
  x <- table(X[z=="A"]); x <- x/sum(x); x <- x[match(topVb, names(x), nomatch=0)]
  y <- table(X[z=="E"]); y <- y/sum(y); y <- y[match(topVb, names(y), nomatch=0)]
  res[,k] <- y/x
}
w <- pnorm(topVb.earat9192, mean=apply(res, 1, mean), sd=apply(res, 1, sd))
w <- pmin(w, 1-w)
w <- w*length(topVb)		#Bonferroni correction
w <- (-log10(w))
w[w<1.4] <- 0
topVb.earat9192 <- data.frame(topVb, as.numeric(round(topVb.earat9192, digits=3)), round(w, digits=2), stringsAsFactors=F)
topVb.earat9192 <- topVb.earat9192[order(topVb.earat9192[,2]),]
colnames(topVb.earat9192) <- c("TRBV", "EtoA ratio", "-log10(padj)")


### Analyze charact of survivor cells with deleted Vbs in E mice
EdelVb <- rownames(Vb.ae9192.freq)[(Vb.ae9192.freq[,"E"]/Vb.ae9192.freq[,"A"])<0.3]
EdelVb.A.cells9192 <- datae$IGT.cellID[is.element(datae$beta.vgene, EdelVb) & samp.type=="A"]
EdelVb.E.cells9192 <- datae$IGT.cellID[is.element(datae$beta.vgene, EdelVb) & samp.type=="E"]

# Plot on mde
b <- sample((1:nrow(GSE297097.mde.incremental.level1.All.data)), 20000)
za <- is.element(GSE297097.mde.incremental.level1.All.data$IGT.cellID, EdelVb.A.cells9192)
ze <- is.element(GSE297097.mde.incremental.level1.All.data$IGT.cellID, EdelVb.E.cells9192)
plot(GSE297097.mde.incremental.level1.All.data$X[b], GSE297097.mde.incremental.level1.All.data$Y[b], pch=".", col=16)
points(GSE297097.mde.incremental.level1.All.data$X[za], GSE297097.mde.incremental.level1.All.data$Y[za], pch=16, cex=0.6, col=2)
title(main="CD4+ cells expressing E+mtv-sensitive Vh in H2-A+ mice 91.92")
plot(GSE297097.mde.incremental.level1.All.data$X[b], GSE297097.mde.incremental.level1.All.data$Y[b], pch=".", col=16)
points(GSE297097.mde.incremental.level1.All.data$X[ze], GSE297097.mde.incremental.level1.All.data$Y[ze], pch=16, cex=0.6, col=5)
title(main="CD4+ cells expressing E+mtv-sensitive Vh in H2-E+ mice 91.92")

#cell-types (91/92 were sorted on CD4+)
za <- is.element(annot$IGT.cellID, EdelVb.A.cells9192)
ze <- is.element(annot$IGT.cellID, EdelVb.E.cells9192)
x <- table(annot$level1[za]); print(round(100*x/sum(x), digits=1))
x <- table(annot$level1[ze]); print(round(100*x/sum(x), digits=1))
x <- table(annot$level2[ze]); print(round(100*x/sum(x), digits=1))

#compare Va and Jb
EdelVb.A.Va9192 <- EdelVb.E.Va9192 <- rep(0, length(Va))
x <- datae$alpha.vgene[is.element(datae$beta.vgene, EdelVb) & samp.type=="A"]
y <- datae$alpha.vgene[is.element(datae$beta.vgene, EdelVb) & samp.type=="E"]
for (i in 1:length(Va)) {EdelVb.A.Va9192[i] <- sum(x==Va[i]); EdelVb.E.Va9192[i] <- sum(y==Va[i])}
plot(EdelVb.E.Va9192, EdelVb.A.Va9192, xlab="Va frequency in E+ Vb deletion survivors", ylab="Va frequency in A +"); title(main="IGT91.92")
plot(EdelVb.E.Va9192, EdelVb.A.Va9192, pch=".", xlab="Va frequency in E+ Vb deletion survivors", ylab="Va frequency in A +"); title(main="IGT91.92")
text(EdelVb.E.Va9192+0.1, EdelVb.A.Va9192, Va, cex=0.7, adj=0)

EdelVb.A.Jb9192 <- EdelVb.E.Jb9192 <- rep(0, length(Jb))
x <- datae$beta.jgene[is.element(datae$beta.vgene, EdelVb) & samp.type=="A"]
y <- datae$beta.jgene[is.element(datae$beta.vgene, EdelVb) & samp.type=="E"]
for (i in 1:length(Jb)) {EdelVb.A.Jb9192[i] <- sum(x==Jb[i]); EdelVb.E.Jb9192[i] <- sum(y==Jb[i])}
plot(EdelVb.E.Jb9192, EdelVb.A.Jb9192, xlab="Jb numbers in E+ Vb deletion survivors", ylab="Jb numbers in A +"); title(main="IGT91.92")
plot((100*EdelVb.E.Jb9192/sum(EdelVb.E.Jb9192)), (100*EdelVb.A.Jb9192/sum(EdelVb.A.Jb9192)), pch=".", xlab="Jb pct in E+ Vb deletion survivors", ylab="Jb pct in A +")
text((100*EdelVb.E.Jb9192/sum(EdelVb.E.Jb9192)), (100*EdelVb.A.Jb9192/sum(EdelVb.A.Jb9192)), Jb, cex=0.7, adj=0)
title(main="IGT91.92")

dev.off()
message("Plots saved to: ", file.path(script_dir, "8. I-A vs I-E IGT91.92.pdf"))
