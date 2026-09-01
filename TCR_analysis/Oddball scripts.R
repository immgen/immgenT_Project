# Oddballs scripts for R conversion

######### Load dataprep outputs #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

load(file.path(script_dir, "1. Dataprep and basics 0226.RData"))
# datc, datntr, va.duplic.simplif (if saved) ...

random.TRA.seqs <- read.csv(
  file.path(script_dir, "random_TRA_seqs.csv"),
  stringsAsFactors = FALSE
)
if (!exists("va.duplic.simplif")) {
  va.duplic.simplif <- read.csv(
    file.path(script_dir, "va.duplic.simplif.csv"),
    stringsAsFactors = FALSE
  )
}

pdf(file.path(script_dir, "Oddball scripts.pdf"), width = 10, height = 10)

# TCRa Oddballs analysis
# use datc and datntr from 1.Dataprep

## quantitate N across VaJa pairs - only use the primary alpha
apair <- paste(datc$alpha.vgene, datc$alpha.jgene, sep=".")
apairs.cnt <- table(apair)
apairs.cnt[apairs.cnt>1000] <- 1000
apairs <- names(apairs.cnt)
apairs.cnt <- as.numeric(apairs.cnt)
names(apairs.cnt) <- apairs
apairs.main3p <- apairs.maintonext3p <- apairs.n3p <- rep(0, length(apairs))
for(i in 1:length(apairs)) {
  y <- rev(sort(table(datc$alpha.junction[apair==apairs[i]])))
  if (length(y) == 0) next
  apairs.n3p[i] <- length(y)
  apairs.main3p[i] <- round(y[1]/sum(y), digits=3)
  apairs.maintonext3p[i] <- if (length(y) > 1) round(y[1]/y[2], digits=3) else NA
}
#plots
plot(as.numeric(apairs.cnt), type="h", xlab="VaJa pair index", ylab="Count")
plot(apairs.cnt, apairs.main3p, pch=".", xlim=c(0, 1200), ylim=c(0,1), xlab="VaJa pair occurence", ylab="Frequency of most abundant CDR3")
plot(apairs.n3p, apairs.main3p, pch=".", ylim=c(0,1), xlab="Nb of junctions for that VaJa pair", ylab="Frequency of most abundant CDR3")
plot(apairs.cnt, apairs.n3p, pch=".", xlim=c(0, 1200), xlab="VaJa pair occurence", ylab="Nb of junctions for that VaJa pair")

z <- (apairs.cnt>100 & apairs.main3p>0.15) | (apairs.cnt>50 & apairs.main3p>0.3)
plot(apairs.cnt, apairs.main3p, xlim=c(0, 1200), ylim=c(0,1), pch=".", xlab="VJ pair occurence", ylab="Frequency of most abundant CDR3 for VJ pair")
text(apairs.cnt[z]+5, apairs.main3p[z], apairs[z], adj=0, cex=0.5)
plot(apairs.cnt, apairs.n3p, pch=".", xlim=c(0, 1200), xlab="VaJa pair occurence", ylab="Nb of junctions for that VaJa pair")
text(apairs.cnt[z]+5, apairs.n3p[z], apairs[z], adj=0, cex=0.5)

plot(apairs.cnt, apairs.main3p, pch=".", xlim=c(0, 1200), ylim=c(0,1), xlab="VaJa pair occurence", ylab="Frequency of most abundant CDR3")
fit <- supsmu(apairs.cnt, apairs.main3p)
lines(fit, lwd=2, col=3)
lines(c(0,1200), c(0.09, 0.09))

plot(apairs.maintonext3p, apairs.main3p, ylim=c(0,1), pch=".", xlab="VaJa First/Second CDR3 ratio", ylab="Frequency of most abundant CDR3")
plot(apairs.cnt, apairs.maintonext3p, xlim=c(0, 1200),  xlab="VaJa pair occurence", ylab="VaJa First/Second CDR3 ratio")
plot(apairs.cnt, apairs.maintonext3p, pch=".", xlim=c(0, 1200),  xlab="VaJa pair occurence", ylab="VaJa First/Second CDR3 ratio")
text(apairs.cnt+5, apairs.maintonext3p, apairs, cex=0.5, adj=0)

## Calculate p.values on  main3p
apairs.main3pval <- rep(0, length(apairs))
sel <- apairs.cnt>50
z <- quantile(apairs.cnt[sel], 0:5/5); z[1] <- z[1]-1; z[6] <- z[6]+1
for (i in 1:(length(z)-1)) {
  comp <- apairs.cnt>z[i] & apairs.cnt<(z[i+1]+1)
  w <- apairs.main3p[comp]
  med <- median(w)
  MAD <- median(abs(w - med))
  zsc <- 0.6745 * (w - med) / MAD
  p <- 2 * pnorm(-abs(zsc))
  apairs.main3pval[comp] <- -log10(p.adjust(p, method = "BH"))
}

z <- apairs.main3pval>1.3
plot(apairs.cnt, apairs.main3p, pch=".", xlim=c(0, 1200), ylim=c(0,1), xlab="VaJa pair occurence", ylab="Frequency of most abundant CDR3")
points(apairs.cnt[z], apairs.main3p[z], pch=16, col=3)

#Make a table of oddballs
oddballs <- apairs[(apairs.cnt>70 & apairs.main3p>0.15) | (apairs.cnt>50 & apairs.main3p>0.3)]
odb.char0226 <- data.frame(
  matrix(0, length(oddballs),7),
  matrix("---", length(oddballs),3),
  stringsAsFactors=F
)
colnames(odb.char0226) <- c("TotPairOccr", "FreqMainAjunct", "FirstToSecondMain3aRatio", "CD4", "Treg", "CD8", "NonConv", "MainAjunct", "MostFreqVB", "MainOrgan")
for (i in 1:length(oddballs)) {
  odb.char0226[i,1] <- sum(apair==oddballs[i])
  odb.char0226[i,2] <- apairs.main3p[apairs==oddballs[i]]
  odb.char0226[i,3] <- apairs.maintonext3p[apairs==oddballs[i]]
  z <- apair==oddballs[i]
  y <- table(datc$alpha.junction[z])
  w <- names(y)[y==max(y)][1]
  odb.char0226[i,8] <- w
  z <- z & datc$alpha.junction==w
  y <- table(datc$beta.vgene[z]); y <- y[y==max(y)][1]
  odb.char0226[i,9] <- paste(names(y), round(as.numeric(y)/sum(apair==oddballs[i]), digits=2), sep="/")
  y <- table(datc$organ[z]); y <- y[y==max(y)][1]
  odb.char0226[i,10] <- paste(names(y), as.numeric(y), sep="/")
  odb.char0226[i,4] <- sum(datc$level1[z]=="CD4")
  odb.char0226[i,5] <- sum(datc$level1[z]=="Treg")
  odb.char0226[i,6] <- sum(datc$level1[z]=="CD8")
  odb.char0226[i,7] <- sum(is.element(datc$level1[z], c("Tz", "CD8aa", "DN", "DP")))
}
Ndiv <- NbNucj <- FreqMainAnuc <- rep(0, length(oddballs))
NInMainAnuc <- MainAnuc <- rep(".", length(oddballs))
for (i in 1:length(oddballs)) {
  z <- is.element(apair, oddballs[i]) & is.element(datc$alpha.junction, odb.char0226$MainAjunct[i])
  x <- datc$alpha.junction.nt[z]
  y <- rev(sort(table(x)))
  NbNucj[i] <- length(y); FreqMainAnuc[i] <- round(max(y)/sum(y), digits=2)
  z <- z & datc$alpha.junction.nt==names(y)[1]
  NInMainAnuc[i] <- datc$alpha.vj.nregion[z][1]
  MainAnuc[i] <- names(y)[1]
}

odb.char0226 <- data.frame(1:length(oddballs), oddballs, odb.char0226, NbNucj, FreqMainAnuc, NInMainAnuc, MainAnuc, stringsAsFactors=F)
colnames(odb.char0226)[1:2] <- c("odbID", "VaJa")
odb.char0226 <- odb.char0226[rev(order(odb.char0226$FreqMainAjunct)),]

z <- match(odb.char0226[,2], apairs, nomatch=0)
Padj <- round(apairs.main3pval[z], digits=2)
odb.char0226 <- data.frame(odb.char0226, Padj, stringsAsFactors=F)
z <- match(odb.char0226[,2], apair, nomatch=0)
Pgen <- datc$alpha.Pgen.aa[z]
odb.char0226 <- data.frame(odb.char0226, Pgen, stringsAsFactors=F)

### Oddballs in TCRlo DPs ? version 0226 >> Table S6A
x <- datntr[is.element(datntr$IGT, c("IGT71","IGT72")) & regexpr("TRAV", datntr$alpha.vgene)>0 & datntr$alpha.functionality=="productive",]
apairx <- paste(x$alpha.vgene, x$alpha.jgene, sep=".")

odb.freqs.dp <- matrix(0, nrow(odb.char0226), 4)
colnames(odb.freqs.dp) <- c("Nb.VaJain10K.T", "Nb.canonic3In10K.T", "Nb.VaJain10K.DP", "Nb.canonic3in10K.DP")
rownames(odb.freqs.dp) <- odb.char0226$VaJa
for (i in 1:nrow(odb.char0226)) {
  z <- is.element(apair, odb.char0226$VaJa[i])
  odb.freqs.dp[i,1] <- (10^4)*sum(z)/length(apair)
  z <- is.element(apair, odb.char0226$VaJa[i]) & is.element(datc$alpha.junction, odb.char0226$MainAjunct[i])
  odb.freqs.dp[i,2] <- (10^4)*sum(z)/length(apair)
  z <- is.element(apairx, odb.char0226$VaJa[i])
  odb.freqs.dp[i,3] <- (10^4)*sum(z)/length(apairx)
  z <- is.element(apairx, odb.char0226$VaJa[i]) & is.element(x$alpha.junction, odb.char0226$MainAjunct[i])
  odb.freqs.dp[i,4] <- (10^4)*sum(z)/length(apairx)
}
plot(odb.freqs.dp[,1]+0.01, odb.freqs.dp[,3]+0.01, log="xy", xlim=c(0.01,300), ylim=c(0.01,300), xlab="VaJaN/10K in whole datc", ylab="VaJaN/10K in 71-72 DPs")
plot(odb.freqs.dp[,2]+0.01, odb.freqs.dp[,4]+0.01, log="xy", xlim=c(0.01,300), ylim=c(0.01,300), xlab="Canonic N/10K in whole datc", ylab="Canonic N/10K in 71-72 DPs")

plot(odb.freqs.dp[,2]+0.01, odb.freqs.dp[,4]+0.01, type="n", log="xy", xlim=c(0.01,300), ylim=c(0.01,300), xlab="Canonic N/10K in whole datc", ylab="Canonic N/10K in 71-72 DPs")
text(odb.freqs.dp[,2]+0.01, odb.freqs.dp[,4]+0.01, 1:nrow(odb.char0226), cex=0.7)
text((odb.freqs.dp[,1]+0.01)[1:2], (odb.freqs.dp[,2]+0.01)[1:2], c("N", "M"), col=3, cex=0.8)
lines(c(0.01,300), c(0.01,300))

odb.char0226 <- data.frame(odb.char0226, round(odb.freqs.dp, digits=2), stringsAsFactors=F)

#### Same for beta chain
bpair <- paste(datc$beta.vgene, datc$beta.jgene, sep=".")
bpairs.cnt <- table(bpair)
bpairs <- names(bpairs.cnt)
bpairs.cnt <- as.numeric(bpairs.cnt)
names(bpairs.cnt) <- bpairs
bpairs.main3p <- rep(0, length(bpairs))

for(i in 1:length(bpairs)) {
  y <- table(datc$beta.junction[bpair==bpairs[i]])
  bpairs.main3p[i] <- round(max(y)/sum(y), digits=3)
}

plot(bpairs.cnt, bpairs.main3p, ylim=c(0,1), pch=16, cex=0.5, xlab="VbJb pair occurence", ylab="Frequency of most abundant CDR3")
plot(bpairs.cnt, bpairs.main3p, pch=".", ylim=c(0,1), xlab="VbJb pair occurence", ylab="Frequency of most abundant CDR3")

z <- (bpairs.cnt>100 & bpairs.main3p>0.15) | (bpairs.cnt>50 & bpairs.main3p>0.3)
plot(bpairs.cnt, bpairs.main3p, ylim=c(0,1), pch=".", xlab="VbJb pair occurence", ylab="Frequency of most abundant CDR3")
if (any(z)) text(bpairs.cnt[z]+20, bpairs.main3p[z], bpairs[z], cex=0.7, adj=0)

#### Analyze CDR3 distribution in OLGA random data
dt <- random.TRA.seqs
colnames(dt)[colnames(dt) == "apha.junction"] <- "alpha.junction"

x <- match(dt$alpha.vgene, va.duplic.simplif$From, nomatch=0)
dt$alpha.vgene[x>0] <- va.duplic.simplif$To[x]

sim.apair <- paste(dt$alpha.vgene, dt$alpha.jgene, sep=".")
sim.apairs <- unique(sim.apair)

sim.apairs.cnt <- sim.apairs.main3p <- rep(0, length(sim.apairs))

for (i in 1:length(sim.apairs)) {
  z <- sim.apair==sim.apairs[i]
  sim.apairs.cnt[i] <- sum(z)
  x <- table(dt$alpha.junction[z])
  sim.apairs.main3p[i] <- round(max(x)/sum(z), digits=2)
}

z <- (sim.apairs.cnt>70 & sim.apairs.main3p>0.15) | (sim.apairs.cnt>50 & sim.apairs.main3p>0.3)
plot(sim.apairs.cnt,sim.apairs.main3p,pch=".", xlab="Total occurences, all individuals", ylab="Frequency of most frequent junction")
points(sim.apairs.cnt[z], sim.apairs.main3p[z], pch=16, col=2, cex=0.6)

plot(sim.apairs.cnt, sim.apairs.main3p, pch=".", xlim=c(0, 1200), ylim=c(0,1), xlab="VaJa pair occurence", ylab="Frequency of most abundant CDR3")
fit <- supsmu(sim.apairs.cnt, sim.apairs.main3p)
lines(fit, lwd=2, col=2)
title(main="OLGA random joins, dominant CDR3 in VaJa pairs")

z <- sim.apairs.cnt>70 & sim.apairs.main3p>0.2
text(800, 0.98, paste(sum(z), "AJ with occurence>70 & main3p>0.2"), adj=0, cex=0.7)

dev.off()
message("Plots saved to: ", file.path(script_dir, "Oddball scripts.pdf"))
