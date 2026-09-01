# 11. Public specificities
# Use dt from datntr (no transgenics)

######### Load dataprep outputs #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

load(file.path(script_dir, "1. Dataprep and basics 0226.RData"))
# datntr, datbsl, annot, cell.types

paired.igts <- read.csv(file.path(script_dir, "paired.igts.csv"))

if (!exists("paired.igts")) {
  paired.igts <- NULL
  message("NOTE: paired.igts not found; skipping paired-IGT collapse.")
}

map <- data.frame(
  IGT.cellID = annot$IGT.cellID[annot$level1 == "Treg"],
  X = annot$mdeX[annot$level1 == "Treg"],
  Y = annot$mdeY[annot$level1 == "Treg"],
  stringsAsFactors = FALSE
)

pdf(file.path(script_dir, "11.Public Specificities.pdf"), width = 10, height = 10)

##select the dataset to use (then as dt)
dt <- datntr

#shorten to TCRs with junctions in both TRA and TRB
z <- dt$alpha.junction!="" & dt$beta.junction!=""
dt <- dt[z,]

all.Nclonos <- paste(dt$alpha.vgene, dt$alpha.junction.nt, dt$beta.vgene, dt$beta.junction.nt, sep=".")
all.Pclonos <- paste(dt$alpha.vgene, dt$alpha.junction, dt$beta.vgene, dt$beta.junction, sep=".")
dim(all.Pclonos)

#reset IGT and org vectors
IGT <- dt$IGT
IGT[is.element(IGT, c("IGT7", "IGT8", "IGT9"))] <- "IGT789"
if (!is.null(paired.igts)) {
  for(i in 1:nrow(paired.igts)) IGT[is.element(IGT, unlist(paired.igts[i,1:2]))] <- paired.igts[i,3]
}
org <- dt$organ

## Identify between-IGT Pclono repeats
rep.Pclonos <- table(all.Pclonos)
rep.Pclonos <- rep.Pclonos[rep.Pclonos>1]
rep.Pclonos <- rep.Pclonos[regexpr("TRAV", names(rep.Pclonos))>0 & regexpr("TRBV", names(rep.Pclonos))>0]

x <- rep(0, length(rep.Pclonos))
nm <- names(rep.Pclonos)
for (i in 1:length(rep.Pclonos)) {
  z <- all.Pclonos==nm[i]
  x[i] <- length(unique(dt$IGT[z]))
}
rep.Pclonos <- rep.Pclonos[x>1]

# run through rep.Pclonos, get characteristics
shared.Nclono <- b <- y <- n <- nn <- m <- susp <- rep(0, length(rep.Pclonos))
ct <- o <- w <- rep("-", length(rep.Pclonos))
shared.Pclono <- rep(FALSE, length(rep.Pclonos))
nm <- names(rep.Pclonos)
a <- c(10, 100, seq(1000,100000, by=1000))
for (i in 1:length(rep.Pclonos)) {
  x <- all.Pclonos==nm[i]
  y[i] <- length(unique(IGT[x]))
  z <- table(annot$level1[is.element(annot$IGT.cellID, dt$IGT.cellID[x])])
  ct[i] <- paste(paste(names(z)[z>0], z[z>0], sep=":"), collapse="/")
  z <- table(org[x])
  o[i] <- paste(paste(names(z), z, sep=":"), collapse="/")
  m[i] <- sum(all.Pclonos==nm[i])
  z <- substring(dt$IGT.cellID[x], regexpr("\\.", dt$IGT.cellID[x])+1)
  b[i] <- length(unique(z))
  susp[i] <- sum(duplicated(z))
  n[i] <- length(unique(all.Nclonos[x]))
  z <- table(IGT[x])
  w[i] <- paste(paste(names(z), z, sep=":"), collapse="/")
  z <- table(all.Nclonos[x]); cnt <- rep(0, length(z))
  for (j in 1:length(z)) cnt[j] <- length(unique(IGT[all.Nclonos==names(z)[j]]))
  nn[i] <- max(z)
  shared.Nclono[i] <- max(cnt)
  if (is.element(i, a)) print(paste(i, date()))
}
z <- y>1
rep.Pclonos <- data.frame(nm[z], ct[z], m[z], y[z], b[z], susp[z], w[z], o[z], n[z], nn[z], shared.Nclono[z], stringsAsFactors=F)
colnames(rep.Pclonos) <- c("Pclonos", "CellType", "TotCount", "NbIGTs", "NbBcds", "DupCelBcds",  "InIGT","InOrgans", "NbNucSqs", "CntMainNuc","IGTsSharingNclono" )
rep.Pclonos <- rep.Pclonos[order(rep.Pclonos$Pclonos),]

#Identify and remove contamination suspects
z <- rep.Pclonos$DupCelBcds>0
w <- (regexpr("IGT100", rep.Pclonos$InIGT)>0 | regexpr("IGT101", rep.Pclonos$InIGT)>0) & regexpr("IGT9", rep.Pclonos$InIGT)>0
w <- w | ((regexpr("IGT3:", rep.Pclonos$InIGT)>0) & regexpr("IGT44", rep.Pclonos$InIGT)>0)
suspect.rep.Pclonos <- rep.Pclonos[z|w,]
suspct.clono.cells <- dt$IGT.cellID[is.element(all.Pclonos, suspect.rep.Pclonos$Pclonos)]
rep.Pclonos <- rep.Pclonos[z==F & w==F,]

#rep.Pclonos with shared Nclonos
z <- rep.Pclonos$IGTsSharingNclono>1
sharedN.rep.Pclonos <- rep.Pclonos[z,]

#Plots
par(mfrow=c(2,2))
n <- nrow(rep.Pclonos)
plot(rep.Pclonos$TotCount*runif(n,0.92, 1.08), rep.Pclonos$NbNucSqs*runif(n,0.95, 1.05), log="xy", ylim=c(1,20), pch=".")
plot(rep.Pclonos$IGTsSharingNclono, rep.Pclonos$NbNucSqs, log="xy", xlim=c(1,50), ylim=c(1,20), pch=16)
plot(rep.Pclonos$TotCount*runif(n,0.92, 1.08), rep.Pclonos$NbIGTs*runif(n,0.95, 1.05), ylim=c(1,10), log="xy", pch=16, xlab="Total nb of Occurences", ylab="Nb of IGT harboring the clonotype")
y <- rep(0, max(rep.Pclonos$NbNucSqs))
for (i in 1:length(y)) y[i] <- sum(rep.Pclonos$NbNucSqs==i)
plot(1:length(y), y, type="b", xlim=c(0,30), log="y", xlab="Nb of Nucleotides sequences encoding the TRAB clonotype", ylab="Nb of Occurences")
x <- round(100*sum(rep.Pclonos$NbNucSqs==1)/nrow(rep.Pclonos), digits=1)
text(10, 0.9*max(y), paste("Pct of public TRAB clonotypes encoded\nby a single Nuc sequence:", x), cex=0.7, adj=0)

#label TRAV11
z <- regexpr("TRAV11", rep.Pclonos$Pclonos)>0
x <- rep.Pclonos$TotCount*runif(n,0.95, 1.05)
y <- rep.Pclonos$NbIGTs*runif(n,0.95, 1.05)
plot(x, y, ylim=c(1,60), log="xy", pch=16, xlab="Total nb of Occurences", ylab="Nb of IGT harboring the clonotype")
points(x[z], y[z], col=3, pch=16); text(2, 40, "TRAV11", col=3, adj=0)

#Count the IGTs in which same Nucleotide publics occur
table(IGT[is.element(all.Pclonos, sharedN.rep.Pclonos$Pclonos)])

#Different representations (all publics)
z <- is.element(all.Pclonos, rep.Pclonos$Pclonos)
sum(z)
length(unique(dt$IGT[z]))

## Extract MTB- and SFB-specific rep.Pclonos
x <- (regexpr("IGT47", rep.Pclonos$InIGT)>0)+ (regexpr("IGT61", rep.Pclonos$InIGT)>0) + (regexpr("IGT66", rep.Pclonos$InIGT)>0)>1
MTB.rep.Pclonos <- rep.Pclonos[x,]
y <- (as.numeric(regexpr("IGT37", rep.Pclonos$InIGT)>0) + as.numeric(regexpr("IGT39", rep.Pclonos$InIGT)>0) + as.numeric(regexpr("IGT90", rep.Pclonos$InIGT)>0))>1
SFB.rep.Pclonos <- rep.Pclonos[y,]
c.rep.Pclonos <- rep.Pclonos[x==F & y==F,]

# count of iNKT-related and non-repeats
z <- regexpr("TRAV11", c.rep.Pclonos$Pclonos)>0 & regexpr("ALGRLHF", c.rep.Pclonos$Pclonos)>0
w <- regexpr("TRAV1", c.rep.Pclonos$Pclonos)>0 & regexpr("NYQLIW", c.rep.Pclonos$Pclonos)>0
nktm.rep.Pclonos <- c.rep.Pclonos[z | w,]
c.rep.Pclonos <- c.rep.Pclonos[z==F & w==F,]

# Plot different classes of pclonos on basic NxNbIGT plot
par(mfrow=c(1,1))
n <- nrow(rep.Pclonos)
plot(rep.Pclonos$TotCount*runif(n,0.92, 1.08), rep.Pclonos$NbIGTs*runif(n,0.92, 1.08), ylim=c(1,10), log="xy", pch=".", xlab="Total nb of Occurences", ylab="Nb of IGT harboring the clonotype")
title(main="All Public Clonotypes")
n <- nrow(MTB.rep.Pclonos)
points(MTB.rep.Pclonos$TotCount*runif(n,0.92, 1.08), MTB.rep.Pclonos$NbIGTs*runif(n,0.92, 1.08), pch=16, col=3)
title(main="MTB or SFB-linked Public clonos")
n <- nrow(SFB.rep.Pclonos)
points(SFB.rep.Pclonos$TotCount*runif(n,0.92, 1.08), SFB.rep.Pclonos$NbIGTs*runif(n,0.92, 1.08), pch=16, col=5)
n <- nrow(nktm.rep.Pclonos)
points(nktm.rep.Pclonos$TotCount*runif(n,0.92, 1.08), nktm.rep.Pclonos$NbIGTs*runif(n,0.92, 1.08), pch=16, col=4)
title(main="iNKT or MAIT-linked Public clonos")
n <- nrow(c.rep.Pclonos)
points(c.rep.Pclonos$TotCount*runif(n,0.92, 1.08), c.rep.Pclonos$NbIGTs*runif(n,0.92, 1.08), pch=16, col=2)
title(main="Other Public clonos")

#organ distribution
organs <- sort(unique(dt$organ))
y <- rep(0, length(organs)); names(y) <- organs
x <- table(dt$organ[z])
w <- match(organs, names(x), nomatch=0)
y[w>0] <- x[w]
publicAB.organ.freq <- round(100*y/sum(y), digits=1)
y <- rep(0, length(organs)); names(y) <- organs
x <- table(dt$organ)
w <- match(organs, names(x), nomatch=0)
y[w>0] <- x[w]
all.organ.freq <- round(100*y/sum(y), digits=1)
plot(all.organ.freq, publicAB.organ.freq, pch=".")
text(all.organ.freq, publicAB.organ.freq, names(all.organ.freq), cex=0.7, adj=0)

# N region (take out "Possibles and No call)
z <- is.element(all.Pclonos, rep.Pclonos$Pclonos[rep.Pclonos$NbNucSqs==1])
w <- regexpr("ossible", dt$alpha.vj.nregion)<0 & regexpr("call", dt$alpha.vj.nregion)<0
x <- nchar(dt$alpha.vj.nregion[z&w])
y <- sample(nchar(dt$alpha.vj.nregion[w]), size=length(x))
wilcox.test(x, y)$p.value

w <- regexpr("ossible", dt$beta.vd.nregion)<0 & regexpr("call", dt$beta.vd.nregion)<0
summary(nchar(dt$beta.vd.nregion[z&w]))
summary(nchar(dt$beta.vd.nregion[w]))

w <- regexpr("ossible", dt$beta.dj.nregion)<0 & regexpr("call", dt$beta.dj.nregion)<0
summary(nchar(dt$beta.dj.nregion[z&w]))
summary(nchar(dt$beta.dj.nregion[w]))

##cell-type
z <- is.element(all.Pclonos, rep.Pclonos$Pclonos)
x <- table(annot$level1[is.element(annot$IGT.cellID, dt$IGT.cellID)])
y <- table(annot$level1[is.element(annot$IGT.cellID, dt$IGT.cellID[z])])
round(100*x/sum(x), digits=1)
round(100*y/sum(y), digits=1)

z <- is.element(all.Pclonos, rep.Pclonos$Pclonos[rep.Pclonos$NbNucSqs==1])
x <- table(annot$level1[is.element(annot$IGT.cellID, dt$IGT.cellID)])
y <- table(annot$level1[is.element(annot$IGT.cellID, dt$IGT.cellID[z])])
round(100*x/sum(x), digits=1)
round(100*y/sum(y), digits=1)

## Tabulate the sharing of clonotypes between IGTs (all or with single N)
z <- rep.Pclonos$InIGT
w <- unique(dt$IGT)
res <- matrix(0, length(w), length(w))
tot <- rep(0, length(w))
colnames(res) <- rownames(res) <- w
for (i in 1:length(w)) {
  x <- regexpr(paste(w[i], ":", sep=""), z)>0
  tot[i] <- sum(x)
  for (j in 1:length(w)) {
    if (j==i) {res[i,j] <- NA; next}
    y <- regexpr(paste(w[j], ":", sep=""), z)>0
    res[i,j] <- sum(x&y)
  }
}
summary(res)

## Where do different classes of public clonotypes map?
set <- c.rep.Pclonos$Pclonos
x <- is.element(all.Pclonos, set)
nm <- "Public Treg clonotypes"
y <- sample((1:nrow(map)), size=min(20000, nrow(map)))
z <- is.element(map$IGT.cellID, dt$IGT.cellID[x])
plot(map$X[y], map$Y[y], pch=".", col=16)
points(map$X[z], map$Y[z], pch=16, col=4)
title(main=nm)

## Add Pgen probabilities to the rep.Pclonos (from SummaryTable columns)
z <- match(rep.Pclonos$Pclonos, all.Pclonos, nomatch=0)
z1 <- match(dt$IGT.cellID[z], dt$IGT.cellID, nomatch=0)
z2 <- match(dt$IGT.cellID, dt$IGT.cellID, nomatch=0)

test <- data.frame(rep.Pclonos$Pclonos, all.Pclonos[z], dt$IGT.cellID[z], dt[z1,c("alpha.Pgen.nt","beta.Pgen.nt","alpha.Pgen.aa")], stringsAsFactors=F)

summary(dt$alpha.Pgen.nt[z1])
summary(dt$beta.Pgen.nt[z1])
summary(dt$alpha.Pgen.nt[z2], na.rm=T)
summary(dt$beta.Pgen.nt[z2][is.infinite(dt$beta.Pgen.nt[z2])==F], na.rm=T)

x <- dt$alpha.Pgen.nt[z1][is.na(dt$alpha.Pgen.nt[z1])==F]
x <- sample(x, min(500, length(x)))
y <- dt$alpha.Pgen.nt[z2][is.na(dt$alpha.Pgen.nt[z2])==F]
y <- sample(y, min(500, length(y)))
(wilcox.test(x,y))$p.value

x <- dt$beta.Pgen.nt[z1][is.na(dt$beta.Pgen.nt[z1])==F]
x <- sample(x, min(500, length(x)))
y <- dt$beta.Pgen.nt[z2][is.na(dt$beta.Pgen.nt[z2])==F & is.infinite(dt$beta.Pgen.nt[z2])==F]
y <- sample(y, min(500, length(y)))
(wilcox.test(x,y))$p.value

dev.off()
message("Plots saved to: ", file.path(script_dir, "11.Public Specificities.pdf"))
