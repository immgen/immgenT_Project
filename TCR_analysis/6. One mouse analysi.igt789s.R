# Repeat analysis of single-mouse data, IGT7,8,9
#inherits dat, annot and others from 1.Dataprep

######### Load dataprep outputs #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

load(file.path(script_dir, "1. Dataprep and basics 0226.RData"))
# dat, annot

# lineage-specific MDE files not in repo; use annot mde coordinates filtered by level1
make_mde_map <- function(lineage) {
  z <- annot$level1 == lineage
  data.frame(
    IGT.cellID = annot$IGT.cellID[z],
    X = annot$mdeX[z],
    Y = annot$mdeY[z],
    stringsAsFactors = FALSE
  )
}
mde.incremental.CD4 <- make_mde_map("CD4")
mde.incremental.CD8 <- make_mde_map("CD8")
mde.incremental.Treg <- make_mde_map("Treg")
mde.incremental.CD4 <- mde.incremental.CD4[mde.incremental.CD4$X<4,]
mde.incremental.Treg <- mde.incremental.Treg[mde.incremental.Treg$X>(-5),]

pdf(file.path(script_dir, "6. One mouse analysi.igt789s.pdf"), width = 10, height = 10)

dat1 <- dat[is.element(dat$IGT, c("IGT7", "IGT8", "IGT9")),]
#make sample names uniform, remove "-", " "		#####CAUTION: different clean up for 100/101 #########
dat1$sample.name <- gsub("-", ".", dat1$sample.name, fixed = TRUE)
dat1$sample.name <- gsub(" ", "", dat1$sample.name, fixed = TRUE)
dat1$sample.name[dat1$sample.name=="skinback"] <- "skin.back"
dat1$sample.name[dat1$sample.name=="Skin.ear"] <- "skin.ear"
dat1$sample.name[dat1$sample.name=="SILP"] <- "SI.LP"
dat1$sample.name[dat1$sample.name=="Colon.IEL"] <- "colon.IEL"
#hardcode tissue order, here for 789
tiss.ord <- c("BoneMarrow", "LN.axillary", "LN.inguinal", "spleen", "blood", "skin.ear","skin.back", "kidney", "liver", "lung", "colon.IEL", "SI.IEL", "SI.LP", "colon.LP")

#Distribution of activated/resting cell-types in different organs
ActResInTissues789 <- matrix(0,2,length(tiss.ord)); colnames(ActResInTissues789) <- tiss.ord; rownames(ActResInTissues789) <- c("resting", "activated")

for (i in 1:length(tiss.ord)) {
  z <- is.element(annot$IGT.cellID, dat1$IGT.cellID[dat1$sample.name==tiss.ord[i]])
  ActResInTissues789[1,i] <- sum(annot$level2.group[z]=="resting")
  ActResInTissues789[2,i] <- sum(annot$level2.group[z]=="activated")
}

## Identify full nuc-level clonotype repeats
abN1 <- paste(dat1$alpha.vgene,  dat1$alpha.jgene,  dat1$alpha.junction.nt, dat1$beta.vgene, dat1$beta.jgene,  dat1$beta.junction.nt, sep=".")
ab1 <- paste(dat1$alpha.vgene,  dat1$alpha.jgene,  dat1$alpha.junction, dat1$beta.vgene, dat1$beta.jgene,  dat1$beta.junction, sep=".") #protein, when needed
x <- table(abN1[regexpr("TRAV", dat1$alpha.vgene)>0 & regexpr("TRBV", dat1$beta.vgene)>0] )
reps <- rev(sort(x[x>1]))


## make a dataframe of reps characteristics, where repeated clonotypes live in the mouse
N <- paste(dat$alpha.vgene,  dat$alpha.jgene,  dat$alpha.junction.nt, dat$beta.vgene, dat$beta.jgene,  dat$beta.junction.nt, sep=".")
N <- N[is.element(dat$IGT, c("IGT7", "IGT8", "IGT9"))==F]
P <- paste(dat$alpha.vgene,  dat$alpha.jgene,  dat$alpha.junction, dat$beta.vgene, dat$beta.jgene,  dat$beta.junction, sep=".")
P <- P[is.element(dat$IGT, c("IGT7", "IGT8", "IGT9"))==F]
x <- matrix(0, length(reps), length(tiss.ord)); colnames(x) <- tiss.ord
igt <- rep("-", length(reps))
othN <- othP <- nkt <- treg <- cd4 <- cd8 <- cd8aa <- rep(0, length(reps))

for (i in 1:length(reps)) {
  z <- abN1==names(reps)[i]
  igt[i] <- paste(substring(unique(dat1$IGT[z]), 4), collapse=";")
  nkt[i] <- sum(dat1$isNKT[z]=="TRUE")
  w <- table(annot$level1[is.element(annot$IGT.cellID, dat1$IGT.cellID[z])])
  treg[i] <- sum(w[names(w)=="Treg"], na.rm=TRUE)
  cd4[i] <- sum(w[names(w)=="CD4"], na.rm=TRUE)
  cd8[i] <- sum(w[names(w)=="CD8"], na.rm=TRUE)
  cd8aa[i] <- sum(w[names(w)=="CD8aa"], na.rm=TRUE)
  othN[i] <- sum(is.element(names(reps)[i], N))
  othP[i] <- sum(is.element(ab1[z][1], P))
  for (j in 1:length(tiss.ord)) { x[i,j] <- sum(abN1==names(reps)[i] & dat1$sample.name==tiss.ord[j])}
  print(paste(i, date()))
}

replics.789 <- data.frame(igt, as.numeric(reps), nkt, treg, cd4, cd8, cd8aa, othN, othP, x, stringsAsFactors=F)
rownames(replics.789) <- names(reps)

## Plot main replics.789 in the group.specific MDE space
cd4reps <- (1:nrow(replics.789))[replics.789$cd4> replics.789$treg & replics.789$cd4> replics.789$cd8]
cd8reps <- (1:nrow(replics.789))[replics.789$cd8> replics.789$treg & replics.789$cd8> replics.789$cd4]
tregreps <- (1:nrow(replics.789))[replics.789$treg> replics.789$cd4 & replics.789$treg> replics.789$cd8]

tss <- (1:ncol(replics.789))[is.element(colnames(replics.789), tiss.ord)]
for (j in (1:3)) {
  if (j==1) {set <- cd4reps; map <- mde.incremental.CD4; nm <- "CD4"}
  if (j==2) {set <- cd8reps; map <- mde.incremental.CD8; nm <- "CD8"}
  if (j==3) {set <- tregreps; map <- mde.incremental.Treg; nm <- "Treg"}
  par(mfrow=c(2,3), cex.main=0.6)
  y <- sample((1:nrow(map)), size=min(20000, nrow(map)))
  for (i in 1:min(6, length(set))) {
    x <- abN1==names(reps)[set[i]]
    mts <- colnames(replics.789)[tss][replics.789[set[i],tss]==max(replics.789[set[i],tss])]
    w <- paste(dat1$alpha.vgene[x][1], dat1$alpha.jgene[x][1], dat1$beta.vgene[x][1], dat1$beta.jgene[x][1], sep=".")
    z <- is.element(map$IGT.cellID, dat1$IGT.cellID[x])
    plot(map$X[y], map$Y[y], pch=".", col=16)
    points(map$X[z], map$Y[z], pch=16, col=i+1)
    title(main=paste(nm, w, "\nN=", sum(x), ", most in", mts))
  }
}

## More specific look at kidney CD8s (pool all)
set <- (1:nrow(replics.789))[replics.789$cd8> replics.789$treg & replics.789$cd8> replics.789$cd4 & apply(replics.789[,tss, drop=FALSE], 1, max)==replics.789$kidney]
map <- mde.incremental.CD8; nm <- "Kidney CD8 expansions"
x <- is.element(abN1, names(reps)[set])
y <- sample((1:nrow(map)), size=min(20000, nrow(map)))
z <- is.element(map$IGT.cellID, dat1$IGT.cellID[x])
par(mfrow=c(1,1))
plot(map$X[y], map$Y[y], pch=".", col=16)
points(map$X[z], map$Y[z], pch=16, col=3)
title(main=paste("Kidney CD8 expansions", "N=", sum(x)))


## Is the ear/back clonotype sharing reproduced in IGT5 (had both)
dat5 <- dat[is.element(dat$IGT, "IGT5"),]
table(dat5$sample.name)
abN5 <- paste(dat5$alpha.vgene,  dat5$alpha.jgene,  dat5$alpha.junction.nt, dat5$beta.vgene, dat5$beta.jgene,  dat5$beta.junction.nt, sep=".")
reps5 <- table(abN5[regexpr("TRAV", dat5$alpha.vgene)>0 & regexpr("TRBV", dat5$beta.vgene)>0] )
reps5 <- reps5[reps5>1]
reps5 <- rev(sort(reps5))
y <- unique(dat5$sample.name)
x <- matrix(0, length(reps5), length(y)); colnames(x) <- y
tr5 <- rep(FALSE, length(reps5))
for (i in 1:length(reps5)) {
  z <- abN5==names(reps5)[i]
  w <- is.element(dat1$IGT.cellID[z], annot$IGT.cellID[annot$level1=="Treg"])
  tr5[i] <- (sum(w)/sum(z))>0.2
  for (j in 1:length(y)) { x[i,j] <- sum(abN5==names(reps5)[i] & dat5$sample.name==y[j])}
}
replics5 <- data.frame(as.numeric(reps5),tr5, x, stringsAsFactors=F)
rownames(replics5) <- names(reps5)


## For Reviewer figure, permutation demonstration that tissue differences are not from sampling
X <- paste(dat1$alpha.vgene,  dat1$alpha.jgene,  dat1$alpha.junction.nt, dat1$beta.vgene, dat1$beta.jgene,  dat1$beta.junction.nt, sep=".")
zk <- dat1$sample.name=="kidney"
zc <- dat1$sample.name=="colon.LP"
z <- regexpr("TRAV", X)>0 & regexpr("TRBV", X)>0 & dat1$alpha.functionality=="productive" & dat1$beta.functionality=="productive"
X <- X[z]; zk <- zk[z]; zc <- zc[z]
x <- unique(X)
w <- match(X, x, nomatch=0)
dupc <- sum(table(w[zc])>1)
dupk <- sum(table(w[zk])>1)
kc <- length(intersect(w[zc], w[zk]))
iter <- 10
dupc.p <- dupk.p <- kc.p <- rep(0, length(iter))
for (i in 1:iter) {
  w <- sample(w, replace=FALSE)
  dupc.p[i] <- sum(table(w[zc])>1)
  dupk.p[i] <- sum(table(w[zk])>1)
  kc.p[i] <- length(intersect(w[zc], w[zk]))
}

dev.off()
message("Plots saved to: ", file.path(script_dir, "6. One mouse analysi.igt789s.pdf"))
