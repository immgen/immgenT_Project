# Clonotype sharing between lineage clusters, rev for 0226 releases
# uses datntr, annot, cell.types, igts
# matching of N clonotypes is only within a given IGT.sample, so can't calculate frequencies
# only run for true clones (N match)

######### Load dataprep outputs #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

load(file.path(script_dir, "1. Dataprep and basics 0226.RData"))
# datntr, annot, cell.types, igts

# Use the "all productive"  dt subselection from datntr
z <- regexpr("TRBV", datntr$beta.vgene)>0 & regexpr("TRAV", datntr$alpha.vgene)>0 & datntr$alpha.functionality=="productive" & datntr$beta.functionality=="productive" #"true clonotypes with alpha and beta
dt <- datntr[z,]

L.clust.clone.sharingN.jac <- L.clust.clone.sharingN.incid <- vector("list", length=length(cell.types)-2)
for (k in 1:8) {
  print(paste(k, date()))
  X <- annot[annot$level1==cell.types[k],]	#trim to lineage
  w <- table(X$level2)
  if (k<8) z <- w>100 & regexpr("w", names(w))<0 #limit to clusters with reasonable chance of duplication and ditch miniverse and prolif
  if (k==8) z <- w>100  #keep all>100 for DPs
  w <- names(w)[z]
  w <- w[order(substring(w, regexpr("\\.",w)+1))]
  resN <- matrix(0, length(w), length(w))
  colnames(resN) <- rownames(resN) <- w
  resNincid <- resNu <- resN

  for (m in 1:length(igts)) {
    Y <- dt[dt$IGT==igts[m],]
    tcN <- paste(Y$alpha.vgene, Y$alpha.jgene, Y$alpha.junction.nt, Y$beta.vgene, Y$beta.jgene, Y$beta.junction.nt, sep=".")

    for (i in 1:length(w)) {
      x <- is.element(Y$IGT.cellID,X$IGT.cellID[X$level2==w[i]])
      for (j in 1:length(w)) {
        y <- is.element(Y$IGT.cellID,X$IGT.cellID[X$level2==w[j]])
        resN[i,j] <- resN[i,j]+length(intersect(tcN[x], tcN[y]))
        resNu[i,j] <- resNu[i,j]+length(union(tcN[x], tcN[y]))
        resNincid[i,j] <- resNincid[i,j]+ as.numeric(length(intersect(tcN[x], tcN[y]))>0)
      }
    }
  } #end m loop in igts

  for (i in 1:length(w)) {
    for (j in 1:length(w)) {
      resN[i,j] <- resN[i,j]/mean(c(diag(resN)[i], diag(resN)[j]))
    }
  }

  L.clust.clone.sharingN.jac[[k]] <- resN/resNu
  L.clust.clone.sharingN.incid[[k]] <- resNincid
} #end k loop through lineages

for (k in 1:6) {
  write.csv(L.clust.clone.sharingN.jac[[k]], file=file.path(script_dir, paste0("L.clust.clone.sharingN.jac.", cell.types[k], ".csv")), row.names=TRUE)
  write.csv(L.clust.clone.sharingN.incid[[k]], file=file.path(script_dir, paste0("L.clust.overlap.incidence.", cell.types[k], ".csv")), row.names=TRUE)
}



## Output as a matrix, shared clonos in rows, clusters in columns (run one lineage at a time
#Run across all, but select IGTs where duplicated, and only score sharing if in same IGT
k <- 1
X <- annot[annot$level1==cell.types[k],]	#trim to lineage
w <- table(X$level2)
z <- w>100 & regexpr("w", names(w))<0 #limit to clusters with reasonable chance of duplication and ditch miniverse and prolif
w <- names(w)[z]
w <- w[order(substring(w, regexpr("\\.",w)+1))]
Y <- dt[dt$level1==cell.types[k],]
tcN <- paste(Y$alpha.vgene, Y$alpha.jgene, Y$alpha.junction.nt, Y$beta.vgene, Y$beta.jgene, Y$beta.junction.nt, sep=".")
dups <- unique(tcN[duplicated(tcN)])

res <- matrix(0, length(dups), length(w)); shared <- rep(FALSE, length(dups)); inIGTs <- rep("-", length(dups))
colnames(res) <- w

for (i in 1:length(dups)) {
  z <- tcN==dups[i]
  if (length(unique(Y$level2[z]))<2) next
  if (length(unique(Y$level2[z]))>1) shared[i] <- TRUE
  y <- unique(Y$IGT[z])
  inIGTs[i] <- paste(y, collapse=";")
  for (j in 1:length(y)) {
    x <- table(Y$level2[z & Y$IGT==y[j]])
    if (length(x)<2) next
				if(length(x)>1) {
					idx <- match(names(x), w)
					ok <- !is.na(idx)
					res[i, idx[ok]] <- res[i, idx[ok]] + 1
				}
  } #end j loop through IGTs with that duplication
}	#end i loop through dups


## TCR sharing plots between CD4 clusters PAPER ? ######

#rescale the Jaccard matrix
x <- round(100*L.clust.clone.sharingN.jac[[1]], digits=2)
x[x==100] <- 0
x <- x[,colnames(x)!="CD4.P"]; x <- x[rownames(x)!="CD4.P",]

testQ <- data.frame(rep("CD4.Q", 10), rownames(x)[7:16], x[7:16,colnames(x)=="CD4.Q"], stringsAsFactors=F)
testG <- data.frame(rep("CD4.G", 10), rownames(x)[7:16], x[7:16,colnames(x)=="CD4.Q"], stringsAsFactors=F)

testAll <- data.frame(rep(colnames(x), each=ncol(x)),rep(colnames(x), times=ncol(x)),rep(0, (ncol(x))^2), stringsAsFactors=F)
for (i in 1:nrow(testAll)) testAll[i,3] <- as.numeric(x[rownames(x)==testAll[i,1], colnames(x)==testAll[i,2]])


for (i in 1:2) testAll[i,3] <- as.numeric(x[rownames(x)==testAll[i,1], colnames(x)==testAll[i,2]])

pdf(file.path(script_dir, "9d. cross level2 clonotype sharing within IGT.pdf"), width = 10, height = 10)

## Alluvial plots for CD4 paper, clonal sharing between clusters
set <- c("CD4.G", "CD4.H", "CD4.I")

#reset dt, X and Y on CD4 clusters G, H, I
z <- regexpr("TRBV", datntr$beta.vgene)>0 & regexpr("TRAV", datntr$alpha.vgene)>0 & datntr$alpha.functionality=="productive" & datntr$beta.functionality=="productive" #"true clonotypes with alpha and beta
z <- z & is.element(datntr$level2, set)
dt <- datntr[z,]

tcN <- paste(dt$alpha.vgene, dt$alpha.jgene, dt$alpha.junction.nt, dt$beta.vgene, dt$beta.jgene, dt$beta.junction.nt, sep=".")
x <- unique(tcN)

res <- matrix (0, length(x), length(set))
colnames(res) <- set
for (i in 1:length(x)) {
  y <- tcN==x[i]
  for (j in 1:length(set)) { res[i,j] <- sum(y & dt$level2==set[j]) }
}

CD4.GHI.clonesharing <- data.frame(x, res, stringsAsFactors=F)
colnames(CD4.GHI.clonesharing)[1] <- "clono"

#shorter set that only keep the more frequent sharings, butalso keep all the unique ones
z <- (rowSums(res[,-1]>0)>1 & rowSums(res[,-1])>9) | rowSums(res[,-1]>0)==1
restrim <- res[z,]

## Make scatter representation of repeat clonotypes in H and I
z <- CD4.GHI.clonesharing$CD4.H>0 & CD4.GHI.clonesharing$CD4.I >0
x <- CD4.GHI.clonesharing$CD4.H[z] / sum(dt$level2=="CD4.H")
y <- CD4.GHI.clonesharing$CD4.I[z] / sum(dt$level2=="CD4.H")
plot(x+0.0001, y+0.0001,log="xy", xlim=c(0.0001, 0.05)	, pch=16, ylim=c(0.0001, 0.05), xlab="Representation in CD4.H", ylab="Representation in CD4.I")
lines(c(0.0001, 0.05), c(0.0001, 0.05))
title(main="Shared TCR clones, CD4.H vs CD4.I")

test <- res[res[,2]>0 & res[,3]>0 , , drop=FALSE]
plot(1:2, c(test[1,1], test[1,2]), type="b", log="y", xlim=c(0.1,9), ylim=c(0, max(test)))
for (i in 2:nrow(test)) lines(1:2, c(test[i,1], test[i,2]))

dev.off()
message("CSV and PDF outputs saved in: ", script_dir)
