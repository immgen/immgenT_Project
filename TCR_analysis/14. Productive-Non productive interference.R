#14. Suppression by non-productive rearrangement
#Noticed by chance, following up on allelic inclusion preferences
#that a non-productive Vb strongly reduced the same Vb in the productive allele

######### Load dataprep outputs #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

load(file.path(script_dir, "1. Dataprep and basics 0226.RData"))
# dat, datntr, Va, Vb, annot, cell.types, ...

pdf(file.path(script_dir, "14. Productive-Non productive interference.pdf"), width = 10, height = 10)

vb.npppair.mat <- matrix(0, length(Vb), length(Vb))
colnames(vb.npppair.mat) <- rownames(vb.npppair.mat) <- Vb
n <- mtch <- rep(0, length(Vb))

w <- datntr[datntr$beta2.functionality!="productive" & datntr$beta2.vgene!="",]		#work from all the cells with a non productive beta2

for ( i in 1:length(Vb)) {
  x <- w$beta2.vgene==Vb[i]
  n[i] <- sum(x)
  for ( j in 1:length(Vb)) {
    y <- w$beta.vgene[x]==Vb[j]
    vb.npppair.mat[i,j] <- round(100*sum(y)/sum(x), digits=1)
  }
}
for (i in 1:length(Vb)) mtch[i] <- vb.npppair.mat[i,i]
pct <- round(100*n/sum(n),digits=1)
vb.npppair.mat[is.na(vb.npppair.mat)] <- 0
#p.value on the lower frequenies of matching in unproductive and productive, paired t.test on the observed "self" pairing pct vs the expected (colMeans)
testb <- t.test(mtch, colMeans(vb.npppair.mat, na.rm=TRUE), paired=TRUE)

vb.npppair.mat <- data.frame(n, pct, vb.npppair.mat, stringsAsFactors=F)

par(mfrow=c(1,1))
image(seq_len(nrow(vb.npppair.mat)), seq_len(ncol(vb.npppair.mat)-2), as.matrix(vb.npppair.mat[,-(1:2)]), axes=FALSE)
text(0, 1:length(Vb), Vb, adj=1, cex=0.7)
text(1:length(Vb), 0, Vb, adj=1, cex=0.7, srt=90)

#calculate a normalized version of this (divided by the mean usage of each TRBV)
vb.npppair.mat.norm <- as.matrix(vb.npppair.mat[,-c(1:2)])
for (i in 1:length(Vb)) vb.npppair.mat.norm[,i] <- vb.npppair.mat.norm[,i]/colMeans(vb.npppair.mat.norm, na.rm=TRUE)[i]
vb.npppair.mat.norm[is.na(vb.npppair.mat.norm) | is.infinite(vb.npppair.mat.norm)] <- 0
vb.npppair.mat.norm <- round(log2(vb.npppair.mat.norm + 1e-6), digits=2)

#recompute p.value on the normalized matrix
testbn <- t.test(diag(vb.npppair.mat.norm), colMeans(vb.npppair.mat.norm, na.rm=TRUE), paired=TRUE)


###Same for alpha locus
va.npppair.mat <- matrix(0, length(Va), length(Va))
colnames(va.npppair.mat) <- rownames(va.npppair.mat) <- Va
n <- mtch <- rep(0, length(Va))

w <- datntr[datntr$alpha2.functionality!="productive" & datntr$alpha2.vgene!="",]		#work from all the cells with a non productive alpha2
for ( i in 1:length(Va)) {
  x <- w$alpha2.vgene==Va[i]
  n[i] <- sum(x)
  for ( j in 1:length(Va)) {
    y <- w$alpha.vgene[x]==Va[j]
    va.npppair.mat[i,j] <- round(100*sum(y)/sum(x), digits=1)
  }
}
pct <- round(100*n/sum(n),digits=1)
for (i in 1:length(Va)) mtch[i] <- va.npppair.mat[i,i]
va.npppair.mat[is.na(va.npppair.mat)] <- 0
testa <- t.test(mtch, colMeans(va.npppair.mat, na.rm=TRUE), paired=TRUE)

va.npppair.mat <- data.frame(n, pct, va.npppair.mat, stringsAsFactors=F)

image(seq_len(nrow(va.npppair.mat)), seq_len(ncol(va.npppair.mat)-2), as.matrix(va.npppair.mat[,-(1:2)]), axes=FALSE)
text(0, 1:length(Va), Va, adj=1, cex=0.7)
text(1:length(Va), 0, Va, adj=1, cex=0.7, srt=90)

#trim to take only the alphas with non-productive N>1000
z <- va.npppair.mat$n>800 & (colMeans(va.npppair.mat[,-(1:2)]))>1
va.npppair.mat <- va.npppair.mat[z,c(1,2,(3:ncol(va.npppair.mat))[z])]

#calculate a normalized version of this (divided by the mean usage of each TRAV)
va.npppair.mat.norm <- as.matrix(va.npppair.mat[,-c(1:2)]) +0.02
va.npppair.mat.norm <- t(t(va.npppair.mat.norm)/colMeans(va.npppair.mat.norm))
va.npppair.mat.norm <- round(log2(va.npppair.mat.norm), digits=2)



### test in thymic DPs
vb.npppair.mat.dp7172 <- matrix(0, length(Vb), length(Vb))
colnames(vb.npppair.mat.dp7172) <- rownames(vb.npppair.mat.dp7172) <- Vb
n <- rep(0, length(Vb))
w <- datntr[datntr$beta2.functionality!="productive"  & datntr$beta2.vgene!="" & is.element(datntr$IGT, c("IGT71", "IGT72")),]		#thymic DPs, d3 and adult
for ( i in 1:length(Vb)) {
  x <- w$beta2.vgene==Vb[i]
  n[i] <- sum(x)
  for ( j in 1:length(Vb)) {
    y <- w$beta.vgene[x]==Vb[j]
    vb.npppair.mat.dp7172[i,j] <- round(100*sum(y)/sum(x), digits=1)
  }
}
pct <- round(100*n/sum(n),digits=1)
vb.npppair.mat.dp7172 <- data.frame(n, pct, vb.npppair.mat.dp7172, stringsAsFactors=F)
# in d3 DPs only
vb.npppair.mat.dp7172.d3 <- matrix(0, length(Vb), length(Vb))
colnames(vb.npppair.mat.dp7172.d3) <- rownames(vb.npppair.mat.dp7172.d3) <- Vb
n <- rep(0, length(Vb))
w <- datntr[datntr$beta2.functionality!="productive" & datntr$beta2.vgene!="" & is.element(datntr$IGT, c("IGT71", "IGT72")) & is.element(datntr$sample.name, c("3d DP Thymocytes F", "3d DP Thymocytes M")),]		#thymic DPs, d3
for ( i in 1:length(Vb)) {
  x <- w$beta2.vgene==Vb[i]
  n[i] <- sum(x)
  for ( j in 1:length(Vb)) {
    y <- w$beta.vgene[x]==Vb[j]
    vb.npppair.mat.dp7172.d3[i,j] <- round(100*sum(y)/sum(x), digits=1)
  }
}
pct <- round(100*n/sum(n),digits=1)
vb.npppair.mat.dp7172.d3 <- data.frame(n, pct, vb.npppair.mat.dp7172.d3, stringsAsFactors=F)

## Frequency of non productive rearrangements
#Seems pretty different from mature/selected
#all cells with p and np pairs
plot(vb.npppair.mat$pct, colMeans(vb.npppair.mat[,-(1:2)]), pch=".", xlab="Vb frequency in non-productive", ylab="Vb frequency in matched productive")
text(vb.npppair.mat$pct+0.1, colMeans(vb.npppair.mat[,-(1:2)]), Vb, adj=0, cex=0.7)
#np vs overall productive
y <- rep(0, length(Vb))
for (i in 1:length(Vb)) y[i] <- sum(datntr$beta.vgene==Vb[i])
y <- round(100*y/sum(y), digits=1)
plot(vb.npppair.mat$pct, y, pch=".", xlab="Vb frequency in non-productive", ylab="Vb frequency in productive TCRs overall")
text(vb.npppair.mat$pct+0.1, y, Vb, adj=0, cex=0.7)
# in DPs
y <- rep(0, length(Vb))
for (i in 1:length(Vb)) y[i] <- sum(datntr$beta.vgene==Vb[i] & is.element(datntr$IGT, c("IGT71", "IGT72")) )
y <- round(100*y/sum(y), digits=1)
plot(vb.npppair.mat.dp7172$pct, y, pch=".", xlab="Vb frequency in DP non-productive", ylab="Vb frequency in DP productive", xlim=c(0,14), ylim=c(0,14))
text(vb.npppair.mat.dp7172$pct+0.1, y, Vb, adj=0, cex=0.7)

dev.off()
message("Plots saved to: ", file.path(script_dir, "14. Productive-Non productive interference.pdf"))
