# 1. Basic dataprep and counting
# Performs harmonization, simplifies Va nomenclature, generates different data subsets (datc, datntr,etc)
#For 0226 data release, suspect duplicates (clonotype+barcode) cleaned out
#Annotation carried in the SummaryTable, but also keep using annot (compatibility with old code and stay current)

######### Input files #####
# All.IGTs.Enhanced.Summary.Table.2026-02-12.csv
# GSE297097_annotation_table_20260206_IGT1_104_cleaned.csv
# va.duplic.simplif.csv


######### Outputs #####
# dat
# datc
# datntr
# datbsl 
# datnt2
# igts
# Va
# Vb
# Ja
# Jb
# abP
# abN
# abP.ntr
# abN.ntr
# annot
# cell.types
# saved to: 1. Dataprep and basics 0226.RData


######### Load input data #####
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

dat <- read.csv(
  file.path(script_dir, "All.IGTs.Enhanced.Summary.Table.2026-02-12.csv"),
  stringsAsFactors = FALSE
)
annot <- read.csv(
  file.path(script_dir, "GSE297097_annotation_table_20260206_IGT1_104_cleaned.csv"),
  stringsAsFactors = FALSE
)
va.duplic.simplif <- read.csv(
  file.path(script_dir, "va.duplic.simplif.csv"),
  stringsAsFactors = FALSE
)
#edit column names to match previous code
colnames(dat)[colnames(dat)=="IGT.ID"] <- "IGT"
#correct N lowercase
x <- c("alpha.vj.nregion", "alpha2.vj.nregion", "beta.vd.nregion", "beta.dj.nregion" ,"beta2.vd.nregion", "beta2.dj.nregion")
for (i in 1:length(x)) {
	y <- colnames(dat)==x[i]
	z <- regexpr("Possible", dat[,y])<0 & regexpr("call", dat[,y])<0
	dat[z,y] <- toupper(dat[z,y])
}


cell.types <- c("CD4", "CD8", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

## Remove ambiuous and simplify the V region calls in dat (will carry down in all derived data tables, don't need in later scripts)
dat$beta.jgene[dat$beta.jgene=="TRBJ1-1  | TRBJ1-2"] <- "TRBJ1-1"; dat$beta.jgene[dat$beta.jgene=="TRBJ2-4  | TRBJ2-5"] <- "TRBJ2-5";
#drop the ambiguous Ja, Jb and Vb calls (small proportion)
z <- regexpr("\\|", dat$alpha.jgene)
dat$alpha.jgene[z>0] <- substring(dat$alpha.jgene[z>0],1,z[z>0]-1)
dat$alpha.jgene <- gsub(",", "", dat$alpha.jgene, fixed = TRUE)
dat$alpha.jgene <- gsub(" ", "", dat$alpha.jgene, fixed = TRUE); dat$alpha.jgene <- gsub(" ", "", dat$alpha.jgene, fixed = TRUE)	#do twice, otherwise leaves some
z <- regexpr("\\|", dat$beta.jgene)
dat$beta.jgene[z>0] <- substring(dat$beta.jgene[z>0],1,z[z>0]-1)
dat$beta.jgene <- gsub(",", "", dat$beta.jgene, fixed = TRUE)
dat$beta.jgene <- gsub(" ", "", dat$beta.jgene, fixed = TRUE); dat$beta.jgene <- gsub(" ", "", dat$beta.jgene, fixed = TRUE)
z <- regexpr("\\|", dat$beta.vgene)
dat$beta.vgene[z>0] <- substring(dat$beta.vgene[z>0],1,z[z>0]-1)
dat$beta.vgene <- gsub(",", "", dat$beta.vgene, fixed = TRUE)
dat$beta.vgene <- gsub(" ", "", dat$beta.vgene, fixed = TRUE); dat$beta.vgene <- gsub(" ", "", dat$beta.vgene, fixed = TRUE)
z <- regexpr("\\|", dat$alpha.vgene)
dat$alpha.vgene[z>0] <- substring(dat$alpha.vgene[z>0],1,z[z>0]-1)
dat$alpha.vgene <- gsub(",", "", dat$alpha.vgene, fixed = TRUE)
dat$alpha.vgene <- gsub(" ", "", dat$alpha.vgene, fixed = TRUE); dat$alpha.vgene <- gsub(" ", "", dat$alpha.vgene, fixed = TRUE)
#DEV: could have done more subtle for Va (many instances) but doesn't matter because next simplification will harmonize anyway)
#Condense the duplicated TRAV when possible (>95% aa match, usually >98% or complete), from curated va.duplic.simplif table
x <- match(dat$alpha.vgene, va.duplic.simplif$From, nomatch=0)
dat$alpha.vgene[x>0] <- va.duplic.simplif$To[x]

## Do the same in alpha2 and beta2
dat$beta2.jgene[dat$beta2.jgene=="TRBJ1-1  | TRBJ1-2"] <- "TRBJ1-1"; dat$beta2.jgene[dat$beta2.jgene=="TRBJ2-4  | TRBJ2-5"] <- "TRBJ2-5";
z <- regexpr("\\|", dat$alpha2.jgene)
dat$alpha2.jgene[z>0] <- substring(dat$alpha2.jgene[z>0],1,z[z>0]-1)
dat$alpha2.jgene <- gsub(",", "", dat$alpha2.jgene, fixed = TRUE)
dat$alpha2.jgene <- gsub(" ", "", dat$alpha2.jgene, fixed = TRUE); dat$alpha2.jgene <- gsub(" ", "", dat$alpha2.jgene, fixed = TRUE)	#do twice, otherwise leaves some
z <- regexpr("\\|", dat$beta2.jgene)
dat$beta2.jgene[z>0] <- substring(dat$beta2.jgene[z>0],1,z[z>0]-1)
dat$beta2.jgene <- gsub(",", "", dat$beta2.jgene, fixed = TRUE)
dat$beta2.jgene <- gsub(" ", "", dat$beta2.jgene, fixed = TRUE); dat$beta2.jgene <- gsub(" ", "", dat$beta2.jgene, fixed = TRUE)
z <- regexpr("\\|", dat$beta2.vgene)
dat$beta2.vgene[z>0] <- substring(dat$beta2.vgene[z>0],1,z[z>0]-1)
dat$beta2.vgene <- gsub(",", "", dat$beta2.vgene, fixed = TRUE)
dat$beta2.vgene <- gsub(" ", "", dat$beta2.vgene, fixed = TRUE); dat$beta2.vgene <- gsub(" ", "", dat$beta2.vgene, fixed = TRUE)
z <- regexpr("\\|", dat$alpha2.vgene)
dat$alpha2.vgene[z>0] <- substring(dat$alpha2.vgene[z>0],1,z[z>0]-1)
dat$alpha2.vgene <- gsub(",", "", dat$alpha2.vgene, fixed = TRUE)
dat$alpha2.vgene <- gsub(" ", "", dat$alpha2.vgene, fixed = TRUE); dat$alpha2.vgene <- gsub(" ", "", dat$alpha2.vgene, fixed = TRUE)
#DEV: could have done more subtle for Va (many instances) but doesn't matter because next simplification will harmonize anyway)#Condense the duplicated TRAV when possible (>95% aa match, usually >98% or complete), from curated va.duplic.simplif table
x <- match(dat$alpha2.vgene, va.duplic.simplif$From, nomatch=0)
dat$alpha2.vgene[x>0] <- va.duplic.simplif$To[x]

## Cleanup where beta2 is exactly the same as beta, or alpha2==alpha (528 and 329 occurences)
z <- dat$beta.vgene==dat$beta2.vgene & dat$beta.vgene!="" & dat$beta.junction.nt==dat$beta2.junction.nt
dat$beta2.vgene[z] <- ""; dat$beta2.jgene[z] <- ""; dat$beta2.junction[z] <- ""; dat$beta2.junction.nt[z] <- ""; dat$beta2.functionality[z] <- ""
z <- dat$alpha.vgene==dat$alpha2.vgene & dat$alpha.vgene!="" & dat$alpha.junction.nt==dat$alpha2.junction.nt
dat$alpha2.vgene[z] <- ""; dat$alpha2.jgene[z] <- ""; dat$alpha2.junction[z] <- ""; dat$alpha2.junction.nt[z] <- ""; dat$alpha2.functionality[z] <- ""

## Variants of the basic dat matrix
#datntr : take out all transgenics to avoid bias of some frequency analyses
z <- is.element(dat$IGT, c("IGT33", "IGT35","IGT36","IGT38", "IGT40", "IGT45","IGT46","IGT48","IGT55","IGT58", "IGT60", "IGT77","IGT95", "IGT96" ))
datntr <- dat[z==F,]

#datbsl : take out all transgenics, infected, tissue, thymus, tissues or "frontline" datasets, keep SLOs
z <- is.element(dat$IGT, c("IGT2", "IGT3", "IGT6", "IGT7", "IGT8", "IGT9", "IGT10", "IGT11", "IGT12", "IGT13", "IGT14", "IGT19", "IGT33", "IGT44", "IGT77", "IGT91", "IGT92", "IGT100", "IGT101", "IGT105", "IGT106", "IGT108", "IGT109")) #no infection IGTs
w <- dat$organ!="thymus"
x <- regexpr("colon", dat$organ)<0 & regexpr("intestine", dat$organ)<0 & regexpr("lung", dat$organ)<0 & regexpr("skin", dat$organ)<0 & regexpr("kidney", dat$organ)<0 & regexpr("liver", dat$organ)<0
datbsl <- dat[w&z&x,]
dim(datbsl)

#datnt2: no transgenic, no thymus
z <- datntr$organ=="thymus"
datnt2 <- datntr[z==F,]

#"datc", the "clean" data (periphery only, paired productive only, amplifications collapsed)
#Unlike previous, this one removes transgenics completely
datc <- dat
#remove thymus
datc <- datc[datc$organ!="thymus",]
#remove transgenic IGTs
z <- is.element(dat$IGT, c("IGT33", "IGT35","IGT36","IGT38", "IGT40", "IGT45","IGT46","IGT48","IGT55","IGT58", "IGT60", "IGT77","IGT95", "IGT96" ))
datc <- datc[z==F,]
#remove cells with only one chain
datc <- datc[datc$alpha.vgene!="" & datc$alpha.functionality=="productive" & datc$beta.vgene!="" & datc$beta.functionality=="productive",]
#remove pmel beta-only cells
z <- datc$beta.vgene=="TRBV14" & 	datc$beta.jgene=="TRBJ1-6"	& datc$beta.junction.nt=="TGTGCCAGCAGTTTCCACAGGGACTATAATTCGCCCCTCTACTTT"
datc <- datc[z==F,]
#remove amplified full clonotypes
x <- paste(datc$alpha.vgene, datc$alpha.jgene, datc$alpha.junction.nt, datc$beta.vgene, datc$beta.junction.nt, sep=".")
z <- duplicated(x)==F
datc <- datc[z,]

#datcb, as datc but removing the non-B6 IGTs (NOD and Foxp3.p327fs)
z <- is.element(datc$IGT, c("IGT24", "IGT49", "IGT74", "IGT99", "IGT104"))
datcb <- datc[z==F,]

#data with non productive rearangements (on either a or b chains)
datnpa <- dat[dat$alpha.functionality=="unproductive" & dat$alpha2.functionality!="productive",]	#to make sure that the unproductive chain doesn't hide a productive one
datnpb <- dat[dat$beta.functionality=="unproductive" & dat$beta2.functionality!="productive",]


## Define Va and Vb vectors of V regions
Va <- sort(unique(dat$alpha.vgene))
Va <- Va[is.element(Va, c("TRAV20", ""))==F]	#remove, never used
Vb <- sort(unique(dat$beta.vgene))
Vb <- Vb[Vb!="TRBV21" & Vb!="TRBV24" &Vb!=""]		#take out 21 and 24, pseudogenes
Ja <- unique(dat$alpha.jgene)
Ja <- Ja[order(as.numeric(substring(Ja,5)))]
Jb <- c("TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", "TRBJ1-6", "TRBJ2-1", "TRBJ2-2", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", "TRBJ2-7")	#drop 2.6, no instances

# igts vector
x <- unique(dat$IGT)
igts <- x[order(as.numeric(substring(x,4)))]

## Full protein or nucleic clonotypes in dat or variants
abP <- paste(dat$alpha.vgene, dat$alpha.junction, dat$alpha.jgene, dat$beta.vgene, dat$beta.junction, dat$beta.jgene, sep=".")
abN <- paste(dat$alpha.vgene, dat$alpha.junction, dat$alpha.jgene, dat$alpha.junction.nt, dat$beta.vgene, dat$beta.junction, dat$beta.jgene, dat$beta.junction.nt, sep=".")

abP.ntr <- paste(datntr$alpha.vgene, datntr$alpha.jgene, datntr$alpha.junction, datntr$beta.vgene, datntr$beta.jgene, datntr$beta.junction, sep=".")
abN.ntr <- paste(datntr$alpha.vgene, datntr$alpha.jgene, datntr$alpha.junction.nt, datntr$beta.vgene, datntr$beta.jgene, datntr$beta.junction.nt, sep=".")
aN.ntr <- paste(datntr$alpha.vgene, datntr$alpha.jgene, datntr$alpha.junction.nt, sep=".")
bN.ntr <- paste(datntr$beta.vgene, datntr$beta.jgene, datntr$beta.junction.nt, sep=".")


## Table of Tregs from datc for CONGA analysis
datc.Tregs <- datc$IGT.cellID[is.element(datc$IGT.cellID, annot$IGT.cellID[annot$level1=="Treg"])]
length(datc.Tregs)

		
## Table of mixed cells for an "overall conga" run	
datc.CD4 <- datc$IGT.cellID[is.element(datc$IGT.cellID, annot$IGT.cellID[annot$level1=="CD4"])]
datc.CD8 <- datc$IGT.cellID[is.element(datc$IGT.cellID, annot$IGT.cellID[annot$level1=="CD8"])]
datc.CD8aa <- datc$IGT.cellID[is.element(datc$IGT.cellID, annot$IGT.cellID[annot$level1=="CD8aa"])]
# 2026 annotation: legacy "nonconv" is now level1 == "Tz"
datc.nonconv <- datc$IGT.cellID[is.element(datc$IGT.cellID, annot$IGT.cellID[annot$level1=="Tz"])]

sample_cells <- function(ids, n) {
  if (length(ids) == 0) stop("Cannot sample cells: subset is empty")
  sample(ids, min(n, length(ids)), replace = length(ids) < n)
}

x <- c(
  sample_cells(datc.CD4, 10000),
  sample_cells(datc.CD8, 10000),
  sample_cells(datc.Tregs, 4000),
  sample_cells(datc.CD8aa, 2000),
  sample_cells(datc.nonconv, 2000)
)
y <- c(rep("CD4", 10000), rep("CD8", 10000),rep("Treg", 4000),rep("CD8aa", 2000),rep("nonconv", 2000))
datc.mix <- cbind(x,y)

######### Save outputs #####
output_objects <- c(
  "dat", "datc", "datntr", "datbsl", "datnt2",
  "igts", "Va", "Vb", "Ja", "Jb",
  "abP", "abN", "abP.ntr", "abN.ntr",
  "annot", "cell.types"
)
save(
  list = output_objects,
  file = file.path(script_dir, "1. Dataprep and basics 0226.RData")
)
