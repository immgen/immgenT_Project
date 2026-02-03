args = commandArgs(TRUE)
sc_IGT = args[1] 
genelist_file = args[2]


libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

sc_IGT = readRDS(sc_IGT)

sc_IGTM<- as.data.frame(cbind(sc_IGT$nCount_RNA,sc_IGT$nCount_ADT))
colnames(sc_IGTM) <- c('nCount_RNA', 'nCount_ADT')


genelist<-read.csv(genelist_file)

mnp.genelist<- genelist[c(1:27), 2]
granul.genelist<- genelist[c(28:54), 2]
B.genelist<- genelist[c(55:81), 2]
T.genelist<- genelist[c(82:108), 2]
ILC.genelist<- genelist[c(109:135), 2]
stroma.genelist<- genelist[c(136:162), 2]
mast.genelist<- genelist[c(163:189), 2]

subset_keep<-sc_IGT
DefaultAssay(subset_keep)<-"RNA"

subset_keep<-AddModuleScore_UCell(subset_keep, features=list(c(T.genelist)),  name = "T.enriched")
subset_keep<-AddModuleScore_UCell(subset_keep, features=list(c(mnp.genelist)), name = "mnp.enriched")
subset_keep<-AddModuleScore_UCell(subset_keep, features=list(c(B.genelist)),  name = "B.enriched")
subset_keep<-AddModuleScore_UCell(subset_keep, features=list(c(ILC.genelist)),  name = "ILC.enriched")

percentages<- rbind(subset_keep$signature_1T.enriched,
                    subset_keep$signature_1mnp.enriched,
                    subset_keep$signature_1B.enriched,
                    subset_keep$signature_1ILC.enriched)
row.names(percentages) <- c("T Percentages","MNP Percentages","B Percentages","ILC Percentages")

percentages<-as.data.frame(t(percentages))


p1_T<- ggplot()+ 
  geom_point(data=percentages ,aes(`T Percentages`, `B Percentages`),color='blue' )+ 
  xlim(0, max(percentages[,1])) + ylim(0, max(percentages[,3])) +theme(text = element_text(size=10))+
  ggtitle("ADT T cells")+ scale_y_continuous(limits = c(0, max(percentages[,3])),breaks = (seq(0,max(percentages[,3]), by = 0.01)))

p2_T<-ggplot()+ 
  geom_point(data=percentages ,aes(`T Percentages`, `ILC Percentages`),color='blue' )+
  xlim(0, max(percentages[,1])) + ylim(0, max(percentages[,3])) +theme(text = element_text(size=10))+
  ggtitle("ADT T cells")+ scale_y_continuous(limits = c(0, max(percentages[,3])),breaks = (seq(0,max(percentages[,3]), by = 0.01)))

p3_T<-ggplot()+ 
  geom_point(data=percentages ,aes(`T Percentages`, `MNP Percentages`),color='blue' )+
  xlim(0, max(percentages[,1])) +theme(text = element_text(size=10))+
  ggtitle("ADT T cells")+ scale_y_continuous(limits = c(0, max(percentages[,3])),breaks = (seq(0,max(percentages[,3]), by = 0.01)))

pdf(file = paste0('Tcell_cleanup:Step1-RNA_CellScoring.pdf'), width = 19, height = 20)
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(p1_T, vp = define_region(row = 1, col = 1))
print(p2_T, vp = define_region(row = 1, col = 2))
print(p3_T, vp = define_region(row = 1, col = 3))
ggplot(sc_IGTM, aes(nCount_RNA, nCount_ADT)) + geom_point(color='gray')+
scale_y_continuous(breaks = round(seq(min(sc_IGTM[2]),max(sc_IGTM[2]), by = 250),1))


dev.off()

saveRDS(subset_keep, file = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds")
