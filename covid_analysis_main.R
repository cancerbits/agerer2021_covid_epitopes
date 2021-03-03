#Single-cell analysis for Agerer et al, 2021
#Author: Luis Fernando Montaño-Gutierrez
#Affiliation: Developmental Cancer genomics Group, St. Anna Children's Cancer Research Insitute, Vienna, Austria.

#Importing libraries
library(Seurat)
library(patchwork)
library(devtools)
library(SeuratObject)
library(ggrepel)
library(AUCell)
library(RColorBrewer)
library(simpleCache)
library(wesanderson)
attach(wes_palettes)

###############################
#declaring recurring variables
###############################
chd="./rcache"
##############################
# declaring relevant gene sets
##############################
simpleCache("genesets", {
  exhaustion=c("TOX"     ,"HAVCR2"   ,"LAG3"     ,"ENTPD1"   ,"PDCD1"    ,"CTLA4"    ,"CD38"     ,"TIGIT"    ,
               "VCAM1"    ,"CD27"     ,"SNAP47"   ,"IGFLR1"   ,"RAD51"    ,"CCNB1"    ,"BUB1"     ,
               "SIRPG"    ,"SEMA4A"   ,"CXCR6"    ,"FUT8"     ,"HLA-DMA"  ,"ITGAE"    ,"UBE2F"    ,
               "NDFIP2"   ,"CD63"     ,"FKBP1A"   ,"TPI1"     ,"CDCA8"    ,"NCAPG2"   ,"CDKN3"    ,"CCL4L2"   ,
               "RGS2"     ,"NAB1"     ,"ID3"      ,
               "CCR5"     ,"GOLIM4"   ,"ACP5"     ,"HLA-DRA"  ,"FCRL3"    ,"OSBPL3"   ,"ICOS"     ,"FAM3C"    ,
               "PTPN11"   ,"CKS2"     ,"GALM"     ,"SNX9"     ,"IRF4"     ,"STMN1"    ,"PRDM1"    ,"CD2BP2"   ,
               "RAB27A"   ,"DUSP4"    ,"PHLDA1"   ,"ITM2A"    ,"IFI35"    ,"ISG15"    ,"STAT3"    ,"WARS"     ,
               "SYNGR2"   ,"GBP2"     ,"LYST"     ,"BST2"     ,"PARK7")
  
  cytotoxic=c("CTSW","GNLY","GZMA","GZMB","GZMH","IFNG","KLRB1","KLRD1","KLRK1","NKG7","PRF1")
  
  ifn=c("IFI6"    ,"MX1"       ,"ISG20"     ,"ISG15"     ,"LY6E"      ,"IFIT3"     ,"MX2"       ,"XAF1"      ,
        "OAS3"      ,"EIF2AK2"   ,"SAMD9L"    ,"EPSTI1"    ,"IFITM1"    ,"OAS1"      ,"SAMD9"     ,
        "IFIT1"     ,"RNF213"    ,"USP18"     ,"IFI44"     ,"SP100"     ,"IFI35"     ,"BST2"      ,"PLSCR1"    ,"TRIM22"    ,
        "RSAD2"     ,"HERC6"     ,"SP110"     ,"IRF7"      ,"OAS2"      ,"PARP14"    ,"CMPK2"     ,"PARP9"     ,"GBP1"      ,"DDX60"     ,
        "PSMB9"     ,"IFIH1"     ,"DDX58"     ,"NUB1"      ,"IFIT2"     ,"ADAR"      ,"PSME2"     ,"LAMP3"     ,"TNFSF10"   ,"NMI"       ,
        "UBE2L6"    ,"APOL6"     ,"PSME1"     ,"LAP3"      ,"GBP4"      ,"STAT2"     ,"B2M"       ,"IRF1"      ,"IFITM3"    ,
        "PSMB8"     ,"PNPT1"     ,"HLA-C"     ,"MT2A"      ,"GBP5"      ,"IRF9"      ,"IFITM2"    ,"SAMHD1"    ,"HLA-A"     ,"TAP1"      ,"TRIM25"    ,
        "RTP4"      ,"TMEM140"   ,"PML"       ,"TRIM38"    ,"TRIM21"    ,"ZBP1"      ,"DHX58"     ,"TAPBP"     ,"PARP12"    ,"LGALS3BP"  ,
        "GBP3"      ,"ELF1"      ,"ZNFX1"     ,"GBP2"      ,"IRF2"      ,"CD38"      ,"RBCK1"     ,"TRIM14"    ,"MYD88"     ,"NCOA3"     ,
        "CNP"       ,"TRAFD1"    ,"LYSMD2"    ,"MOV10"     ,"CASP1"     ,"CASP4"     ,"HLA-F"     ,"OGFR"      ,"CD274"     ,"SOCS1"     ,
        "NLRC5"     ,"JAK1"      ,"IFI27"     ,"CD74"      ,"GCH1"      ,"IL15RA"    ,"PSMA3"     ,"TDRD7"     ,"TXNIP"     ,"TRIM5"     ,
        "TRIM26"    ,"CASP7"     ,"CD47"      ,"IFI30"     ,"PTPN1"     ,"JAK2"      ,"PSMB10"    ,"HLA-E"     ,"IL15"      ,"MVP"       ,
        "FAS"       ,"HIF1A"     ,"VAMP5"     ,"WARS"      ,"CFH"       ,"RIPK2"     ,"STAT3"     ,"CASP8"     ,"EIF4E3"    ,"VAMP8"     ,
        "UBA7"      ,"IL10RA"    ,"PIAS1"     ,"IRF3"      ,"NCOA7"     ,"NOD1"      ,"CAMK2G"    ,"ARL4A"     ,"SPPL2A"    ,"RNF31"     ,
        "IRF5"      ,"SLAMF7"    ,"HLA-B"     ,"ITGB7"     ,"FAM46A"    ,"TOR1B"     ,"SELL"      ,"CSF1"      ,"TRIM35"    ,"CD69"      ,
        "SOCS3"     ,"PIM1"      ,"RAPGEF6"   ,"IFNAR1"    ,"RNASEL"    ,"RIPK1"     ,"SOD2"      ,"IL18BP"    ,"PTPN2"     ,"IL2RB"     ,"PRKCD"     ,
        "HLA-G"     ,"PDE4B"     ,"AUTS2"     ,"ST3GAL5"   ,"PTPN6"     ,"TRIM68"    ,"CDKN1A"    ,"IP6K2"     ,"LPAR6"     ,"CASP3"     ,"BPGM"      ,"ICAM1"     ,
        "ISOC1"     ,"NUP93"     ,"PSMB2"     ,"GZMA"      ,"ARID5B"    ,"GPR18"     ,"LCP2"      ,"HLA-DMA"   ,"TRIM8"     ,"TYK2"      ,"PFKP"      ,
        "CCL5"      ,"PNP"       ,"ABCE1"     ,"CAMK2D"    ,"IL4R"      ,"STAT4"     ,"PELI1"     ,"IFNG"      ,"PTPN11"    ,"SRI"       ,"IFNGR1"    ,
        "SUMO1"     ,"FGL2"      ,"HLA-DPA1"  ,"HLA-DQA1"  ,"IRF8"      ,"IFNAR2"    ,"HLA-DRA"   ,"TNFAIP3"   ,"UPP1"      ,"HLA-DRB5"  ,"NFKB1"     ,"XCL1"      ,
        "IRF4"      ,"HLA-DPB1"  ,"CD44"      ,"HLA-DQB1"  ,"HLA-DRB1"  ,"ST8SIA4"   ,"MTHFD2"    ,"NFKBIA"    ,"BTG1"      ,"NAMPT"     ,
        "PROCR"     ,"HLA-DQA2"  ,"TRIM34"    ,"VCAM1"     ,"CD40"      ,"CMKLR1-KLRK1")
  
  viral=c("XCL1","XCL2","TNFRSF9","CRTAM","ATP1B3","TAGAP","SLA","NR4A2","IFNG","TNFRSF1B","PTPN6",
          "CCT6A","CCT3","GHITM","EGR2","MAP2K3","PSMA4","RBPJ","HSPA9","FABP5","PAM","SNRPB","MTHFD2",
          "ILF2","RBM3","SHMT2","CD82","HSPD1","SERBP1","CCND2","SEC61B","POMP","CYCS","PGAM1","CCT2",
          "IMPDH2","VDAC1","NOL7","DDX21","NAMPT","BTG3","RPA3","FDPS","MIR155HG","NUTF2","SNRPE","TUBA1B",
          "ZBED2","NHP2","TOMM5")
  
  unhelped=c("LRRC61","KNTC1","ARSB","KLF12","CMKLR1","RBM38","SH3TC1","CLSPN","PRR5L","KIAA1671","STK32C","ZMIZ1",
             "IFITM1","IFITM2","IFITM3","PRR5","CKAP4","STMN1","ITGAM","CCNF","CHN2","ENTPD1","VKORC1L1","P2RX7",
             "OSBPL3","UBE2C","SPECC1","ZEB2","ST3GAL1","KLRG1","SYTL2","KLRB1","CORO2A","GPR55","S1PR5","EMILIN2",
             "RUNX1","HAVCR2","ANXA1","HECTD2","ARHGAP11A","ARHGAP11B","PHLDB2","NCAPG2","ADAM8","SLC43A3","SUFU",
             "CX3CR1","BHLHE40","WIPI1","KLRC4","KLRC4-KLRK1","KLRC3","KLRC2","KLRC1","ESPL1","MICAL3","IRF4","L1CAM",
             "SEMA4A","CCNB1","TSPAN32","TIAM1","BUB1","MYADM","IL18RAP","IL9R","LDLR","KIAA0513","NFIC","APOBR","BMPR2",
             "CDC6","GZMA")
  
  genesets=list(exhaustion= exhaustion, cytotoxic=cytotoxic, unhelped=unhelped, viral=viral, ifn=ifn)
  
}, recreate=TRUE, cacheDir=chd)


#####################################
#declaring recurring atomic functions
#####################################
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))
getvar=function(x, var) x[[var]][[1]]
getmetavars=function(x) colnames(x@meta.data)
getcats=function(x, var) unique(x@meta.data[[var]])
wrtnonames<-function(x, filename) {write.table(x, file=filename, row.names=F, col.names=F, quote=F)}
grepvec<-function(vec,ifpos, ifneg=FALSE) sapply(vec, function(x) if(grepl(ifpos, x)){return(ifpos) }else{return(ifneg)}, USE.NAMES=F)
genotype= function(x){ if ((grepl("wt", x))){
  y="WT"
}else{y="MUT"}
  y}
response= function(x){ if ((grepl("pos", x))){
  y="pos"
}else{y="neg"}
  y}


wtpos="WT YLQ+"
wtneg="WT YLQ-"
mutpos="MUT YLQ+"
mutneg= "MUT YLQ-"
sublabels<-function(x){
  
  y=gsub('HTO-WT_pos', wtpos, x=x)
  y=gsub('HTO-WT_neg', wtneg, x=y)
  y=gsub('HTO-MUT_pos', mutpos, x=y)
  y=gsub('HTO-MUT_neg', mutneg, x=y)
  y
}

###############################
###data preparation and loading
###############################
#1. Download the GEO directory to the folder where all the GIT scripts are located. 
#2. Rename the GEO directory as GEOupload if it is not already called as such.
#3. run the prepdirs.sh script. this will create the following directories:
#sars042- with files ready for analysis 
#sars060- with files ready for analysis
#plots- tom dump all figures

MY_GIT_DIR<-"type in the path to the GIT folder"
setwd(MY_GIT_DIR)

#actual loading of the data.
BSFdemux=TRUE
if (BSFdemux==True){

  sars042dir="./sars042"
  sars060dir="./sars060"
  

  simpleCache("loadseqsBSF", {
  sars042.umi=Read10X(sars042dir, gene.column=2)
  sars042.umi<-sars042.umi$`Gene Expression`
  sars060.umi=Read10X(sars060dir, gene.column=2)
  sars060.umi<-sars060.umi$`Gene Expression`
  
#importing Hashtag oligo barcoding and initial quality control metadata
qcpath042<-paste0(sars042dir, "/QC042rep.csv")
qcpath060<-paste0(sars060dir, "/QC060rep.csv")

sars042.BSF<-read.delim(qcpath042, sep=",", row.names=1)   
sars042.hashtag=CreateSeuratObject(counts = sars042.umi, meta.data=sars042.BSF)


sars060.BSF<-read.delim(qcpath060, sep=",", row.names=1)   
sars060.hashtag=CreateSeuratObject(counts = sars060.umi, meta.data=sars060.BSF)

#############################################################################################
#Appending both patient datasets and basic annotation
#############################################################################################

covid.hashtag <- merge(x= sars042.hashtag, y= sars060.hashtag,  add.cell.ids = c("SARS042", "SARS060"), project= "covid")

}, cacheDir="./rcache", reload=T, assignToVar="covid.hashtags")
}



simpleCache("covidhashtags3", {
  #loading previous block
  simpleCache("covid.hashtags", cacheDir="./rcache", reload=T, assignToVar="covid.hashtags")
   
  mtcutoff=10

labs=covid.hashtag[["hto_demux"]][[1]]
labs2=covid.hashtag[["sample"]][[1]]

covid.hashtag[["genotype"]]<-unname(sapply(labs, genotype))
covid.hashtag[["response"]]<-unname(sapply(labs, response))
covid.hashtag[["percent.mt"]] <- PercentageFeatureSet(covid.hashtag, pattern = "^MT-")
covid.hashtag[["percent.trv"]] <- PercentageFeatureSet(covid.hashtag, pattern = "^TRV")
covid.hashtag[["percent.trbv"]] <- PercentageFeatureSet(covid.hashtag, pattern = "^TRBV")
covid.hashtag[["patientlabel"]]<- grepvec(labs2, "SARS042", "SARS060")
covid.hashtag[["htolabels"]]<- sapply(labs, sublabels, USE.NAMES=F)

###############################
#Standard quality control plots
###############################

#distribution of genes per barcode
qc1=ggplot(covid.hashtag@meta.data, aes(x=nCount_RNA))+
  geom_histogram(binwidth=20, alpha=0.5)+
  geom_histogram(data=covid.hashtag@meta.data[covid.hashtag@meta.data[, "hto_demux"]=="Doublet", ], fill="red", alpha=0.4, binwidth=50)+
  labs(x = "name1",y = "Density")+
  guides(fill=guide_legend())+
  coord_cartesian(xlim=c(0, 15000))+
  labs(title="Reads per barcode: singlet vs doublet")
qc1

#distribution of counts per barcode
qc2=ggplot(covid.hashtag@meta.data, aes(x=nFeature_RNA))+
  geom_histogram(binwidth=20)+
  geom_histogram(data=covid.hashtag@meta.data[covid.hashtag@meta.data[, "hto_demux"]=="Doublet", ], fill="red", alpha=0.4, binwidth=50)+
  labs(title="Genes per barcode: all vs doublet")
  coord_cartesian(xlim=c(0, 1000));qc2



####################
#applying QC filters
####################

#counting cells before filtering  
numcells= function(sc) length(sc[["nFeature_RNA"]][[1]])
cat(sprintf("number of cells before filtering: %d\n", numcells(covid.hashtag)))
  
covid.hashtag <- subset(covid.hashtag, subset =  percent.mt < mtcutoff & pass_QC=="True")
cat(sprintf("applying filter...\nnumber of cells after filtering: %d\n", numcells(covid.hashtag)))

################################################################
#preliminary PCA exploration of the data for further QC curation
################################################################

#normalize count depth
covid.hashtag <- NormalizeData(covid.hashtag)
# Find and scale variable features
nfeatures=1500
covid.hashtag <- FindVariableFeatures(covid.hashtag, selection.method = "vst", nfeatures=nfeatures)
covid.hashtag <- ScaleData(covid.hashtag, features = VariableFeatures(covid.hashtag))
# Run PCA
covid.hashtag <- RunPCA(covid.hashtag, features = VariableFeatures(covid.hashtag))
covid.hashtag=ProjectDim(covid.hashtag, reduction = "pca")


dms=10 

covid.hashtag <- FindNeighbors(covid.hashtag, reduction = "pca", dims = 1:dms )
covid.hashtag <- FindClusters(covid.hashtag, resolution = 0.6, verbose = FALSE)
covid.hashtag <- RunTSNE(covid.hashtag, reduction = "pca", dims = 1:dms)


#visualising principal components to see if any specific cluster is dominating the variance
ppc1=FeaturePlot(object = subset(covid.hashtag, subset = percent.mt<10), features = "PC_1", cols= c("grey", "blue"), reduction = "tsne");

#filtering  seurat clusters based on our findings after inspection
covid.hashtag[["QC_clusters"]]=sapply(covid.hashtag[["seurat_clusters"]], function(x) paste0("QC", x))



qcnames<-sapply(0:13, function(x) paste0("QC", x))
qccols<-mycolors[1:14]
getqccol= function(x) qccols[qcnames==x]
covid.hashtag[["qc_cols"]]= sapply(covid.hashtag@meta.data[,"QC_clusters"], getqccol, USE.NAMES=F)
Idents(covid.hashtag)<-"QC_clusters"
covid.hashtag@active.ident <- factor(x = covid.hashtag@active.ident, levels = qcnames)



tsne1=DimPlot(covid.hashtag, group.by = "QC_clusters", reduction="tsne", label=T, label.size=6, label.color="black")#+scale_color_manual( values=qccols);tsne1 


qc3 <- FeatureScatter(covid.hashtag, feature1 = "nFeature_RNA", feature2 = "percent.mt")#+scale_color_manual(values=qccols)
qc4 <- FeatureScatter(covid.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")#+scale_color_manual(values=qccols)
qc5 <- FeatureScatter(covid.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")#+scale_color_manual(values=qccols)  

covid.hashtag <- subset(covid.hashtag, idents = c("QC8", "QC12"), invert=TRUE)

#after removing the bad clusters, we re run the PCA and embeddings as those previous clusters may have been biasing the variance.
tsne2<-DimPlot(covid.hashtag, group.by = "htolabels", reduction="tsne");

Idents(covid.hashtag)="QC_clusters"
covid.hashtag@active.ident <- factor(x = covid.hashtag@active.ident, levels = qcnames)


qc6 <- FeatureScatter(covid.hashtag, feature1 = "nFeature_RNA", feature2 = "percent.mt")#+scale_color_manual( values=qccols[!grepl("QC8|QC12", qcnames)])
qc7 <- FeatureScatter(covid.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")#+scale_color_manual( values=qccols[!grepl("QC8|QC12", qcnames)])
qc8 <- FeatureScatter(covid.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")#  +scale_color_manual(values=qccols[!grepl("QC8|QC12", qcnames)])

####################################################
###Plotting the Quality control supplemental figure.
###################################################

qc_assembly=(tsne1+ppc1)/(qc3+qc4+qc5)/(qc6+qc7+qc8)

tsnecl=DimPlot(covid.hashtag, group.by = "hto_demux", reduction="umap");
tsnecl

covid.hashtag
}, reload=T, cacheDir=chd, assignToVariable="covid.hashtag3")


pdf("./plots/qc_outliercluster.pdf", width=12, height=12)
qc_assembly
dev.off()

###################################################
#Batch correction of Patient datasets (dataset integration) by Seurat Anchoring
##################################################
mergedims=30
simpleCache("covidmerged", {
  
covid.list<- SplitObject(covid.hashtag, split.by = "patientlabel")
for (i in 1:length(covid.list)) {
  covid.list[[i]] <- NormalizeData(covid.list[[i]], verbose = FALSE)
  covid.list[[i]] <- FindVariableFeatures(covid.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

covid.anchors <- FindIntegrationAnchors(object.list = covid.list, dims = 1:50)
covid.merged<- IntegrateData(anchorset = covid.anchors, dims = 1:30)

}, cacheDir=chd, reload=T, assignToVariable="covid.merged" )

DefaultAssay(covid.merged)<-"integrated"

################################################
#Run PCA, TSNE and UMAP for integrated dataset.
################################################
covid.merged <- ScaleData(covid.merged, features = VariableFeatures(covid.merged))
covid.merged <- RunPCA(covid.merged, features = VariableFeatures(covid.merged))
covid.merged <-ProjectDim(covid.merged, reduction = "pca")
covid.merged <- FindNeighbors(covid.merged, reduction = "pca", dims = 1:dms )
covid.merged <- FindClusters(covid.merged, resolution = 0.6, verbose = FALSE)
covid.merged <- RunTSNE(covid.merged, reduction = "pca", dims = 1:dms)
covid.merged <- RunUMAP(covid.merged, reduction = "pca", dims = 1:dms)


##############################################################################################
#Replacing covid.hashtag with the contents of covid.merged integrated dataset, inverting the first UMAP dimension
##############################################################################################
simpleCache("covidAUCmerged_invumap",{
  simpleCache("covidAUCmerged", cacheDir=chd ,reload=T, assignToVariable="covid.merged")
  covid.hashtag=covid.merged
  covid.hashtag@reductions$umap@cell.embeddings[,"UMAP_1"]<- -covid.hashtag@reductions$umap@cell.embeddings[,"UMAP_1"]
  covid.hashtag
}, cacheDir= chd, reload=TRUE, assignToVariable="covid.hashtag")




########################
##### ANALYSIS
########################


########################
###AUC workflow
########################


simpleCache("covidAUC",{
  
  library(AUCell)  
  cells_rankings <- AUCell_buildRankings(covid.hashtag@assays$RNA@counts, nCores=1, plotStats=TRUE)
  #make sure to import the gene lists at the bottom
  simpleCache("genesets", cacheDir=chd,  reload=T)
  
  cells_AUC <- AUCell_calcAUC(genesets, cells_rankings)
  
  #AUCell_plotHist(cells_AUC["unhelped",], aucThr=0.1)
  
  covid.hashtag=AddMetaData(covid.hashtag, metadata=as.data.frame(t(getAUC(cells_AUC))))
}, cacheDir= chd, reload=TRUE, assignToVariable="covid.hashtag")

simpleCache("covidAUCmerged",{
  
  library(AUCell)  
  cells_rankings <- AUCell_buildRankings(covid.merged@assays$RNA@counts, nCores=1, plotStats=TRUE)
  simpleCache("genesets", cacheDir="./rcache",  reload=T)
  cells_AUC <- AUCell_calcAUC(genesets, cells_rankings)
  
  #histograms
  #AUCell_plotHist(cells_AUC["unhelped",], aucThr=0.1)
  
  covid.merged=AddMetaData(covid.merged, metadata=as.data.frame(t(getAUC(cells_AUC))))
}, cacheDir= chd, reload=TRUE, assignToVariable="covid.merged")


#UMAP plotting of gene signature AUCs.

col="dark blue"; col2="yellow"
auc1<-FeaturePlot(object = covid.hashtag, features = "cytotoxic", cols= c(col, col2), reduction = "umap", combine=FALSE)[[1]]; auc1
auc2<-FeaturePlot(object = covid.hashtag, features = "viral",  cols= c(col, col2), reduction = "umap", combine=FALSE)[[1]]; auc2
auc3<-FeaturePlot(object = covid.hashtag, features = "unhelped",  cols= c(col, col2), reduction = "umap", combine=FALSE)[[1]]; auc3
auc4<-FeaturePlot(object = covid.hashtag, features = "ifn",  cols= c(col, col2), reduction = "umap", combine=FALSE)[[1]]; auc4
auc5<-FeaturePlot(object = covid.hashtag, features = "exhaustion", cols= c(col, col2), reduction = "umap", combine=FALSE)[[1]]; auc5


###########
#UMAP plots
###########

tsneHTO=DimPlot(covid.hashtag, group.by = "htolabels", reduction= "umap")+scale_color_brewer(palette="Set2")+labs(title="Treatment"); tsneHTO
tsnept=DimPlot(covid.hashtag, group.by = "patientlabel", reduction= "umap")+labs(title="Patient")+scale_color_manual(values = c(FantasticFox1[3], FantasticFox1[4]));tsnept
tsnecl=DimPlot(covid.hashtag, group.by = "seurat_clusters", reduction= "umap", label=TRUE)+labs(title="Clusters")+scale_color_manual(values = mycolors);tsnecl
tsneresp=DimPlot(covid.hashtag, group.by = "response", reduction= "umap")+labs(title="Response")+scale_color_manual(values = c(Darjeeling1[2], Darjeeling2[3]));tsneresp

pdf("./plots/umaps_overview.pdf", width=15, height=7)
tsnecl+tsneresp+tsnept+tsneHTO
dev.off()

pdf("./plots/auc_umaps.pdf", width=7, height=10)
allaucs<- FeaturePlot(object = covid.hashtag, features = c("cytotoxic","exhaustion", "viral", "unhelped", "ifn"), cols= c(col, col2), reduction = "umap"); allaucs
dev.off()

pdf("./plots/PC_umaps.pdf", width=10, height=10)
FeaturePlot(object = covid.hashtag, features = c("PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6", "PC_7", "PC_8", "PC_9"), cols= c(col, col2), reduction = "umap")#Genes along the pc1 and pc2
dev.off()

###########################
#Plotting interesting genes
##########################
col="dark blue"; col2="yellow"; 
tgene1<-FeaturePlot(object = covid.hashtag, features = "GNLY", cols= c(col, col2), reduction = "umap", combine=FALSE)[[1]]; tgene1
tgene2<-FeaturePlot(object = covid.hashtag, features = "GZMK", cols= c(col, col2), reduction = "umap", combine=FALSE)[[1]]; tgene2


pdf("./plots/gene_features.pdf", width=15, height=15)
tgenes<-FeaturePlot(object = covid.hashtag, features =allgenes, cols= c(col, col2), reduction = "umap"); tgenes
dev.off()


##############################################
###Finding markers, positive wt vs negative wt
###############################################

Idents(covid.hashtag)<- "response"

markerspos2=FindMarkers(covid.hashtag, ident.1="pos", ident.2="neg", test.use="MAST")

#enrichment in negative response
Idents(covid.hashtag)<- "htolabels"
mutpos="MUT YLQ+"
wtpos="WT YLQ+"

markersmw2=FindMarkers(covid.hashtag, ident.2=wtpos, ident.1=mutpos, test.use="MAST")


#################################################
###VIOLIN PLOTS OF SIGNIFICANT GENES. POS vs NEG
#################################################
respcolors=c(Darjeeling1[2], Darjeeling2[3])

pdf("./plots/violinplots_neg_pos_mast_negativeFC.pdf",width=20, height=20 )
allmarkertable=markerspos2
sign="avg_log2FC"
thresh=200
tops=20
allmarkertable=allmarkertable[order(allmarkertable[,sign]),]
VlnPlot(covid.hashtag, rownames(head(allmarkertable[order(allmarkertable[sharedp,sign]),], tops)), pt.size=0, col=respcolors )
dev.off() 



###
###VIOLIN PLOTS OF SIGNIFICANT GENES. MUT vs WT
###
mwcolors=brewer.pal(name="Set2", n=4)[c(2,4)]

pdf("./plots/violinplots_mt_wt_mast_negativeFC.pdf",width=20, height=20 )
allmarkertable=markersmw2
sign="avg_log2FC"
thresh=200
tops=20
allmarkertable=allmarkertable[order(allmarkertable[,sign]),]
VlnPlot(covid.hashtag, rownames(head(allmarkertable[order(allmarkertable[,sign]),], tops)), idents= c(wtpos, mutpos), flip=T,  pt.size=0, col=mwcolors )
dev.off() 

pdf("./plots/violinplots_mt_wt_mast_positiveFCtop20.pdf",width=20, height=20 )
allmarkertable=markersmw2
sign="avg_log2FC"
thresh=200
tops=20
allmarkertable=allmarkertable[order(allmarkertable[,sign], decreasing=T),]
VlnPlot(covid.hashtag, rownames(head(allmarkertable[order(allmarkertable[,sign]),], tops)), idents= c(wtpos, mutpos), flip=T,  pt.size=0, col=mwcolors )
dev.off() 



###
###VIOLIN PLOTS OF SIGNIFICANT GENES. ALL FOUR CONDITIONS
###
all4colors=brewer.pal(name="Set2", n=4)[c(3,1,4,2)]
covid.hashtag@active.ident <- factor(x = covid.hashtag@active.ident, levels = c("WT YLQ-", "MUT YLQ-", "WT YLQ+", "MUT YLQ+"))

pdf("./plots/violinplots_all4_mast_negativeFC_top30.pdf",width=20, height=20 )
allmarkertable=markersmw2
sign="avg_log2FC"
thresh=200
tops=30
allmarkertable=allmarkertable[order(allmarkertable[,sign]),]
VlnPlot(covid.hashtag, rownames(head(allmarkertable[order(allmarkertable[,sign]),], tops)),  pt.size=0, col=all4colors )+labs(ylab="expr.level \n lognorm RC")
dev.off() 

pdf("./plots/violinplots_all4_mast_positiveFCtop30.pdf",width=20, height=20 )
allmarkertable=markersmw2
sign="avg_log2FC"
thresh=200
tops=30
allmarkertable=allmarkertable[order(allmarkertable[,sign], decreasing=T),]
VlnPlot(covid.hashtag, rownames(head(allmarkertable[order(allmarkertable[,sign]),], tops)),  pt.size=0, col=all4colors )+labs(ylab="expr.level \n lognorm RC")
dev.off() 


###
###DE significance thresholding
###
#we make two tests: one for neg vs pos and one for mut vs wt
ntests=2   

capvalue=350
pdf("./plots/volcanoplot_wt_mut_mast.pdf", width=7, height=7)
colhigh=Zissou1[5]
collow=Zissou1[1]
allmarkertable=markersmw2
sign="p_val_adj"
thresh=120; upperthresh<-0; lowerthresh<- 0;

capoutliers<-function(x) if (x>350){return(capvalue)}else{return(x)}
#multiply by 2 for bonferroni correction and cap outliers
allmarkertable[,"p_val_adj"]=allmarkertable[,"p_val_adj"]*ntests #bonferroni correction
allmarkertable=allmarkertable[order(allmarkertable[,sign]),]
allmarkertable[["minuslogval"]]= -log10(allmarkertable[, sign])
allmarkertable[["minuslogval"]]= sapply(allmarkertable[["minuslogval"]], capoutliers)
##selection of markers above thresh for plotting
seltab=allmarkertable[allmarkertable[,"minuslogval"]>thresh,];
seltab[["name"]]=rownames(seltab); 
seltabminus=seltab[seltab[,"avg_log2FC" ]<lowerthresh,];
seltabplus=seltab[seltab[,"avg_log2FC" ]>upperthresh,];
gv= ggplot(allmarkertable, aes(x= avg_log2FC, y=minuslogval))+
  geom_point()+
  geom_hline(yintercept = capvalue, lty=3, col="grey", lwd=1.5)+
  geom_point(data=seltabplus, col=colhigh)+
  geom_text_repel(data=seltabplus, aes(label=name), col=colhigh, max.overlaps=100)+
  geom_point(data=seltabminus, col=collow)+
  geom_text_repel(data=seltabminus, aes(label=name), col=collow, max.overlaps=100)+
  geom_hline(yintercept = thresh, lty=2)+
  coord_cartesian(ylim=c(0, capvalue+10))
gv
mastmw=allmarkertable
dev.off()

#####positive vs negative volcano plot, MAST

pdf("./plots/volcanoplot_pos_neg_mast.pdf", width=7, height=7)
colhigh=Zissou1[5]
collow=Zissou1[1]
allmarkertable=markerspos2
sign="p_val_adj"
thresh=250; upperthresh<- 1; lowerthresh<- -1;
allmarkertable[,"p_val_adj"]=allmarkertable[,"p_val_adj"]*ntests #bonferroni correction
allmarkertable=allmarkertable[order(allmarkertable[,sign]),]
allmarkertable[["minuslogval"]]= -log10(allmarkertable[, sign])
#multiply by to for bonferroni correction and cap outliers
allmarkertable[["minuslogval"]]= sapply(allmarkertable[["minuslogval"]], capoutliers)
##selection of markers above thresh for plotting
seltab=allmarkertable[allmarkertable[,"minuslogval"]>thresh,];
seltab[["name"]]=rownames(seltab); 
seltabminus=seltab[seltab[,"avg_log2FC" ]<lowerthresh,];
seltabplus=seltab[seltab[,"avg_log2FC" ]>upperthresh,];
gvpn= ggplot(allmarkertable, aes(x= avg_log2FC, y=minuslogval))+
  geom_point()+
  geom_hline(yintercept = capvalue, lty=3, col="grey", lwd=1.5)+
  geom_point(data=seltabplus, col=colhigh)+
  geom_text_repel(data=seltabplus, aes(label=name), col=colhigh, max.overlaps=100)+
  geom_point(data=seltabminus, col=collow)+
  geom_text_repel(data=seltabminus, aes(label=name), col=collow, max.overlaps=100)+
  geom_hline(yintercept = thresh, lty=2)+
  geom_vline(xintercept = c(lowerthresh,upperthresh), lty=2)+
  coord_cartesian(ylim=c(0, capvalue+10))
gvpn
mastpos=allmarkertable
dev.off()


###
###exporting MAST table
###
colnames(mastpos)= sapply(colnames(mastpos), function(x) paste0("pos_vs_neg_", x), USE.NAMES=F)
colnames(mastmw)= sapply(colnames(mastmw), function(x) paste0("mut_vs_wt_" ,x ), USE.NAMES=F)


write.table(merge(mastpos, mastmw, by="row.names", incomparables=NA), quote=FALSE, file="./pvalues_table.csv", sep=",")










