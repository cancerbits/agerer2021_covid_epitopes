#Single-cell analysis for Agerer et al, 2021
#-- Analysis of T cell receptor sequences --
#Author: Luis Fernando Montaño-Gutierrez, Florian Halbritter
#Affiliation: Developmental Cancer Genomics Group, St. Anna Children's Cancer Research Insitute, Vienna, Austria.

#Importing libraries
library(Seurat)
library(data.table)
library(devtools)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(simpleCache)
library(ggalluvial)
library(viridis)
library(scRepertoire)


###############################
#declaring recurring variables
###############################
chd="./rcache"
options("RCACHE.DIR"=chd)
colNames <- list(hto="hto_demux", sample="sample", group="htolabels", patient="patient", genotype="genotype", cluster="RNA_snn_res.0.6")
sampleIds <- c("sars042","sars060")
groupOrder <- c("WT YLQ-","MUT YLQ-","WT YLQ+","MUT YLQ+")
freqThresh <- 10
topN <- 5
fixedCols <- c("n/a"="#EEEEEE", "other"="#DDDDDD")
clonoCols <- c("clonotypes_top","response","genotype")
cloneCallTypes <- c("aa") # c("aa","nt","gene","gene+nt")
prefCC <- "aa"
ctUse <- paste0("CT", prefCC)

##############################
#declaring recurring functions
##############################
# add clonotype info to Seurat object:
# from https://raw.githubusercontent.com/j-andrews7/EZscRNA/52d762e59ee9ce59f2d221a9f837354d2ba89a03/R/vdj.R
# (extracted to avoid additional dependencies)
AddClonotype <- function(vdj.dir, scrna) {
	vdj <- read.csv(sprintf("%s/filtered_contig_annotations.csv", vdj.dir))
	vdj <- vdj[!duplicated(vdj$barcode), ]
	vdj <- vdj[, c("barcode", "raw_clonotype_id")]
	names(vdj)[names(vdj) == "raw_clonotype_id"] <- "clonotype_id"
	clono <- read.csv(sprintf("%s/clonotypes.csv", vdj.dir))
	vdj <- merge(vdj, clono[, c("clonotype_id", "cdr3s_aa")])
	vdj <- vdj[, c(2,1,3)]
	rownames(vdj) <- vdj[, 1]
	vdj[, 1] <- NULL
	clono.seurat <- AddMetaData(object = scrna, metadata = vdj)
	return(clono.seurat)
}
pdfPlot <- function(name,width,height) {
	pdf(file=paste0(name,".pdf"), height=height, width=width, pointsize=12, useDingbats=FALSE)
}
dtToDf <- function(dt, rownameCol=1) {
	df <- as.data.frame(dt)
	rownames(df) <- df[,rownameCol]
	if(is.numeric(rownameCol)) {
		 df <- df[,-rownameCol,drop=FALSE] 
	}
	else {
		 df <- df[,setdiff(colnames(df),rownameCol),drop=FALSE] 
	}
	df
}

# load in pre-processed Seurat object from gene expression analysis
simpleCache("covidAUCmerged", cacheDir= chd, reload=TRUE, assignToVariable="sc")

# plot UMAP with metadata to confirm consistency with previous analyses:
pdfPlot("umap_combined", 16, 8)
print(DimPlot(sc, group.by = unlist(colNames), label = T, combine = T))
dev.off()

# load results of CellRanger VDJ analysis:
contiqList <- unlist(sapply(sampleIds, function(s) {	
	tbl <- read.table(baseResultsDir(sprintf("cellranger/bsf/COUNT/%s_VDJ/filtered_contig_annotations.csv", s)), sep=",", header=T)	
	tbl$cell_id <- as.character(paste_(gsub("Covid_","",s),tbl$barcode))	
	tbl$sample_name <- gsub("-\\d+$","",sc@meta.data[tbl$cell_id, colNames$sample])
	return(split(tbl, tbl$sample_name))
}, simplify=F), recursive=F)
	
ids <- gsub("^(.+)\\.(.+)$","\\2",names(contiqList))
samps <- gsub("^Covid_(.+)\\.(.+)$","\\1",names(contiqList))
fullNames <- paste_(samps, ids)
combinedContiqList <- combineTCR(contiqList, ID=ids, samples=samps, cells ="T-AB") #
	
if(grepl("Covid_",colnames(sc)[1])) sc@meta.data$orig_rn <- colnames(sc)

sc@meta.data$bio_id <- gsub("^(.+)_([ATCG]+(-1)?)$", "\\1", colnames(sc))
sc@meta.data$bc <- gsub("^(.+)_([ATCG]+(-1)?)$", "\\2", colnames(sc))

# use scRepertoire's preferred naming style:
newId <- paste_(sc@meta.data[,"bio_id"], sc@meta.data[,colNames$sample], sc@meta.data$bc) 
scX <- RenameCells(sc, new.names=newId)	
scX <- combineExpression(combinedContiqList, scX, cloneTypes=c(None = 0, Single = 1, Small = 10, Medium = 50, Large = 250, Hyperexpanded = 500), cloneCall=prefCC, groupBy = "sample")
	
pdfPlot("screp_umap", 20, 4)
pUmap <- DimPlot(scX, split.by=colNames$sample, group.by = "cloneType", combine = T, cols=rev(viridis(5)))
print(pUmap)
print(occupiedscRepertoire(scX, x.axis = colNames$sample))
dev.off()
pdfPlot("screp_umap_simple", 7, 4)
pUmap <- DimPlot(scX, group.by = "cloneType", combine = T, cols=rev(viridis(5)))
print(pUmap)
dev.off()

scX@reductions$umap@cell.embeddings[,1] <- -1 * scX@reductions$umap@cell.embeddings[,1]
plotCols <- sapply(clonoCols, function(curCol) {
	structure(brewer.pal("Dark2",topN),names=paste0("top-",1:topN))
}, simplify=F)
plotCols[["clonotypes_top"]][names(fixedCols)] <- fixedCols
plotCols[["genotype"]] <- c(WT="#DDDDDD", MUT=plotCols[["clonotypes_top"]][[1]])
plotCols[["response"]] <- c(neg="#DDDDDD", pos=plotCols[["clonotypes_top"]][[1]])
plotCols[["group"]] <- c("MUT YLQ-"="#66C2A5", "MUT YLQ+"="#FC8D62", "WT YLQ-"="#8DA0CB", "WT YLQ+"="#E78AC3") 
plotColsFlat <- unlist(plotCols)
names(plotColsFlat) <- gsub("^.+\\.([^\\.]+)$", "\\1", names(plotColsFlat))


# determine and label top clonotypes (after QC):
dtClonos <- as.data.table(scX@meta.data, keep.rownames=T)
allClonotypes <- dtClonos[, .N, by=.(ct=get(ctUse), bio_id)][!is.na(ct),][,.(N, rnk=rank(-N, ties="random"), ct),by=bio_id][,.(bio_id, ct, N, rnk, ct_id=paste_(bio_id, rnk))]
dtClonos <- merge(dtClonos, allClonotypes, by.x=c("bio_id",ctUse), by.y=c("bio_id","ct"), all=T)
dtClonos[, clonotypes_top:=ifelse(rnk<=topN, paste0("top-",rnk), "other")]
dtClonos[is.na(clonotypes_top), clonotypes_top:="n/a"]
allClonotypes[, clonotypes_top:=ifelse(rnk<=topN, paste0("top-",rnk), "other")]
allClonotypes[is.na(clonotypes_top), clonotypes_top:="n/a"]
scX@meta.data <- dtToDf(dtClonos, "rn")


pdfPlot("umap_clonicity", 6.25, 3.5)
print(DimPlot(scX, group.by = "cloneType", raster=T, combine = T, cols=rev(viridis(5))))
dev.off()
pdfPlot("umap_split_top_clonotypes", 6, 3.5)
print(DimPlot(scX, split.by="bio_id", group.by = "clonotypes_top", cols=plotCols[["clonotypes_top"]], combine = T, raster=T))
dev.off()

pData <- dtClonos[!is.na(rnk),.N,by=c("bio_id", colNames$group, "clonotypes_top")][order(bio_id, get(colNames$group), -N),]
pData[, grp:=factor(get(colNames$group),levels=groupOrder)]
pBars <- ggplot(pData, aes(x=grp, y=N, fill=clonotypes_top)) + scale_fill_manual(values=plotCols[["clonotypes_top"]]) + xlab(NULL) + ylab("# Cells with clonotype") + facet_wrap(~bio_id) + geom_bar(position="fill", stat="identity")
pBars <- pBars + defTheme(flipX=T)  #+NoLegend() 
gg(pBars, "bars_clonotypes", 6.5, 3, type="pdf")


pdfPlot("streams_clonotypes", 6.25, 4)
p <- ggplot(pData, aes(x=grp, fill=clonotypes_top, group=clonotypes_top, stratum=clonotypes_top, alluvium=clonotypes_top,
        y=N, label=clonotypes_top)) + defTheme(flipX=T, topLegend=T) + geom_stratum() + geom_flow(stat = "alluvium", data=pData[grepl("top", clonotypes_top),]) + scale_fill_manual(values=plotCols[["clonotypes_top"]]) + facet_wrap(~bio_id) + xlab(NULL) + ylab("Number of cells")
print(p)
dev.off()

fwrite(dtClonos, file=resultsDir("meta.csv"))
fwrite(allClonotypes, file=resultsDir("clonos.csv"))


pData <- rbindlist(combinedContiqList, idcol="sample_name")
pData[, tcra:=stringr::str_split(stringr::str_split(pData[, CTgene], "_", simplify = TRUE)[,1], "[.]", simplify = TRUE)[, 1]]
pData[, tcrb:=stringr::str_split(stringr::str_split(pData[, CTgene], "_", simplify = TRUE)[,2], "[.]", simplify = TRUE)[, 1]]

pDataM <- melt(pData, measure.vars=c("tcra","tcrb"))[!is.na(value),]
pDataM[, value:=gsub(";?NA;?","",value)]
pDataM <- pDataM[value!="",]
pDataM[, value:=gsub("tr[ab]","",value,ignore.case=T)]

pDataX <- merge(pDataM[,.N, by=.(ID, variable, value)], unique(as.data.table(scX@meta.data)[,.(group=get(colNames$group),ID=get(colNames$sample),patient=get(colNames$patient))]), by="ID", all.x=T)
pDataX[, perc:=N/sum(N), by=ID]
pDataX[, group:=factor(group, levels=groupOrder)]

pDataY <- pDataX[,.(lower=min(perc), upper=max(perc), perc=mean(perc)),by=.(group, variable, value)]

pdfPlot("vgenes_mean", 8, 8)
p <- ggplot(pDataY, aes(x=reorder(value,perc), y=perc, fill=group)) + geom_bar(stat="sum", position="stack") + defTheme() + coord_flip() + facet_wrap(~variable, scales="free") + xlab(NULL) + ylab("Relative frequency (mean across patients)") + scale_fill_manual(values=plotCols[["group"]])
print(p)
dev.off()
pdfPlot("vgenes", 8, 8)
p <- ggplot(pDataX, aes(x=reorder(value,perc), y=perc, color=group, group=group)) + geom_point(aes(shape=patient), position=position_dodge(width=1)) + geom_pointrange(aes(ymax=upper, ymin=lower), shape="|", position=position_dodge(width=1), data=pDataY) + defTheme() + coord_flip() + facet_wrap(~variable, scales="free") + xlab(NULL) + ylab("Relative frequency per sample") + scale_color_manual(values=plotCols[["group"]])
print(p)
dev.off()
pdfPlot("vgenes_larger", 8, 12)
print(p)
dev.off()
	
fwrite(pData, file=resultsDir("clonos2.csv"))	

