library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

#Pre-processing workflow
counts <- Read10X_h5(filename = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(file = "atac_v1_pbmc_10k_singlecell.csv", header = TRUE, row.names = 1)

chrom_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), fragments = 'atac_v1_pbmc_10k_fragments.tsv.gz', min.cells = 10, min.features = 200)

pbmc <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg19"

# add the gene information to the object
Annotation(pbmc) <- annotations

#Computing QC Metrics
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

png(filename="DensityScatter_Plot.png")
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

png(filename="TSS_Plot.png")
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
dev.off()

png(filename="FragmentHistogram_Plot.png")
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
dev.off()

png(filename="Violin_Plot.png")
VlnPlot(object = pbmc, features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'), pt.size = 0.1, ncol = 5)
dev.off()

#Removing the cells which are outliers for the QC metrics
pbmc <- subset(x = pbmc, subset = nCount_peaks > 3000 &  nCount_peaks < 30000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 4 & TSS.enrichment > 3)

#Normalization and linear dimensional reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

png(filename="DepthCorrelation_Plot.png")
DepthCor(pbmc)
dev.off()

#Non-linear dimension reduction and clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

png(filename="NonLinearDimesnionReduction_DimPlot.png")
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()

#Create a gene activity matrix
gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(object = pbmc, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(pbmc$nCount_RNA))

DefaultAssay(pbmc) <- 'RNA'
png(filename="FeaturePlot.png")
FeaturePlot(object = pbmc, features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
dev.off()

#Integrating with scRNA-seq data
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

transfer.anchors <- FindTransferAnchors(reference = pbmc_rna, query = pbmc, reduction = 'cca')
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = pbmc_rna$celltype, weight.reduction = pbmc[['lsi']], dims = 2:30)
pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

png(filename="DimPlot_scRNA_scATAC.png")
plot1 <- DimPlot(object = pbmc_rna, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(object = pbmc, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2
dev.off()

# replace each label with its most likely prediction
for(i in levels(pbmc)) {
  cells_to_reid <- WhichCells(pbmc, idents = i)
  newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
  Idents(pbmc, cells = cells_to_reid) <- newid
}

#Find differentially accessible peaks between cell types
# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(object = pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14+ Monocytes", test.use = 'LR', latent.vars = 'nCount_peaks')

png(filename="Violin_Plot_Differentially_Accesible_Peaks.png")
plot1 <- VlnPlot(object = pbmc, features = rownames(da_peaks)[1], pt.size = 0.1, idents = c("CD4 Naive","CD14+ Monocytes"))
plot2 <- FeaturePlot(object = pbmc, features = rownames(da_peaks)[1], pt.size = 0.1)

plot1 | plot2
dev.off()

fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14+ Monocytes")
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]

open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

#finding closest gene to each of these peaks
closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

#Plotting genomic regions
# set plotting order

levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 effector","Double negative T cell","NK dim", "NK bright", "pre-B cell",'B cell progenitor',"pDC","CD14+ Monocytes",'CD16+ Monocytes')

png(filename="Genomic_Regions_CoveragePlot.png")
CoveragePlot(object = pbmc, region = rownames(da_peaks)[1], extend.upstream = 40000, extend.downstream = 20000)
dev.off()
