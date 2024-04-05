#Code to make figures from single cell RNAseq analysis
library(Seurat)
library(data.table)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggbeeswarm)
library(ComplexHeatmap)
library(ggpubr)
library(SCPA)
library(msigdbr)
library(ggrepel)
library(dplyr)
library(patchwork)

a1 <- data.table(patient_response_paper = rep(c("Patient 1-Resistant",
                                                "Patient 2-Resistant",
                                                "Patient 3-Responding",
                                                "Patient 4-Partial Responding"),
                                              c(2,3,4,5)),
                 biopsy = c(0,2,0,4,48,0,5,14,60,0,5,11,39,62),
                 time = c(2,2,48,48,48,60,60,60,60,62,62,62,62,62))


cell_col <- data.table(Cell = c("CD8+ T-cells",
                                "CD4+ T-cells",
                                "NK cells",
                                "B-cells",
                                "Malignant",
                                "Macrophages"),
                       Color = c("#D62728FF",
                                 "#2CA02CFF",
                                 "#9467BDFF",
                                 "#E377C2FF",
                                 "#7F7F7FFF",
                                 "#FF7F0EFF"))

cell_col2 <- cell_col$Color
names(cell_col2) <- cell_col$Cell

cell_col3 <- data.table(Cell = c("CD8+ T-cells",
                                 "CD4+ T-cells",
                                 "NK cells",
                                 "B-cells",
                                 "Malignant",
                                 "TTR macrophages", 
                                 "TTS macrophages"),
                        Color = c("#D62728FF",
                                  "#2CA02CFF",
                                  "#9467BDFF",
                                  "#E377C2FF",
                                  "#7F7F7FFF",
                                  "#0073C2FF",
                                  "#EFC000FF"))

cell_col4 <- cell_col3$Color
names(cell_col4) <- cell_col3$Cell

pat_col <- data.table(Patient = c("Patient 1-Resistant",
                                  "Patient 2-Resistant",
                                  "Patient 3-Responding",
                                  "Patient 4-Partial Responding"),
                      Color = c("#374E55FF",
                                "#DF8F44FF",
                                "#00A1D5FF",
                                "#B24745FF"))
pat_col2 <- pat_col$Color
names(pat_col2) <- pat_col$Patient

mac_col <- c( "#0073C2FF",
              "#EFC000FF")

names(mac_col) <- c("TTR macrophages", 
                    "TTS macrophages")

a2 <- data.table(patient_response_paper = c("Patient 1-Resistant",
                                            "Patient 2-Resistant",
                                            "Patient 3-Responding",
                                            "Patient 4-Partial Responding",
                                            "Patient 4-Partial Responding"),
                 biopsy = c(2,48,60,11,62),
                 Response = c("Non-responder",
                              "Non-responder",
                              "Responder",
                              "Responder",
                              "Non-responder"))

j1 <- readRDS("FNA_Seurat.rds")
j1$cell_simp2 <- factor(j1$cell_simp, levels = c("Macrophages", "CD4+ T-cells", "CD8+ T-cells", "NK cells", "B-cells", "Malignant"))

j1t <- j1@meta.data

#Figure 1B
ggplot(a1, aes(x = biopsy, y = patient_response_paper)) +
  geom_col(data = a2, aes(fill = Response)) +
  geom_point(size = 3) +
  scale_fill_npg() +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 16, family = "URWHelvetica"),
        axis.title = element_text(size = 20, family = "URWHelvetica"),
        strip.text = element_text(size = 16, family = "URWHelvetica"),
        legend.text = element_text(size = 16, family = "URWHelvetica"),
        legend.title = element_text(size = 20, family = "URWHelvetica")) +
  ylab("Patient") +
  xlab("Days on treatment")
ggsave("Figure 1B", width = 12, height= 8)

#Figure 1C
pdf("Figure 1C.pdf", width = 12, height = 8)
DimPlot(j1, group.by = "cell_simp", cols = cell_col2)
dev.off()

#Figure 1D

ggplot(j1t, aes(x = biopsy, fill = cell_simp2)) +
  geom_bar(position = "fill") +
  theme_bw() +
  facet_wrap(~patient_response_paper, nrow = 1) +
  scale_fill_manual(values = cell_col2, name = "Cell types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20)) +
  ylab("Proportion") +
  xlab("Biopsy")
ggsave("Figure 1D", width = 12, height= 8)

#Supplemental Figure 1A
pdf("Supplemental Figure 1A.pdf", width = 12, height = 8)

DimPlot(j1, group.by = "seurat_clusters")
dev.off()

#Supplemental Figure 1B
j1.count <- as.matrix(GetAssayData(j1, slot="counts"))
j1.annot <- data.frame(V1 = j1$cell_type)

j1.cnv.cd4 <- CreateInfercnvObject(raw_counts_matrix=j1.count,
                                   annotations_file= j1.annot,
                                   delim = "\t",
                                   gene_order_file= "gene_pos.txt",
                                   ref_group_names="CD4+ T-cells")

j1.cnv.cd4.run <- infercnv::run(j1.cnv.cd4,
                                analysis_mode = "subclusters",
                                #tumor_subcluster_partition_method = "random_trees",
                                tumor_subcluster_pval = 0.1,
                                k_nn= 30,
                                leiden_resolution = 0.15 ,
                                num_ref_groups = 1,
                                cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir="infercnv/",  # dir is auto-created for storing outputs
                                cluster_by_groups=FALSE,   # cluster
                                denoise=TRUE,
                                HMM=TRUE,
                                num_threads = 20,
                                no_plot = TRUE)

cgroup <- data.table(Sample_barcode = j1$SampBar,
                     Cell_type = j12$cell_type,
                     Paper_name = j1$patient_response_paper)
setkey(cgroup, Sample_Barcode)
gexp <- data.table(Gene = dimnames(j1.cnv.cd4.run@expr.data)[[1]],
                   j1.cnv.cd4.run@expr.data)
gch <- data.table(Gene  = dimnames(j1.cnv.cd4.run@count.data)[[1]],
                  j1.cnv.cd4.run@gene_order)

gch[, chr := factor(chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                                    "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))]

gech <- merge(gch, gexp, by = "Gene")

gech[, chrn := as.numeric(chr)]

setkey(gech, chr, start)

gech.mt <- data.table::transpose(gech[, setdiff(colnames(gech), c("chr", "chrn", "start", "stop")), with = FALSE],
                                 keep.names = "Cell", make.names = "Gene")

cnv1 <- as.matrix(gech.mt[, setdiff(colnames(gech.mt), "Cell"), with = FALSE])
rownames(cnv1) <- gech.mt$Cell

row_ha <- rowAnnotation(df = data.frame(
  Cell_type = cgroup[rownames(cnv1), Cell_type],
  Sample = cgroup[rownames(cnv1), Paper_name]),
  col = list(Cell_type = cc_color,
             Sample = sample_color))
png("Supplemental Figure 1B.png", width = 1200, height = 800)

Heatmap(cnv1, 
        column_split = gech$chrn, 
        use_raster = TRUE,
        raster_resize_mat = mean,
        cluster_columns = FALSE,
        left_annotation = row_ha,
        cluster_rows = j1.cnv.cd4.run@tumor_subclusters$hc[["all_observations"]],
        show_row_names = FALSE,
        show_column_names = FALSE)
dev.off()

#Figure S1C
ggplot(j1t[j1t$cell_type %in% "Malignant",], aes(x = Melanoma_UMAP_1, y = Melanoma_UMAP_2, color = patient_response_paper)) +
  geom_point() +
  scale_color_manual(values = pat_col2, name = "Patient") +
  theme_bw() +
  xlab("UMAP.1") +
  ylab("UMAP.2")
ggsave("Supplemental Figure 1C", width = 12, height= 8)

#Supplemental Figure 1D
s1c <- as.matrix(j1t[j1t$cell_type %in% "Malignant", c("HOEK_Proliferative6", "HOEK_Invasive7", "VERFAILLIE_Proliferative8", 
                                                       "VERFAILLIE_INVASIVE9", "TIROSH_MITF10", "TIROSH_AXL11", 
                                                       "TSOI_Undifferentiated12", "TSOI_Undifferetiated_Neural_crest_like13", 
                                                       "TSOI_Neural_crest_like14", "TSOI_Neural_crest_like_transitory15", 
                                                       "TSOI_Transitory16", "TSOI_Transitory_melanocytic17", "TSOI_Melanocytic18",
                                                       "RAMBOW_Mitosis19", "RAMBOW_Neuro20", "RAMBOW_MITF_targets21", 
                                                       "RAMBOW_Pigmentation22", "RAMBOW_Immune23", "RAMBOW_Invasion24", 
                                                       "RAMBOW_Hypometabolic25", "RAMBOW_MSC26")])

rownames(s1c) <- j1t[j1t$cell_type %in% "Malignant", "SampBar"]

s1c2 <- t(s1c)

annot_s1c <- HeatmapAnnotation(df = data.frame(Patient = j1t[j1t$cell_type %in% "Malignant", "patient_response_paper"]),
                               col = list(Patient = pat_col2))


png("Supplemental Figure 1D.png", width = 1200, height = 800, units = "px")
Heatmap(s1c2, use_raster = FALSE,
        show_column_names = FALSE,
        row_labels = c("HOEK Proliferative", "HOEK Invasive", "VERFAILLIE Proliferative", 
                       "VERFAILLIE INVASIVE", "TIROSH MITF", "TIROSH AXL", 
                       "TSOI Undifferentiated", "TSOI Undifferetiated Neural crest like", 
                       "TSOI Neural crest like", "TSOI Neural crest like transitory", 
                       "TSOI Transitory", "TSOI Transitory melanocytic", "TSOI Melanocytic",
                       "RAMBOW Mitosis", "RAMBOW Neuro", "RAMBOW MITF targets", 
                       "RAMBOW Pigmentation", "RAMBOW Immune", "RAMBOW Invasion", 
                       "RAMBOW Hypometabolic", "RAMBOW MSC"),
        column_split = j1t[j1t$cell_type %in% "Malignant", "patient_response_paper"],
        top_annotation = annot_s1c)
dev.off()

#Figure 2A
sec <- fread("Full_Secretome_with _Mphage_Factors.txt", header = FALSE)
Idents(j1.edit2) <- j1.edit2$cell_type

malig.markers.deseq <- FindMarkers(j1, test.uses = "DESeq2", min.pct = 0.01,
                                   subset.ident = "Malignant", 
                                   group.by = "response", 
                                   ident.1 = "Resistant",
                                   ident.2 = "Responding")
malig.markers.deseq$label <- rownames(malig.markers.deseq)
setDT(malig.markers.deseq)
malig.markers.deseq[, Secreted := "No"]
malig.markers.deseq[label %in% sec$V1, Secreted := "Yes"]
malig.markers.deseq[is.na(Secreted), Secreted := "No"]
#assign minimum value to p-values equal to zero
malig.markers.deseq$logFDR <- -log10(malig.markers.deseq$p_val_adj + 2e-308)

malig.markers.deseq[, label2 := label]
malig.markers.deseq[p_val_adj > 0, label2 := NA]
malig.markers.deseq[abs(avg_log2FC) < 5, label2 := NA]
malig.markers.deseq[abs(avg_log2FC) < 1, label2 := NA]
malig.markers.deseq[abs(avg_log2FC) < 1, Significance := "None"]
malig.markers.deseq[avg_log2FC > 1, Significance := "Upregulated"]
malig.markers.deseq[avg_log2FC < -1, Significance := "Down-regulated"]


ggplot(malig.markers.deseq, aes(x = avg_log2FC, y = logFDR, color = Significance, label = label2)) +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = c(-1,1)) +
  geom_point(size = 3) +
  
  geom_label_repel(aes(fill = Secreted), ylim = c(-5, 350), 
                   color = "black", force = 1.5, max.overlaps = 30, size = 6) +
  scale_color_jama()+
  scale_fill_lancet(alpha = 0.4) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  xlab("Fold change (log2)") +
  ylab("FDR (-log10)") +
  scale_y_continuous(limits = c(-5, 350))
ggsave("Figure 2A.pdf", width = 12, height = 8)

#Figure 2B
fig2b <- FetchData(object = j1, vars = c("patient_paper", "cell_simp", "POSTN"), layer = "data")
fig2b$cell_simp <- factor(fig2b$cell_simp, levels = c("Macrophages", "CD4+ T-cells", "CD8+ T-cells",
                                                      "NK cells", "B-cells", "Malignant"))
ggplot(fig2b, aes(x = cell_simp, y = POSTN, color = cell_simp)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col2, name = "Cell type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~patient_paper) +
  xlab("")
ggsave("Figure 2B.pdf", width = 12, height = 8)

#Supplemental Figure 2A
datas2a <- FetchData(object = j1, vars = c("patient_paper", "cell_simp", "biopsy", "HIST1H1E", "MGP", "POSTN", "SPP1", "SEMA3C", "APOC1", "COL1A2", 
                                                 "DCT", "PLCG2", "PSAT1", "MDK", "SCIN", "DMD", "SYNE2", "BSG", 
                                                 "CENPF", "HSP90AB1", "NREP", "TOMM7", "SH3BGRL3", "APP", "HIST1H4C", 
                                                 "MBP", "TENM3", "GSTP1", "CDC42", "GRK3", "MT2A", "KCNAB2", "CD81", 
                                                 "CTSB", "KRTAP19-1", "S100A10", "BHLHE41", "DLL3", "S100B", "L1CAM", 
                                                 "SGK1", "FXYD5", "IFI27", "SEMA6A", "SCRN1", "CD63", "CAPG", 
                                                 "LGALS3", "MSRB2", "TSTD1", "CAV1", "QPCT", "IFI6", "MT1E", "MFSD12", 
                                                 "RUNX3", "KIT", "ARMCX3"))

datas2a.mal <- datas2a[datas2a$cell_simp == "Malignant",]

figs2a.scale <- t(apply(datas2a.mal[, c("HIST1H1E", "MGP", "POSTN", "SPP1", "SEMA3C", "APOC1", "COL1A2", 
                                        "DCT", "PLCG2", "PSAT1", "MDK", "SCIN", "DMD", "SYNE2", "BSG", 
                                        "CENPF", "HSP90AB1", "NREP", "TOMM7", "SH3BGRL3", "APP", "HIST1H4C", 
                                        "MBP", "TENM3", "GSTP1", "CDC42", "GRK3", "MT2A", "KCNAB2", "CD81", 
                                        "CTSB", "KRTAP19-1", "S100A10", "BHLHE41", "DLL3", "S100B", "L1CAM", 
                                        "SGK1", "FXYD5", "IFI27", "SEMA6A", "SCRN1", "CD63", "CAPG", 
                                        "LGALS3", "MSRB2", "TSTD1", "CAV1", "QPCT", "IFI6", "MT1E", "MFSD12", 
                                        "RUNX3", "KIT", "ARMCX3")], 2, scale))
colnames(figs2a.scale) <- rownames(datas2a.mal)

figs2a.annot <- HeatmapAnnotation(df = data.frame(Patient = datas2a.mal$patient_paper,
                                                  Biopsy = datas2a.mal$biopsy),
                                  col = list(Patient = c("Patient 1" = "#374E55FF",
                                                         "Patient 2" = "#DF8F44FF",
                                                         "Patient 3" = "#00A1D5FF",
                                                         "Patient 4" = "#B24745FF"),
                                             Biopsy = c("Baseline" =  "#1F77B4FF",
                                                        "PostTreatment1" = "#FF7F0EFF",
                                                        "PostTreatment2" = "#2CA02CFF",
                                                        "PostTreatment3" = "#D62728FF",
                                                        "PostTreatment4" = "#9467BDFF")))
pdf("Supplemental Figure 2A.pdf", height = 12, width= 8)
Heatmap(figs2a.scale,
        top_annotation = figs2a.annot,
        circlize::colorRamp2(c(-4, -2, 0, 2, 4), c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020")),
        cluster_columns = FALSE,
        raster_resize_mat = mean,
        cluster_rows = TRUE,
        column_order = order(datas2a.mal$patient_paper, decreasing = FALSE),
        show_column_names = FALSE,
        use_raster = TRUE)
dev.off()

#Supplemental Figure 2B
figs2b <- FetchData(object = j1, vars = c("patient_paper", "cell_simp", "MDK", "SPP1"), layer = "data")
figs2b$cell_simp <- factor(figs2b$cell_simp, levels = c("Macrophages", "CD4+ T-cells", "CD8+ T-cells",
                                                        "NK cells", "B-cells", "Malignant"))
fig2b_top <- ggplot(figs2b, aes(x = cell_simp, y = MDK, color = cell_simp)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col2, name = "Cell types") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~patient_paper) +
  xlab("")

fig2b_bot <-ggplot(figs2b, aes(x = cell_simp, y = SPP1, color = cell_simp)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col2, name = "Cell types") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~patient_paper) +
  xlab("")

fig2b_top / fig2b_bot + plot_layout(guides = "collect") &
  theme(legend.position='bottom')
ggsave("Supplemental Figure 2B.pdf", height = 12, width = 8)

#Supplemental Figure 2C
figs2c <- FetchData(object = j1, vars = c("patient_paper", "cell_simp", "biopsy", "POSTN", "MDK", "SPP1"), layer = "data")
figs2c2 <- figs2c[figs2c$cell_simp == "Malignant",]
figs2c3 <- reshape2::melt(figs2c2, id.vars = c("patient_paper", "cell_simp", "biopsy"))

ggplot(figs2c3, aes(x = biopsy, y = value, color = biopsy)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_d3(name = "Biopsy") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_grid(variable~patient_paper)
ggsave("Supplemental Figure 2C.pdf", height = 12, width = 8)

#Figure 3A
fig3a <- FetchData(object = j1, vars = c("patient_paper", "cell_simp", "ITGB5"), layer = "data")
fig3a$cell_simp <- factor(fig3a$cell_simp, levels = c("Macrophages", "CD4+ T-cells", "CD8+ T-cells",
                                                      "NK cells", "B-cells", "Malignant"))
ggplot(fig3a, aes(x = cell_simp, y = ITGB5, color = cell_simp)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col2, name = "Cell types") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~patient_paper, nrow = 1) +
  xlab("")
ggsave("Figure 3A.pdf", height = 12, width = 8)

#Figure 3B
fig3b_top<- DimPlot(j1, group.by = "cell_simp", cols = cell_col2)

fig3b_bot <- ggplot(j1t[j1t$cell_simp %in% "Macrophages",], aes(x = Macrophage_UMAP_1, y = Macrophage_UMAP_2, color = Macrophage_cluster)) +
  geom_point() +
  scale_color_manual(values = c("Cluster 1" = "#D62728FF",
                                "Cluster 2" = "#2CA02CFF",
                                "Cluster 3" = "#FF7F0EFF",
                                "Cluster 4" = "#1F77B4FF"), 
                     name = "Cluster") +
  theme_bw() +
  xlab("UMAP.1") +
  ylab("UMAP.2")

fig3b_top / fig3b_bot
ggsave("Figure 3B.pdf", height = 12, width = 8)

#Figure 3C
data3c <- FetchData(object = j1, vars = c("patient_paper", "cell_type", "Macrophage_cluster", "C1QC", "SPP1","FCN1", "MNDA", "CD163", "GPNMB", "LIPA",
                                                "MSR1", "MRC1", "CCL2", "MARCO", "CD63",
                                                "LGALS3", "CD48"), layer = "data")
data3c.mac <- data3c[!is.na(data3c$Macrophage_cluster),]

fig3.scale <- t(apply(data3c.mac[, c("C1QC", "SPP1","FCN1", "MNDA", "CD163", "GPNMB", "LIPA",
                                     "MSR1", "MRC1", "CCL2", "MARCO", "CD63",
                                     "LGALS3", "CD48")], 2, scale))
colnames(fig3.scale) <- rownames(data3c.mac)

fig3.annot <- HeatmapAnnotation(df = data.frame(Patient = data3c.mac$patient_paper,
                                                Macrophage = data3c.mac$cell_type,
                                                Cluster = data3c.mac$Macrophage_cluster),
                                col = list(Patient = c("Patient 1" = "#374E55FF",
                                                       "Patient 2" = "#DF8F44FF",
                                                       "Patient 3" = "#00A1D5FF",
                                                       "Patient 4" = "#B24745FF"),
                                           Macrophage = mac_col,
                                           Cluster = c("Cluster 1" = "#D62728FF",
                                                       "Cluster 2" = "#2CA02CFF",
                                                       "Cluster 3" = "#FF7F0EFF",
                                                       "Cluster 4" = "#1F77B4FF")))
pdf("Figure 3C.pdf", height = 12, width = 8)
Heatmap(fig3.scale,
        top_annotation = fig3.annot,
        circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")),
        cluster_columns = FALSE,
        column_order = order(data3c.mac$Macrophage_cluster, decreasing = TRUE),
        show_column_names = FALSE)
dev.off()

#Figure 3D
fig3d <- reshape2::melt(data3c.mac[, c("cell_type", "CD63", "CD163","LIPA", 
                                       "MSR1", "MRC1", "GPNMB", "LGALS3",
                                       "MARCO", "CCL2", "SPP1")], id.vars = "cell_type")

ggplot(fig3d, aes(x = cell_type, y = value, fill = cell_type)) +
  geom_violin() +
  geom_quasirandom(size = 2, pch = 21) +
  scale_fill_manual(values = mac_col, name = "Macrophage") +
  theme_bw() +
  facet_wrap(~variable, nrow = 4) +
  stat_compare_means(method = "wilcox", comparisons = list(c("TTR macrophages", "TTS macrophages")), 
                     label = "p.signif", label.y = 11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  ylab("Expression level") +
  xlab("Macrophage")
ggsave("Figure 3D.pdf", height = 12, width = 8)

#Figure 3E
ggplot(j1t, aes(x = biopsy, fill = cell_type)) +
  geom_bar(position = "fill") +
  theme_bw() +
  facet_wrap(~patient_response_paper, nrow = 1) +
  scale_fill_manual(values = cell_col4, name = "Cell types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20)) +
  ylab("Proportion") +
  xlab("Biopsy")
ggsave("Figure 3E.pdf", width = 12, height= 8)

#Supplemental Figure 4A
figs4a <- FetchData(object = j1, vars = c("patient_paper", "cell_simp", "ITGAV", "ITGB3"), layer = "data")

figs4a$cell_simp <- factor(figs4a$cell_simp, levels = c("Macrophages", "CD4+ T-cells", "CD8+ T-cells",
                                                        "NK cells", "B-cells", "Malignant"))

figs4a2 <- reshape2::melt(figs4a, id.vars = c("patient_paper", "cell_simp"))

ggplot(figs4a2, aes(x = cell_simp, y = value, color = cell_simp)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col2, name = "Cell types") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_grid(variable~patient_paper) +
  xlab("")
ggsave("Supplemental Figure 4A.pdf", width = 12, height = 8)

#Supplemental Figure 4B
figs4b <- FetchData(object = j1, vars = c("patient_paper", "cell_simp", "LRP1", "LRP2",
                                                "ALK", "NCAM1", "NOTCH2", "PTPRB", "PTPRZ1"), layer = "data")

figs4b$cell_simp <- factor(figs4b$cell_simp, levels = c("Macrophages", "CD4+ T-cells", "CD8+ T-cells",
                                                        "NK cells", "B-cells", "Malignant"))

figs4b2 <- reshape2::melt(figs4b, id.vars = c("patient_paper", "cell_simp"))

ggplot(figs4b2, aes(x = cell_simp, y = value, color = cell_simp)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col2, name = "Cell types") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_grid(variable~patient_paper) +
  xlab("")
ggsave("Supplemental Figure 4B.pdf", width = 12, height = 16)

#Supplemental Figure 4C
figs4c <- FetchData(object = j1, vars = c("patient_paper", "cell_type", "M1_literature4", "M2_literature5"), layer = "data")
figs4c2 <- figs4c[figs4c$cell_type %in% c("TTR macrophages", "TTS macrophages"),]

ggplot(figs4c2, aes(x = M1_literature4, y = M2_literature5, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = mac_col, name = "Macrophage") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  xlab("Martinez M1 enrichment") +
  ylab("Martinez M2 enrichment")
ggsave("Supplemental Figure 4C.pdf", width = 12, height = 8)

#Supplemental Figure 4D
figs4d <- FetchData(object = j1, vars = c("patient_paper", "cell_type", "M1_Newman2", "M2_Newman3"), layer = "data")
figs4d2 <- figs4d[figs4d$cell_type %in% c("TTR macrophages", "TTS macrophages"),]

ggplot(figs4d2, aes(x = M1_Newman2, y = M2_Newman3, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = mac_col, name = "Macrophage") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  xlab("Newman M1 enrichment") +
  ylab("Newman M2 enrichment")
ggsave("Supplemental Figure 4D.pdf", width = 12, height = 8)

#Supplemental Figure 5A
figs5a <- FetchData(object = j1, vars = c("patient_paper", "Macrophage_cluster", "cell_type", "ITGAV", "ITGB3", "ITGB5"), layer = "data")

figs5a2 <- figs5a[figs5a$cell_type %in% c("TTR macrophages", "TTS macrophages"),]

figs5a3 <- reshape2::melt(figs5a2, id.vars = c("patient_paper", "cell_type", "Macrophage_cluster"))

ggplot(figs5a3, aes(x = cell_type, y = value, color = cell_type)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = mac_col, name = "Macrophage") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_grid(variable~patient_paper) +
  stat_compare_means(method = "wilcox", comparisons = list(c("TTR macrophages", "TTS macrophages")), 
                     label = "p.signif", label.y = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  xlab("")
ggsave("Supplemental Figure 5A.pdf", width = 12, height = 8)

#Supplemental Figure 5B
figs5b <- FetchData(object = j1, vars = c("patient_paper", "cell_type", "LRP1", "LRP2",
                                                "ALK", "NCAM1", "NOTCH2", "PTPRB", "PTPRZ1"), layer = "data")

figs5b2 <- figs5b[figs5b$cell_type %in% c("TTR macrophages", "TTS macrophages"),]

figs5b3 <- reshape2::melt(figs5b2, id.vars = c("patient_paper", "cell_type"))

ggplot(figs5b3, aes(x = cell_type, y = value, color = cell_type)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = mac_col, name = "Macrophage") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_grid(variable~patient_paper) +
  stat_compare_means(method = "wilcox", comparisons = list(c("TTR macrophages", "TTS macrophages")), 
                     label = "p.signif", label.y = 4.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  xlab("")
ggsave("Supplemental Figure 5A.pdf", width = 12, height = 16)

#Supplemental Figure 5C
cc1 <- createCellChat(object = GetAssayData(object = j1, assay = "RNA", slot = "counts"),
                      meta = j1t,
                      group.by = "cell_type")

ccdb <- CellChatDB.human

cc1@DB <- ccdb

cc1 <- subsetData(cc1)

future::plan("multisession", workers = 60)

cc1 <- identifyOverExpressedGenes(cc1)
cc1 <- identifyOverExpressedInteractions(cc1)
# project gene expression data onto PPI network (optional)
cc1 <- projectData(cc1, PPI.human)
cc1 <- computeCommunProb(cc1,
                         raw.use = FALSE,
                         population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cc1 <- filterCommunication(cc1, min.cells = 10)

cc1 <- computeCommunProbPathway(cc1)

cc1 <- aggregateNet(cc1)

vertex.receiver <- 1:7
pdf("Supplemental Figure 5C top.pdf", width = 12, height = 8)
netVisual_aggregate(cc1, signaling = "MK",  vertex.receiver = vertex.receiver)
dev.off()

pdf("Supplemental Figure 5C bottom.pdf", width = 12, height = 8)
netVisual_aggregate(cc1, signaling = "PERIOSTIN",  vertex.receiver = vertex.receiver)
dev.off()

#Supplemental Figure 6A
figs6a <- FetchData(object = j1, vars = c("cell_type", "VEGFA", "TNF"), layer = "data")

figs6a2 <- figs6a[figs6a$cell_type %in% c("TTR macrophages", "TTS macrophages"),]

figs6a3 <- reshape2::melt(figs6a2, id.vars = c("cell_type"))

ggplot(figs6a3, aes(x = cell_type, y = value, color = cell_type)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = mac_col, name = "Macrophage") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~variable, nrow = 2) +
  xlab("") +
  ylab("Expression level")
ggsave("Supplemental Figure 6A.pdf", width = 12, height = 16)

#Supplemental Figure 6B
cc1.net <- subsetCommunication(cc1, slot.name = "net")
setDT(cc1.net)

tts.mel <- cc1.net[source %in% c("TTS macrophages") & target %in% c("Malignant"), interaction_name]
ttr.mel <- cc1.net[source %in% c("TTR macrophages") & target %in% c("Malignant"), interaction_name]

pdf("Supplemental Figure 6B.pdf", width = 8, height = 16)
netVisual_bubble(cc1, 
                 sources.use = c(6,7), 
                 targets.use = c(4), 
                 pairLR.use = data.frame(interaction_name = union(tts.mel, ttr.mel)),
                 remove.isolate = FALSE)
dev.off()

#Figure 7C
top10 <- cc1.net[source %in% c("TTR macrophages", "TTR macrophages") & target %in% c("Malignant"), ][order(prob,decreasing = TRUE), interaction_name]

pdf("Figure 7C.pdf", width = 12, height = 8)
netVisual_bubble(cc1, 
                 sources.use = c(6,7), 
                 targets.use = c(4), 
                 pairLR.use = data.frame(interaction_name = top10[1:10]),
                 remove.isolate = FALSE)
dev.off()

#Supplemental Figure 7A
figs7a <- FetchData(object = j1, vars = c("patient_paper", "cell_type", "SPP1"), layer = "data")

figs7a$cell_type <- factor(figs7a$cell_type, levels = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", 
                                                        "Malignant",  "NK cells", "TTR macrophages",
                                                        "TTS macrophages"))

ggplot(figs7a, aes(x = cell_type, y = SPP1, color = cell_type)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col4, name = "Cell type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~patient_paper) +
  xlab("")

#Supplemental Figure 7B
figs7b <- FetchData(object = j1, vars = c("patient_paper", "cell_type", "CD44"), layer = "data")

figs7b$cell_type <- factor(figs7b$cell_type, levels = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", 
                                                        "Malignant",  "NK cells", "TTR macrophages",
                                                        "TTS macrophages"))

ggplot(figs7b, aes(x = cell_type, y = CD44, color = cell_type)) +
  geom_violin() +
  geom_quasirandom(size = 2)+
  theme_bw()+
  scale_color_manual(values = cell_col4, name = "Cell type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12, family = "URWHelvetica"),
        legend.text = element_text(size = 12, family = "URWHelvetica"),
        legend.title = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 14, family = "URWHelvetica"),
        title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~patient_paper) +
  xlab("")
