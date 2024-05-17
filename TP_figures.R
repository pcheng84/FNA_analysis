#Figures for Vasilevska et al. Tumor profiler analysis
library(scater)
library(scran)
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

cell_tp <- data.table(Cell = c("B.cells",
                               "CAFs",
                               "Dendritic.cells",
                               "Endothelial.cells",
                               "Macrophages",
                               "Melanoma.melanocytic",
                               "Melanoma.mesenchymal",
                               "Plasma.cells",
                               "Plasmcytoid.dendritic.cell",
                               "T.cells",
                               "uncertain",
                               "unknown"),
                      Color = c("#E377C2FF",
                                "#1F77B4FF",
                                "#2CA02CFF",
                                "#D62728FF",
                                "#FF7F0EFF",
                                "#7F7F7FFF",
                                "#7F7F7FFF",
                                "#9467BDFF",
                                "#8C564BFF",
                                "#BCBD22FF",
                                "#17BECFFF",
                                "#AEC7E8FF"))
cell_tp2 <- cell_tp$Color
names(cell_tp2) <- cell_tp$Cell

cell_tp3 <- data.table(Cell = c("B.cells",
                                "CAFs",
                                "Dendritic.cells",
                                "Endothelial.cells",
                                "Macrophages",
                                
                                "TTR macrophages", 
                                "TTS macrophages",
                                "Malignant",
                                "Plasma.cells",
                                "Plasmacytoid.dendritic.cell",
                                "T.cells",
                                "uncertain",
                                "unknown"),
                       Color = c("#E377C2FF",
                                 "#92C5DE",
                                 "#2CA02CFF",
                                 "#D62728FF",
                                 "#FF7F0EFF",
                                 
                                 "#0073C2FF",
                                 "#EFC000FF",
                                 
                                 "#7F7F7FFF",
                                 
                                 "#9467BDFF",
                                 "#8C564BFF",
                                 "#BCBD22FF",
                                 "#17BECFFF",
                                 "#E0E0E0"))

cell_tp4 <- cell_tp3$Color
names(cell_tp4) <- cell_tp3$Cell

cell_mac <- c("#FF7F0EFF", "#EFC000FF", "#0073C2FF")
names(cell_mac) <- c("Macrophages",  "TTS macrophages", "TTR macrophages")


m2b.time <- fread("TuPro_previous_treatment.txt")
m2b.time[, Biopsy := 0]
m2b.time[Previous_treatment_class == "", Previous_treatment_class := NA]
m2b.time[During_treatment_class == "", During_treatment_class := NA]

m2b.time[Next_treatment_class == "", Next_treatment_class := NA]

tp.braf2 <- readRDS("TP_SCE.rds")

tpt <- makePerCellDF(tp.braf2, features = c("POSTN", "MLANA", "MSR1", "GPNMB", "CD63", "CD68", "CD163", "MDK", "SPP1",
                                             "AXL", "SOX9", "SOX10", "LIPA", "TNF", "VEGFA", "CCL2", 
                                             "LGALS3", "CXCL16", "MARCO", "TGFB1", "LRP1", "IL1B", 
                                             "IL15", "CD274", "IDO1", "PTN", "MURF1", "WNT5A", "CAPN3",
                                             "LRP1", "ITGAV", "ITGB3", "ITGB5"))
setDT(tpt)
#Figure 2H
ggplot(tpt, aes(x = celltype_paper, y = POSTN, color = celltype_paper)) +
  geom_quasirandom(bandwidth = 2) +
  scale_color_manual(values = cell_tp2, name = "Cell type")+
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
        axis.text.y = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 16, family = "URWHelvetica"),
        legend.text = element_text(size = 14, family = "URWHelvetica"),
        legend.title = element_text(size = 16, family = "URWHelvetica"))
ggsave("Figure 2H.pdf", width = 12, height = 8)

#Figure 2I
ggplot(tpt, aes(x = celltype_paper, y = MDK, color = celltype_paper)) +
  geom_quasirandom(bandwidth = 2) +
  scale_color_manual(values = cell_tp2, name = "Cell type")+
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
        axis.text.y = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 16, family = "URWHelvetica"),
        legend.text = element_text(size = 14, family = "URWHelvetica"),
        legend.title = element_text(size = 16, family = "URWHelvetica"))
ggsave("Figure 2I.pdf", width = 12, height = 8)

#Figure 6B
ggplot(m2b.time) +
  geom_segment(aes(x = Prev_start, xend = Prev_end, y = Patient4, yend = Patient4, color = Previous_treatment_class), size = 3) +
  geom_segment(aes(x = During_start, xend = During_end, y = Patient4, yend = Patient4, color = During_treatment_class), size = 3) +
  geom_segment(aes(x = Next_start, xend = Next_end, y = Patient4, yend = Patient4, color = Next_treatment_class), size = 3) +
  geom_point(aes(x = Biopsy, y = Patient4), size = 2) +
  theme_bw() +
  scale_color_d3("category20", name = "Previous treatment class") +
  scale_x_continuous(breaks = seq(-24,24, by = 12)) +
  facet_wrap(~Prev_Dur_BRAF, scales = "free_y") +
  xlab("Timeline (Months)") +
  ylab("Patient") +
  theme(axis.text = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 16, family = "URWHelvetica"),
        legend.text = element_text(size = 14, family = "URWHelvetica"),
        legend.title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica"))
ggsave("Figure 6B.pdf", width = 12, height = 8)  

#Figure 6C
ggplot(tpt, aes(x = Patient, fill = Macrophage_phenotype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cell_tp4, name = "Cell type")+
  theme_bw() +
  xlab("") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 16, family = "URWHelvetica"),
        legend.text = element_text(size = 14, family = "URWHelvetica"),
        legend.title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica")) +
  facet_wrap(~Previous_Treatment, scales = "free_x")
ggsave("Figure 6C.pdf", width = 12, height = 8)

#Figure 6D
ggplot(tpt[celltype_paper == "Macrophages"], aes(x = UMAP.1, y = UMAP.2, color = Macrophage_phenotype)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cell_mac, name = "Macrophage")+
  theme_bw() +
  facet_wrap(~Previous_Treatment) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  coord_cartesian(xlim = c(-8, -2), ylim = c(0, 4)) +
  theme(axis.text = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 16, family = "URWHelvetica"),
        legend.text = element_text(size = 14, family = "URWHelvetica"),
        legend.title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica"))
ggsave("Figure 6D.pdf", width = 12, height = 8)

#Figure 6E
tpt.cells <- tpt[, .N, by = c("Previous_Treatment", "Patient", "Macrophage_phenotype")]

tpt.cells[, Proportion := N/sum(N) * 100, by = Patient]

ggplot(tpt.cells[Macrophage_phenotype == "TTR macrophages"], aes(x = Previous_Treatment, y = Proportion)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  stat_compare_means(comparisons = list(c("BRAF + MEK treated", "Other treated")), size = 6) +
  ylab("% of TTR M\u03D5") +
  xlab("")  +
  theme(axis.text = element_text(size = 14, family = "URWHelvetica"),
        axis.title = element_text(size = 16, family = "URWHelvetica"),
        legend.text = element_text(size = 14, family = "URWHelvetica"),
        legend.title = element_text(size = 16, family = "URWHelvetica"),
        strip.text = element_text(size = 14, family = "URWHelvetica"))
ggsave("Figure 6E.pdf", width = 12, height = 8, device = cairo_pdf)
