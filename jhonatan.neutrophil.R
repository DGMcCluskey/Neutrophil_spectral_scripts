#packages
library(CATALYST)
library(harmony)
library(cowplot)
library(flowCore)
library(scater)
library(SingleCellExperiment)
library(viridis)
library(RColorBrewer)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(data.table)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggsignif)
library(rstatix)
library(multcomp)
library(scales)
library(ggthemes)
library(slingshot)
library(ggbeeswarm)
library(ggbreak)
library(DescTools)
library(pak)
library(gghighlight)
library(scattermore)
library(clustree)
library(gridExtra)
library(patchwork)
library(ggbump)


#for july 2025 data the files are in separate files for each repeat, so load them in separately and add a prefix then go from there
setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/Jhonatan neutrophil data/CD4+PMN_coculture_experiment/CD4+PMN_coculture_experiment/")
#load in fcs files (files exported from sony software to fcs files and then renamed for easier use prior)
fcs <- list.files(pattern = ".fcs$")
fcs

fcs <- fcs[3:22]
fs <- read.flowSet(fcs, transformation = F, truncate_max_range = F)

#removing FCS, SSC, CTV and the extra CD62L
all_cols <- colnames(fs)
all_cols
fs <- fs[,-c(1:4,7,20)]
colnames(fs)
#export remaining colnames to csv to create first of 2 metadata files needed for analysis (as I've done this analysis before I already have this excel so can skip)
cols <- colnames(fs)
cols <- as.data.frame(cols)
#add marker details in excel as second column
#write.csv(cols, file = "Marker.metadata.csv")
markers <- read.csv("Marker.metadata.csv")
markers

ld <- data.frame(fluorophore = "Zombie-NIR-A", antigen = "Live_dead")

markers <- rbind(markers,ld)
markers

#make sample metadata
#save list of file names and then edit in excel and reload in
names <- data.frame(fcs)
#write.csv(names, file = "sample.metadata.csv")
samples <- read.csv("sample.metadata.csv")

#samples <- samples %>% filter(!(grepl("1", time)))
#samples <- samples %>% filter(!(grepl("2", time)))
#samples <- samples %>% filter(!(grepl("DiHOME_0.1", treatment)))
#samples <- samples %>% filter(!(grepl("EpOME_0.1", treatment)))


colnames(samples)

#now all data is prepped, the singlecellexperiment object is made
#try a variety of cofactors (5,10,50,100,150) to see which transforms the data the best
#default cofactor is 5, which is used for cytof data, whilst it's often recommended to use 150 for flow cytometry data
setdiff(colnames(fs), markers$fluorophore)
colnames(fs)
fs
setdiff(fcs, samples$file_name)

sce.custom <- prepData(x = fs, panel = markers, md = samples, FACS = T, transform = T, 
  cofactor = c("AF700-A" = 1000, "BUV563-A" = 1000, "BV421-A" = 3000,"BV510-A" = 2000, 
  "BV605-A" = 1000, "BV650-A" = 1000, "BV711-A" = 1000, "FITC-A" = 1000, "PE-A" = 1000, 
  "PE-CF594-A" = 5000, "PE-Fire640-A" = 1000, "PerCP-Cy5.5-A" = 1000, "SparkUV-387-A" = 3000,
"Zombie-NIR-A" = 500), 
  panel_cols = list(channel = "fluorophore", antigen = "antigen"),
  md_cols = list(file = "file_name", id = "sample_id",
  factors = c("donor_id", "neut", "ratio")))

table(sce.custom$sample_id)

#downsample here if required
# Extract metadata (colData) to work with sample IDs
coldata <- colData(sce.custom)
# Create a function to downsample the cells for each sample
downsample_cells <- function(sce.custom, percent = 0.10) {
  # Create an empty list to store indices of the downsampled cells
  downsampled_cells <- list()
  # Loop through each sample_id and downsample cells
  for (sample in unique(coldata$sample_id)) {
    # Get the cells that belong to this sample
    sample_cells <- which(coldata$sample_id == sample)
    # Determine the number of cells to sample (20% of total cells in the sample)
    n_cells <- length(sample_cells)
    n_to_sample <- floor(n_cells * percent)
    # Randomly select 20% of the cells
    downsampled_cells[[sample]] <- sample(sample_cells, n_to_sample)
  }
  # Combine the downsampled cells into a single vector
  downsampled_indices <- unlist(downsampled_cells)
  # Subset the SingleCellExperiment object to keep only the downsampled cells
  sce_downsampled <- sce.custom[, downsampled_indices]
  return(sce_downsampled)
}
# Downsample 5% of the cells from each sample
sce_downsampled <- downsample_cells(sce.custom, percent = 0.10)

sum(table(sce_downsampled$time))

cocustom <- plotExprs(sce.custom, color_by = "donor_id")
cocustom
#using custom cofactor setup
sce <- sce.custom
rm(sce.custom)
rm(sce_downsampled)

#colours for plots
palette <- brewer.pal(n = 12, name = "Paired")

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


#plotting number of cells in each sample
CATALYST::plotCounts(sce, group_by = "neut", color_by = "ratio")+scale_fill_manual(values = col_vector)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  theme(axis.text = element_text(colour = "black", size = 12), 
        axis.title = element_text(size = 18, colour = "black"))+
  xlab(NULL)+ylab("Number of cells")

#PCA AND MDS PLOTS
pbMDS(sce, by = "sample_id", color_by = "donor_id", shape_by = "ratio", label_by = NULL)+
  scale_colour_manual(values = col_vector)+
  theme_classic()+geom_point(size = 5)+scale_colour_manual(values = c("indianred2", "steelblue2",
"forestgreen", "darkorchid2"))

pbMDS(sce, by = "sample_id", color_by = "patient_id", shape_by = "time", label_by = NULL)+
  theme_bw()+scale_colour_manual(values = col_vector)+geom_point(size = 5)

table(sce$sample_id)

#HARMONY BATCH CORRECTION - don't do until more than one patient
expdata <- sce@assays@data$exprs
meta <- sce@colData$donor_id
expdata <- t(expdata)
nrow(expdata)
nrow(meta)
harmony_output <- HarmonyMatrix(expdata, meta_data = meta, do_pca = F)
sce@assays@data$uncorrected <- sce@assays@data$exprs
sce@assays@data$exprs <- t(harmony_output)

saveRDS(sce, file = "neuts.pre-clustering.rds", compress = F)

#clustering
set.seed(61)
sce <- CATALYST::cluster(sce, features = NULL, xdim = 10, ydim = 10, maxK = 30, verbose = T, seed = 61)

delta <- delta_area(sce)
delta <- delta$data

ggplot(delta, aes(x = k, y = y))+geom_point(size = 5)+geom_line(linewidth = 1.5)+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 18), 
        axis.title = element_text(size = 18, colour = "black"))+
  xlab("K (number of clusters)")+ylab("Relative change in \n area under CDF curve")



sce <- runDR(sce, "UMAP", features = NULL)

saveRDS(sce, file = "neuts.clustered.rds", compress = F)
table(sce$neut, sce$donor_id)
getwd()
sce <- readRDS(file = "neuts.clustered.rds")
format(object.size(data), units = "auto")
#order groups so that easier to interpret


plotDR(sce, "UMAP", color_by = "meta10")+
  geom_scattermore()+theme_bw()

res <- c("meta10", "meta12", "meta15", "meta20")
umaps <- list()

for (i in res) {
umaps[[i]] <- plotDR(sce, "UMAP", color_by = i)+theme_bw()+
    geom_scattermore()
} 

p <- plot_grid(umaps$meta10, umaps$meta12, umaps$meta15, umaps$meta20)

ggsave(plot = p, filename = "umaps.k.comparison.png", dpi = 500,
       height = 10, width = 14)

nrs <- plotNRS(sce, color_by = "donor_id")
nrs <- nrs$data

p <- ggplot(nrs, aes(x = antigen, y = NRS))+geom_boxplot(outlier.shape = NA, fill = "darkseagreen")+theme_bw()+
  theme(axis.text = element_text(size = 12, colour = "black"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(plot = p, filename = "NRS.png", dpi = 500,
       height = 6, width = 10)


markers.input <- markers$antigen

plotDR(sce, color_by = markers.input)
plotDR(sce, color_by = "neut")+facet_wrap(~ratio)

plotExprs(sce, features = c("Live_dead", "CD4"), color_by = "donor_id")


#SECOND CLUSTERING=========================================================
#there are some small island clusters i believe represent dead cells, so will remove these

#whilst here I'll also exclude debris clusters if identifiable
plotPbExprs(sce, k = "meta10", color_by = "cluster_id", group_by = "cluster_id", facet_by = "antigen",
            features = NULL)

sce <- filterSCE(sce, k = "meta10", cluster_id %in% c(3,4,7,8,9))

plotDR(sce, "UMAP", color_by = "meta10")+
  geom_scattermore()+theme_bw()

sce <- CATALYST::cluster(sce, features = c("CD49d", "CD62L", "CD69", "CXCR4", "CD28",
"CCR7", "HLA-DR", "CD25", "CD45RA", "CD27", "CD45RO"),
xdim = 10, ydim = 10, maxK = 30, verbose = T, seed = 61)

delta <- delta_area(sce)
delta <- delta$data

ggplot(delta, aes(x = k, y = y))+geom_point(size = 5)+geom_line(linewidth = 1.5)+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 18), 
        axis.title = element_text(size = 18, colour = "black"))+
  xlab("K (number of clusters)")+ylab("Relative change in \n area under CDF curve")



sce <- runDR(sce, "UMAP", features = c("CD49d", "CD62L", "CD69", "CXCR4", "CD28",
"CCR7", "HLA-DR", "CD25", "CD45RA", "CD27", "CD45RO"))


res <- c("meta3", "meta5", "meta7", "meta10")
umaps <- list()

for (i in res) {
  umaps[[i]] <- plotDR(sce, "UMAP", color_by = i)+theme_bw()+
    geom_scattermore()
} 

p <- plot_grid(umaps$meta3, umaps$meta5, umaps$meta7, umaps$meta10)

ggsave(plot = p, filename = "umaps.k.comparison.2ndclustering.png", dpi = 500,
       height = 10, width = 14)


markers.input <- markers$antigen

p3 <- plotDR(sce, "UMAP", color_by = markers.input)+theme_bw()+
  geom_scattermore()+scale_colour_viridis(option = "B")+theme(strip.text = element_text(size = 16, colour = "black",
                                                                                        face = "bold"))
p3

ggsave(plot = p3, filename = "UMAP.marker.expression.2ndclustering.png", dpi = 500,
       height = 14, width = 18)

saveRDS(sce, file ="neuts_2nd_clustering.rds", compress = F)



#run clustree to view when certain clusters split at each k
sce$k5 <- cluster_ids(sce, "meta5")
sce$k6 <- cluster_ids(sce, "meta6")
sce$k7 <- cluster_ids(sce, "meta7")
sce$k8 <- cluster_ids(sce, "meta8")
sce$k9 <- cluster_ids(sce, "meta9")
sce$k10 <- cluster_ids(sce, "meta10")
sce$k11 <- cluster_ids(sce, "meta11")
sce$k12 <- cluster_ids(sce, "meta12")
sce$k13 <- cluster_ids(sce, "meta13")
sce$k14 <- cluster_ids(sce, "meta14")
sce$k15 <- cluster_ids(sce, "meta15")


clustree::clustree(sce, prefix = "k")


#check number of cells in each cluster
counts <- table(cluster_ids(sce, "meta5"))
counts


markers.input <- markers$antigen

hm <- plotExprHeatmap(sce, features = markers.input, 
                      by = "cluster_id", k = "meta7", bars = T, perc = T, scale = "last",
                      col_clust = F, hm_pal = rev(hcl.colors(10, "Reds")))

hm <- hm@matrix

hm2 <- pheatmap::pheatmap(hm, scale = "none", color = rev(hcl.colors(100, "RdYlBu")), border_color = "black",
                          cluster_cols = T, cluster_rows = T, treeheight_col = 0, cellheight = 20, cellwidth = 20)


#ANNOTATION========================================
#annotation - only done retrospectively once downstream plots inspected
merging.table <- data.frame("original" = c(1,2,3,4,5,6,7,8,9,10),
                            labels = c("CD16+CD11b+", "Early apoptotic", "Apoptotic (CD38+)", "Dead", "Dead",
                                       "Apoptotic (CD38+)", "CD16-", "Apoptotic (CD38+)", 
                                       "CD16+CD11b-", "Apoptotic (CD38-)"))

sce <- mergeClusters(sce, k = "meta10", table = merging.table, id = "labels", overwrite = T)

sce$labels <- cluster_ids(sce, k = "labels")

saveRDS(sce, file = "neuts.3rdclustering.annotated.rds", compress = F)

sce <- readRDS("neuts.3rdclustering.annotated.rds")

plotPbExprs(sce, features = "type", group_by = "cluster_id", k = "labels")

table(sce$)
#============================================================================================================================
#pulling out data to use for plotting
#first get data from catalyst object and make it into a dataframe
#will also add different resolutions (narrowed down based on clustree)
rm(neut.data)
ex <- assay(sce, "exprs")
neut.data <- data.frame(t(ex), sample_id = sce$sample_id, donor = sce$donor_id, 
                        neut = sce$neut, ratio = sce$ratio, 
                        k5 = sce$k5,
                        k6 = sce$k6, k7 = sce$k7,
                        k8 = sce$k8, k9 = sce$k9,
                        k10 = sce$k10, k11 = sce$k11,
                        k12 = sce$k12, k13 = sce$k13,
                        k14 = sce$k14, k15 = sce$k15)
#add umap dimensions
umap.coords <- data.frame(reducedDim(sce, "UMAP"))
neut.data$UMAP_1 <- umap.coords$UMAP1
neut.data$UMAP_2 <- umap.coords$UMAP2

write.csv(neut.data, file = "t.cell.df.csv")


old.cluster.palette <- c("indianred2", "cyan4", "steelblue2", "mediumvioletred", "blue", "mediumseagreen", "darkorchid2",
                     "darkred", "goldenrod2", "goldenrod4", "grey30", "forestgreen", "black", "red", "grey60", "hotpink")

cluster.palette <- c("orange", "indianred2", "mediumseagreen", "cyan3", "steelblue2", "blue", "darkred", "black")

label_data <- neut.data %>%
  group_by(k10) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p <- ggplot(neut.data, aes(x = UMAP_1, y = UMAP_2, colour=k10, fill = k10))+geom_scattermore(pointsize = 1)+
  scale_colour_manual(values = c(old.cluster.palette))+
  scale_fill_manual(values = c(old.cluster.palette))+
  geom_label_repel(data = label_data, aes(label = k10), colour = "white", force_pull = 0.1, force = 5)+
  theme_void()#+facet_wrap(~patient)
p

ggsave(plot = p, filename = "UMAP.annotated.k10.png", dpi = 500,
       height = 6, width = 10)


p <- ggplot(neut.data, aes(x = UMAP_1, y = UMAP_2, colour=k10, fill = k10))+geom_scattermore(pointsize = 2)+
  scale_colour_manual(values = c(old.cluster.palette))+
  scale_fill_manual(values = c(old.cluster.palette))+
  theme_bw()+facet_wrap(~donor+neut, ncol = 4)+theme(strip.text = element_text(size = 12, colour = "black", face = "bold"))
p

ggsave(plot = p, filename = "UMAP.annotated.patient.png", dpi = 500,
       height = 8, width = 10)


hm <- plotExprHeatmap(sce, by = "cluster_id", k = "meta10", bars = T, perc = T, row_dend = F, col_dend = F,
                      k_pal = cluster.palette)

hm

hm <- hm@matrix

hm2 <- pheatmap::pheatmap(hm, scale = "column", color = rev(hcl.colors(100, "RdBu")), border_color = "black",
                          cluster_cols = T, cluster_rows = T, treeheight_col = 0, 
                          cellheight = 30, cellwidth = 30)

ggsave(plot = p2, filename = "heatmap.annotated.png", dpi = 500,
       height = 8, width = 16)



neut.data2 <- pivot_longer(neut.data, names_to = "marker", values_to = "expression", 1:22)
ggplot(neut.data2, aes(x = k16, y = expression, fill = k16))+geom_violin(scale = "width")+
  facet_wrap(~marker, scales = "free_y", ncol = 2, strip.position = "right")+
  theme_bw()+scale_fill_manual(values = cluster.palette)



p <- ggplot(neut.data, aes(x = UMAP_1, y = UMAP_2, colour=k5))+geom_scattermore(pointsize = 2)+
  scale_colour_manual(values = c("indianred2", "steelblue2", "mediumseagreen", "orange", "hotpink"))+facet_wrap(~donor, ncol = 1)+
  theme_bw()+
  theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))
p

ggsave(plot = p, filename = "T_cells_per_donor_umap.png", dpi = 500,
       height = 9, width = 6)


#marker boxplot expression

k5_cols <- c("indianred2", "steelblue2", "mediumseagreen", "orange", "hotpink")

neut.data2 <- pivot_longer(neut.data, names_to = "marker", values_to = "expression", 1:14)

p <- ggplot(neut.data2, aes(x = k5, y = expression, fill = k5))+geom_boxplot(outlier.shape = NA)+
  theme_bw()+scale_fill_manual(values = k5_cols)+facet_wrap(~marker, scales = "free_y", ncol = 7)+
  theme(strip.text = element_text(size = 12, colour = "black", face = "bold"))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 14))
p

ggsave(plot = p, filename = "cluster.boxplot.expression.png", dpi = 500,
       height = 4, width = 12)

#violin
neut.data3 <- filter(neut.data2, !marker %in% c("Live_dead", "CD4", "CD15"))
p.violin <- ggplot(neut.data3, aes(x = k5, y = expression, fill = k5))+geom_violin(scale = "width", trim = T)+
  facet_grid(rows = vars(marker), scales = "free", switch = "y")+theme_cowplot()+
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], ""))+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0, size = 14))+
  scale_fill_manual(values = k5_cols)+
  theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1))
p.violin

ggsave(plot = p.violin, filename = "k5.violin.png", dpi = 500,
       height = 14, width = 8, bg = "white")

#marker expression over time

ggplot(filter(neut.data2, marker %in% c("CD62L", "CD11b")), 
       aes(x = time, y = expression, fill = time))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~labels+marker)+theme_bw()

#counts <- table(cluster_ids(sce_downsampled_filtered, "meta12"))
#counts
#the baseline untreated sample makes plotting difficult, so removing and plotting separately
baseline <- filter(neut.data, condition %in% c("A. Baseline (0h)"))
no.baseline <- filter(neut.data, !condition %in% c("A. Baseline (0h)"))

p1 <- ggplot(baseline, aes(x = UMAP_1, y = UMAP_2, colour=time))+geom_scattermore(pointsize = 2)+
  facet_wrap(~time+treatment, ncol = 4)+scale_colour_manual(values = c("goldenrod2"))+
  theme_bw()+theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "none")+xlab(NULL)
p1

p2 <- ggplot(no.baseline, aes(x = UMAP_1, y = UMAP_2, colour=time))+geom_scattermore(pointsize = 2)+
  facet_wrap(~time+treatment, ncol = 8)+scale_colour_manual(values = c("darkorchid2", "mediumseagreen"))+
  theme_bw()+theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "none")+ylab(NULL)+xlab(NULL)
p2

p3 <- ggplot(no.baseline, aes(x = UMAP_1, y = UMAP_2, colour=labels))+geom_scattermore(pointsize = 2)+
  scale_colour_manual(values = cluster.palette)+
  theme_bw()+theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "right")+guides(color = guide_legend(ncol = 1))+
  ggtitle("Clusters (k12)")
p3

p4 <- ggplot(no.baseline, aes(x = UMAP_1, y = UMAP_2, colour=time))+geom_scattermore(pointsize = 2)+
  facet_wrap(~treatment, ncol = 8)+scale_colour_manual(values = c("darkorchid2", "mediumseagreen"))+
  theme_bw()+theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "none")+ylab(NULL)
p4

combo <- plot_grid(p1,p2,p3,p4, ncol = 2, rel_widths = c(1,3), rel_heights = c(1.5,1))
combo
ggsave(plot = combo, filename = "UMAP.facet.time.treatment.png", dpi = 500,
       height = 8, width = 20)

#before annotating populations, need to find out what they are
#violin
violin.input <- pivot_longer(neut.data, names_to = "marker", values_to = "expression", 1:20)

violin.input <- filter(violin.input, marker %in% c("CD54", "CD101", "CD33", "CD49d",
                                                    "CD62L", "CD47", "CD38", "AnnexinV", "CD16",
                                                    "CD64", "CD11a", "CD182", "CD15", "CD11b", "CD45RA"))

p.violin <- ggplot(violin.input, aes(x = labels, y = expression, fill = marker))+geom_violin(scale = "width", trim = T)+
  facet_grid(rows = vars(marker), scales = "free", switch = "y")+theme_cowplot()+
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], ""))+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0))+
  scale_fill_manual(values = col_vector)
p.violin

ggsave(plot = p.violin, filename = "violin.plot.annotated.png", dpi = 500, bg = "white",
       height = 8, width = 10)



#proportions total
proportions <- neut.data %>% group_by(k5, neut, ratio, donor) %>% dplyr::count(k5, neut, ratio, donor)
proportions <- pivot_wider(proportions, names_from = k5, values_from = n)
proportions[is.na(proportions)] <- 0
proportions <- pivot_longer(proportions, names_to = "cluster", values_to = "number", 4:8)
proportions <- proportions

proportions$neut <- factor(proportions$neut, levels = c("none", "young", "mature", "aged"))

proportions$condition <- paste0(proportions$neut, "_", proportions$ratio)

proportions$condition <- factor(proportions$condition, levels = c("none_none",
"young_1", "young_2", "young_3", 
"mature_1", "mature_2", "mature_3",
"aged_1", "aged_2", "aged_3"))

p1 <- ggplot(proportions, aes(x = ratio, y = number, fill = cluster))+geom_col(position = "dodge")+theme_bw()+
  scale_fill_manual(values = k5_cols)+facet_wrap(~donor+neut, scales = "free_x", ncol = 4)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100000))+
  theme(axis.text = element_text(size = 12, colour = "black"))+xlab("Patient")+ylab("absolute number of cells")+
  theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1))+
  theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  theme(axis.title = element_text(size = 16, colour = "black", face = "bold"))+
  xlab("neutrophil:T cell ratio")+ylim(0,75000)
p1

p2 <- ggplot(proportions, aes(x = ratio, y = number, fill = cluster))+geom_col(position = "fill")+theme_bw()+
  scale_fill_manual(values = k5_cols)+facet_wrap(~donor+neut, scales = "free_x", ncol = 4)+
  theme(axis.text = element_text(size = 12, colour = "black"))+xlab(NULL)+ylab("% of cells")+
  theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1))+
  theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  theme(axis.title = element_text(size = 16, colour = "black", face = "bold"))+
  xlab("neutrophil:T cell ratio")
p2

combo <- p1+theme(legend.position = "none")+p2+plot_layout(ncol = 2)
combo
ggsave(plot = combo, filename = "t_cell_numbers_proportions.png", dpi = 500,
       height = 6, width = 12)
#stacked bar plot of proportions
#neut.1 <- filter(neut.data, patient %in% c("3"))
proportions <- neut.data %>% group_by(labels, time, treatment) %>% dplyr::count(labels)
proportions <- pivot_wider(proportions, names_from = labels, values_from = n)
proportions[is.na(proportions)] <- 0
proportions$total <- rowSums(proportions[, 3:9])
proportions <- pivot_longer(proportions, names_to = "cluster", values_to = "number", 3:9)
proportions$prop <- (proportions$number/proportions$total)*100

#proportions$cluster <- factor(proportions$cluster, levels = c("0", "1", "2", "3", "4", "5", "6", "7",
                                                              #"8", "9", "10", "11", "12"))


proportions$time <- factor(proportions$time, levels = c("0hr", "6hr", "22hr"))


ggplot(proportions, aes(x = time, y = prop, colour = treatment, group = interaction(treatment, patient)))+geom_bump()+
  facet_wrap(~cluster, scales = "free")+scale_colour_manual(values = cluster.palette)+theme_bw()

proportions$time_numeric <- as.numeric(gsub("hr", "", proportions$time))

ggplot(proportions, aes(x = time_numeric, y = prop, colour = treatment)) +
  geom_smooth(se = F, aes(group = treatment)) +  # one smoothed line per treatment
  facet_wrap(~cluster, scales = "free") +
  theme_bw() +
  scale_colour_manual(values = cluster.palette) +
  xlab("Time (hours)") +
  ylab("Proportion (%)")

summary_df <- proportions %>%
  group_by(time, treatment, cluster) %>%
  summarise(
    mean_prop = mean(prop),
    sd_prop = sd(prop),
    .groups = "drop"
  )

ggplot(summary_df, aes(x = time, y = mean_prop, colour = treatment, group = treatment)) +
  geom_smooth(se=F) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop), width = 0.2) +
  facet_wrap(~cluster, scales = "free") +
  theme_bw() +
  scale_colour_manual(values = cluster.palette) +
  ylab("Proportion (%)") +
  xlab("Time")



stacked <- ggplot(proportions, aes(x = treatment, y = number, fill = cluster))+geom_col(position = "fill")+
  theme_bw()+scale_fill_manual(values = cluster.palette)+ylab("% of neutrophils")+xlab("Time (hours)")+
  scale_y_continuous(expand = c(0,0))+theme(legend.position = 'right', 
                                            legend.key.spacing.y = unit(0.5, 'cm'))+
  facet_grid(~time, scales = "free", space = "free")+
  theme(axis.text.x = element_text(size = 11, colour = "black", face = "bold", angle = 45, hjust=1))+
  theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  xlab("Treatment")
stacked

ggsave(plot = stacked, filename = "annotated.proportions.time.stacked.plot.png", dpi = 500,
       height = 6, width = 12)


#doing delta proportions (difference of each treatment at each timepoint compared to untreated control)
#i'm going to perform this after calculting % change in proportion first
proportions <- as.data.frame(neut.data %>% group_by(labels, time, treatment, patient) %>% dplyr::count(labels))
proportions <- pivot_wider(proportions, names_from = labels, values_from = n)
proportions[is.na(proportions)] <- 0
proportions$total <- rowSums(proportions[4:10])
proportions <- pivot_longer(proportions, names_to = "cluster", values_to = "number", 4:10)
proportions$prop <- (proportions$number/proportions$total)*100
proportions <- proportions[, -c(4,6)]
proportions <- pivot_wider(proportions, names_from = "treatment", values_from = "prop")

proportions <- proportions %>% mutate(across( .cols = 4:11, .fns = ~ . - Non_treated,
                                      .names = "delta_{.col}"))

proportions <- proportions[, -c(4:11)]
proportions <- pivot_longer(proportions, names_to = "treatment", values_to = "delta.proportion", 4:11)

proportions <- proportions %>%
  mutate(delta.proportion = if_else(time == "0hr" & is.na(delta.proportion), 0, delta.proportion))

proportions$treatment <- factor(proportions$treatment, levels = c("delta_Non_treated", "delta_DMSO", "delta_MA", "delta_DMSO_MA",
                                                                  "delta_DiHOME", "delta_EpOME", "delta_GSK", "delta_EpOME_GSK"), 
                                labels = c("Untreated", "DMSO", "MA", "DMSO_MA", "DiHOME", "EpOME", "GSK", "EpOME_GSK"))
proportions$time <- factor(proportions$time, levels = c("0hr", "6hr", "22hr"))

#treatment.palette <- c("mediumseagreen", "indianred2", "royalblue2", "darkorchid2", "grey50", "black", "goldenrod4")

treatment.cols <- c("black", "grey50", "goldenrod", "orange", "mediumseagreen", "indianred2", "steelblue", "darkorchid2")


#line plot
delta.plot <- ggplot(proportions, aes(x = time, y = delta.proportion, colour = cluster, group = cluster))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  stat_summary(fun = mean, geom = "point", size = 3) +  # mean points
  #stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2) + 
  stat_summary(fun = mean, geom = "line", size = 2)+
  theme_bw()+scale_colour_manual(values = cluster.palette)+
  ylab("Proportion relative to untreated")+xlab("Time (hours)")+
  scale_y_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0.1,0.1))+
  theme(strip.text = element_text(size = 12, colour = "black", face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16))+
  facet_wrap(~treatment, scales = "fixed", ncol = 4)
delta.plot

ggsave(plot = delta.plot, filename = "delta.proportions.no.errorbar.png", dpi = 500,
       height = 8, width = 12)



#marker expression relative to non_treated over time (delta expression)
neut.data2 <- pivot_longer(neut.data, names_to = "marker", values_to = "expression", 1:22)
avg_expression <- neut.data2 %>%
  group_by(time, treatment, marker) %>%
  summarise(avg = mean(expression), .groups = "drop") %>%
  complete(time, treatment, marker, fill = list(avg = NA))

avg_expression <- pivot_wider(avg_expression, names_from = "treatment", values_from = "avg")

avg_expression <- avg_expression %>% mutate(across( .cols = 3:10, .fns = ~ . - Non_treated,
                                              .names = "delta_{.col}"))

avg_expression <- avg_expression[, -c(3:10)]
avg_expression <- pivot_longer(avg_expression, names_to = "treatment", values_to = "delta.expression", 3:10)

avg_expression <- avg_expression %>%
  mutate(delta.expression = if_else(time == "0hr" & is.na(delta.expression), 0, delta.expression))

unique(avg_expression$marker)

input <- filter(avg_expression, marker %in% c("CD62L", "CD11b", "AnnexinV", "Live_dead", "CD38", "CD184"))

delta.plot <- ggplot(input, 
                     aes(x = time, y = delta.expression, colour = treatment, group = treatment))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  stat_summary(fun = mean, geom = "point", size = 3) +  # mean points
  #stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2) + 
  stat_summary(fun = mean, geom = "line", size = 2)+
  theme_bw()+scale_colour_manual(values = cluster.palette)+
  ylab("nMFI relative to untreated")+xlab("Time (hours)")+
  scale_y_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0.1,0.1))+
  theme(strip.text = element_text(size = 12, colour = "black", face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16))+
  facet_wrap(~marker, scales = "free", ncol = 3)
delta.plot

#heatmap of marker expression
markers <- read.csv("Marker.metadata.csv")
markers <- markers$antigen

#remove cd45 as uneeded
markers <- markers[-12]
#markers <- markers[-6]

hm <- plotExprHeatmap(sce, features = markers, 
                      by = "cluster_id", k = "labels", bars = T, perc = T, scale = "never",
                      col_clust = F, hm_pal = rev(hcl.colors(10, "Reds")))


hm <- hm@matrix

hm2 <- pheatmap::pheatmap(hm, scale = "column", color = rev(hcl.colors(11, "Spectral")), border_color = "black",
                          cluster_cols = T, cluster_rows = T, treeheight_col = 0, 
                          cellheight = 20, cellwidth = 20, clustering_distance_rows = "correlation")


ggsave(plot = hm2, filename = "annotated.cluster.heatmap.png", dpi = 500,
       height = 4, width = 6)

#annotating clusters
neut.data$labels <- factor(neut.data$k10, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                           labels = c("Blood_Immature", "Blood_Immature", "Blood_Mature", "Blood_Mature", "Blood_Aged", 
                                      "Blister_Immature", "Blister_Mature", "Blood_Aged", "Blister_Aged", "Blister_Aged"))

table(neut.data$labels)

#annotated umaps
centroids <- neut.data %>%
  group_by(labels) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

palette <- brewer.pal(n = 8, name = "Set2")

p.umap <-ggplot(neut.data, aes(x = UMAP_1, y = UMAP_2, colour = labels))+
  geom_scattermore(pointsize = 1)+scale_colour_manual(values = palette)+
  geom_label_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = labels, fill = labels), 
                   size = 5, colour = "white", box.padding = 2, min.segment.length = 10)+
  theme_bw()+scale_fill_manual(values = palette)+theme(legend.position = "none")
p.umap

ggsave(plot = p.umap, filename = "annotated.UMAP.png", dpi = 500,
       height = 4, width = 5)

#umap tissue by time

p.umap.tissue <- ggplot(neut.data, aes(x = UMAP_1, y = UMAP_2, colour = tissue))+
  geom_scattermore()+facet_wrap(~condition, ncol = 4)+theme_bw()+
  scale_colour_manual(values = c("mediumseagreen", "darkorchid"))
p.umap.tissue


p.umap.density <- ggplot(neut.data, aes(x = UMAP_1, y = UMAP_2))+
  stat_density_2d_filled(aes(fill = after_stat(level)), geom = "polygon",
                         colour = "ivory", linewidth = 0.01,
                         contour_var = "ndensity")+
  theme_bw()+facet_wrap(~condition, ncol = 4)+scale_fill_viridis_d(option = "F", direction = -1)+
  scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
  theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))+
  labs(fill = "Density")

p.umap.density

ggsave(plot = p.umap.density, filename = "UMAP.density.png", dpi = 500,
       height = 5, width = 12)


#plotting umap, umap expression and stacked proportions
top.row <- plot_grid(p.umap, stacked, p.umap.density, ncol = 3, rel_widths = c(0.7,0.3,1))
bottom.row <- plot_grid(p.umap.expression, p.boxplots, ncol = 2)

combo <- plot_grid(top.row, bottom.row, ncol = 1)
combo

ggsave(plot = combo, filename = "umaps.and.proportions.png", dpi = 500,
       height = 12, width = 24)

#checking if the blister cells that are called as blood clusters are represented in all blister or
#are a one off artifact
#this is because it could just be blood contamination when collecting the blister
proportions <- neut.data %>% group_by(labels, time, patient, tissue, sample_id) %>% dplyr::count(labels)
proportions <- pivot_wider(proportions, names_from = labels, values_from = n)
proportions[is.na(proportions)] <- 0
proportions$total <- rowSums(proportions[, c(5:12)])
proportions <- pivot_longer(proportions, names_to = "labels", values_to = "number", cols = 5:12)
proportions$percent <- (proportions$number/proportions$total)*100


proportions <- filter(proportions, tissue %in% "Blister")
proportions <- filter(proportions, labels %in% c("Blood_1", "Blood_2", "Blood_3", "Blood_4", "Blood_5"))


ggplot(proportions, aes(x = time, y = percent, fill = labels))+
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75))+
  geom_point(shape = 21, position = position_dodge(width = 0.75))+
  theme_bw()+facet_wrap(~tissue, scales = "free")+scale_fill_manual(values = palette)

#proportions over time

proportions <- neut.data %>% group_by(labels, condition, patient, tissue, sample_id) %>% dplyr::count(labels)
proportions <- pivot_wider(proportions, names_from = labels, values_from = n)
proportions[is.na(proportions)] <- 0
proportions$total <- rowSums(proportions[, c(5:10)])
proportions <- pivot_longer(proportions, names_to = "labels", values_to = "number", cols = 5:10)
proportions$percent <- (proportions$number/proportions$total)*100

proportions$patient_tissue <- paste0(proportions$patient, "_", proportions$tissue)

proportions$labels <- factor(proportions$labels, levels = c("Blister_Immature", "Blister_Mature", "Blister_Aged",
                                                            "Blood_Immature", "Blood_Mature", "Blood_Aged"))

proportions.boxplots <- ggplot(proportions, aes(x = condition, y = percent, fill = condition))+geom_boxplot(outlier.shape = NA)+
  theme_bw()+facet_wrap(~labels, scales = "free_y", ncol = 3)+
  scale_fill_manual(values = c("indianred2", "darkred", "mediumvioletred", "hotpink",
                               "royalblue2", "blue", "cyan4", "mediumseagreen"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("% of neutrophils")+
  xlab(NULL)+theme(strip.text = element_text(size = 14, colour = "black", face = "bold"))
proportions.boxplots

ggsave(plot = proportions.boxplots, filename = "proportions.boxplots.png", dpi = 500,
       height = 6, width = 10)


#plotting proportion boxplots alongside density umaps

combo <- proportions.boxplots+p.umap.density+plot_layout(ncol = 1)
combo

ggsave(plot = combo, filename = "proportions.boxplot.density.umap.png", dpi = 500,
       height = 10, width = 18)



proportions$time <- as.character(proportions$time)
proportions$time <- as.numeric(proportions$time)


ggplot(proportions, aes(x = time, y = percent, colour = tissue, group = tissue))+geom_smooth(se=F)+
  theme_bw()+facet_wrap(~labels, scales = "free")+scale_colour_manual(values = palette)


p <- ggplot(proportions, aes(x = time, y = percent, group = tissue, colour = tissue, linetype = tissue))+geom_point()+
  geom_smooth(se=F, method = "loess")+
  theme_bw()+facet_wrap(~labels, scales = "free")+scale_colour_manual(values = c("seagreen3", "darkorchid"))+
  ylab("Percentage of total neutrophils")+xlab("Time (hours)")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  theme(axis.title = element_text(size = 14, colour = "black"))+
  theme(strip.text = element_text(face = "bold"))
p

ggsave(plot = p, filename = "netrophil.clusters.by.time.png", dpi = 500,
       height = 5, width = 10)

ggplot(neut.small, aes(x = UMAP_1, y = UMAP_2, colour=labels))+geom_point(pct = ".", size = 0.5)+
  scale_colour_manual(values = palette)+theme_bw()+
  facet_wrap(~labels)

table(neut.data$labels)

proportions$time_tissue <- paste(proportions$tissue,"_",proportions$time)

proportions$time_tissue <- factor(proportions$time_tissue, levels = c("Blister _ 3", "Blister _ 5", "Blister _ 7", "Blister _ 9", "Blister _ 24",
                                                                      "Blood _ 0", "Blood _ 6", "Blood _ 10", "Blood _ 25"),
                                  labels = c("Blister_3", "Blister_5", "Blister_7", "Blister_9", "Blister_24",
                                             "Blood_0", "Blood_6", "Blood10", "Blood_25"))

p <- ggplot(proportions, aes(x = sample_id, y = percent, fill = labels))+geom_col()+
  theme_bw()+scale_fill_manual(values = palette)+scale_y_continuous(expand = c(0,0))+
  facet_wrap(~time_tissue, scales = "free", ncol = 5)+xlab("Patient")+ylab("Percentage of neutrophils")+
  theme(axis.text.x = element_blank())+
  theme(legend.position = c(1,0), legend.justification = c(1,-0.1))+
  theme(legend.title = element_blank(), legend.text = element_text(size = 8))+
  theme(legend.key.size = unit(1,"line"))
p

ggsave(plot = p, filename = "proportions.stacked.plot.png", dpi = 500,
       height = 5, width = 14)


#pca on clusters to see if any from blister and blood are near each other
matrix <- neut.data[, c(1:14,16,28,17)]

#optional: remove cd62l (and CD45?) as this is such a strong determinant of blister and blood
#if you do this check index numbers as theres one column less
matrix <- matrix[, -3]
matrix <- pivot_longer(matrix, names_to = "marker", values_to = "expression", 1:13)
matrix <- matrix %>% group_by(labels, patient, marker, condition) %>% summarise(mean = mean(expression))
matrix <- pivot_wider(matrix, names_from = "marker", values_from = "mean")
#optional - remove data where blister or blood cells are called as clusters from the opposite tissue
# Function to extract source name
extract_source <- function(x) {
  sub("_.*|\\s.*", "", x)  # keep only part before "_" or before space
}

# Apply extraction
label_source <- extract_source(matrix$labels)
condition_source <- extract_source(matrix$condition)

# Keep rows where they match
matrix <- matrix[label_source == condition_source, ]


matrix$label_time <- paste0(matrix$labels, "-", matrix$condition, "-", matrix$patient)
matrix.names <- matrix$label_time
matrix <- matrix[, -c(1,2,3,17)]
rownames(matrix) <- matrix.names
pca.res <- prcomp(matrix)

data <- as.data.frame(pca.res$x)
data$cluster_time <- rownames(data)

data <- data %>% separate(cluster_time, into = c("cluster", "time", "patient"), sep = "-")

data <- data %>% mutate(tissue = ifelse(grepl("Blood", cluster), "Blood", "Blister"))

#data$time <- as.character(data$time)
#data$time <- as.numeric(data$time)

#data$time <- factor(data$time, levels = c(0,3,5,6,7,9,10,25))

p1 <- ggplot(data, aes(x = PC1, y = PC2, colour = cluster, shape = time))+geom_point(size = 3)+theme_bw()+
  scale_colour_manual(values = c("indianred2", "darkred", "mediumvioletred",
                                 "cyan2","cornflowerblue", "blue", "cyan4", "mediumseagreen"))+
  scale_shape_manual(values = c(20,15,17,18,21,22,24,23))
p1

ggsave(plot = p1, filename = "PCA.png", dpi = 500,
       height = 5, width = 6)

p2 <- ggplot(data, aes(x = PC2, colour = cluster))+geom_density()+theme_bw()+
  scale_colour_manual(values = c("indianred2", "darkred", "mediumvioletred",
                                 "cyan2","cornflowerblue", "blue", "cyan4", "mediumseagreen"))
p2

#baseline pca




#pca at patient level (all clusters)
matrix <- neut.data[, c(1:14,16,17)]
matrix <- pivot_longer(matrix, names_to = "marker", values_to = "expression", 1:14)
matrix <- matrix %>% group_by(patient, marker, condition) %>% summarise(mean = mean(expression))
matrix <- pivot_wider(matrix, names_from = "marker", values_from = "mean")
matrix$label_time <- paste0(matrix$condition, "-", matrix$patient)
matrix.names <- matrix$label_time
matrix <- matrix[, -c(1,2,17)]
rownames(matrix) <- matrix.names
pca.res <- prcomp(matrix)

data <- as.data.frame(pca.res$x)
data$cluster_time <- rownames(data)

data <- data %>% separate(cluster_time, into = c("time", "patient"), sep = "-")

data <- data %>% mutate(tissue = ifelse(grepl("Blood", time), "Blood", "Blister"))

p3 <- ggplot(data, aes(x = PC1, y = PC2, colour = tissue, shape = time))+geom_point(size = 3)+theme_bw()+
  scale_colour_manual(values = c("indianred2", "royalblue"))+
  scale_shape_manual(values = c(20,15,17,18,21,22,24,23))
p3


combo <- plot_grid(p3,p1,p2, ncol = 3, rel_widths = c(1,1,1))
combo

ggsave(plot = combo, filename = "pca.clusters.png", dpi = 500,
       height = 5, width = 18)

#note for wednesday 9th april:
#run just at cluster level not time
#some kind of correlation between clusters
#heatmap hierachical clustering?
#pseudotime on whole object just to see what happens

#heatmap with clustering
#use matrix and data made above for expression and metadata
meta <- data[, c(15:16)]

meta$time

meta$time <- factor(meta$time, levels = c("Blister (03h)", "Blister (05h)","Blister (07h)",
                                          "Blister (09h)", "Blood (00h)", "Blood (06h)",
                                          "Blood (10h)", "Blood (25h)"))

colnames(meta)[2] <- "condition"

#meta$time <- as.character(meta$time)
#meta$time <- as.numeric(meta$time)

time.col <- brewer.pal(n = 8, name = "RdBu")
meta.colours <- list(cluster =  c(Blister_1 = "royalblue2", Blister_2 = "blue", Blister_3 = "cyan4",
                                  Blood_1 = "indianred2", Blood_2 = "red", Blood_3 = "darkred", 
                                  Blood_4 = "hotpink", Blood_5 = "mediumvioletred"),
                     condition = c( `Blister (03h)` = time.col[5], `Blister (05h)` = time.col[6], 
                                    `Blister (07h)` = time.col[7], `Blister (09h)` = time.col[8],
                                    `Blood (00h)` = time.col[4], `Blood (06h)` = time.col[3], 
                                    `Blood (10h)` = time.col[2], `Blood (25h)` = time.col[1])) 

#heatmap 1 - clusters by time
pheatmap::pheatmap(t(matrix), scale = "row", cluster_cols = T, annotation_col = meta, 
                   border_color = "black",
                   color = brewer.pal(n = 11, name = "PuOr"), annotation_colors = meta.colours,
                   fontsize_col = 3, clustering_distance_cols = "euclidean")


#heatmap 2- just clusters
matrix <- neut.data[, c(1:14,28)]
matrix <- pivot_longer(matrix, names_to = "marker", values_to = "expression", 1:14)
matrix <- matrix %>% group_by(labels, marker) %>% summarise(mean = mean(expression))
matrix <- pivot_wider(matrix, names_from = "marker", values_from = "mean")
matrix.names <- matrix$labels
matrix <- matrix[, -c(1)]
rownames(matrix) <- matrix.names


pheatmap::pheatmap(t(matrix), scale = "row", cluster_cols = T, 
                   border_color = "black",
                   color = brewer.pal(n = 11, name = "PuOr"),
                   fontsize_col = 10, clustering_distance_cols = "correlation")

#heatmap 3 - clusters by time, but each patient represented
matrix <- neut.data[, c(1:14,28,19,16)]
matrix <- pivot_longer(matrix, names_to = "marker", values_to = "expression", 1:14)
matrix <- matrix %>% group_by(labels, marker, time, patient) %>% summarise(mean = mean(expression))
matrix <- pivot_wider(matrix, names_from = "marker", values_from = "mean")
matrix$label_time <- paste0(matrix$labels, "-", matrix$time)
matrix.names <- matrix$label_time
matrix <- matrix[, -c(1,2,17)]
rownames(matrix) <- matrix.names
pca.res <- prcomp(matrix)

data <- as.data.frame(pca.res$x)
data$cluster_time <- rownames(data)

data <- data %>% separate(cluster_time, into = c("cluster", "time"), sep = "-")

data <- data %>% mutate(tissue = ifelse(grepl("Blood", cluster), "Blood", "Blister"))


#marker expression over time
input <- pivot_longer(neut.data, names_to = "marker", values_to = "expression", 1:14)
input <- input %>% group_by(time, marker, tissue, patient) %>% summarise(mean = mean(expression))

p1 <- ggplot(input, aes(x = time, y = mean, colour = tissue, group = tissue))+geom_smooth(se=F, method = "glm")+
  theme_bw()+facet_wrap(~ marker, scales = "free", ncol = 7)+
  scale_colour_manual(values = c("mediumseagreen", "darkorchid"))+
  xlab("Time")+ylab("Expression")

p1

ggsave(plot = p1, filename = "marker.expression.over.time.png", dpi = 500,
       height = 5, width = 12)

#marker expression over time but boxplot
input <- pivot_longer(neut.data, names_to = "marker", values_to = "expression", 1:14)
input <- input %>% group_by(marker, tissue, patient, condition) %>% summarise(mean = mean(expression))

p1 <- ggplot(input, aes(x = condition, y = mean, colour = condition))+geom_boxplot()+
  theme_bw()+facet_wrap(~ marker, scales = "free", ncol = 7)+
  #scale_colour_manual(values = c("mediumseagreen", "darkorchid"))+
  xlab("Time")+ylab("Expression")

p1


#================================================================================================================================================
#pseudotime analysis
#subset 1000 cells from each sample
idx <- seq_len(ncol(sce))
idx <- split(idx, sce$sample_id)
n <- 1000
idx <- lapply(idx, \(.) {
  n <- min(n, length(.))
  sample(.,n)
})
sub <- sce[, unlist(idx)]
table(sub$sample_id)

#organise data
sling <- sub
expdata <- sling@assays@data$exprs
expdata <- t(expdata)
nrow(expdata)
clusterdata <- colData(sling)$labels
length(clusterdata)
lin <- getLineages(expdata, clusterdata, reducedDim = "PCA", omega = F, 
                   start.clus = c("CD16+CD11b-"))
curve <- getCurves(lin, approx_points = 150, 
                   allow.breaks = T, extend = "n", stretch =0.01)


pseudotime.lineages <- as.data.frame(slingPseudotime(curve))

average.pseudotime <- slingAvgPseudotime(curve)

average.pseudotime <- as.data.frame(average.pseudotime)

sling$avgpseudo <- average.pseudotime

SlingshotDataSet(curve)

sling$lin1 <- pseudotime.lineages$Lineage1
sling$lin2 <- pseudotime.lineages$Lineage2
#sling$lin3 <- pseudotime.lineages$Lineage3
#sling$lin4 <- pseudotime.lineages$Lineage4

data <- plotDR(sling, "UMAP", color_by = "labels")

data <- data$data
slingshot.data <- data

p1 <- ggplot(data, aes(x = x, y = y, colour = average.pseudotime))+
  geom_point(size = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  scale_colour_viridis(option = "B")+
  ggtitle("Average pseudotime")+xlab("UMAP_1")+ylab("UMAP_2")
p1

ggsave(plot = p1, filename = "average.pseudotime.umap.png", dpi = 500,
       height = 4, width = 6)

umap <- reducedDims(sling)$UMAP

umap_curve_embedding <- embedCurves(curve, umap) 

umap_curve_embedding

line1 <- as.data.frame(SlingshotDataSet(umap_curve_embedding)@curves[[1]]$s)
line2 <- as.data.frame(SlingshotDataSet(umap_curve_embedding)@curves[[2]]$s)
#line3 <- as.data.frame(SlingshotDataSet(umap_curve_embedding)@curves[[3]]$s)
#line4 <- as.data.frame(SlingshotDataSet(umap_curve_embedding)@curves[[4]]$s)
#umap with all lineage curves plotted
p2 <- ggplot(slingshot.data, aes(x = x, y = y, colour = labels))+
  geom_point(size = 0.5)+theme_bw()+
  geom_path(data = line1, aes(x = UMAP1, y = UMAP2), inherit.aes = F, linetype = "solid",
            size = 0.7, arrow = arrow(type = "closed", angle = 20, length = unit(0.4, "cm")))+
  geom_path(data = line2, aes(x = UMAP1, y = UMAP2), inherit.aes = F, linetype = "dashed",
            size = 0.7, arrow = arrow(type = "closed", angle = 20, length = unit(0.4, "cm")))+
  #geom_path(data = line3, aes(x = UMAP1, y = UMAP2), inherit.aes = F, linetype = "twodash", 
  #size = 0.7, arrow = arrow(type = "closed", angle = 20, length = unit(0.4, "cm")))+
  #geom_path(data = line4, aes(x = UMAP1, y = UMAP2), inherit.aes = F, linetype = "twodash", 
  #size = 1.2, arrow = arrow(type = "closed", angle = 20, length = unit(0.4, "cm")))+
  scale_colour_manual(values = cluster.palette)+
  guides(colour = guide_legend(override.aes = list(size = 3)))+
  theme(legend.text = element_text(size = 9))+xlab("UMAP_1")+ylab("UMAP_2")
p2

combo <- p1+p2+plot_layout(ncol = 2)

ggsave(plot = combo, filename = "pseudotime.umaps.png", dpi = 500,
       height = 4, width = 14)

#getting the mst backbone

centroids <- curve@metadata$mst
centroids <- as.data.frame(centroids)

#adding expression data from beginning to the slingshot.data object so that i can
#plot marker expression along pseudotime
slingshot <- cbind(slingshot.data, expdata)

slingshot$cellid <- rownames(slingshot)
slingshot <- slingshot[, -29]
slingshot <- gather(slingshot, key = "lineage", value = "pseudo", 27:28)
slingshot <- gather(slingshot, key = "marker", value = "expression", 29:50)

p1 <- ggplot(slingshot, aes(x = pseudo, y = expression, colour = lineage))+geom_smooth(se=F)+
  theme_bw()+facet_wrap(~ marker, scales = "free")+scale_colour_manual(values = c("indianred2", "orange"))+
  xlab("Pseudotime")+ylab("Expression")

p1

ggsave(plot = p1, filename = "pseudotime.marker.expression.png", dpi = 500,
       height = 5, width = 12)


p <- ggplot(slingshot, aes(x = pseudo, colour = treatment))+geom_density(alpha = 0.3)+
  theme_bw()+facet_wrap(~lineage, scales = "free", ncol = 1)+scale_colour_manual(values = cluster.palette)
p

ggsave(plot = p, filename = "pseudotime.population.desities.png", dpi = 500,
       height = 6, width = 14)

p <- ggplot(slingshot, aes(x = pseudo, colour = time))+geom_density(alpha = 0.3, size = 3)+
  theme_bw()+facet_wrap(~treatment+lineage, scales = "fixed", ncol = 2)+scale_colour_manual(values = cluster.palette)+
  theme(strip.text = element_text(size = 12, colour = "black", face = "bold"))
p

ggsave(plot = p, filename = "pseudotime.densities.time.treatment.png", dpi = 500,
       height = 12, width = 10)
