library(ComplexHeatmap)
library(viridisLite)
library(colorRamp2)
library(BrepPhylo)
library(Platypus)
library(optparse)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
source("utils.R")

#read in user input:
option_list <- list(
  make_option(
    c("-r", "--reports"),
    type = "character",
    default = "./MiXCR",
    help = "Path to MiXCR output (IGH) [default %default]"
  ),
  make_option(
    c("-s", "--species"),
    type = "character",
    default = "hsa",
    help = "Species of samples [default %default]"
  ),
  make_option(
    c("-m", "--metadata"),
    type = "character",
    default = file.path("./samplelist.csv"),
    help = "Path to sample metadata table [default %default]"
  ),
  make_option(
    c("-o", "--output_dir"),
    type = "character",
    default = file.path("."),
    help = "Output directory [default %default]"
  )
)

arguments <- parse_args(OptionParser(option_list = option_list))
#arguments$metadata <- file.path("/ix/drajasundaram/drajasundaram/shared_bts76_dhr11/BCR_pipeline/BCR_nextflow/sampleslist.csv")
#arguments$reports  <- file.path("/ix/drajasundaram/drajasundaram/shared_bts76_dhr11/BCR_pipeline/bcr-test/test_results/MiXCR")
#arguments$output_dir <- file.path("/ix/drajasundaram/drajasundaram/shared_bts76_dhr11/BCR_pipeline/bcr-test/test_results/platypus")

#argument handling:
species <- arguments$species
valid_species <- c("hsa", "mmu")
stopifnot(species %in% valid_species)

outdir <- file.path(arguments$output_dir)
IGH_reports <-
  list.files(file.path(arguments$reports), full.names = T, pattern = "_IGH.tsv",recursive = T)

metadata <- read.csv(arguments$metadata)

stopifnot(nrow(metadata) > 0)

rownames(metadata) <- metadata$SampleID

if (!dir.exists(outdir)) {
  dir.create(outdir,showWarnings = T)
} else {
  print("Output directory already exists...")
}
setwd(outdir)

if (!dir.exists(file.path("./reformatted"))) {
  dir.create(file.path("./reformatted"))
}

#Step 1: read in MiXCR output & reformat
for (file in IGH_reports) {
  table <-
    read.table(
      file,
      header = T,
      sep = "\t",
      quote = "",
      encoding = "UTF-8",
      fill = TRUE
    )
  colnames(table) <- str_replace_all(colnames(table), "Imputed", "")
  write.table(
    table,
    file =
      paste0("./reformatted/",
            basename(file)),
    col.names = T,
    sep = "\t"
  )
}
IGH_reports <- list.files("./reformatted", full.names = T)

print("Metadata:")
print(metadata[basename(IGH_reports) %>% str_replace_all(., "_IGH.tsv", ""), "Group"])

#Step 2: build Platypus VGM matrix from reformatted MiXCR output
vgm <- custom_VDJ_bulk_to_vgm(
  IGH_reports[1:length(IGH_reports)],
  input.type = "MIXCR",
  integrate.MIXCR.output = TRUE,
  vgm.expanded = TRUE,
  clone.strategy = "cdr3.aa",
  group.id = metadata[basename(IGH_reports) %>% str_replace_all(., "_IGH.tsv", ""), "Group"],
  cell.type = "B.cell",
  batches = rep("1", length(IGH_reports)),
  best.match.only = TRUE
)

#create column of vgene + allele for CSR calculation
vgm[[1]]$VDJ_vgene_allele <- vgm[[1]]$VDJ_vgene
vgm[[1]]$VDJ_vgene <- vgm[[1]]$VDJ_vgene %>% gsub("\\*.*", "", .)
vgm_1 <- vgm[[1]]

#relabel sample ID based on metadata
vgm_1$sample_id <- metadata$SampleID[word(vgm_1$sample_id,2,sep="s") %>% as.numeric()]

#logo plot of CDR3 region across samples:
logoplot <- VDJ_logoplot_vector(
  cdr3.vector = vgm_1$VDJ_cdr3s_aa,
  length_cdr3 = "auto",
  seq_type = "auto"
)

logoplot+ggtitle("CDR3 Amino Acid Composition")
png("CDR3_logoplot.png",width = 700,height=500)
logoplot+ggtitle("CDR3 Amino Acid Composition")
dev.off()

#Step 3: quantifying Somatic Hyper-Mutation (SHM)
SHM_plots <-
  VDJ_plot_SHM.custom(VDJ = vgm_1,
                      group.by = "group_id",
                      quantile.label = 0.95)
shm.plot <- SHM_plots[[1]]+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("SHM per Group (VDJ)")
shm.data <- shm.plot$data[shm.plot$data$name =="VDJ_SHM",]

#calculate SHM rate:
shm.data$barcode <- NULL
shm.data$rate <-
  shm.data$value / str_length(vgm_1$VDJ_sequence_nt_trimmed)
shm.data$VDJ_length <- str_length(vgm_1$VDJ_sequence_nt_trimmed)
shm.data$CDR3_length <- str_length(vgm_1$VDJ_aaSeqCDR3)
shm.data$isotype <- vgm_1$VDJ_cgene
shm.data$isotype[shm.data$isotype == ""] <- NA
shm.data <- shm.data[complete.cases(shm.data),]

#plot SHM rate across all groups & BCR classes:
shm.plot <- shm.data %>% ggplot(aes(x = group, y = rate, color = group)) + ylab("Rate") +xlab(NULL)+ geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap( ~ isotype, scales = "free") + scale_y_continuous(labels = scales::percent) + ggtitle("Somatic Hypermutation")

png("SHM_plot.png",width=750,height=550)
shm.plot
dev.off()

#Step 4: CDR3 length distribution & overlap:
cdr3.histogram <-
  shm.data %>% ggplot(aes(x = CDR3_length,color=group)) +
  geom_density(position = "jitter")+ scale_y_continuous(labels = scales::percent)
cdr3.violin <-
  shm.data %>% ggplot(aes(x = group, y = CDR3_length, fill = group)) +
  geom_violin()
cdr3.plots <- cdr3.violin + ggtitle("CDR3 Length") | cdr3.histogram

png("CDR3_length_plots.png",width=750,height=550)
cdr3.plots
dev.off()

cdr3_overlap <- VDJ_overlap_heatmap(
  VDJ = vgm_1,
  feature.columns = c("VDJ_cdr3s_aa"),
  grouping.column = "group_id",
  axis.label.size = 20,
  pvalues.label.size = 12,
  add.barcode.table = T,
  plot.type = "ggplot"
)
cdr3_overlap

#Step 5: repertoire diversity estimation
#Shannon Evenness:
diversity_plt_Shannon <- VDJ_diversity(
  VDJ = vgm_1,
  feature.columns = c("VDJ_sequence_nt_trimmed"),
  grouping.column = "group_id",
  metric = c("shannonevenness"),
  subsample.to.same.n = T
)
#Gini-Simpson Index:
diversity_plt_Gini.Simpson <- VDJ_diversity(
  VDJ = vgm_1,
  feature.columns = c("VDJ_sequence_nt_trimmed"),
  grouping.column = "group_id",
  metric = c("ginisimpson"),
  subsample.to.same.n = T
)
#Jaccard similarity (matrix):
diversity_plt_Jaccard <- VDJ_diversity(
  VDJ = vgm_1,
  feature.columns = c("VDJ_vgene"),
  grouping.column = "group_id",
  metric = c("jaccard"),
  subsample.to.same.n = T
)

#save plots:
diversity_plt_Gini.Simpson | diversity_plt_Shannon |
  diversity_plt_Jaccard + ggtitle("Jaccard Similarity Index", subtitle = "VDJ V-gene") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle =
          element_text(hjust = 0.5))

png("repertoire_diversity_plots.png",width=750,height=500)
diversity_plt_Gini.Simpson | diversity_plt_Shannon |
diversity_plt_Jaccard + ggtitle("Jaccard Similarity Index", subtitle = "VDJ V-gene") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle =
          element_text(hjust = 0.5))
dev.off()

#Step 6: IGH V-gene usage heatmap
##top 100 variable genes across repertoires
vgm_gene_usage <-
  VGene_usage(vgm_1, 0, 100, FALSE) %>% select(c(Vgene, Percentage, Sample)) %>% pivot_wider(
    names_from = Sample,
    values_from = c(Vgene, Percentage),
    values_fn = list
  ) %>% unnest(cols = everything()) %>% data.frame

rownames(vgm_gene_usage) <- vgm_gene_usage[,1]
vgm_gene_usage <- vgm_gene_usage[,-1:-length(unique(metadata$SampleID))]
vgm_gene_usage[1:length(unique(metadata$SampleID))] %>% unique()
colnames(vgm_gene_usage) <- str_remove_all(colnames(vgm_gene_usage),"Percentage_")

vgm_annotation_df <-
  data.frame("Sample" = colnames(vgm_gene_usage),"Group" = metadata[colnames(vgm_gene_usage),"Group"])

#sample colors:
sample_cols <- custom_colors$discrete[1:length(vgm_annotation_df$Sample)]
names(sample_cols) <- vgm_annotation_df$Sample

#group colors:
group_cols <- colors_spanish[1:length(unique(vgm_annotation_df$Group))]
names(group_cols) <- vgm_annotation_df$Group %>% unique()

#scale gene util ratios:
vgm_gene_usage <- log(vgm_gene_usage+1)
col_fun = colorRamp2(c(0,round(max(vgm_gene_usage)/2) , max(vgm_gene_usage)), c("white","blue", "red"))

vgm_annotation_df <- vgm_annotation_df[order(vgm_annotation_df$Group),]
vgm_gene_usage <- vgm_gene_usage[,vgm_annotation_df$Sample]

#ComplexHeatmap:
vgm_gene_usage_heatmap <- Heatmap(
  column_title = "IGH-V Gene Usage", 
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  cluster_columns =  F,
  cluster_rows = T,
  vgm_gene_usage,
  show_row_dend = T,
  show_column_dend = F,
  name = "log(Usage %)",
  col = col_fun,
  top_annotation = columnAnnotation(
    df = vgm_annotation_df,
    col = list(Sample = sample_cols, Group=group_cols),
    show_annotation_name = TRUE
  ),
  column_labels = NULL,
  row_names_gp = gpar(fontsize = 10, lwd = 2),
  column_names_gp = gpar(fontsize = 10, hjust = 0.5),
  show_column_names = F,
  column_dend_side = "bottom",
  row_names_side = "left",
  row_names_centered = F
)
vgm_gene_usage_heatmap

png("IGHV_gene_usage_heatmap.png",width=500,height=650)
vgm_gene_usage_heatmap
dev.off()

#Step 7: Quantifying class-switch recombination (CSR) events w/ BrepPhylo
clones <- vgm_1
clones$CloneID <- str_remove_all(clones$clonotype_id, "clonotype")
clones$CloneID <- paste(clones$group_id, clones$CloneID, sep = "_")

clones <-
  clones[, c("CloneID",
             "sample_id",
             "group_id",
             "VDJ_cgene",
             "VDJ_vgene",
             "VDJ_vgene_allele",
             "VDJ_cdr3s_nt")]
clones$VDJ_cgene[clones$VDJ_cgene == ""] <- NA
clones <- clones[complete.cases(clones),]

clones$Subclass <- str_replace(clones$VDJ_cgene, "IGH", "")
clones$Class <- clones$Subclass %>% gsub('[0-9.]', '', .)
clones <- clones %>% data.frame()

#note: dnapars is available w/ the BrepPhylo package, but only works with Linux!
dnapars_executable <-
  system.file("exe/dnapars", package = "BrepPhylo")
outputFolder <- path.expand("./CSR_batchAnalysis")
dir.create(outputFolder, showWarnings = TRUE)

csr_species <- "Homo_sapiens"
if(species=="mmu"){
  csr_species <- "Mus_musculus"
}

#clonal lineage analysis:
batch_results <- doBatchCloneAnalysis(
  clones,
  outputFolder = outputFolder,
  species = csr_species,
  sequence_column = "VDJ_cdr3s_nt",
  IGHVgeneandallele_column = "VDJ_vgene_allele",
  plotFormat = "pdf",
  cloneID_column = "CloneID",
  label_column = "VDJ_cgene",
  phyloTreeType = "dnapars",
  phyloTreeOptions = list("executable" = dnapars_executable),
  useTempDir = FALSE,
  minCloneSize = 3
)
batch_summary <- getSummaryFromBatch(batch_results)
batch_summary.csr <- batch_summary$csr_events
batch_summary.csr$Group <- gsub('_[^_]+$','',batch_summary.csr$CloneID)

#summarize CSR events:
csr_summary <- summariseCSR(
  batch_summary.csr,
  dist_column = "distFromGermline",
  cloneID_column = "CloneID",
  summarise_variables = "Group"
)
#plot:
CSR_plot <- plotCSRsummary(csr_summary) + facet_wrap(~ Group) +
  scale_fill_viridis_c(name = "mean distance\nfrom germline") + ggtitle(label = "Class Switch Recombination", subtitle =
                                                                          "BrepPhylo")
png("CSR_events_plot.png",width=700,height=500)
CSR_plot
dev.off()
