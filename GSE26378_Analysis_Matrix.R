
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery"))


library(GEOquery)
library(limma)
library("tidyverse")

# Download and load the GEO dataset
geo_data <- getGEO("GSE26378", GSEMatrix = TRUE)
expression_set <- geo_data[[1]]  # Assuming the dataset has one object

# Extract expression matrix and sample data
expr_matrix <- exprs(expression_set)
sample_data <- pData(expression_set)


# Check available columns for group information
print(colnames(sample_data))

# Check which columns for identification of subclasses or control info (e.g. "health status:ch1")
print(unique(sample_data$"outcome:ch1"))  # This should give you "A", "B", "C", "Control" or similar

# Create new identification
sample_data$condition <- ifelse(("normal control"==sample_data$"disease state:ch1"),"Control","Sepsis")

sample_data$condition <- ifelse(("normal control"==sample_data$"disease state:ch1"),"Control",
                                ifelse(("septic shock patient"==sample_data$"disease state:ch1") & ("Survivor"==sample_data$"outcome:ch1"),
                                       "Sep_Surv","Seps_Nonsurv"))

# Create groups based on the subclass or control column
groups <- factor(sample_data$condition)

# Rename the levels to be syntactically valid (replace spaces with underscores)
#levels(groups) <- c("normal_control", "septic_shock_patient")

# Check the distribution of the groups to ensure it has more than one level (replicate)
print(table(groups))

#  Redefine the design matrix based on the subclasses
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
design

# Check if data are normalized
pdf(file="Boxplot_data_matrix.pdf")
boxplot(expr_matrix) #col=c(rep(2,10),rep(3,20),rep(4,11)))
dev.off()

# Normalize data using the control samples if matrix not norm
#expr_matrix_norm <- normalizeBetweenArrays(expr_matrix, method = "quantile")
#pdf(file="Boxplot_quantileNorm_data_matrix.pdf")
#boxplot(expr_matrix_norm, col=c(rep(2,10),rep(3,20),rep(4,11)))
#dev.off()

# Fit linear model and perform differential expression analysis
fit <- lmFit(expr_matrix, design)

# Using the new valid group names
contrast_matrix <- makeContrasts(
  sepS_vs_control = Sep_Surv - Control,
  sepNS_vs_control = Seps_Nonsurv - Control,
  levels = design
)

# Fit the contrast model
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract differential expressed genes
res_sepS_vs_control <- topTable(fit2, coef = "sepS_vs_control", number = Inf)
res_sepS_vs_control <- rownames_to_column(res_sepS_vs_control,var="ID")
head(res_sepS_vs_control)
colnames(res_sepS_vs_control) <- c("ID","LogFC_Sep_S","Mean","t_Sep_S","pval_Sep_S","padj_Sep_S","B_Sep_S")
head(res_sepS_vs_control)

res_sepNS_vs_control <- topTable(fit2, coef = "sepNS_vs_control", number = Inf)
res_sepNS_vs_control <- rownames_to_column(res_sepNS_vs_control,var="ID")
head(res_sepNS_vs_control)
colnames(res_sepNS_vs_control) <- c("ID","LogFC_Sep_NS","Mean","t_Sep_NS","pval_Sep_NS","padj_Sep_NS","B_Sep_NS")
head(res_sepNS_vs_control)

results_GSE26378 <- merge(res_sepS_vs_control[,c(1,3,2,4,6)], res_sepNS_vs_control[,c(1,2,4,6)], by = "ID")
head(results_GSE26378)

## Annotation via biomart
library(biomaRt)
mart.mm <- useMart(biomart="ensembl",dataset ="hsapiens_gene_ensembl")
#listDatasets(mart.mm)
# get probeset IDs
ID <- results_GSE26378$ID
# get gene symbols
listAttributes(mart.mm)
genes <- getBM(attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id","description","chromosome_name",
                              "start_position","end_position","entrezgene_id","external_gene_name","gene_biotype"),
               filters = "affy_hg_u133_plus_2", values = ID, mart = mart.mm)
head(genes)
#genes_uniq <- genes[!duplicated(genes$external_gene_name),]
genes_uniqID <- genes[!duplicated(genes$affy_hg_u133_plus_2),]

results_GSE26378 <- full_join (genes_uniqID,results_GSE26378,by=c("affy_hg_u133_plus_2"="ID"))
head(results_GSE26378)

write.csv(results_GSE26378,file="results_GSE26378.csv")


### For heatmap
library("pheatmap")

colnames(expr_matrix)
rownames(expr_matrix)
list_genes <- read.table(file="List_genes.txt",header = T)
list_genes

ID <- list_genes$ID
heatmap_matrix <- expr_matrix[(c(ID)),]
heatmap_matrix
sample_data$"condition"

heatmap_matrix <- as.data.frame(heatmap_matrix)

heatmap_matrix$Healthy <- rowMeans(heatmap_matrix[,c(7,8,10:13,40,49,63,64,78:82,95:100)])
heatmap_matrix$Sepsis <- rowMeans(heatmap_matrix[,c(1:6,14:18,20:22,25:33,35,36,38,39,41,42,44:47,50,52:55,57:62,
                                                    66:77,83:90,92:94,101:103)])
heatmap_matrix$Sepsis_Non_Surv <- rowMeans(heatmap_matrix[,c(9,19,23,24,34,37,43,48,51,56,65,91)])
head(heatmap_matrix,n=5)

heatmap_matrix <- heatmap_matrix[,c(104:106)]
heatmap_matrix
heatmap_matrix <- rownames_to_column(heatmap_matrix,var="ID")
heatmap_matrix <- left_join(list_genes,heatmap_matrix,by="ID")
heatmap_matrix <- heatmap_matrix[,c(2:5)]
heatmap_matrix <- column_to_rownames(heatmap_matrix,var="Symbol")
heatmap_matrix

breaklist=seq(-1.5,1.5,by=0.1)
pdf("Heatmap_GSE26378.pdf",width = 2, height = 3.5)
pheatmap(heatmap_matrix,scale="row",main = "GSE26378_Sepsis_Ctr_Children",
         show_rownames = T, show_colnames = T, fontsize = 4, fontsize_row = 4, fontsize_col = 4,
         cluster_rows = F,cluster_cols=F,clustering_method = "complete",
         clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",
         cutree_rows = NA, cutree_cols = NA, gaps_col= NULL,gaps_row =NULL, breaks = breaklist,
         color = colorRampPalette(c("deepskyblue","gray20","yellow")) (length(breaklist)),
         annotation_row = NA,cellwidth = 4, cellheight = 4,legend = T,treeheight_row= NA,border_color = NA,
         labels_col = NULL)
dev.off()











## Annotation via GPL platform
#  Download the GPL platform file from GEO to get probe-to-gene mapping
gpl_data <- getGEO("GPL570", AnnotGPL = TRUE)
gpl_table <- Table(gpl_data)

#  Check the column names in the platform data to find where gene symbols are stored
print(colnames(gpl_table))  # Find the column that contains gene symbols

#  Merge the differential expression data with the platform annotation
# Assuming the column 'ID' is for ILMN probe IDs and 'Gene symbol' is for gene names
top_table_septic_vs_control$ID <- rownames(top_table_septic_vs_control)

# Merge to add gene symbols for septic vs control
results_with_genes_septic_vs_control <- merge(res_sepsis_vs_control, gpl_table[, c("ID", "Gene symbol")], by = "ID", all.x = TRUE)

write.csv(results_with_genes_septic_vs_control, file = "septic_vs_control_with_genes.csv")
