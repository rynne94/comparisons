# set the working directoyy
setwd('~/Documents/Jen/Sequencing/OneDrive_1_01-03-2021')

# installs required
install.packages('tidyverse')
library(tidyverse)

# read the tables of fold changes and p-values
dex_total <- read_csv('results_total_control_DEXcomparison.csv')
head(dex_total)
dex_total <- rename(dex_total, gene_ID = ...1)
dex_80S <-  read_csv('results_control_DEXcomparison_outlier_80S.csv')
dex_80S <- rename(dex_80S, 'gene_ID' = ...1) 
dex_poly <- read_csv('results_control_DEXcomparison_poly.csv')
dex_poly <- rename(dex_poly, 'gene_ID' = ...1)

kd_total <- read_csv('results_total_control_KD_noDEXcomparison.csv')
kd_total <- rename(kd_total, 'gene_ID' = ...1)
kd_80S <- read_csv('results_noDEXcomparison_outlier_80S.csv')
kd_80S <- rename(kd_80S, 'gene_ID' = ...1)
kd_poly <- read_csv('results_noDEX_CKDcomparison_poly.csv')
kd_poly <- rename(kd_poly, 'gene_ID' = ...1)

kd_total_dex <- read_csv('results_jen_tables/results_knockdown_DEXcomparison_total_annotated.csv')
head(kd_total_dex)
kd_total_dex <- rename(kd_total_dex, 'hgnc_symbol' = 'gene_name')
kd_80S_dex <- read_csv('results_jen_tables/results_knockdown_DEXcomparison_80S_annotated.csv')
head(kd_80S_dex)
kd_80S_dex <- rename(kd_80S_dex, 'hgnc_symbol' = 'gene_name')
kd_poly_dex <- read_csv('results_jen_tables/results_knockdown_DEXcomparison_poly_annotated.csv')
head(kd_poly_dex)
kd_poly_dex <- rename(kd_poly_dex, 'hgnc_symbol' = 'gene_name')

# read the table of annotation of genes
annotation <- read_csv('~/Documents/Jen/Sequencing/OneDrive_1_01-03-2021/annotation.csv')
head(annotation)
annotation <- rename(annotation, 'gene_ID' = ensembl_gene_id)

dex_total <- full_join(dex_total, annotation, by = 'gene_ID')
head(dex_total)


# select genes that are statistically significant
dex_total_p <- dex_total %>% 
  filter(padj <0.1) %>% 
  inner_join(annotation, by = 'gene_ID') %>% 
  dplyr::select(gene_ID,hgnc_symbol, log2FoldChange )

dex_80S_p <- dex_80S %>% 
  filter(padj <0.1)%>% 
  inner_join(annotation, by = 'gene_ID') %>% 
  dplyr::select(gene_ID,hgnc_symbol, log2FoldChange)

dex_poly_p <- dex_poly %>% 
  filter(padj<0.1)%>% 
  inner_join(annotation, by = 'gene_ID') %>% 
  dplyr::select(gene_ID,hgnc_symbol, log2FoldChange)

kd_total_p <- kd_total %>% 
  filter(padj <0.1) %>% 
  inner_join(annotation, by = 'gene_ID') %>% 
  dplyr::select(gene_ID,hgnc_symbol, log2FoldChange )

kd_80S_p <- kd_80S %>% 
  filter(padj <0.1)%>% 
  inner_join(annotation, by = 'gene_ID') %>% 
  dplyr::select(gene_ID,hgnc_symbol, log2FoldChange)

kd_poly_p <- kd_poly %>% 
  filter(padj<0.1)%>% 
  inner_join(annotation, by = 'gene_ID') %>% 
  dplyr::select(gene_ID,hgnc_symbol, log2FoldChange)

kd_total_dex_p <- kd_total_dex %>% 
  filter(padj <0.1) %>% 
  dplyr::select(hgnc_symbol, log2FoldChange )

kd_80S_dex_p <- kd_80S_dex %>% 
  filter(padj <0.1)%>% 
  dplyr::select(hgnc_symbol, log2FoldChange)

kd_poly_dex_p <- kd_poly_dex %>% 
  filter(padj<0.1)%>% 
  dplyr::select(hgnc_symbol, log2FoldChange)

head(kd_poly_dex_p)

#write for biomart analysis
write_csv(dex_total_p, 'dexa_total_significant.csv')
write_csv(dex_80S_p, 'dexa_80S_significant.csv')
write_csv(dex_poly_p, 'dexa_poly_significant.csv')
write_csv(kd_total_p, 'kd_total_significant.csv')
write_csv(kd_80S_p, 'kd_80S_significant.csv')
write_csv(kd_poly_p, 'kd_poly_significant.csv')

# COMPARISONS 

# mostafa et al

mostafa <- read_csv('most_data.csv')
mostafa <- rename(mostafa, hgnc_symbol = 'Gene')

#compare with GMD candidates
gmd <-  read_csv('gmd_candidates.csv')
head(gmd) # we retain the column gene
gmd <- gmd %>% 
  dplyr::select(gene)
head(gmd)

gmd_total <- inner_join(gmd, dex_total_p, by = "hgnc_symbol")
gmd_80S <- inner_join(gmd, dex_80S_p, by = "hgnc_symbol")
gmd_poly <- inner_join(gmd, dex_poly_p, by = "hgnc_symbol")

# gr-ip comparison
gr_ip <- read_csv('gr_ip.csv')
gr_ip_total <-  inner_join(gr_ip, dex_total_p, by = "hgnc_symbol")
gr_ip_80S <- inner_join(gr_ip, dex_80S_p, by = "hgnc_symbol")
gr_ip_poly <- inner_join(gr_ip, dex_poly_p, by = "hgnc_symbol")

# top mRNAs
top <- read_csv('top_RNAs.csv')
top_poly <-inner_join(top, dex_poly_p, by = "hgnc_symbol")

# compare with ARE
are <- read_csv('~/Documents/Jen/Sequencing/ARE_3UTR.csv')
head(are)
are <- are %>% 
  dplyr::select(gene_name) %>% 
  rename('hgnc_symbol' = 'gene_name')
dim(are)

are_total_kd <- inner_join(are, kd_total_p, by = "hgnc_symbol")
are_80S_kd <- inner_join(are, kd_80S_p, by = "hgnc_symbol")
are_poly_kd <- inner_join(are, kd_poly_p, by = "hgnc_symbol")

dim(are_poly_dex_kd %>% 
     filter(log2FoldChange > 0))

are_total_dex_kd <- inner_join(are, kd_total_dex_p, by = "hgnc_symbol")
are_80S_dex_kd <- inner_join(are, kd_80S_dex_p, by = "hgnc_symbol")
are_poly_dex_kd <- inner_join(are, kd_poly_dex_p, by = "hgnc_symbol")

are_total <- inner_join(are, dex_total_p, by = "hgnc_symbol")
are_80S <- inner_join(are, dex_80S_p, by = "hgnc_symbol")
are_poly <- inner_join(are, dex_poly_p, by = "hgnc_symbol")

# compare zhang
zhang <- read_csv('zhang.csv')
head(zhang)

zhang_total <- inner_join(zhang, kd_total_p, by = "hgnc_symbol")
zhang_80S <- inner_join(zhang, kd_80S_p, by = "hgnc_symbol")
zhang_poly <- inner_join(zhang, kd_poly_p, by = "hgnc_symbol")

head(kd_total_dex_p)

zhang_total_kd <- inner_join(zhang, kd_total_dex_p, by = "hgnc_symbol")
head(zhang_total_kd)
zhang_80S_kd <- inner_join(zhang, kd_80S_dex_p, by = "hgnc_symbol")
zhang_poly_kd <- inner_join(zhang, kd_poly_dex_p, by = "hgnc_symbol")
head(zhang_80S_kd)
dim(zhang_80S)
dim(zhang_80S_kd %>% 
  filter(log2FoldChange > 0))

# epithelial

epi <- read_csv('epithelial_genes.csv')
head(epi)

epi_total_kd <- inner_join(epi, kd_total_p, by = "hgnc_symbol")
epi_80S_kd <- inner_join(epi, kd_80S_p, by = "hgnc_symbol")
epi_poly_kd <- inner_join(epi, kd_poly_p, by = "hgnc_symbol")

epi_total_dex_kd <- inner_join(epi, kd_total_dex_p, by = "hgnc_symbol")
epi_80S_dex_kd <- inner_join(epi, kd_80S_dex_p, by = "hgnc_symbol")
epi_poly_dex_kd <- inner_join(epi, kd_poly_dex_p, by = "hgnc_symbol")

epi_total <- inner_join(epi, dex_total_p, by = "hgnc_symbol")
epi_80S <- inner_join(epi, dex_80S_p, by = "hgnc_symbol")
epi_poly <- inner_join(epi, dex_poly_p, by = "hgnc_symbol")

# venn diagrams

install.packages('VennDiagram')
library(VennDiagram)
library(RColorBrewer)


set1 <- dex_total_p$gene_ID
set2 <- dex_80S_p$gene_ID
set3 <- dex_poly_p$gene_ID
myCol <- brewer.pal(3, "PiYG")
myCol2 <- c('')
png('venn_dex.png')
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Total", "80S", "Poly"), 
  filename = 'venn_dex.png',
  output = TRUE,
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  col=c("#3446eb", '#eb34e8', '#fde725ff'),
  fill = c(alpha("#3446eb",0.1), alpha('#eb34e8',0.1), alpha('#fde725ff',0.1)),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 180),
  cat.dist = c(-0.02, -0.02, -0.02),
  cat.fontfamily = "sans",
  rotation = 1
)

dev.off()
