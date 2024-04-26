setwd("/fast/AG_Forslund/shared/DZIF_CKD_Holle/DZIF_Lotus2SLV138/")
#install.packages("devtools")
#devtools::install_github("grunwaldlab/metacoder")
#install.packages("BiocManager")
#BiocManager::install("phyloseq")

library(readr)
library(dplyr)
library(phyloseq)
library(metacoder)
library(svglite)
otu_data <- read_tsv("OTU.txt")
tax_data <- read_tsv("hiera_BLAST.txt")
biom=import_biom("OTU.biom")
obj=parse_qiime_biom("OTU.biom", class_regex = "(.*)", class_key = "taxon_name")
otu_data <- left_join(otu_data, tax_data,
                      by = c("OTU" = "OTU")) # identifies cols with shared IDs
DZIF_rej=read_tsv("../Rosa/tables/heattree_SampleID_rejectionstatus.tsv",col_types = "cc")
DZIF_rej$rej <- ifelse(DZIF_rej$rej == 1, "rejection", "non rejection")
#gsub 1: rejection, 0: non rejection
#obj$data$otu_table <- calc_obs_props(obj, "otu_table")
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table",
                                       cols = DZIF_rej$SampleID)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = DZIF_rej$rej, cols = DZIF_rej$SampleID)


obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = DZIF_rej$SampleID, # What columns of sample data to use
                                      groups = DZIF_rej$rej) # What category each sample is assigned to
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value,
                                               method = "fdr")
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.1] <- 0



## calculate cliff delta 
indices <- match(names(data.frame(obj$data$tax_abund))[-1], DZIF_rej$SampleID)

# Reorder the data frame based on the indices
DZIF_rej<- DZIF_rej[indices, ]
rejvec<- factor(DZIF_rej$rej)
names(rejvec) <- DZIF_rej$SampleID

cliff_deltas<- t(apply( data.frame(obj$data$tax_abund) ,1, function(x) { c( x[1], cliff.delta(x[-1], rejvec )$estimate)  }))

obj$data$diff_table$cliff_delta <- cliff_deltas[,2]

set.seed(1) # This makes the plot appear the same each time it is run 
#pdf("heat_tree_matrix.pdf", width = 10, height = 10)
heat_tree_matrix( filter_taxa( obj,!startsWith(taxon_names, "s" )) ,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = gsub(pattern = ".*_", replacement = "", taxon_names),
                 node_color = cliff_delta, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-0.55, 0.55), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-0.55, 0.55), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions")
#######         
#dev.off()

obj_genus<-   filter_taxa( obj,!startsWith(taxon_names, "s" ))

prep_heattree<- function(obj){
  obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table",
                                         cols = DZIF_rej$SampleID)
  
  obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = DZIF_rej$rej, cols = DZIF_rej$SampleID)
  
  obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                        cols = DZIF_rej$SampleID, # What columns of sample data to use
                                        groups = DZIF_rej$rej) # What category each sample is assigned to
  obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value,
                                                 method = "fdr")
  obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.1] <- 0
  
  
  cliff_deltas<- t(apply( data.frame(obj$data$tax_abund) ,1, function(x) { c( x[1], cliff.delta(x[-1], rejvec)$estimate)  }))
  
  obj$data$diff_table$cliff_delta <- cliff_deltas[,2]
  
  return(obj)
  
}


obj_genus <- prep_heattree(obj_genus)

obj_genus$data$diff_table$cliff_delta <- as.numeric(( obj_genus$data$diff_table$cliff_delta ))*(-1)

obj_genus$data$diff_table$cliff_delta[obj_genus$data$diff_table$wilcox_p_value > 0.1] <- 0


pdf("../Rosa/plots/heat_tree_matrix_genus_og_layout.pdf", width = 12, height = 12)


heat_tree_plot <-heat_tree_matrix(  obj_genus,
                   data = "diff_table",
                   node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                   node_label = gsub(pattern = ".*_", replacement = "", taxon_names),
                   node_color = cliff_delta, # A column from `obj$data$diff_table`
                   node_color_range = rev(brewer.pal(10,"PiYG")), # The built-in palette for diverging data
                   node_color_trans = "linear", # The default is scaled by circle area
                   node_color_interval = c(-0.5, 0.5), # The range of `log2_median_ratio` to display
                   edge_color_interval = c(-0.5, 0.5), # The range of `log2_median_ratio` to display
                   node_size_axis_label = "Number of OTUs",
                   node_color_axis_label = "Cliff Delta",
                   initial_layout = "re", layout = "da",
) + theme(legend.text = element_text(size=6), legend.title =  element_text(size=8))



ggsave("../Rosa/plots/heat_tree_matrix_genus_re_da_layout.pdf",heat_tree_matrix(  obj_genus,
                  data = "diff_table",
                  node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                  node_label = gsub(pattern = ".*_", replacement = "", taxon_names),
                  node_color = cliff_delta, # A column from `obj$data$diff_table`
                  node_color_range = rev(brewer.pal(10,"PiYG")), # The built-in palette for diverging data
                  node_color_trans = "linear", # The default is scaled by circle area
                  node_color_interval = c(-0.5, 0.5), # The range of `log2_median_ratio` to display
                  edge_color_interval = c(-0.5, 0.5), # The range of `log2_median_ratio` to display
                  node_size_axis_label = "Number of OTUs",
                  node_color_axis_label = "Cliff Delta",
                 initial_layout = "re", layout = "da",
                ) + theme(legend.text = element_text(size=6), legend.title =  element_text(size=8)),   width  =10, height=12)            

#node_color_range=brewer.pal(4,"PiYG"))
#######         
dev.off()

          
