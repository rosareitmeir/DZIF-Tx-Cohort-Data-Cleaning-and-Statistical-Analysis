### LongDat run for complete Cohort 

library(LongDat)

###################### functions


merge_taxalevels <- function(path, meta.data, num_meta=5){
  genus <- generate_LongDat_MasterTable("Genus", path, meta.data)
  family <- generate_LongDat_MasterTable("Family", path, meta.data)
  order<- generate_LongDat_MasterTable("Order", path, meta.data)
  class <- generate_LongDat_MasterTable("Class", path, meta.data)
  phylum <-generate_LongDat_MasterTable("Phylum", path, meta.data)
  
  alltaxalevels<- cbind(genus, family[,num_meta:ncol(family)], order[,num_meta:ncol(order)], class[,num_meta:ncol(class)], phylum[,num_meta:ncol(phylum)])
  
  return(alltaxalevels)
}


generate_LongDat_MasterTable <- function(rank, path, meta.input){
  
  #reading in 
  
  file <- list.files(path= path, pattern=paste("^", rank, ".*_n_0.tsv", sep=""))
  
  feature.input <- read.table(paste(path, file[1], sep=""), row.names=NULL, header=T, sep="\t")
  
  feature.input <- feature.input[,-2]
  
  # changing to row names 
  feature.input[,1] <- make.unique(feature.input[,1], sep = "_")
  
  row.names(feature.input) <- feature.input[,1]
  
  # convert and transpose
  
  feature.input[,2:ncol(feature.input)] <- lapply(feature.input[,2:ncol(feature.input)], as.integer )
  
  feature.input <- t(feature.input[, -1])
  
  # getting rid of OTUs not from the kingdom bacteria 
  feature.input<- as.data.frame(feature.input) %>% select( starts_with("Bacteria"))
  
  
  # subsetting and sorting it 
  feature.input <- feature.input[rownames(feature.input) %in% rownames(meta.input), ] 
  
  feature.input <- as.data.frame(feature.input[order(rownames(feature.input)), ])
  
  feature.input$SampleID <- rownames(feature.input)
  
  
  ## make master table 
  master.table <- make_master_table(metadata_table = meta.input,
                                    feature_table = feature.input , sample_ID = "SampleID", individual = "PatientID")
  
  
  return(master.table)
  
}






### taking the plotting function from longdat and modify it to add the taxa levels to the plot 
cuneiform_plot2<- function (result_table, x_axis_order = NULL, covariate_panel = TRUE, 
                            pos_color = "red", neg_color = "blue", panel_width = 4, title = "LongDat result cuneiform plot", 
                            title_size = 20, covariate_text_size = 5, x_label_size = 12, ylabels,
                            y_label_size = 12, legend_title_size = 10, legend_text_size = 10,taxa_size=10, direction=1) 
{
  if (missing(result_table)) {
    stop("Error! Necessary argument \"result_table\" is missing.")
  }
  result_table <- result_table
  sig_result <- result_table %>% dplyr::filter(Signal != "NS")
  if (nrow(sig_result) == 0) {
    stop("All results are non-significant! There is nothing to plot.")
  }
  sig_wide <- sig_result %>% dplyr::select(c(Feature, Signal, taxa,
                                             stringr::str_which(string = colnames(sig_result), pattern = "Effect")))
  Effect_wide <- sig_wide %>% dplyr::select(stringr::str_which(string = colnames(sig_wide), 
                                                               pattern = "EffectSize", negate = TRUE))
  EffectSize_wide <- sig_wide %>% dplyr::select(c(Feature, 
                                                  Signal, stringr::str_which(string = colnames(sig_wide), 
                                                                             pattern = "EffectSize")))
  Effect_long <- Effect_wide %>% tidyr::pivot_longer(stringr::str_which(string = colnames(Effect_wide), 
                                                                        pattern = "Effect"), names_to = "Effect_name", values_to = "Effect")
  EffectSize_long <- EffectSize_wide %>% tidyr::pivot_longer(stringr::str_which(string = colnames(EffectSize_wide), 
                                                                                pattern = "EffectSize"), names_to = "EffectSize_name", 
                                                             values_to = "EffectSize")
  All_long <- cbind(Effect_long, EffectSize_long[, c(3:4)])
  All_long$Alpha <- ifelse(All_long$Effect != "NS", yes = "Significant", 
                           no = "Non-significant")
  All_long$Shape <- ifelse(All_long$EffectSize > 0, yes = "24", 
                           no = ifelse(All_long$EffectSize < 0, yes = "25", no = "1"))
  if (is.null(x_axis_order)) {
    All_long$Effect_name <- All_long$Effect_name
  }
  else {
    All_long$Effect_name <- factor(All_long$Effect_name, 
                                   levels = x_axis_order)
  }
  All_long$Alpha <- factor(All_long$Alpha, levels = c("Significant", 
                                                      "Non-significant"))
  All_long$Shape <- factor(All_long$Shape, levels = c("1", 
                                                      "24", "25"))
  g1 <- ggplot2::ggplot(All_long, aes(x = Effect_name, y = reorder(Feature, EffectSize))) + 
    geom_point(aes(shape = Shape, fill = EffectSize, alpha = Alpha), 
               size = 3.5) + scale_shape_manual(values = c(1, 24, 
                                                           25), breaks = c("1", "24", "25"), labels = c("No change", 
                                                                                                        "Enriched", "Decreased"), name = "Effect", drop = FALSE) + 
    scale_fill_viridis_c(direction= direction) + scale_alpha_manual(breaks = c("Non-significant", 
                                                                               "Significant"), values = c(0.4, 1), drop = FALSE) + ggtitle(title) + 
    labs(fill = "Effect size", alpha = "Significance") +  
    facet_grid(taxa~., scales = "free", space = "free_y", switch = "y") +
    scale_y_discrete(breaks=ylabels$Feature, labels=ylabels$label)+
    
    theme_light() + theme(title = element_text(size = title_size),  strip.text = element_text(size=taxa_size),
                          axis.text.y.left = element_markdown(size = y_label_size), axis.title.x = element_blank(), 
                          axis.text.x = element_text(size = x_label_size), axis.title.y = element_blank(), 
                          legend.title = element_text(size = legend_title_size), 
                          legend.text = element_text(size = legend_text_size)) + 
    guides(fill = guide_colorbar(raster = FALSE, nbin = 30))
  if (covariate_panel == TRUE) {
    g2 <- ggplot2::ggplot(Effect_wide %>%mutate(Feature=factor(Feature, levels=rev(ylabels$Feature))), aes(x = "Covariate status", 
                                                                                                           y = Feature)) + geom_text(aes(label = Signal), size = covariate_text_size, 
                                                                                                                                     color = "gray30") + theme_light() + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.text.x = element_text(size = x_label_size), 
                                                                                                                                                                               axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
                                                                                                                                                                               axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                                                                                                                                                                               panel.grid.major.x = element_blank())
    final_plot <- (g1 | g2) + patchwork::plot_layout(guides = "collect", 
                                                     widths = c(panel_width, 1))
  }
  else {
    final_plot <- g1
  }
  print("Finished plotting successfully!")
  return(final_plot)
}

######################## 


### prepare meta input

##basic.info is the table containing the basic variables for the patient and samples: patientId, sampleid, dates of KT,Event,Sample, rejection status...
##basic.info is given as table in the supplement

#basic.info <- read_delim("../basic_info.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)


CC.meta.data <- basic.info %>% filter(origin == "Empfänger" & Transplantation_day_diffs > 0) %>% select(SampleID, PatientID, Transplantation_day_diffs, Event_day_diffs, rej) 
rownames(CC.meta.data) <- CC.meta.data$SampleID

# path to rarefied counts 

### prepare feature input 
CC.input <- merge_taxalevels(path, CC.meta.data, num_meta = 5)


## run longdat
cc.longdat.res <- longdat_cont(input = CC.input, data_type = "count",
                               test_var = "Transplantation_day_diffs", variable_col = 5, fac_var = c(1,4))

## longdat result
cc.res <- cc.longdat.res$Result_table %>% mutate(taxa= factor(case_when(Feature %in% names(genus) ~ "Genus",
                                                                        Feature %in% names(family) ~ "Family",
                                                                        Feature %in% names(class) ~ "Class",
                                                                        Feature %in% names(order) ~ "Order",  
                                                                        Feature %in% names(phylum) ~ "Phylum"), levels = c("Genus", "Family", "Order", "Class","Phylum"))) 


## plotting results as cuneiform plot 
feature.levels <-cc.res %>%  mutate(Feature=  make.unique( sub(".*\\;", "",Feature))) %>% arrange(taxa, Feature) %>% pull(Feature)


cuneiform_plot2(result_table = cc.longdat.res$Result_table %>% mutate(Feature= make.unique(sub(".*\\;", "",Feature)))) 

longdat_plt_total_cohort<-  cuneiform_plot2(result_table = cc.res %>% mutate(Feature= factor(make.unique(sub(".*\\;", "", gsub("_", " ", Feature))), levels= feature.levels)),  title="",covariate_panel=F,x_label_size = 9, 
                                            y_label_size = 9, legend_title_size = 9, legend_text_size = 9, taxa_size = 8) + 
  theme( strip.placement = "outside", 
         strip.background =element_rect(fill="white", color="grey"), 
         strip.text.y.left  = element_text(angle=0, color="black", size=10),
         axis.text.y = element_text(color="black"),
         # panel.border = element_rect(colour = "black", fill = NA)
         
  )
