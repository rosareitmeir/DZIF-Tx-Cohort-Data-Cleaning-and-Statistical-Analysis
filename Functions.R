### Functions

#### generate meta input for metadeconfoundR, choosing meta variables to test for 
generate_meta_input <- function(sample.data, split_var, morevar=F){
  if( morevar){
    sample.data<- data.frame(sample.data) %>% select( {{split_var}},  age, origin,  Event_day_diffs,  Transplantation_day_diffs,  
                                                      Geschlecht, sex_donor,
                                                      starts_with("age"), 
                                                      distance_trans_rej,
                                                      Kreaanstieg, Prednianderung, MinCreatinine,
                                                      BMI, BMI.donor, height.cm, height.donor.cm, weightAtTx.kg, weightAtTx.donor.kg,  
                                                      primary.kidney.condition, ICDcode1, bloodAB0, typeOfTx, 
                                                      viralinfec_before_sample, viralinfec_before_30d_sample,
                                                      IV_before_sample, IV_30d_before_sample,  total_antibioclass, total_antibioclass_30d, distinct_antibioclass, distinct_antibioclass_30d,  
                                                      antibac_before_sample, antibac_before_30d_sample,
                                                      ebv, ebv.donor, cmv, cmv.donor, 
                                                      Lab_Sample_Distance,  Albumin,  CRPFinal, CreatinineFinal, HbA1c,  HemoglobinFinal,  LeukocytesFinal,  LymphocytesAbsFinal2,   
                                                      MonocytesAbsFinal2,  TotalneutroAbsFinal2,  PhosphateFinal,  PTHFinal,  Urea,  UricacidFinal,  Urine.Erythrocytes, Urine.Leukocytes, Urine.Nitrite,Urine.ProteinFinal)
    
  }
  else{sample.data<- data.frame(sample.data) %>% select( {{split_var}}, Event_day_diffs, Transplantation_day_diffs, Geschlecht, age, sex_donor, starts_with("age_donor"))}
  
  sample.data<- sample.data[order(rownames(sample.data)), ]
  return(sample.data)
}


####### run metadeconfoundR 
run_metadeconfound <- function(rank, path, meta.input, abundance_filter=0.05, prevalence_filter=0.1){
  
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
  feature.input<- as.data.frame(feature.input) %>% dplyr::select( starts_with("Bacteria"))
  
  
  # subsetting and sorting it 
  feature.input <- feature.input[row.names(feature.input) %in% rownames(meta.input), ] 
  
  
  ## prevalence and abundance filter 
  feature.input <- feature.input[, colSums(feature.input, na.rm=T)/nrow(feature.input) >= abundance_filter]
  
  prevalences <- c()
  for ( i in 1: ncol(feature.input)){
    prevalences[i] <- sum(feature.input[,i] >0, na.rm = T)/nrow(meta.input)
    
  }
  
  
  feature.input <- feature.input[, prevalences >= prevalence_filter]
  
  
  ## sorting it 
  feature.input <- as.data.frame(feature.input[order(rownames(feature.input)), ])
  
  
  ## run  metadeconfound 
  
  # checking 
  if (all(rownames(feature.input) == rownames(meta.input)) & all(order(rownames(feature.input)) == order(rownames(meta.input)) )){
    print(meta.input)
    output<- MetaDeconfound(featureMat = feature.input, metaMat = as.data.frame(meta.input), nnodes=16, returnLong = T)
    
    return(output)
    
  }
  
  else{
    return(all(rownames(feature.input) == rownames(meta.input)))
  }
  
}

##### workflow metadeconfoundR
#- input: taxonomic level, meta.input generate by the above function, path to rarefaction files, splitting variable is optional to generate a volcano plot 
#- executes metadeconfoundR, saves results as tsv, heatmap (with short names) and if desired also a volcano plot 
#- returns: metadeconfoundR output (table with results)

workflow_metadeconfound <- function(rank,  meta.input, path, outpath, split.var = "", abundance_filter=0.05, prevalence_filter=0.1){
  
  output <- run_metadeconfound(rank, path, meta.input, abundance_filter,prevalence_filter)
  signature <- output %>% filter( Qs <= 0.1 & Ds != 0.01) %>% mutate_at( vars(Ps,Qs,Ds),
                                                                         function(x){format(as.numeric(as.character(x)), scientific=F)} ) 
  
  outfile <- paste(outpath, paste(rank, "significant.tsv", sep="_"), sep= "/")
  
  outheatmap <- paste(outpath, paste(rank, "heatmap", sep="_"), sep= "/")
  
  if (split.var != ""){
    
    volcano  <- create_volcanoplot(output, split.var)
    
    outvolcano <- paste(outpath, paste(rank, split.var,  "volcano.png", sep="_"), sep= "/")
    
    ggsave(outvolcano, volcano)
  }
  
  
  write.table( signature, outfile, quote = F, col.names = T, row.names = F, sep="\t" )
  
  
  ggsave(paste(outheatmap, "short_names.png", sep="_") , BuildHeatmap(output, featureNames = data.frame(output$feature, make.unique(sub(".*\\;", "", output$feature)))))
  
  ggsave(paste(outheatmap, ".png", sep=""), BuildHeatmap(output), width=10, height=14 )
  
  
  
  return(output)
}


### get count table for a taxa level 
get_count_table <- function(path, rank, filter=T, abundance_filter=0.05, prevalence_filter=0.1, samplelist){
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
  feature.input<- as.data.frame(feature.input) %>% dplyr::select( starts_with("Bacteria"))
  
  
  # subsetting and sorting it 
  feature.input <- feature.input[row.names(feature.input) %in% samplelist, ] 
  
  if(filter){
    ## prevalence and abundance filter 
    feature.input <- feature.input[, colSums(feature.input, na.rm=T)/nrow(feature.input) >= abundance_filter]
    
    prevalences <- c()
    for ( i in 1: ncol(feature.input)){
      prevalences[i] <- sum(feature.input[,i] >0, na.rm = T)/length(samplelist)
      
    }
    
    
    feature.input <- feature.input[, prevalences >= prevalence_filter]
    
  }
  
  
  ## sorting it 
  feature.input <- as.data.frame(feature.input[order(rownames(feature.input)), ])
  
  return(feature.input)
}
### volcano genus
generate_volcano_plot <- function(metadeconfound.output,annotation=T, highfeatures=NULL, high_n=6, absolute=T, ds=0.05, qs=0.1 , colors=c("blue", "red","black"),text_size=2){
  volcano_df <- metadeconfound.output %>% filter(metaVariable== "rej") %>% mutate(OTU_type = case_when(Ds > ds & Qs < qs ~ "increased",
                                                                                                       Ds < (-1*ds) & Qs < qs  ~ "depleted",
                                                                                                       TRUE ~ "not significant"), feature= sub(".*\\;", "", feature))   
  
  if (is.null(highfeatures) & absolute){
    highfeatures <- volcano_df %>% filter(Qs <0.1) %>% top_n(high_n,abs(Ds)) %>% pull(feature)}
  if (is.null(highfeatures) & !absolute){
    highfeatures <- volcano_df %>% filter(Qs <0.1) %>% top_n(high_n/2,Ds) %>% pull(feature)
    highfeatures <- c(highfeatures, volcano_df %>% filter(Qs <0.1) %>% slice_min(Ds, n=high_n/2) %>% pull(feature) )
  }
  
  volcano_df <- volcano_df %>% mutate(label = case_when( feature %in% highfeatures ~ feature, T ~ NA)) 
  
  if (annotation){
    volcano_plot <- ggplot(data=volcano_df, aes(x=Ds, y=-log10(Qs) , label=label) ) + geom_point(aes(col= OTU_type), size=2.5) + theme(legend.title= element_blank())+ xlab("Cliff's Delta") + ylab("-log10(adjusted p-value)") + geom_text_repel(size=text_size) +scale_color_manual(values=colors) +
      geom_vline(xintercept=c((-1*ds), ds), col="black", linetype="dashed") +
      geom_hline(yintercept=-log10(qs), col="black",linetype="dashed") 
  }
  else{
    volcano_plot <- ggplot(data=volcano_df, aes(x=Ds, y=-log10(Qs)) ) + geom_point(aes(col= OTU_type), size=2.5) + theme(legend.title= element_blank())+ xlab("Cliff's Delta") + ylab("-log10(adjusted p-value)") +scale_color_manual(values=colors) +
      geom_vline(xintercept=c((-1*ds), ds), col="black", linetype="dashed") +
      geom_hline(yintercept=-log10(qs), col="black",linetype="dashed")
  }
  return(volcano_plot +  theme_bw() + 
           theme(
             axis.line = element_line(colour = "darkgrey"),
             axis.ticks.y  = element_line(color = "black"),
             axis.text.y = element_text(color="black", size=11),
             
             axis.text.x = element_text(color="black", size=11),
             axis.ticks.x = element_line(color = "black"),
             
             axis.title.x  =element_text(size=13),
             axis.title.y  =element_text(size=13),
             
             legend.title = element_blank(),
             legend.text =  element_text(color="black", size=11),
             legend.position = "top",
             
             panel.border = element_blank(),
             panel.background = element_blank()
           ))
}

### violin alpha plot 
create_alpha_violin <- function(input, split.var, metric, ylabel=metric){
  
  sample_sizes <- input %>% group_by(!!!syms(split.var)) %>% summarise(num=n())
  
  plot <- ggplot(input %>%  left_join(sample_sizes) %>%  mutate(myaxis = paste0(!!!syms(split.var), "\n", "n=", num)), aes_string(x= "myaxis", y= metric) ) + 
    geom_violin(aes_string(fill = split.var))  +
    geom_boxplot(color="black",  alpha=0, width=0.45) +
    geom_jitter() + theme(legend.position = "none") +
    xlab("") + ylab(ylabel) +
    stat_compare_means(size=5,label.x.npc = "center") +
    theme(
      axis.text.x = element_text(color="black", size=13),
      axis.text.y = element_text(color="black", size=12),
      axis.ticks = element_line(color = "black"),
      axis.title=element_text(size=13)
    )
  return(plot)
}



############ Prettify Heatmap

prettify_heatmap <- function(output, out.path,  mvs= metavariables_plot, stars_size=5, q_cutoff=0.1, keepMeta=NULL, change_order=c()){
  
  heatmap <- BuildHeatmap2(output, 
                           metaVariableNames = mvs,
                           keepMeta=keepMeta,
                           featureNames = data.frame(output$feature, make.unique(sub(".*\\;", "", output$feature))),
                           stars_size= stars_size,  q_cutoff= q_cutoff) + theme(
                             plot.subtitle = element_blank(),
                             plot.title = element_blank(),
                             axis.text.x = element_text(color="black", size=11),
                             axis.text.y = element_text(color="black", size=11),
                             axis.ticks = element_line(color = "black"),
                             axis.title=element_text(size=14)) + xlab("") + ylab("")
  
  if (length(change_order) >1){
    
    heatmap <- heatmap + scale_x_discrete(limits = change_order)
  }   
  
  return(heatmap)
  
}



BuildHeatmap2 <- function (metaDeconfOutput, q_cutoff = 0.1, d_cutoff = 0.01, 
                           cuneiform = FALSE, coloring = 0, showConfounded = TRUE, intermedData = FALSE, stars_size=2,
                           featureNames = NULL, metaVariableNames = NULL, d_range = "fit", 
                           d_col = c("blue", "white", "red"), keepMeta = NULL, keepFeature = NULL, 
                           trusted = c("OK_sd", "OK_nc", "OK_d", "AD")) 
{
  if (length(d_col) != 3) {
    stop("wrong number of colors in d_col!\nSupply colors for c(min, middle, max)!")
  }
  if (!(d_range %in% c("fit", "full"))) {
    stop("d_range must be either \"fit\" or \"full\"!")
  }
  if (length(trusted) == 0) {
    stop("\"trusted\" must contain at least one trusted status label")
  }
  allLables <- c("OK_sd", "OK_nc", "OK_d", "AD", "NS")
  notTrusted <- allLables[!(allLables %in% trusted)]
  if (class(metaDeconfOutput) == "list") {
    effectSize <- reshape2::melt(data = metaDeconfOutput$Ds, 
                                 varnames = c("feature", "metaVariable"), value.name = "Ds")
    effectSize$Ds[effectSize$Ds == Inf] <- 0
    effectSize$Ds[is.na(effectSize$Ds)] <- 0
    fdr <- reshape2::melt(data = metaDeconfOutput$Qs, varnames = c("feature", 
                                                                   "metaVariable"), value.name = "Qs")
    fdr$Qs[is.na(fdr$Qs)] <- 1
    status <- reshape2::melt(data = metaDeconfOutput$status, 
                             varnames = c("feature", "metaVariable"), value.name = "status")
  }
  else {
    effectSize <- metaDeconfOutput[, c("feature", "metaVariable", 
                                       "Ds")]
    fdr <- metaDeconfOutput[, c("feature", "metaVariable", 
                                "Qs")]
    status <- metaDeconfOutput[, c("feature", "metaVariable", 
                                   "status")]
    effectSize$Ds[effectSize$Ds == Inf] <- 0
    effectSize$Ds[is.na(effectSize$Ds)] <- 0
    fdr$Qs[is.na(fdr$Qs)] <- 1
  }
  insignificant <- unlist(lapply(strsplit(as.character(status$status), 
                                          split = ", "), function(l) any(l %in% notTrusted)))
  trueDeconf <- unlist(lapply(strsplit(as.character(status$status), 
                                       split = ", "), function(l) any(l %in% trusted)))
  effectSize$stars <- cut(fdr$Qs, breaks = c(-Inf, 0.001, 0.01, 
                                             0.1, Inf), label = c("***", "**", "*", ""))
  effectSize$stars[insignificant] <- ""
  effectSize$status <- trueDeconf
  effectSize$stars <- as.character(effectSize$stars)
  for (m in seq_along(effectSize$stars)) {
    if (!effectSize$status[m] && length(effectSize$stars[m]) > 
        0) {
      if (showConfounded) {
        effectSize$stars[m] <- gsub("*", "°", effectSize$stars[m], 
                                    fixed = TRUE)
      }
      else {
        effectSize$stars[m] <- ""
      }
    }
  }
  effectSize$insignificant <- insignificant
  effectSize$trueDeconf <- !trueDeconf
  if (coloring == 1) {
    effectSize$Ds[effectSize$insignificant] <- 1e-06
  }
  if (coloring == 2) {
    effectSize$Ds[effectSize$trueDeconf] <- 1e-06
  }
  remove_metavariables <- vector()
  for (i in unique(effectSize$metaVariable)) {
    aMetaVariable <- fdr[fdr$metaVariable == i, ]
    aMetaVariableD <- effectSize[effectSize$metaVariable == 
                                   i, ]
    if (sum(na.exclude(abs(aMetaVariable$Qs)) > q_cutoff) == 
        length(na.exclude(aMetaVariable$Qs)) || sum(na.exclude(abs(aMetaVariableD$Ds)) < 
                                                    d_cutoff) == length(na.exclude(aMetaVariableD$Ds)) || 
        all(aMetaVariableD$stars == "")) {
      remove_metavariables <- c(remove_metavariables, i)
    }
  }
  if (!is.null(keepMeta)) {
    remove_metavariables <- remove_metavariables[!(remove_metavariables %in% 
                                                     keepMeta)]
  }
  effectSize <- effectSize[!(effectSize$metaVariable %in% remove_metavariables), 
  ]
  remove <- vector()
  for (i in unique(effectSize$feature)) {
    aGenus <- fdr[fdr$feature == i, ]
    aGenusD <- effectSize[effectSize$feature == i, ]
    if (sum(na.exclude(abs(aGenus$Qs)) > q_cutoff) == length(na.exclude(aGenus$Qs)) | 
        sum(na.exclude(abs(aGenusD$Ds)) < d_cutoff) == length(na.exclude(aGenusD$Ds)) | 
        all(aGenusD$stars == "")) {
      remove <- c(remove, i)
    }
  }
  remove <- as.vector(unique(remove))
  if (length(remove) > 0) {
    if (!is.null(keepFeature)) {
      remove <- remove[!(remove %in% keepFeature)]
    }
    if (length(remove) > 0) {
      effectSize <- effectSize[!(effectSize$feature %in% 
                                   remove), ]
    }
  }
  effectSize <- droplevels(effectSize)
  if (length(unique(effectSize$metaVariable)) == 0) {
    stop("No associations pass current q_cutoff and/or d_cutoff filters!")
  }
  eff_cast <- reshape2::dcast(effectSize, effectSize[[1]] ~ 
                                metaVariable, value.var = "Ds")
  rownames(eff_cast) <- eff_cast[[1]]
  eff_cast[[1]] <- NULL
  if (nrow(eff_cast) > 1) {
    ord <- hclust(dist(eff_cast, method = "euclidean"), method = "ward.D")$order
    effectSize$feature <- factor(as.factor(effectSize$feature), 
                                 levels = levels(as.factor(effectSize$feature))[ord])
  }
  if (ncol(eff_cast) > 1) {
    eff_cast <- scale(t(eff_cast))
    ord2 <- hclust(dist(eff_cast, method = "euclidean"), 
                   method = "ward.D")$order
    effectSize$metaVariable <- droplevels(effectSize$metaVariable)
    effectSize$metaVariable <- factor(as.factor(effectSize$metaVariable), 
                                      levels = levels(as.factor(effectSize$metaVariable))[ord2])
  }
  effectSize$featureNames <- effectSize$feature
  effectSize$metaVariableNames <- effectSize$metaVariable
  if (!is.null(featureNames)) {
    if (class(featureNames)[[1]] != "data.frame") {
      warning("class(featureNames) was coerced to \"data.frame\"")
      featureNames <- as.data.frame(featureNames)
    }
    if (length(unique(featureNames[[2]])) != length(featureNames[[2]])) {
      featureNames[[2]] <- make.unique(featureNames[[2]])
      warning("non-unique human-readable feature names where made unique using base::make.unique")
    }
    map = stats::setNames(featureNames[[2]], featureNames[[1]])
    effectSize$featureNames <- map[as.vector(effectSize$feature)]
    effectSize$featureNames <- factor(as.factor(effectSize$featureNames), 
                                      levels = map[levels(effectSize$feature)])
  }
  if (!is.null(metaVariableNames)) {
    if (class(metaVariableNames)[[1]] != "data.frame") {
      warning("class(metaVariableNames) was coerced to \"data.frame\"")
      metaVariableNames <- as.data.frame(metaVariableNames)
    }
    if (length(unique(metaVariableNames[[2]])) != length(metaVariableNames[[2]])) {
      metaVariableNames[[2]] <- make.unique(metaVariableNames[[2]])
      warning("non-unique human-readable metaVariable names where made unique using base::make.unique")
    }
    map = stats::setNames(metaVariableNames[[2]], metaVariableNames[[1]])
    effectSize$metaVariableNames <- map[as.vector(effectSize$metaVariable)]
    effectSize$metaVariableNames <- factor(as.factor(effectSize$metaVariableNames), 
                                           levels = map[levels(effectSize$metaVariable)])
  }
  if (intermedData == TRUE) {
    return(effectSize)
  }
  lowerLim <- min(effectSize$Ds)
  upperLim <- max(effectSize$Ds)
  if (d_range == "full") {
    lowerLim <- -1
    upperLim <- 1
  }
  signifCol <- c("gray45", "black")
  signifMeaning <- c("confounded", "deconfounded")
  if (all(effectSize$status)) {
    signifCol <- c("black")
    signifMeaning <- c("deconfounded")
  }
  if (cuneiform) {
    divShapes <- c()
    divShapesMeaning <- c()
    signs <- unique(sign(effectSize$Ds))
    if (-1 %in% signs) {
      divShapes <- c(divShapes, 25)
      divShapesMeaning <- c(divShapesMeaning, "negative association")
    }
    if (0 %in% signs) {
      divShapes <- c(divShapes, 23)
      divShapesMeaning <- c(divShapesMeaning, "no association/no data")
    }
    if (1 %in% signs) {
      divShapes <- c(divShapes, 24)
      divShapesMeaning <- c(divShapesMeaning, "positive association")
    }
    heatmapGGplot <- ggplot(effectSize, aes(x = metaVariable, 
                                            y = feature)) + geom_point(aes(fill = Ds, shape = as.factor(sign(Ds)), 
                                                                           color = status)) + scale_shape_manual(name = "Direction", 
                                                                                                                 values = divShapes, labels = divShapesMeaning) + 
      scale_fill_gradient2(low = d_col[1], mid = d_col[2], 
                           high = d_col[3], midpoint = 0, guide = guide_colorbar(raster = F), 
                           limits = c(lowerLim, upperLim)) + scale_color_manual(name = "Confounding status", 
                                                                                values = signifCol, labels = signifMeaning) + guides(color = guide_legend(override.aes = list(shape = 24))) + 
      theme_classic() + theme(axis.text.x = element_text(size = 7, 
                                                         angle = 90, hjust = 1, vjust = 0.3), axis.text.y = element_text(size = 7, 
                                                                                                                         angle = 0, hjust = 1, vjust = 0.35), plot.title.position = "plot", 
                              plot.title = element_text(hjust = 0)) + labs(title = "Summarizing cuneiform plot", 
                                                                           x = "Metadata variables", y = "Omics features")
  }
  else {
    heatmapGGplot <- ggplot(effectSize, aes(x = metaVariableNames, 
                                            y = featureNames)) + geom_tile(aes(fill = Ds)) + 
      scale_fill_gradient2(low = d_col[1], mid = d_col[2], 
                           high = d_col[3], midpoint = 0, guide = guide_colorbar(raster = F), 
                           limits = c(lowerLim, upperLim)) + geom_text(aes(label = stars, 
                                                                           colour = status), size = stars_size, key_glyph = "point") + 
      scale_color_manual(name = "Confounding status", values = signifCol, 
                         labels = signifMeaning) + guides(color = guide_legend(override.aes = list(shape = c(1, 
                                                                                                             8)))) + theme_classic() + theme(axis.text.x = element_text(size = 7, 
                                                                                                                                                                        angle = 90, hjust = 1, vjust = 0.3), axis.text.y = element_text(size = 7, 
                                                                                                                                                                                                                                        angle = 0, hjust = 1, vjust = 0.35), plot.title.position = "plot", 
                                                                                                                                             plot.title = element_text(hjust = 0), plot.subtitle = element_text(size = 8)) + 
      labs(title = "Summarizing heatmap", subtitle = "FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ", 
           x = "Metadata variables", y = "Omics features")
  }
  return(heatmapGGplot)
}

### other plots for metadeconfoundR results
create_volcanoplot<- function(df, split.var, q_cutoff=0.1, d_cutoff=0.01){
  
  
  df<- df %>%  filter(metaVariable== split.var) %>% mutate(OTU_type = case_when(Ds > d_cutoff & Qs < q_cutoff ~ "up",
                                                                                Ds < -d_cutoff & Qs < q_cutoff  ~ "down",
                                                                                TRUE ~ "not significant"), feature= sub(".*\\;", "", feature))   
  
  
  return(ggplot(data=df, aes(x=Ds, y=-log10(Qs) , col= OTU_type) ) + geom_point() + theme(legend.title= element_blank())+ xlab("Cliff's Delta") + ylab("-log10(p_adj)")
  )
  
  
}


create_enhance_volcano <- function(output, split.var, label=F, q_cutoff=0.1, d_cutoff=0.01, xlim=c(-1,1), ylim=c(0,2)){
  output<- filter(output, metaVariable == split.var)
  
  
  rownames(output)<- make.unique(sub(".*\\;", "", output$feature))
  volcano <- EnhancedVolcano(output,
                             lab = ifelse(label, rownames(output), NA) ,
                             x = "Ds",
                             y = 'Qs',
                             pCutoff = q_cutoff,
                             FCcutoff = d_cutoff,
                             xlim = xlim,
                             ylim=ylim,
                             xlab="Cliff's Delta",
                             legendLabels = c("NS", "Cliff's Delta", "p-value", "p-value and Cliff's Delta"),
                             title="", subtitle="")
  return(volcano)
  
  
  
}

create_files <- function(output, outpath, rank ){
  
  signature <- output %>% filter( Qs <= 0.1 & Ds != 0.01)  %>% mutate_at( vars(Ps,Qs,Ds),
                                                                          function(x){format(as.numeric(as.character(x)), scientific=F)} ) 
  
  outfile <- paste(outpath, paste(rank, "significant.tsv", sep="_"), sep= "/")
  write.table( signature, outfile, quote = F, col.names = T, row.names = F, sep="\t" )
  
  
}

create_heatmap_pdf <- function(output, outpath, rank){
  outfile <- paste(outpath, paste(rank, "_heatmap.pdf", sep="_"), sep= "/")
  
  hm<- BuildHeatmap(output, featureNames = data.frame(output$feature, make.unique(sub(".*\\;", "", output$feature))))
  ggsave(outfile, hm, width=10, height=14)
}


