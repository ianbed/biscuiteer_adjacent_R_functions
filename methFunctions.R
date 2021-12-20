
# EXAMPLE USAGE OF HeatmapFromBed Function
# myCaption <- paste('UQCRH','chr1_46303335-46304488')
# 
# ordered <- dplyr::arrange(meta,histotype,cellType)$sample
# new_order <- match(ordered,meta$sample)
# 
# heatmapSampleAnno <- rowAnnotation(
#   histotype = meta$histotype[new_order],
#   cellType = meta$cellType[new_order],
#   MIR200c_avg_beta = meta$MIR200c_avg_beta[new_order],
#   col = col_arg_in_vivo
# )
# 
# myMat <- HeatmapFromBed(myChr = 'chr1',
#                         myRanges = IRanges(46303335,width = 1153),
#                         nameArg = myCaption,
#                         new_order = new_order, # sample order
#                         hm_sample_anno = heatmapSampleAnno # row anno in this case
# )

HeatmapFromBed <- function(myChr,myRanges,nameArg,new_order,hm_sample_anno){
  
  bs_list <- generate_or_load_bsseq_list(
    path='../Biscuit_Snakemake_Workflow/analysis/pileup/',
    vcf_header_path = '../Biscuit_Snakemake_Workflow/analysis/pileup/',
    name=nameArg,
    sampNames=meta$sample,
    readBiscuitMergedFlag=TRUE,
    whichFlag = GRanges(
      seqnames=myChr,
      ranges=myRanges,
      strand='*'
    )
  )
  
  names(bs_list) <-  meta$sample
  
  tictoc::tic("Unionize")
  myBS <- bs_list[[1]]
  for (i in (2:length(bs_list))){
    myBS <- biscuiteer::unionize(myBS,bs_list[[i]])
  }
  tictoc::toc()
  
  betaMatrix.tidy <- do.call('rbind',
                             lapply(
                               seq_along(bs_list), FUN = function(
                                 x,n,i
                               ){
                                 # print(paste(x[i],n[i],i))
                                 cbind(
                                   data.frame(rowRanges(x[[i]])),
                                   beta = as.numeric(round(x[[i]]@assays@data$M / x[[i]]@assays@data$Cov, 2)),
                                   depth = as.numeric(x[[i]]@assays@data$Cov),
                                   sample = rep(n[[i]],nrow(x[[i]]))
                                 )
                               },
                               x=bs_list,
                               n=names(bs_list)
                             )
  )
  
  betaMatrix.tidy$pos <- paste0(betaMatrix.tidy$seqnames,':',betaMatrix.tidy$start)
  
  
  ## GET THE BETA MATRIX
  betaMatrix <- betaMatrix.tidy %>% 
    dplyr::select(one_of(c('sample','beta','pos'))) %>%
    tidyr::pivot_wider(values_from = beta,names_from = sample)
  dim(betaMatrix)
  betaMatrix <- dplyr::arrange(betaMatrix,pos)
  .rownames <- betaMatrix$pos
  betaMatrix <- dplyr::select(betaMatrix,-pos)
  dim(betaMatrix)
  rownames(betaMatrix) <- .rownames
  
  # Specific to project
  
  
  betaMatrix <- t(betaMatrix)
  
  
  betaMatrix2 <- betaMatrix[new_order,]
  myHM <- ComplexHeatmap::Heatmap(
    # t(betaMatrix),
    betaMatrix2,
    show_column_names = FALSE,
    col=heatmapColorPal,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title='Beta'),
    row_title_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_names_gp = gpar(
      fontsize = 8
    ),
    column_title = nameArg, 
    column_title_gp = gpar(fontsize = 22),
    right_annotation = hm_sample_anno,
    heatmap_width = unit(9, "in"),  
    heatmap_height = unit(7, "in"),  
    row_title_rot = 0
  )
  
  print(
    myHM
  )
  
  return(
    betaMatrix
  )
}

densPlot = function(x,Covfilter,xlabels=TRUE,adjust_lvl=4){
    # color_scheme = colorRampPalette(c('Yellow','Blue','Black'))(100)
    # color_scheme = colorRampPalette(c('Blue','Red'))(100)
    color_scheme = matlab::jet.colors(256)
    nCGs <- length(which(x$Cov>=Covfilter))
    medianC=round(median(x$M),1)*100
    if(medianC==0){medianC=1}
    mySample <- unique(x$sample)
    # stopifnot
    plot <- dplyr::filter(x,Cov>=Covfilter) %>% ggplot(aes(x=sample,y=M)) + geom_violin(fill=color_scheme[medianC],adjust=adjust_lvl) + theme_bw() +
    # plot <- dplyr::filter(x,Cov>=Covfilter) %>% ggplot(aes(x=sample,y=M)) + geom_violin(fill=color_scheme[medianC]) + geom_boxplot(aes(x=sample,y=M),fill=NA) + theme_bw() +
      #ggtitle(paste0(nCGs," Solo WCGW In PMDs >=",Covfilter," depth\n")) +
      # ggtitle(medianC) + 
      theme(
        plot.title = element_text(size = 12),
        # axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=if(xlabels){element_text(angle = 90,vjust=0.5)}else{element_blank()},
        plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), "null"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
      ) +
      ylim(c(0,1))
    print(
      plot + ggtitle(paste0(mySample,' ',nCGs," CpGs >=",Covfilter,"\n"))
    )
    return(plot)
}

densPlotLegend = function(){
  require(ggplot2)
  # color_scheme = colorRampPalette(c('Yellow','Blue','Black'))(100)
  # color_scheme = colorRampPalette(c('Blue','Red'))(100)
  color_scheme = matlab::jet.colors(256)
  ggplot(data.frame(legend=1:256),aes(legend)) + geom_histogram(binwidth = 1,fill=color_scheme) + scale_x_continuous(breaks=c(0,100),labels=c('0','1')) + coord_flip() +  theme_classic() + theme(
    plot.title = element_text(size = 12),
    plot.background = element_blank(),
    axis.text.x=element_blank(),
    axis.line.x = element_blank(),
    # axis.line.y = element_blank(),
    axis.text.y= element_text(size=12),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.1, 1, 0.1, 1), "null"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  ) + ggtitle('Median Methylation') 
}

# readBiscuit_as_list wrapper around biscuiteer::readBiscuit
# designed to work well with Biscuit_Snakemake
readBiscuit_as_list <- function(directory, sampNames, referenceGenome, readBiscuitMergedFlag, vcf_path, whichFlag){
  require(tictoc)
  bsseq_list = list()
  for (i in 1:length(sampNames)){#length(sampNames)){
    bsseq_list[i] = biscuiteer::readBiscuit(BEDfile = paste0(directory,sampNames[i],"_mergecg.bed.gz"),
                                            VCFfile = paste0(vcf_path,sampNames[i],".vcf.gz"),
                                            genome = referenceGenome,
                                            merged = readBiscuitMergedFlag,
                                            which = whichFlag)
    }
  return(bsseq_list)
}

# file to load and save the biscuit bsseq obejcts as an Rds
# generally a wrapper around readBiscuit_as_list()
# if the Rds exists, then this will be loaded instead of readBiscuit_as_list FUNCTION
generate_or_load_bsseq_list = function(name,vcf_header_path,path,sampNames,readBiscuitMergedFlag,whichFlag){
  myRDS <- paste0("biscuiteer_list_",name,".Rds")
  if(!file.exists(myRDS)){
  
    # import as list of bsseq objects with biscuiteer
    # tic("readBiscuit_as_list()")
    BSseq_list <- readBiscuit_as_list(directory = paste0(path,"/"),
                                      vcf_path = paste0(vcf_header_path,"/"),
                                      sampNames = sampNames,
                                      referenceGenome = config$ref$fasta,
                                      readBiscuitMergedFlag = readBiscuitMergedFlag,
                                      which = whichFlag)
    # toc()
  
    cat('Loaded',length(BSseq_list),'BSseq objects\n')
    
    saveRDS(BSseq_list, file = myRDS)
    return(BSseq_list)
  }else{
    BSseq_list <- readRDS(file = myRDS)
    return(BSseq_list)
  }
  
}

getMult_densPlots = function(myBSlist,name,ylimit,CovFilter){
  myDensPlots <- list()
  require(patchwork)
  require(bsseq)
  require(tidyverse)
  for(i in 1:length(myBSlist)){
  
    sample <- gsub('.sorted.markdup','',pData(myBSlist[[i]])$sampleNames)
  
    pData(myBSlist[[i]])$group <- meta$group[i]
    
    # get the Methylation values
    x<-data.frame(
      M = bsseq::getCoverage(myBSlist[[i]],type="M")[,1]/bsseq::getCoverage(myBSlist[[i]],type="Cov")[,1],
      Cov = bsseq::getCoverage(myBSlist[[i]],type="Cov")[,1]
    )
    x$sample = rep(sample,nrow(x))
    plot <- densPlot(x,CovFilter)
    plot <- plot + coord_cartesian(ylim=c(0,ylimit))
    myDensPlots <- rlist::list.append(myDensPlots,plot)
  } 
  
  names(myDensPlots) <- meta$sample
  
  pw <- (
        myDensPlots[["RD_29"]] + theme(axis.text.y=element_text(face="bold")) |
        myDensPlots[["RD_30"]] |
        myDensPlots[["ND_27"]] | 
        myDensPlots[["ND_22"]] |
        myDensPlots[["D0_24"]] | 
        myDensPlots[["D0_25"]] | 
        myDensPlots[["D0_26"]] | 
        myDensPlots[["D2_3"]] |
      # ) /
      # (
        myDensPlots[["D2_4"]] |
        myDensPlots[["D4_5"]] |
        myDensPlots[["D4_6"]] |
        myDensPlots[["D5_12"]] |
        myDensPlots[["D5_13"]] |
        myDensPlots[["D8_14"]] |
        myDensPlots[["D8_15"]] |
        myDensPlots[["D8_17"]] |
        myDensPlots[["D8_18"]]
      ) + theme(legend.position="left") +
        patchwork::plot_annotation(tag_levels = 'A',title=paste0('Coverage >=',myCovFilter))
      # myDensPlots[[3]]  + patchwork::plot_annotation(tag_levels = 'A',title='')
  
  myReturn <- list(
      densPlots = myDensPlots,
      patchwork = pw
  )
  saveRDS(myReturn,file=paste0(name,'.Rds'))
  return(
    myReturn  
  )
}

# deprecated, pass instead a granges to readBiscuit
intersectBSseqListGenomicRanges = function(subject.gr,queryList,maxgap){ # returns a subset bsseq object
  returnList=list()
  for(i in 1:length(queryList)){
    myBSSEQ <- queryList[[i]]
    # isec <- GenomicRanges::findOverlaps(
    #     query = myBSSEQ, #fg.gr
    #     subject = subject.gr
    # )
    # isec <- GenomicRanges::intersect(myBSSEQ,subject.gr)
    
    
    isec <- subsetByOverlaps(myBSSEQ,subject.gr,maxgap = maxgap)

    pData(isec)$sampleNames <- pData(myBSSEQ)$sampleNames
    # class(isec)
    # head(getCoverage(isec))
    returnList=rlist::list.append(returnList,isec)
    
  }
  print(length(returnList))
  return(returnList)
}

# deprecated, use snakemake pipeline
slidingWindow = function(dir){
  files <- list.files(dir)[grep('.tsv',list.files(dir))]
  # if(length(files)){
  for(f in files){
    name <- gsub('.tsv','',f)
    data <- read.delim(paste0(dir,f),sep = "\t",check.names = FALSE)
    rownames(data) <- paste(data$chr,data$windowPos,sep="_")
    colnames(data) <- gsub('_mergecg','',colnames(data))
    dim(data)
    data <- data[,-which(colnames(data)%in%c('chr','windowPos'))]
    if(any(is.na(rowSums(data)))){
      data <- data[-which(is.na(rowSums(data))),]
    }
    dim(data)
    if(! exists("masterList")){
      masterList <- list(data)
      names(masterList) <- name
    }else{
      masterList <- rlist::list.append(masterList,data)
      names(masterList)[length(masterList)] <- name
    }
  }
  return(masterList)
  # }else{
  #   return(paste("No files found in",dir))
  # }
}

# deprecated, use snakemake pipeline
BinnedAvgPCA = function(binnedDF,meta = meta, name = NULL){
  pca <- prcomp(t(binnedDF))
  pr_comps <- data.frame(pca$x)
  pr_comps$sample <- rownames(pr_comps)
  pr_comps <- dplyr::left_join(meta,pr_comps,by='sample')
  
  
  pca_plot <- ggplot(pr_comps, aes(x=PC1, y=PC2)) +
    geom_point(size=7,aes(color = group)) +
    theme_bw() +
    viridis::scale_color_viridis(discrete = TRUE)
  # scale_color_brewer(palette = 'Set1') 
  pca_plot_label <- pca_plot +
    geom_label_repel(aes(label=sample),size=2,show.legend = FALSE) #labels for individual datapoints
  
  # Plot percent variation explained
  prop_var <- data.frame(t(summary(pca)$importance))
  names(prop_var) = c('sd', 'prop', 'cum')
  prop_var$num = 1:nrow(prop_var)
  
  var_plot <- ggplot(prop_var, aes(x=num, y=prop)) +
    geom_point(size=1.5) +
    geom_line() +
    scale_x_continuous(limits = c(1, 12), breaks = 1:12) +
    xlab("Principal Component") +
    ylab("Proportion of\n Variance") +
    ggtitle("") +
    theme_bw() +
    theme(
      axis.title.y = element_text(vjust=1),
      plot.margin = unit(c(0,0,0,6), "mm")
    )
  
  
  library(patchwork)
  layout <- '
  A
  A
  A
  A
  B
  '
  pca_plot +
    var_plot + plot_annotation(tag_levels = 'A',title=paste0('Methylation average - ',name,' CpG bins')) + plot_layout(design = layout)
}


## Function from DMRcate
## modified here for mostly colors

dmr.plot.modified <- function (ranges, dmr, CpGs, main_title, what = c("Beta", "M"), arraytype = c("EPIC","450K"), phen.col, genome = c("hg19", "hg38", "mm10"), ...){
  require(Gviz)
  
  eh = ExperimentHub()
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  genome <- match.arg(genome)
  stopifnot(class(CpGs)[1] %in% c("matrix", "BSseq", "GenomicRatioSet"))
  stopifnot(dmr %in% 1:length(ranges))
  group <- unique(names(phen.col))
  if (is(CpGs, "matrix") | is(CpGs, "GenomicRatioSet")) {
    if (is(CpGs, "matrix")) {
      if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
      if (arraytype == "EPIC") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
    }
    else {
      grset <- CpGs
    }
    CpGs <- getBeta(grset)
    RSanno <- getAnnotation(grset)
    RSanno <- RSanno[order(RSanno$chr, RSanno$pos), ]
    CpGs <- CpGs[rownames(RSanno), ]
    cpgs.ranges <- GRanges(RSanno$chr, IRanges(RSanno$pos, 
                                               RSanno$pos))
    values(cpgs.ranges) <- CpGs
    isbsseq <- FALSE
  }
  else {
    if (any(width(CpGs) > 1)) {
      stop("Error: all ranges in the BSseq object must be single nucleotides with width 1.")
    }
    if (is.null(rownames(colData(CpGs)))) {
      stop("Error: BSseq object must be annotated with colData with sample IDs as rownames of the data.frame.")
    }
    stopifnot(ncol(CpGs) == length(phen.col))
    cpgs.ranges <- CpGs
    isbsseq <- TRUE
  }
  ranges$ID <- paste0("DMR_", 1:length(ranges))
  ranges.reduce <- GenomicRanges::reduce(ranges + 5000)
  dmrs.inplot <- ranges[ranges %over% ranges.reduce[subjectHits(findOverlaps(ranges[dmr], 
                                                                             ranges.reduce))]]
  ranges.inplot <- ranges.reduce[ranges.reduce %over% dmrs.inplot]
  cpgs.ranges <- subsetByOverlaps(cpgs.ranges, ranges.inplot)
  if (isbsseq) {
    methRatios <- GRanges(seqnames(cpgs.ranges), ranges(cpgs.ranges), 
                          mcols = as.matrix(getCoverage(cpgs.ranges, type = "M"))/as.matrix(getCoverage(cpgs.ranges, 
                                                                                                        type = "Cov")))
  }
  else {
    methRatios <- cpgs.ranges
  }
  values(methRatios) <- as.matrix(values(methRatios))
  colnames(values(methRatios)) <- gsub("mcols.", "", colnames(values(methRatios)))
  dt.group <- lapply(unique(names(phen.col)), function(i) DataTrack(cex.axis=0.5,methRatios[, 
                                                                               names(phen.col) %in% i], name = i, background.title = phen.col[i], 
                                                                    type = "heatmap", showSampleNames = TRUE, ylim = c(0, 
                                                                                                                       1), genome = genome, gradient = c("blue", "yellow")))
  dt.group <- c(dt.group, list(DataTrack(methRatios, groups = names(phen.col), 
                                         type = "smooth", aggregateGroups = TRUE, aggregation = function(x) mean(x, 
                                                                                                                 na.rm = TRUE), col = phen.col[sort(group)], ylim = c(0, 
                                                                                                                                                                      1), name = "Smoothed\n group means", na.rm = TRUE)))
  switch(genome, hg19 = {
    grt = eh[["EH3133"]]
  }, hg38 = {
    grt = eh[["EH3135"]]
  }, mm10 = {
    grt = eh[["EH3137"]]
  })
  chromosome(grt) <- as.character(seqnames(methRatios)[1])
  extras <- list(AnnotationTrack(dmrs.inplot, name = "DMRs", 
                                 showFeatureId = TRUE, col = NULL, fill = "purple", id = dmrs.inplot$ID, 
                                 fontcolor = "black"))
  values(cpgs.ranges) <- NULL
  basetracks <- list(IdeogramTrack(genome = genome, chromosome = as.character(seqnames(ranges.inplot))), 
                     GenomeAxisTrack(), grt, AnnotationTrack(GRanges(seqnames(cpgs.ranges), 
                                                                     ranges(cpgs.ranges)), name = "CpGs", fill = "green", 
                                                             col = NULL, stacking = "dense"))
  suppressWarnings(plotTracks(c(basetracks, extras, dt.group), 
                              from = start(ranges.inplot), to = end(ranges.inplot),
                              main=main_title,cex.main=0.8,
                              ...))
}

annotate_dmr.ranges = function(dmr.ranges,t2g){
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(org.Hs.eg.db)
  require(GenomicFeatures)
  g <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # remove large "genes" from g
  length(g)
  g <- g[-which(data.frame(g)$width>1e6)]
  length(g)
  
  res <- findOverlaps(
    query = dmr.ranges,
    subject = g,
    maxgap = 1e4# a granges from the transcript database
  )
  
  # bind together the genes from txdb & the dmrcate information
  x <- cbind(
    as.data.frame(g[res@to]) %>% dplyr::rename('gene_start'='start','gene_end'='end','gene_width'='width'),
    as.data.frame(dmr.ranges[res@from])
  )
  # remove duplicated colnames
  remove <- which(duplicated(colnames(x)))
  cat("removing", colnames(x)[remove],"\n")
  x <- x[,-remove]
  
  resDF <- dplyr::left_join(x,t2g,by=c('gene_id'='entrezgene_id'))
  dim(resDF)
  
  # Add hyper v. hypo annotation
  # hyper <- dmr.ranges[mcols(dmr.ranges)$meandiff>0]
  # hypo <- dmr.ranges[mcols(dmr.ranges)$meandiff<0]
  ## HYPER (methylated more in X when X_v_Z)
  resDF$direction <- ifelse(resDF$meandiff>0,'hyper','hypo')
  
  return(resDF)
}
