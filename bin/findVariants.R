findVariants <- function(){
  
  # Load genome file name
  genome_file <- "E:\\zon_lab\\Sequencing\\MGH_SEQUENCING\\custom_genome\\Jon-20.fa"
  # genome_file <- choose.files(caption = "Select indexed genome file...", multi = FALSE,
  #                             filters = c("Fasta file (*.fa)","*.fa"))
  genome_file_no_ext <- tools::file_path_sans_ext(genome_file)
  
  #Load GTF file for custom genome
  gtf_fname <- "E:\\zon_lab\\Sequencing\\MGH_SEQUENCING\\custom_genome\\Jon-20.gtf"
  # gtf_fname <- choose.files(caption = "Select genome gtf file...",
  #                           filters = c("GTF file (*.gtf)","*.gtf"))
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname,format = "gtf")
  
  #Import BED file information of CRISPR target sites
  gd_fname <- "E:\\zon_lab\\Sequencing\\MGH_SEQUENCING\\R_custom_genome_target_locations.bed"
  # gd_fname <- choose.files(caption = "Select bed file with CRISPR targets...",
  #                          filters = c("BED file (*.bed)","*.bed"))
  gd <- rtracklayer::import(gd_fname)
  gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center") #add 5 bp to each side for visualization
  
  
  #Get reference sequences for each gene
  genome_seq <- readDNAStringSet(genome_file)
  reference <- BSgenome::getSeq(genome_seq,gdl)
  
  
  
  #Load metadata table
  md_fname <- choose.files(caption = "Select bam file meta data...",
                           filters = c("Excel file (*.xlsx)","*.xlsx"))
  md <- gdata::read.xls(md_fname,1)
  
  
  
  #Convert metadata columns to characters
  md$bamfile <- file.path(md$bamfile)
  md$short.name <- as.character(md$short.name)
  md$genes <- as.character(md$genes)
  md$group <- as.character(md$group)
  
  
  #Get targeted regions and delete empty spaces
  injected_gRNA <- md$genes
  md <- subset(md, select = -genes)
  md[md == ""] <- NA
  md <- na.omit(md)
  injected_gRNA[injected_gRNA == ""] <- NA
  injected_gRNA <- na.omit(injected_gRNA)
  
  
  #Get total number of aligned reads per sample
  total_aligned_reads <- lapply(md$bamfile, function(x) countBam(x,param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))))
  md$total_reads <- sapply(total_aligned_reads, function(x) x$records)
  md <- subset(md,total_reads>0)
  
  #Find indices of targeted regions in the gdl GRanges
  gene_idx <- pmatch(injected_gRNA, gdl$name)
  
  #Rearrange groups to individual fish
  samples <- split(md, md$short.name)
  num_of_samples <- length(samples)
  
  
  #Get or make directory for excel outputs
  bam_dir <- dirname(dirname(md$bamfile[1]))
  
  if(dir.exists(file.path(bam_dir,"excel_output"))) {
    excel_output_dir <- file.path(bam_dir,"excel_output")
  } else {
    dir.create(file.path(bam_dir,"excel_output"))
    excel_output_dir <- file.path(bam_dir,"excel_output")
  }
  
  
  # generate excel outputs
  crispr_set <- list(length(gene_idx))
  wb <- list()
  excel_fname <- list()
  
  min_count <- 25
  min_freq <- 0.5
  
  for(i in 1:num_of_samples) {
    num_of_groups <- length(samples[[i]]$group)
    excel_fname <- file.path(excel_output_dir,paste(unique(samples[[i]]$short.name),".xlsx",sep=""),fsep="\\")
    crispr_set <- readsToTargets(file.path(samples[[i]]$bamfile), targets = gdl[gene_idx], references = reference[gene_idx],
                                 names = samples[[i]]$group, target.loc = 22, verbose=TRUE)
    
    wb <- createWorkbook(type="xlsx")
    sheet <- createSheet(wb,sheetName="variants")
    
    r_count <- 1
    c_count <- 1
    
    header <- matrix(NA,nrow=1,ncol=num_of_groups,dimnames=list(c("total reads mapped"), c(paste(samples[[i]]$short.name,samples[[i]]$group,sep="_"))))
    
    header[1,] <- samples[[i]]$total_reads
    
    addDataFrame(header,sheet,startRow=1,startColumn=1)
    addDataFrame("sample name",sheet,startRow=1,startColumn=1,row.names=FALSE,col.names=FALSE)
    
    r_count <- r_count+1
    
    ###########################  UNFILTERED  #########################
    
    for(j in 1:length(crispr_set)){
      #get total variant counts
      var_data_reads <- variantCounts(crispr_set[[j]])
      var_data_percent <- variantCounts(crispr_set[[j]],result="proportions")
      var_data_percent <- round(var_data_percent,3)
      
      #get total indel counts
      indel_reads <- var_data_reads[-grep("SNV|no var|Other",rownames(var_data_reads)),,drop = FALSE]
      # indel_percent <- var_data_percent[-grep("SNV|no var|Other",rownames(var_data_percent)),]
      
      #test to make sure that there are any indels
      if(length(indel_reads)>0){
        
        
        indel_cigar <- rownames(indel_reads)
        
        
        #Get rid of locations, colons, and commas from modified CIGAR strings
        #Isolate numbers from strings to detect if frameshift or in-frame
        indel_cigar <- gsub("^[^:]*:|,.*:","",indel_cigar)
        indel_cigar <- strsplit(indel_cigar,"[DI]")
        indel_cigar <- lapply(indel_cigar,function(x) sum(as.numeric(x)))
        indel_cigar <- unlist(indel_cigar)
        fs_test <- indel_cigar%%3>0
      }
      
      #test if all groups of a sample have aligned reads to the amplicon
      amplicon_test <- lapply(samples[[i]]$group, function(x) grep(x, lapply(crispr_set[[j]]$crispr_runs,
                                                                             function(x) x$name)))
      
      amplicon_test <- amplicon_test >= 1
      amplicon_test[is.na(amplicon_test)] <- FALSE
      group_names <- samples[[i]]$group
      
      
      #intialize stored variables
      reads_mapped_to_target <- integer(num_of_groups)
      reads_mapped_to_target_percent <- integer(num_of_groups)
      total_no_var_reads <- integer(num_of_groups)
      total_no_var_percent <- integer(num_of_groups)
      total_indel_reads <- integer(num_of_groups)
      total_indel_percent <- integer(num_of_groups)
      total_fs_reads <- integer(num_of_groups)
      total_fs_percent <-integer(num_of_groups)
      
      for(k in 1:num_of_groups){
        
        
        if(amplicon_test[[k]]){ #if there are aligned reads
          
          #find column index of group
          group_idx <- grep(group_names[k],colnames(var_data_reads))
          
          group_var_data_reads <- var_data_reads[,group_idx, drop = FALSE]
          
          group_indel_reads <- indel_reads[,group_idx, drop = FALSE]
          
          
          #total reads and % mapped to target
          reads_mapped_to_target[k] <- length(crispr_set[[j]]$crispr_runs[[group_idx]]$alns)
          reads_mapped_to_target_percent[k] <- round(100*reads_mapped_to_target[group_idx]/samples[[i]]$total_reads[[group_idx]],3)
          
          #get no variant counts and %
          no_var_reads <- sum(group_var_data_reads[grep("no var",rownames(group_var_data_reads))])
          no_var_percent <- 100*no_var_reads/reads_mapped_to_target[k]
          
          #get SNV reads and %
          snv_reads <- sum(group_var_data_reads[grep("SNV",rownames(group_var_data_reads))])
          snv_percent <- 100*snv_reads/reads_mapped_to_target[k]
          
          #sum the SNV and no var to get total "wildtype" counts and %
          total_no_var_reads[k] <- no_var_reads + snv_reads
          total_no_var_percent[k] <- no_var_percent + snv_percent
          
          #sum total indels and frameshifts
          #test if there are indel reads
          
          if(length(group_indel_reads)>0){
            total_indel_reads[k] <- sum(group_indel_reads)
            total_indel_percent[k] <- 100*total_indel_reads[k]/reads_mapped_to_target[k]
            
            total_fs_reads[k] <- sum(group_indel_reads[fs_test,,drop = FALSE])
            total_fs_percent[k] <- 100*total_fs_reads[k]/reads_mapped_to_target[k]
            
            
            
          } else {
            total_indel_reads[k] <- 0
            total_indel_percent[k] <- 0
            total_fs_reads[k] <- 0
            total_fs_percent[k] <- 0
          }
          
        } else { #if tissue is missing for amplicon, then output all 0's
          reads_mapped_to_target[k] <- 0
          reads_mapped_to_target_percent[k] <- 0
          total_no_var_reads[k] <- 0
          total_no_var_percent[k] <- 0
          total_indel_reads[k] <- 0
          total_indel_percent[k] <- 0
          total_fs_reads[k] <- 0
          total_fs_percent[k] <-0
          
        }
      }
      
      #add target name
      addDataFrame(crispr_set[[j]]$target$name,sheet,startRow=4,startColumn=c_count,row.names=FALSE,col.names=FALSE)
      
      #construct matrix
      total_info <- matrix(NA,nrow=5,ncol=2*num_of_groups,dimnames=list(c("reads mapped to target","--","no variant","indels","frameshifts"),
                                                                        c(paste(samples[[i]]$group,"reads",sep="_"),
                                                                          paste(samples[[i]]$group,"percent",sep="_"))))
      
      total_info[1,] <- c(reads_mapped_to_target,reads_mapped_to_target_percent)
      total_info[3,] <- c(total_no_var_reads,total_no_var_percent)
      total_info[4,] <- c(total_indel_reads,total_indel_percent)
      total_info[5,] <- c(total_fs_reads,total_fs_percent)
      
      addDataFrame(total_info,sheet,startRow=5,startColumn=c_count)
      
      
      #variant information
      #Variant data
      var_data <- cbind(var_data_reads,var_data_percent)
      addDataFrame(var_data,sheet,startRow=12, startColumn = c_count)
      
      c_count <- c_count+2*num_of_groups+2
      
    }
    
    
    ########################### FILTERED  ############################
    sheet <- createSheet(wb,sheetName="filtered variants")
    r_count <- 2
    c_count <- 1
    
    for(j in 1:length(crispr_set)){
      #get total variant counts
      var_data_reads <- variantCounts(crispr_set[[j]],
                                      min.count = min_count,
                                      min.freq = min_freq,
                                      include.chimeras = FALSE)
      var_data_percent <- variantCounts(crispr_set[[j]],
                                        result="proportions",
                                        min.count = min_count,
                                        min.freq = min_freq,
                                        include.chimeras = FALSE)
      var_data_percent <- round(var_data_percent,3)
      
      #get total indel counts
      indel_reads <- var_data_reads[-grep("SNV|no var|Other",rownames(var_data_reads)),,drop = FALSE]
      # indel_percent <- var_data_percent[-grep("SNV|no var|Other",rownames(var_data_percent)),]
      
      #test to make sure that there are any indels
      if(length(indel_reads)>0){
        
        
        indel_cigar <- rownames(indel_reads)
        
        
        #Get rid of locations, colons, and commas from modified CIGAR strings
        #Isolate numbers from strings to detect if frameshift or in-frame
        indel_cigar <- gsub("^[^:]*:|,.*:","",indel_cigar)
        indel_cigar <- strsplit(indel_cigar,"[DI]")
        indel_cigar <- lapply(indel_cigar,function(x) sum(as.numeric(x)))
        indel_cigar <- unlist(indel_cigar)
        fs_test <- indel_cigar%%3>0
      }
      
      #test if all groups of a sample have aligned reads to the amplicon
      amplicon_test <- lapply(samples[[i]]$group, function(x) grep(x, lapply(crispr_set[[j]]$crispr_runs,
                                                                             function(x) x$name)))
      
      amplicon_test <- amplicon_test >= 1
      amplicon_test[is.na(amplicon_test)] <- FALSE
      group_names <- samples[[i]]$group
      
      
      #intialize stored variables
      reads_mapped_to_target <- integer(num_of_groups)
      reads_mapped_to_target_percent <- integer(num_of_groups)
      total_no_var_reads <- integer(num_of_groups)
      total_no_var_percent <- integer(num_of_groups)
      total_indel_reads <- integer(num_of_groups)
      total_indel_percent <- integer(num_of_groups)
      total_fs_reads <- integer(num_of_groups)
      total_fs_percent <-integer(num_of_groups)
      
      for(k in 1:num_of_groups){
        
        
        if(amplicon_test[[k]]){ #if there are aligned reads
          
          #find column index of group
          group_idx <- grep(group_names[k],colnames(var_data_reads))
          
          group_var_data_reads <- var_data_reads[,group_idx, drop = FALSE]
          
          group_indel_reads <- indel_reads[,group_idx, drop = FALSE]
          
          
          #total reads and % mapped to target
          reads_mapped_to_target[k] <- length(crispr_set[[j]]$crispr_runs[[group_idx]]$alns)
          reads_mapped_to_target_percent[k] <- round(100*reads_mapped_to_target[group_idx]/samples[[i]]$total_reads[[group_idx]],3)
          
          #get no variant counts and %
          no_var_reads <- sum(group_var_data_reads[grep("no var",rownames(group_var_data_reads))])
          no_var_percent <- 100*no_var_reads/reads_mapped_to_target[k]
          
          #get SNV reads and %
          snv_reads <- sum(group_var_data_reads[grep("SNV",rownames(group_var_data_reads))])
          snv_percent <- 100*snv_reads/reads_mapped_to_target[k]
          
          #sum the SNV and no var to get total "wildtype" counts and %
          total_no_var_reads[k] <- no_var_reads + snv_reads
          total_no_var_percent[k] <- no_var_percent + snv_percent
          
          #sum total indels and frameshifts
          #test if there are indel reads
          
          if(length(group_indel_reads)>0){
            total_indel_reads[k] <- sum(group_indel_reads)
            total_indel_percent[k] <- 100*total_indel_reads[k]/reads_mapped_to_target[k]
            
            total_fs_reads[k] <- sum(group_indel_reads[fs_test,,drop = FALSE])
            total_fs_percent[k] <- 100*total_fs_reads[k]/reads_mapped_to_target[k]
            
            
            
          } else {
            total_indel_reads[k] <- 0
            total_indel_percent[k] <- 0
            total_fs_reads[k] <- 0
            total_fs_percent[k] <- 0
          }
          
        } else { #if tissue is missing for amplicon, then output all 0's
          reads_mapped_to_target[k] <- 0
          reads_mapped_to_target_percent[k] <- 0
          total_no_var_reads[k] <- 0
          total_no_var_percent[k] <- 0
          total_indel_reads[k] <- 0
          total_indel_percent[k] <- 0
          total_fs_reads[k] <- 0
          total_fs_percent[k] <-0
          
        }
      }
      
      #add target name
      addDataFrame(crispr_set[[j]]$target$name,sheet,startRow=4,startColumn=c_count,row.names=FALSE,col.names=FALSE)
      
      #construct matrix
      total_info <- matrix(NA,nrow=5,ncol=2*num_of_groups,dimnames=list(c("reads mapped to target","--","no variant","indels","frameshifts"),
                                                                        c(paste(samples[[i]]$group,"reads",sep="_"),
                                                                          paste(samples[[i]]$group,"percent",sep="_"))))
      
      total_info[1,] <- c(reads_mapped_to_target,reads_mapped_to_target_percent)
      total_info[3,] <- c(total_no_var_reads,total_no_var_percent)
      total_info[4,] <- c(total_indel_reads,total_indel_percent)
      total_info[5,] <- c(total_fs_reads,total_fs_percent)
      
      addDataFrame(total_info,sheet,startRow=5,startColumn=c_count)
      
      
      #variant information
      #Variant data
      var_data <- cbind(var_data_reads,var_data_percent)
      addDataFrame(var_data,sheet,startRow=12, startColumn = c_count)
      
      c_count <- c_count+2*num_of_groups+2
      
    }
    
    
    ##################### FILTERED FRAMESHIFT  ######################
    sheet <- createSheet(wb,sheetName="filtered frameshift variants")
    r_count <- 2
    c_count <- 1
    
    for(j in 1:length(crispr_set)){
      #get total variant counts
      var_data_reads <- variantCounts(crispr_set[[j]],
                                      min.count = min_count,
                                      min.freq = min_freq,
                                      include.chimeras = FALSE)
      var_data_percent <- variantCounts(crispr_set[[j]],
                                        result="proportions",
                                        min.count = min_count,
                                        min.freq = min_freq,
                                        include.chimeras = FALSE)
      var_data_percent <- round(var_data_percent,3)
      
      #get total indel counts
      indel_reads <- var_data_reads[-grep("SNV|no var|Other",rownames(var_data_reads)),,drop = FALSE]
      # indel_percent <- var_data_percent[-grep("SNV|no var|Other",rownames(var_data_percent)),]
      
      #test to make sure that there are any indels
      if(length(indel_reads)>0){
        
        
        indel_cigar <- rownames(indel_reads)
        
        
        #Get rid of locations, colons, and commas from modified CIGAR strings
        #Isolate numbers from strings to detect if frameshift or in-frame
        indel_cigar <- gsub("^[^:]*:|,.*:","",indel_cigar)
        indel_cigar <- strsplit(indel_cigar,"[DI]")
        indel_cigar <- lapply(indel_cigar,function(x) sum(as.numeric(x)))
        indel_cigar <- unlist(indel_cigar)
        fs_test <- indel_cigar%%3>0
        
        #to filter fs reads, find number of rows that are 'no var' or 'SNV'
        wtNumRows <- nrow(var_data_reads[grep("SNV|no var",rownames(var_data_reads)),,drop = FALSE])
        fs_test_padded <- !logical(length = wtNumRows)
        fs_test_padded <- c(fs_test_padded,fs_test)
        
        var_data_reads <- var_data_reads[fs_test_padded,,drop = FALSE]
        var_data_percent <- var_data_percent[fs_test_padded,,drop = FALSE]
      }
      
      #test if all groups of a sample have aligned reads to the amplicon
      amplicon_test <- lapply(samples[[i]]$group, function(x) grep(x, lapply(crispr_set[[j]]$crispr_runs,
                                                                             function(x) x$name)))
      
      amplicon_test <- amplicon_test >= 1
      amplicon_test[is.na(amplicon_test)] <- FALSE
      group_names <- samples[[i]]$group
      
      
      #intialize stored variables
      reads_mapped_to_target <- integer(num_of_groups)
      reads_mapped_to_target_percent <- integer(num_of_groups)
      total_no_var_reads <- integer(num_of_groups)
      total_no_var_percent <- integer(num_of_groups)
      total_indel_reads <- integer(num_of_groups)
      total_indel_percent <- integer(num_of_groups)
      total_fs_reads <- integer(num_of_groups)
      total_fs_percent <-integer(num_of_groups)
      
      for(k in 1:num_of_groups){
        
        
        if(amplicon_test[[k]]){ #if there are aligned reads
          
          #find column index of group
          group_idx <- grep(group_names[k],colnames(var_data_reads))
          
          group_var_data_reads <- var_data_reads[,group_idx, drop = FALSE]
          
          group_indel_reads <- indel_reads[,group_idx, drop = FALSE]
          
          
          #total reads and % mapped to target
          reads_mapped_to_target[k] <- length(crispr_set[[j]]$crispr_runs[[group_idx]]$alns)
          reads_mapped_to_target_percent[k] <- round(100*reads_mapped_to_target[group_idx]/samples[[i]]$total_reads[[group_idx]],3)
          
          #get no variant counts and %
          no_var_reads <- sum(group_var_data_reads[grep("no var",rownames(group_var_data_reads))])
          no_var_percent <- 100*no_var_reads/reads_mapped_to_target[k]
          
          #get SNV reads and %
          snv_reads <- sum(group_var_data_reads[grep("SNV",rownames(group_var_data_reads))])
          snv_percent <- 100*snv_reads/reads_mapped_to_target[k]
          
          #sum the SNV and no var to get total "wildtype" counts and %
          total_no_var_reads[k] <- no_var_reads + snv_reads
          total_no_var_percent[k] <- no_var_percent + snv_percent
          
          #sum total indels and frameshifts
          #test if there are indel reads
          
          if(length(group_indel_reads)>0){
            total_indel_reads[k] <- sum(group_indel_reads)
            total_indel_percent[k] <- 100*total_indel_reads[k]/reads_mapped_to_target[k]
            
            total_fs_reads[k] <- sum(group_indel_reads[fs_test,,drop = FALSE])
            total_fs_percent[k] <- 100*total_fs_reads[k]/reads_mapped_to_target[k]
            
            
            
          } else {
            total_indel_reads[k] <- 0
            total_indel_percent[k] <- 0
            total_fs_reads[k] <- 0
            total_fs_percent[k] <- 0
          }
          
        } else { #if tissue is missing for amplicon, then output all 0's
          reads_mapped_to_target[k] <- 0
          reads_mapped_to_target_percent[k] <- 0
          total_no_var_reads[k] <- 0
          total_no_var_percent[k] <- 0
          total_indel_reads[k] <- 0
          total_indel_percent[k] <- 0
          total_fs_reads[k] <- 0
          total_fs_percent[k] <-0
          
        }
      }
      
      #add target name
      addDataFrame(crispr_set[[j]]$target$name,sheet,startRow=4,startColumn=c_count,row.names=FALSE,col.names=FALSE)
      
      #construct matrix
      total_info <- matrix(NA,nrow=5,ncol=2*num_of_groups,dimnames=list(c("reads mapped to target","--","no variant","indels","frameshifts"),
                                                                        c(paste(samples[[i]]$group,"reads",sep="_"),
                                                                          paste(samples[[i]]$group,"percent",sep="_"))))
      
      total_info[1,] <- c(reads_mapped_to_target,reads_mapped_to_target_percent)
      total_info[3,] <- c(total_no_var_reads,total_no_var_percent)
      total_info[4,] <- c(total_indel_reads,total_indel_percent)
      total_info[5,] <- c(total_fs_reads,total_fs_percent)
      
      addDataFrame(total_info,sheet,startRow=5,startColumn=c_count)
      
      
      #variant information
      #Variant data
      var_data <- cbind(var_data_reads,var_data_percent)
      addDataFrame(var_data,sheet,startRow = 12, startColumn = c_count)
      
      c_count <- c_count+2*num_of_groups+2
      
    }
    
    ##################### FILTERED NONFRAMESHIFT  ######################
    sheet <- createSheet(wb,sheetName="filtered nonframeshift variants")
    r_count <- 2
    c_count <- 1
    
    for(j in 1:length(crispr_set)){
      #get total variant counts
      var_data_reads <- variantCounts(crispr_set[[j]],
                                      min.count = min_count,
                                      min.freq = min_freq,
                                      include.chimeras = FALSE)
      var_data_percent <- variantCounts(crispr_set[[j]],
                                        result="proportions",
                                        min.count = min_count,
                                        min.freq = min_freq,
                                        include.chimeras = FALSE)
      var_data_percent <- round(var_data_percent,3)
      
      #get total indel counts
      indel_reads <- var_data_reads[-grep("SNV|no var|Other",rownames(var_data_reads)),,drop = FALSE]
      # indel_percent <- var_data_percent[-grep("SNV|no var|Other",rownames(var_data_percent)),]
      
      #test to make sure that there are any indels
      if(length(indel_reads)>0){
        
        
        indel_cigar <- rownames(indel_reads)
        
        
        #Get rid of locations, colons, and commas from modified CIGAR strings
        #Isolate numbers from strings to detect if frameshift or in-frame
        indel_cigar <- gsub("^[^:]*:|,.*:","",indel_cigar)
        indel_cigar <- strsplit(indel_cigar,"[DI]")
        indel_cigar <- lapply(indel_cigar,function(x) sum(as.numeric(x)))
        indel_cigar <- unlist(indel_cigar)
        fs_test <- indel_cigar%%3>0
        
        #to filter fs reads, find number of rows that are 'no var' or 'SNV'
        wtNumRows <- nrow(var_data_reads[grep("SNV|no var",rownames(var_data_reads)),,drop = FALSE])
        fs_test_padded <- logical(length = wtNumRows) #initialize these as false so that they are true later on
        fs_test_padded <- c(fs_test_padded,fs_test)
        
        var_data_reads <- var_data_reads[!fs_test_padded,,drop = FALSE]
        var_data_percent <- var_data_percent[!fs_test_padded,,drop = FALSE]
      }
      
      #test if all groups of a sample have aligned reads to the amplicon
      amplicon_test <- lapply(samples[[i]]$group, function(x) grep(x, lapply(crispr_set[[j]]$crispr_runs,
                                                                             function(x) x$name)))
      
      amplicon_test <- amplicon_test >= 1
      amplicon_test[is.na(amplicon_test)] <- FALSE
      group_names <- samples[[i]]$group
      
      
      #intialize stored variables
      reads_mapped_to_target <- integer(num_of_groups)
      reads_mapped_to_target_percent <- integer(num_of_groups)
      total_no_var_reads <- integer(num_of_groups)
      total_no_var_percent <- integer(num_of_groups)
      total_indel_reads <- integer(num_of_groups)
      total_indel_percent <- integer(num_of_groups)
      total_fs_reads <- integer(num_of_groups)
      total_fs_percent <-integer(num_of_groups)
      
      for(k in 1:num_of_groups){
        
        
        if(amplicon_test[[k]]){ #if there are aligned reads
          
          #find column index of group
          group_idx <- grep(group_names[k],colnames(var_data_reads))
          
          group_var_data_reads <- var_data_reads[,group_idx, drop = FALSE]
          
          group_indel_reads <- indel_reads[,group_idx, drop = FALSE]
          
          
          #total reads and % mapped to target
          reads_mapped_to_target[k] <- length(crispr_set[[j]]$crispr_runs[[group_idx]]$alns)
          reads_mapped_to_target_percent[k] <- round(100*reads_mapped_to_target[group_idx]/samples[[i]]$total_reads[[group_idx]],3)
          
          #get no variant counts and %
          no_var_reads <- sum(group_var_data_reads[grep("no var",rownames(group_var_data_reads))])
          no_var_percent <- 100*no_var_reads/reads_mapped_to_target[k]
          
          #get SNV reads and %
          snv_reads <- sum(group_var_data_reads[grep("SNV",rownames(group_var_data_reads))])
          snv_percent <- 100*snv_reads/reads_mapped_to_target[k]
          
          #sum the SNV and no var to get total "wildtype" counts and %
          total_no_var_reads[k] <- no_var_reads + snv_reads
          total_no_var_percent[k] <- no_var_percent + snv_percent
          
          #sum total indels and frameshifts
          #test if there are indel reads
          
          if(length(group_indel_reads)>0){
            total_indel_reads[k] <- sum(group_indel_reads)
            total_indel_percent[k] <- 100*total_indel_reads[k]/reads_mapped_to_target[k]
            
            total_fs_reads[k] <- sum(group_indel_reads[fs_test,,drop = FALSE])
            total_fs_percent[k] <- 100*total_fs_reads[k]/reads_mapped_to_target[k]
            
            
            
          } else {
            total_indel_reads[k] <- 0
            total_indel_percent[k] <- 0
            total_fs_reads[k] <- 0
            total_fs_percent[k] <- 0
          }
          
        } else { #if tissue is missing for amplicon, then output all 0's
          reads_mapped_to_target[k] <- 0
          reads_mapped_to_target_percent[k] <- 0
          total_no_var_reads[k] <- 0
          total_no_var_percent[k] <- 0
          total_indel_reads[k] <- 0
          total_indel_percent[k] <- 0
          total_fs_reads[k] <- 0
          total_fs_percent[k] <-0
          
        }
      }
      
      #add target name
      addDataFrame(crispr_set[[j]]$target$name,sheet,startRow=4,startColumn=c_count,row.names=FALSE,col.names=FALSE)
      
      #construct matrix
      total_info <- matrix(NA,nrow=5,ncol=2*num_of_groups,dimnames=list(c("reads mapped to target","--","no variant","indels","frameshifts"),
                                                                        c(paste(samples[[i]]$group,"reads",sep="_"),
                                                                          paste(samples[[i]]$group,"percent",sep="_"))))
      
      total_info[1,] <- c(reads_mapped_to_target,reads_mapped_to_target_percent)
      total_info[3,] <- c(total_no_var_reads,total_no_var_percent)
      total_info[4,] <- c(total_indel_reads,total_indel_percent)
      total_info[5,] <- c(total_fs_reads,total_fs_percent)
      
      addDataFrame(total_info,sheet,startRow=5,startColumn=c_count)
      
      
      #variant information
      #Variant data
      var_data <- cbind(var_data_reads,var_data_percent)
      addDataFrame(var_data,sheet,startRow=12, startColumn = c_count)
      
      c_count <- c_count+2*num_of_groups+2
      
    }
    
    saveWorkbook(wb,excel_fname)
    
  }
}