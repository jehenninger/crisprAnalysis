trimReadsAndAlign <- function() {
  
  # Load genome file name
  genome_file <- "E:\\zon_lab\\Sequencing\\MGH_SEQUENCING\\custom_genome\\Jon-20.fa"
  # genome_file <- choose.files(caption = "Select indexed genome file...", multi = FALSE,
  #                             filters = c("Fasta file (*.fa)","*.fa"))
  genome_file_no_ext <- tools::file_path_sans_ext(genome_file)
  
  #Load sample fastq file names and get directory
  fq_fnames <- choose.files(caption = "Select sample fastq files...",
                            filters = c("Fastq file (*.fastq)","*.fastq"))
  fq_dir <- gsub("/","\\\\",dirname(fq_fnames[1]))
  fq_basename <- basename(fq_dir)
  fq_dir_parent <- gsub("/","\\\\",dirname(fq_dir))
  fq_fnames <- basename(fq_fnames)
  
  
  if(dir.exists(file.path(fq_dir_parent,"trimmed_fastq"))) {
    trimmed_fastq_dir <- file.path(fq_dir_parent,"trimmed_fastq",fsep="\\")
  } else {
    dir.create(file.path(fq_dir_parent,"trimmed_fastq"))
    trimmed_fastq_dir <- file.path(fq_dir_parent,"trimmed_fastq",fsep="\\")
  }
  
  #Generate file names for trimmed fastq files
  fq_fnames1 <- gsub(".fastq",".trimmed.paired1.fastq",fq_fnames)
  fq_fnames1 <- file.path(trimmed_fastq_dir,fq_fnames1,fsep="\\")
  fq_fnames2 <- gsub(".fastq",".trimmed.paired2.fastq",fq_fnames)
  fq_fnames2 <- file.path(trimmed_fastq_dir,fq_fnames2,fsep="\\")
  
  fq_fnames_full <- file.path(fq_dir_parent,fq_basename,fq_fnames,fsep="\\")
  
  #Trim fastq by quality using bbduk method
  percent_quality_reads = NULL
  for(i in 1:length(fq_fnames)) {
    cmd <- paste0("java -ea -Xmx100m -cp C:\\Users\\Jon\\bioinformatics_tools\\bbmap\\current jgi.BBDukF", " ",
                  "in=", fq_fnames_full[i], " ",
                  "out1=", fq_fnames1[i], " ",
                  "out2=", fq_fnames2[i], " ",
                  "qtrim=rl", " ",
                  "trimq=30"," ",
                  "minlen=50"," ",
                  "overwrite=t")
    message(cmd, "\n"); system(cmd)
    
    #Find number of reads/percentage before and after trimming 
    origFastq <- fastqq(fq_fnames_full[i])
    origNReads <- nReads(origFastq)/2 #divide by 2 to account for paired ends
    
    newFastq <- fastqq(fq_fnames1[i])
    newNReads <- nReads(newFastq) #since the paired reads are split now, don't need to divide by 2
    
    percent_quality_reads[i] <- 100*newNReads/origNReads
  }
  
  if(dir.exists(file.path(fq_dir_parent,"sorted_bam"))) {
    sorted_bam_dir <- file.path(fq_dir_parent,"sorted_bam",fsep="\\")
  } else {
    dir.create(file.path(fq_dir_parent,"sorted_bam"))
    sorted_bam_dir <- file.path(fq_dir_parent,"sorted_bam",fsep="\\")
  }
  
  #Map, sort, and index the bam files
  sm_fnames <- gsub(".fastq",".paired.sam",fq_fnames)
  sm_fnames <- file.path(sorted_bam_dir,sm_fnames,fsep="\\")
  bm_fnames <- gsub(".fastq",".paired.bam",fq_fnames)
  bm_fnames <- file.path(sorted_bam_dir,bm_fnames,fsep="\\")
  srt_bm_fnames <- gsub(".paired.bam",".sorted.paired.bam",bm_fnames)
  
  for(i in 1:length(fq_fnames)) {
    cmd <- paste0("bowtie2 -x"," ", genome_file_no_ext,
                  " -1 ", fq_fnames1[i],
                  " -2 ", fq_fnames2[i],
                  " -S ", sm_fnames[i], " ","--very-sensitive",
                  " && samtools view -h -b ", sm_fnames[i],
                  " -o ", bm_fnames[i],
                  " && samtools sort -o ", srt_bm_fnames[i], " ",
                  bm_fnames[i],
                  " && samtools index ", srt_bm_fnames[i])
    message(cmd, "\n"); system(cmd)
    unlink(bm_fnames[i])
    unlink(sm_fnames[i])
  }
}