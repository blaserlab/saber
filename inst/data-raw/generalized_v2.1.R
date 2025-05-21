####libraries####
library(CrispRVariants) 
library("Rsamtools")
library(rtracklayer)
library(GenomicFeatures)
library(gridExtra)
library(GenomicRanges)
library(RColorBrewer)
library(rPython)
library(sqldf)
library(rmarkdown)
library(foreach)
library(doParallel)
library(stringr)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(ggdendro)
library(cowplot)
library(scales)
library(viridis)
library(ggExtra)
library(tidyverse)
library(ggrepel)
library(boot)

####version notes####
# version 3.0 uses paired illumina reads from genewiz and needleall aligner\
# version 3.1 changes thresholding method for stacked bar plots and csvs to proportion of all reads
# version 3.2 parallelizes needleall
# no umis
# You need a genome file in a folder named "genome" and your fastqs in a folder named "fastq" for this version to work
# version 4.0 combines main pipeline with clone comparer
# version 4.1 spits out fraction informative data too
# version 4.2_normal_samples is designed to eval normal samples
# analysis_v1.0 derived from pipeline v4.2
# analysis_v2.0 cleans up analysis_v1.0
# generalized_v2.0 derives from analysis_v3.0, reducing the extraneous stuff, and doing two sample classes in parallel
# control samples go in one input folder, experimental in the other.

####set output folder and analysis variables####
outs<-"output_training_validation_20190822_v2"
tvaf<-1#expressed as a percent
tpaf<-0.05# expressed as a proportion

####general functions####
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x))

data_summary<-function(x){
  m<-median(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary1<-function(x){
  m<-median(x)
  ymin<-m-(IQR(x)/2)
  ymax<-m+(IQR(x)/2)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary2<-function(x){
  m<-mean(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary3<-function(x){
  m<-mean(x)
  ymin<-m-sd(x)
  ymax<-m+sd(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

# diversity index functions
shannon_func<-function(x){
  prop<-x/sum(x)
  shannon_entropy<--1*sum(prop*log2(prop))
  return(shannon_entropy)
}

simpson_func<-function(x){
  prop<-x/sum(x)
  simpson_diversity<-1/sum(prop^2)
  return(simpson_diversity)
}

general_diversity_func<-function(x,q){
  prop<-x/sum(x)
  general_diversity<-(sum(prop^q))^(1/(1-q))
  return(general_diversity)
}


# ####input data####
# #Be sure this is correct before running!!!                                                
# #enter the experiment data
# 
# genome<-"gestalt_pipeline4.fasta" # include .fasta.  Genome file has to be in genome folder.
# genome_fp<-paste0(getwd(),"/references/",genome)
# output_file<-"training_set"#for the cripsrvariants plot and prefix on sample output csv and vaf plot
# group_name<-"training_set"
# bed<-"references/gestalt2.bed"
# reference<-"references/gestalt2_ref.txt"
# 
# threshold<-200# number of reads below which they don't appear on the big crispr plot
# 
# use_UMI<-FALSE#works how you think it should work


####core pipeline####
# set path
Sys.setenv(PATH = "/home/OSUMC.EDU/blas02/miniconda3/bin:/opt/bcl2fastq:/opt/cellranger-3.0.2:/home/OSUMC.EDU/blas02/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin")



#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

generalized_func<-function(output_folder,use_UMI, genome_fp, ref_seq, control_name, experimental_name, bed,tvaf,tpaf){
  # make working subdirectories and create variables 
  system(paste0("mkdir ",output_folder))
  
  system("mkdir temp"); temp_folder<-paste0(getwd(),"/temp")
  
  system("mkdir temp/con_fastq_split_R1"); con_fastq_split_R1_folder<-paste0(getwd(),"/temp/con_fastq_split_R1") 
  system("mkdir temp/con_fastq_split_R2"); con_fastq_split_R2_folder<-paste0(getwd(),"/temp/con_fastq_split_R2") 
  system("mkdir temp/exp_fastq_split_R1"); exp_fastq_split_R1_folder<-paste0(getwd(),"/temp/exp_fastq_split_R1") 
  system("mkdir temp/exp_fastq_split_R2"); exp_fastq_split_R2_folder<-paste0(getwd(),"/temp/exp_fastq_split_R2") 
  
  
  system("mkdir temp/con_bam_temp"); con_bam_temp<-paste0(getwd(),"/temp/con_bam_temp")
  system("mkdir temp/exp_bam_temp"); exp_bam_temp<-paste0(getwd(),"/temp/exp_bam_temp")
  
  
  system(paste0("mkdir ",output_folder,"/fastqc_files")); fastqc_files_folder<-paste0(output_folder,"/fastqc_files")
  
  system(paste0("mkdir ",output_folder,"/con_bam_output")); con_bam_output_folder<-paste0(output_folder,"/con_bam_output")
  system(paste0("mkdir ",output_folder,"/exp_bam_output")); exp_bam_output_folder<-paste0(output_folder,"/exp_bam_output")
  
  system(paste0("mkdir ",output_folder,"/PEAR_output")); PEAR_output_folder<-paste0(output_folder,"/PEAR_output")
  
  
  
  system("mkdir temp/con_fastq_trimmed"); con_fastq_trimmed<-paste0(getwd(),"/temp/con_fastq_trimmed")
  system("mkdir temp/exp_fastq_trimmed"); exp_fastq_trimmed<-paste0(getwd(),"/temp/exp_fastq_trimmed")
  
  system("mkdir temp/con_fastq_merged"); con_fastq_merged<-paste0(getwd(),"/temp/con_fastq_merged")
  system("mkdir temp/exp_fastq_merged"); exp_fastq_merged<-paste0(getwd(),"/temp/exp_fastq_merged")
  
  system("mkdir temp/con_sam_temp"); con_sam_temp<-paste0(getwd(),"/temp/con_sam_temp")
  system("mkdir temp/exp_sam_temp"); exp_sam_temp<-paste0(getwd(),"/temp/exp_sam_temp")
  
  system("mkdir temp/con_cutadapt1"); con_cutadapt1_folder<-paste0(getwd(),"/temp/con_cutadapt1")
  system("mkdir temp/con_cutadapt2"); con_cutadapt2_folder<-paste0(getwd(),"/temp/con_cutadapt2")
  system("mkdir temp/exp_cutadapt1"); exp_cutadapt1_folder<-paste0(getwd(),"/temp/exp_cutadapt1")
  system("mkdir temp/exp_cutadapt2"); exp_cutadapt2_folder<-paste0(getwd(),"/temp/exp_cutadapt2")
  
  
  
  
  system("mkdir temp/con_fastq_umi_out"); con_fastq_umi_out<-paste0(getwd(),"/temp/con_fastq_umi_out")
  system("mkdir temp/exp_fastq_umi_out"); exp_fastq_umi_out<-paste0(getwd(),"/temp/exp_fastq_umi_out")
  
  system("mkdir temp/con_bam_dedup"); con_bam_dedup<-paste0(getwd(),"/temp/con_bam_dedup")
  system("mkdir temp/exp_bam_dedup"); exp_bam_dedup<-paste0(getwd(),"/temp/exp_bam_dedup")
  
  system(paste0("mkdir ",output_folder,"/con_bam_dedup_output")); con_bam_dedup_output_folder<-paste0(output_folder,"/con_bam_dedup_output")
  system(paste0("mkdir ",output_folder,"/exp_bam_dedup_output")); exp_bam_dedup_output_folder<-paste0(output_folder,"/exp_bam_dedup_output")
  
  
  
  #load samples from fastq folder
  con_samples<-list.files(paste0(getwd(),"/fastq_control")) # loads all files in the fastq folder
  con_samples_fp<-paste0(getwd(),"/fastq_control/",con_samples)
  exp_samples<-list.files(paste0(getwd(), "/fastq_experimental"))
  exp_samples_fp<-paste0(getwd(), "/fastq_experimental/",exp_samples)
  
  
  #unzip
  foreach(i=1:length(con_samples_fp)) %dopar% {
    cmd<-paste0("gzip -d ",con_samples_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  con_samples<-list.files(paste0(getwd(),"/fastq_control")) 
  con_samples_fp<-paste0(getwd(),"/fastq_control/",con_samples)
  
  foreach(i=1:length(exp_samples_fp)) %dopar% {
    cmd<-paste0("gzip -d ",exp_samples_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  exp_samples<-list.files(paste0(getwd(),"/fastq_experimental")) 
  exp_samples_fp<-paste0(getwd(),"/fastq_experimental/",exp_samples)
  
  
  #put into fastq_split_R1 and R2 folders
  foreach(i=1:length(con_samples_fp)) %dopar% {
    if (grepl("_R1_",con_samples_fp[i])==TRUE){
      cmd<-paste0("cp ",con_samples_fp[i]," ",con_fastq_split_R1_folder)
      message(cmd, "\n"); system(cmd)
    } else{
      cmd<-paste0("cp ",con_samples_fp[i]," ",con_fastq_split_R2_folder)
      message(cmd, "\n"); system(cmd)
    }
  }
  
  con_fastq_split_R1_files<-list.files(con_fastq_split_R1_folder); con_fastq_split_R1_fp<-paste0(con_fastq_split_R1_folder,"/",con_fastq_split_R1_files)
  con_fastq_split_R2_files<-list.files(con_fastq_split_R2_folder); con_fastq_split_R2_fp<-paste0(con_fastq_split_R2_folder,"/",con_fastq_split_R2_files)
  
  foreach(i=1:length(exp_samples_fp)) %dopar% {
    if (grepl("_R1_",exp_samples_fp[i])==TRUE){
      cmd<-paste0("cp ",exp_samples_fp[i]," ",exp_fastq_split_R1_folder)
      message(cmd, "\n"); system(cmd)
    } else{
      cmd<-paste0("cp ",exp_samples_fp[i]," ",exp_fastq_split_R2_folder)
      message(cmd, "\n"); system(cmd)
    }
  }
  
  exp_fastq_split_R1_files<-list.files(exp_fastq_split_R1_folder); exp_fastq_split_R1_fp<-paste0(exp_fastq_split_R1_folder,"/",exp_fastq_split_R1_files)
  exp_fastq_split_R2_files<-list.files(exp_fastq_split_R2_folder); exp_fastq_split_R2_fp<-paste0(exp_fastq_split_R2_folder,"/",exp_fastq_split_R2_files)
  
  
  # merge read pairs with PEAR
  con_mid_names<-substr(con_fastq_split_R1_files, 6,13)
  for(i in 1:length(con_fastq_split_R1_fp)){
    cmd<-paste0("pear -f ",con_fastq_split_R1_fp[i]," -r ",con_fastq_split_R2_fp[i]," -o ",con_fastq_merged,"/",con_mid_names[i], " -j 39 > ",PEAR_output_folder,"/",con_mid_names[i],".PEARreport.txt")
    message(cmd, "\n"); system(cmd)
  }
  unlink(paste0(con_fastq_merged,"/*unassembled*"),recursive = TRUE)
  unlink(paste0(con_fastq_merged,"/*discarded*"),recursive = TRUE)
  con_fastq_merged_files<-list.files(con_fastq_merged); con_fastq_merged_fp<-paste0(con_fastq_merged,"/",con_fastq_merged_files)
  
  exp_mid_names<-substr(exp_fastq_split_R1_files, 6,13)
  for(i in 1:length(exp_fastq_split_R1_fp)){
    cmd<-paste0("pear -f ",exp_fastq_split_R1_fp[i]," -r ",exp_fastq_split_R2_fp[i]," -o ",exp_fastq_merged,"/",exp_mid_names[i], " -j 39 > ",PEAR_output_folder,"/",exp_mid_names[i],".PEARreport.txt")
    message(cmd, "\n"); system(cmd)
  }
  unlink(paste0(exp_fastq_merged,"/*unassembled*"),recursive = TRUE)
  unlink(paste0(exp_fastq_merged,"/*discarded*"),recursive = TRUE)
  exp_fastq_merged_files<-list.files(exp_fastq_merged); exp_fastq_merged_fp<-paste0(exp_fastq_merged,"/",exp_fastq_merged_files)
  
  
  #run trimmomatic on split and ordered files
  con_mid_names<-substr(con_fastq_merged_files, 1,5)
  foreach(i=1:length(con_fastq_merged_fp)) %dopar% {
    cmd<-paste0("java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar SE ", #invoke trimmomatic
                con_fastq_merged_fp[i]," ",# merged read input
                con_fastq_trimmed,"/",con_mid_names[i],".trimmed.fastq ",# merged trimmed output
                "SLIDINGWINDOW:4:15 MINLEN:100")# trim parameters.  
    message(cmd, "\n"); system(cmd)
  }
  con_fastq_trimmed_files<-list.files(con_fastq_trimmed); con_fastq_trimmed_fp<-paste0(con_fastq_trimmed,"/",con_fastq_trimmed_files)
  
  exp_mid_names<-substr(exp_fastq_merged_files, 1,5)
  foreach(i=1:length(exp_fastq_merged_fp)) %dopar% {
    cmd<-paste0("java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar SE ", #invoke trimmomatic
                exp_fastq_merged_fp[i]," ",# merged read input
                exp_fastq_trimmed,"/",exp_mid_names[i],".trimmed.fastq ",# merged trimmed output
                "SLIDINGWINDOW:4:15 MINLEN:100")# trim parameters.  
    message(cmd, "\n"); system(cmd)
  }
  exp_fastq_trimmed_files<-list.files(exp_fastq_trimmed); exp_fastq_trimmed_fp<-paste0(exp_fastq_trimmed,"/",exp_fastq_trimmed_files)
  
  
  #run cutadapt to keep only fastqs that have the full flanking primer sequences
  #5 prime adapter = V6/7F:  TCGAGCTCAAGCTTCGG
  con_mid_names<-substr(con_fastq_trimmed_files, 1,5)
  foreach(i=1:length(con_fastq_trimmed_fp)) %dopar% {
    cmd<-paste0("cutadapt -g ^TCGAGCTCAAGCTTCGG --discard-untrimmed -e 0.01 --action=none -o ",con_cutadapt1_folder,"/",con_mid_names[i],".cutadapt1.fastq ",con_fastq_trimmed_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  con_cutadapt1_files<-list.files(con_cutadapt1_folder); con_cutadapt1_fp<-paste0(con_cutadapt1_folder,"/",con_cutadapt1_files)
  
  #3 prime adapter = V6/7R:  GACCTCGAGACAAATGGCAG (reverse complement of the primer sequence 5'-3')
  con_mid_names<-substr(con_cutadapt1_files, 1,5)
  foreach(i=1:length(con_cutadapt1_fp)) %dopar% {
    cmd<-paste0("cutadapt -a GACCTCGAGACAAATGGCAG$ --discard-untrimmed -e 0.01 --action=none -o ",con_cutadapt2_folder,"/",con_mid_names[i],".cutadapt2.fastq ",con_cutadapt1_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  con_cutadapt2_files<-list.files(con_cutadapt2_folder); con_cutadapt2_fp<-paste0(con_cutadapt2_folder,"/",con_cutadapt2_files)
  
  #run cutadapt to keep only fastqs that have the full flanking primer sequences
  #5 prime adapter = V6/7F:  TCGAGCTCAAGCTTCGG
  exp_mid_names<-substr(exp_fastq_trimmed_files, 1,5)
  foreach(i=1:length(exp_fastq_trimmed_fp)) %dopar% {
    cmd<-paste0("cutadapt -g XTCGAGCTCAAGCTTCGG --discard-untrimmed -e 0.01 --action=none -o ",exp_cutadapt1_folder,"/",exp_mid_names[i],".cutadapt1.fastq ",exp_fastq_trimmed_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  exp_cutadapt1_files<-list.files(exp_cutadapt1_folder); exp_cutadapt1_fp<-paste0(exp_cutadapt1_folder,"/",exp_cutadapt1_files)
  
  #3 prime adapter = V6/7R:  GACCTCGAGACAAATGGCAG (reverse complement of the primer sequence 5'-3')
  exp_mid_names<-substr(exp_cutadapt1_files, 1,5)
  foreach(i=1:length(exp_cutadapt1_fp)) %dopar% {
    cmd<-paste0("cutadapt -a GACCTCGAGACAAATGGCAG$ --discard-untrimmed -e 0.01 --action=none -o ",exp_cutadapt2_folder,"/",exp_mid_names[i],".cutadapt2.fastq ",exp_cutadapt1_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  exp_cutadapt2_files<-list.files(exp_cutadapt2_folder); exp_cutadapt2_fp<-paste0(exp_cutadapt2_folder,"/",exp_cutadapt2_files)
  
  
  # run fastqc on trimmed and demultiplexed input files
  foreach(i=1:length(con_cutadapt2_fp)) %dopar% {
    cmd<-paste0("fastqc -o ",fastqc_files_folder," ",con_cutadapt2_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  foreach(i=1:length(exp_cutadapt2_fp)) %dopar% {
    cmd<-paste0("fastqc -o ",fastqc_files_folder," ",exp_cutadapt2_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  # run multiqc to summarize the qc files
  system(paste0('multiqc -d ',output_folder," -o ",output_folder))
  
  
  # extract UMI tags (optional)
  con_mid_names<-substr(con_cutadapt2_files, 1,5)
  if (use_UMI==TRUE){
    foreach(i=1:length(con_cutadatpt2_fp)) %dopar% {
      cmd<-paste0("umi_tools extract --stdin=",con_cutadapt2_fp[i]," --bc-pattern=NNNNNNNNNN --log=",
                  output_folder,"/con_UMI.log --stdout=",con_fastq_umi_out,"/",con_mid_names[i],".processed.fastq")
      message(cmd, "\n"); system(cmd)
    }
    con_fastq_umi_out_files<-list.files(con_fastq_umi_out); con_fastq_umi_out_fp<-paste0(con_fastq_umi_out,"/",con_fastq_umi_out_files)
  }
  
  exp_mid_names<-substr(exp_cutadapt2_files, 1,5)
  if (use_UMI==TRUE){
    foreach(i=1:length(exp_cutadatpt2_fp)) %dopar% {
      cmd<-paste0("umi_tools extract --stdin=",exp_cutadapt2_fp[i]," --bc-pattern=NNNNNNNNNN --log=",
                  output_folder,"/exp_UMI.log --stdout=",exp_fastq_umi_out,"/",exp_mid_names[i],".processed.fastq")
      message(cmd, "\n"); system(cmd)
    }
    exp_fastq_umi_out_files<-list.files(exp_fastq_umi_out); exp_fastq_umi_out_fp<-paste0(exp_fastq_umi_out,"/",exp_fastq_umi_out_files)
  }
  
  # align with needleall
  if (use_UMI==TRUE){
    foreach(i=1:length(con_fastq_umi_out_files)) %dopar% {
      cmd<-paste0("needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence ",genome_fp," -bsequence ",con_fastq_umi_out_fp[i]," -outfile ",con_sam_temp,"/",con_mid_names[i],".sam -errfile ",output_folder,"/con_needleall.error")
      message(cmd, "\n"); system(cmd)
    }
  } else {
    foreach(i=1:length(con_cutadapt2_files)) %dopar% {
      cmd<-paste0("needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence ",genome_fp," -bsequence ",con_cutadapt2_fp[i]," -outfile ",con_sam_temp,"/",con_mid_names[i],".sam -errfile ",output_folder,"/con_needleall.error")
      message(cmd, "\n"); system(cmd)
    }
  }
  
  con_sam_temp_files<-list.files(con_sam_temp); con_sam_temp_fp<-paste0(con_sam_temp,"/",con_sam_temp_files)
  
  if (use_UMI==TRUE){
    foreach(i=1:length(exp_fastq_umi_out_files)) %dopar% {
      cmd<-paste0("needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence ",genome_fp," -bsequence ",exp_fastq_umi_out_fp[i]," -outfile ",exp_sam_temp,"/",exp_mid_names[i],".sam -errfile ",output_folder,"/exp_needleall.error")
      message(cmd, "\n"); system(cmd)
    }
  } else {
    foreach(i=1:length(exp_cutadapt2_files)) %dopar% {
      cmd<-paste0("needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence ",genome_fp," -bsequence ",exp_cutadapt2_fp[i]," -outfile ",exp_sam_temp,"/",exp_mid_names[i],".sam -errfile ",output_folder,"/exp_needleall.error")
      message(cmd, "\n"); system(cmd)
    }
  }
  
  exp_sam_temp_files<-list.files(exp_sam_temp); exp_sam_temp_fp<-paste0(exp_sam_temp,"/",exp_sam_temp_files)
  
  
  # fix sam file header and select only reads matching at position 1
  for (i in 1:length(con_sam_temp_files)) {
    samdf<-read.delim(con_sam_temp_fp[i], sep="\t", header = FALSE)#read in sam file
    samdf<-samdf[-(1:2),]#chop off the old header
    samdf<-samdf[which(samdf$V4=="1"),]#select only reads mapping to coordinate 1
    sam_header<-read.table("references/sam_header.csv", fill=TRUE, header=FALSE, sep=",", colClasses=(rep("character",13)))# read in standard sam header
    names(sam_header)<-paste("V", 1:13, sep="")
    samdf<-rbind(sam_header,samdf)
    write_tsv(samdf, na = "", path = con_sam_temp_fp[i],col_names = FALSE, append=FALSE)
  }
  con_sam_temp_files<-list.files(con_sam_temp); con_sam_temp_fp<-paste0(con_sam_temp,"/",con_sam_temp_files)
  
  for (i in 1:length(exp_sam_temp_files)) {
    samdf<-read.delim(exp_sam_temp_fp[i], sep="\t", header = FALSE)#read in sam file
    samdf<-samdf[-(1:2),]#chop off the old header
    samdf<-samdf[which(samdf$V4=="1"),]#select only reads mapping to coordinate 1
    sam_header<-read.table("references/sam_header.csv", fill=TRUE, header=FALSE, sep=",", colClasses=(rep("character",13)))# read in standard sam header
    names(sam_header)<-paste("V", 1:13, sep="")
    samdf<-rbind(sam_header,samdf)
    write_tsv(samdf, na = "", path = exp_sam_temp_fp[i],col_names = FALSE, append=FALSE)
  }
  exp_sam_temp_files<-list.files(exp_sam_temp); exp_sam_temp_fp<-paste0(exp_sam_temp,"/",exp_sam_temp_files)
  
  
  
  # convert sam to bam
  foreach(i=1:length(con_sam_temp_files)) %dopar% {
    cmd<-paste0("samtools view -S -b ",con_sam_temp_fp[i]," > ",con_bam_temp,"/",con_mid_names[i],".bam")
    message(cmd,"\n"); system(cmd)
  }
  
  con_bam_temp_files<-list.files(path = con_bam_temp); con_bam_temp_fp<-paste0(con_bam_temp,"/",con_bam_temp_files)
  
  foreach(i=1:length(exp_sam_temp_files)) %dopar% {
    cmd<-paste0("samtools view -S -b ",exp_sam_temp_fp[i]," > ",exp_bam_temp,"/",exp_mid_names[i],".bam")
    message(cmd,"\n"); system(cmd)
  }
  
  exp_bam_temp_files<-list.files(path = exp_bam_temp); exp_bam_temp_fp<-paste0(exp_bam_temp,"/",exp_bam_temp_files)
  
  
  # sort and index bam
  foreach(i=1:length(con_bam_temp_fp)) %dopar% {
    cmd<-paste0("samtools sort ",con_bam_temp_fp[i]," -o ",con_bam_temp_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  foreach(i=1:length(con_bam_temp_fp)) %dopar% {
    cmd<-paste0("samtools index ",con_bam_temp_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  system(paste0("cp -r temp/con_bam_temp ",con_bam_output_folder))#move final bam files and indices to output
  
  foreach(i=1:length(exp_bam_temp_fp)) %dopar% {
    cmd<-paste0("samtools sort ",exp_bam_temp_fp[i]," -o ",exp_bam_temp_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  foreach(i=1:length(exp_bam_temp_fp)) %dopar% {
    cmd<-paste0("samtools index ",exp_bam_temp_fp[i])
    message(cmd, "\n"); system(cmd)
  }
  
  system(paste0("cp -r temp/exp_bam_temp ",exp_bam_output_folder))#move final bam files and indices to output
  
  
  
  
  #deduplicate, sort and index UMI tagged reads (optional)
  if (use_UMI==TRUE){
    foreach(i=1:length(con_fastq_umi_out_files)) %dopar% {
      cmd<-paste0("umi_tools dedup --method=unique -I ",con_bam_temp_fp[i]," --output-stats=",con_bam_dedup_output_folder,"/",con_mid_names[i]," -S ",con_bam_dedup,"/",con_mid_names[i],".dedup.bam") 
      message(cmd, "\n"); system(cmd)
    }
    con_bam_dedup_files<-list.files(path = con_bam_dedup); con_bam_dedup_fp<-paste0(con_bam_dedup,"/",con_bam_dedup_files)
    
    foreach(i=1:length(con_bam_dedup_fp)) %dopar% {
      cmd<-paste0("samtools sort ",con_bam_dedup_fp[i]," -o ",con_bam_dedup_fp[i])
      message(cmd, "\n"); system(cmd)
    }
    
    foreach(i=1:length(con_bam_dedup_fp)) %dopar% {
      cmd<-paste0("samtools index ",con_bam_dedup_fp[i])
      message(cmd, "\n"); system(cmd)
    }
    system(paste0("cp -r temp/con_bam_dedup ", con_bam_dedup_output_folder))
  }
  
  if (use_UMI==TRUE){
    foreach(i=1:length(exp_fastq_umi_out_files)) %dopar% {
      cmd<-paste0("umi_tools dedup --method=unique -I ",exp_bam_temp_fp[i]," --output-stats=",exp_bam_dedup_output_folder,"/",exp_mid_names[i]," -S ",exp_bam_dedup,"/",exp_mid_names[i],".dedup.bam") 
      message(cmd, "\n"); system(cmd)
    }
    exp_bam_dedup_files<-list.files(path = exp_bam_dedup); exp_bam_dedup_fp<-paste0(exp_bam_dedup,"/",exp_bam_dedup_files)
    
    foreach(i=1:length(exp_bam_dedup_fp)) %dopar% {
      cmd<-paste0("samtools sort ",exp_bam_dedup_fp[i]," -o ",exp_bam_dedup_fp[i])
      message(cmd, "\n"); system(cmd)
    }
    
    foreach(i=1:length(exp_bam_dedup_fp)) %dopar% {
      cmd<-paste0("samtools index ",exp_bam_dedup_fp[i])
      message(cmd, "\n"); system(cmd)
    }
    system(paste0("cp -r temp/exp_bam_dedup ", exp_bam_dedup_output_folder))
  }
  
  
  #build metadata for experiment
  if (use_UMI==TRUE){
    con_bam_fnames<-con_bam_dedup_fp
  } else {
    con_bam_fnames <- con_bam_temp_fp
  }
  
  if (use_UMI==TRUE){
    exp_bam_fnames<-exp_bam_dedup_fp
  } else {
    exp_bam_fnames <- exp_bam_temp_fp
  }
  
  
  con_group_desig<-rep(control_name, times=length(con_bam_fnames))
  con_md<-read.csv("references/blank_metadata.csv", header = TRUE)
  newrow<-data.frame(bamfile=con_bam_fnames, directory=getwd(),Short.name=con_mid_names,Targeting.type="",sgRNA1="",sgRNA2="",Group=con_group_desig)
  con_md<-rbind(con_md,newrow)
  
  exp_group_desig<-rep(experimental_name, times=length(exp_bam_fnames))
  exp_md<-read.csv("references/blank_metadata.csv", header = TRUE)
  newrow<-data.frame(bamfile=exp_bam_fnames, directory=getwd(),Short.name=exp_mid_names,Targeting.type="",sgRNA1="",sgRNA2="",Group=exp_group_desig)
  exp_md<-rbind(exp_md,newrow)
  
  #create target region
  gd <- rtracklayer::import(bed)
  gdl <- GenomicRanges::resize(gd, width(gd) + 0, fix = "center") #resize region for analysis
  reference0<-read_file(ref_seq)
  reference1<-substr(reference0,1,310)#this has to be here to remove \n.  Can generalize in future
  reference<-Biostrings::DNAString(reference1)

  
  # make the crispr set
  con_crispr_set <- readsToTarget(con_bam_fnames, 
                              target = gd, 
                              reference = reference, 
                              names = con_md$Short.name, 
                              target.loc = 16, #generalize with distance from primer to first sequence of barcode?
                              collapse.pairs = FALSE,
                              split.snv=FALSE)#split.snv=FALSE adds SNVs into no variant count
  
  exp_crispr_set <- readsToTarget(exp_bam_fnames, 
                                  target = gd, 
                                  reference = reference, 
                                  names = exp_md$Short.name, 
                                  target.loc = 16, #generalize with distance from primer to first sequence of barcode?
                                  collapse.pairs = FALSE,
                                  split.snv=FALSE)#split.snv=FALSE adds SNVs into no variant count
  # plot the variants
  ps<-37#need to generalize
  pam_seq<-seq(ps,280,27)#need to generalize
  
  while (!is.null(dev.list())) dev.off()
  con_p <- plotVariants(con_crispr_set, 
                    col.wdth.ratio = c(1,1),
                    plotAlignments.args = list(pam.start = pam_seq, #c(37,64), #draws a 3-nt box starting including the position noted
                                               target.loc = pam_seq-3, #draws a vertical line after the position noted
                                               guide.loc = IRanges::IRanges(pam_seq-20,pam_seq+2), #first parameter - beginning of target sequence, second - end of target sequence
                                               min.count = 200,
                                               tile.height = 0.9),
                    plotFreqHeatmap.args = list(min.count = 200,
                                                plot.text.size = 3, 
                                                x.size = 8, 
                                                group = con_group_desig,
                                                legend.text.size = 8,
                                                legend.key.height = grid::unit(0.5, "lines")))
  dev.copy2pdf(file=paste0(output_folder,"/control_crispr_plot.pdf"), width = 36, height = 36)  #for 10 samples use 24 x 24, for 30-40 samples use 48 x 48
  
  dev.off()
  
  exp_p <- plotVariants(exp_crispr_set, 
                        col.wdth.ratio = c(1,1),
                        plotAlignments.args = list(pam.start = pam_seq, #c(37,64), #draws a 3-nt box starting including the position noted
                                                   target.loc = pam_seq-3, #draws a vertical line after the position noted
                                                   guide.loc = IRanges::IRanges(pam_seq-20,pam_seq+2), #first parameter - beginning of target sequence, second - end of target sequence
                                                   min.count = 200,
                                                   tile.height = 0.9),
                        plotFreqHeatmap.args = list(min.count = 200,
                                                    plot.text.size = 3, 
                                                    x.size = 8, 
                                                    group = exp_group_desig,
                                                    legend.text.size = 8,
                                                    legend.key.height = grid::unit(0.5, "lines")))
  dev.copy2pdf(file=paste0(output_folder,"/experimental_crispr_plot.pdf"), width = 36, height = 36)  #for 10 samples use 24 x 24, for 30-40 samples use 48 x 48
  
  
  #####generate informative, no variant and common variant tables with vaf and paf thresholds####
  generate_v2<-function(x,class,vaf_thresh,paf_thresh,vc_prop,output_folder){
    thresh<-0
    df<-x 
    rnames<-row.names(df)#duplicate row names as new variable
    df$rnames<-rnames#add row name variabe to df
    ##important:  CrispRVariants proportion generates PERCENT values, not FRACTION values.  So TVAF = 1 => 1% vaf, not 100% VAF.
    common_vars<-names(which(rowSums(vc_prop>vaf_thresh)/ncol(vc_prop)>(paf_thresh)))#returns the names of rows with a propotion greater than vaf_thresh occuring in more than paf_thresh
    ##end important
    remove<-c("no variant")# removes only no variant.  Leaves "Other" in.  Maybe need to pull it out.
    common_vars<-setdiff(common_vars,remove)
    for (i in 1:length(common_vars)) df$rnames<-gsub(paste0("^",common_vars[i],"$"),"common variant",df$rnames)#replaces all common variants as "common variant"
    #df<-df[which(df[1]/sum(df[1])>thresh | df$rnames=="no variant"),]# removes reads below threshold and keeps no variant##20190509 changed to use proportional threshold
    common_variant_sum<-sum(df[1][which(df$rnames=="common variant"),])#sums all common variant calls
    df$is_common_variant<-NA#adds a column that will hold the common variant marker
    df[nrow(df) + 1,] = list(common_variant_sum,"common variant sum","yes")#adds a new row to df with summed novariant calls.  Yes is a no variant marker
    df$is_no_variant<-NA#adds a column that will hold the no_variant marker
    df[which(df$rnames == "no variant"),4]<-"yes"
    df<-df[which(df$rnames!="common variant"),,drop=FALSE]#drop all constituents of no variant sum
    df<-df[which(df$rnames!="Other"),,drop=FALSE]#drop other reads
    df<-df[order(df$is_no_variant, df$is_common_variant,-df[1]),]# sort on no variant marker then on read counts
    df<-df[which(df[1]/sum(df[1])>thresh | df$is_no_variant=="yes" | df$is_common_variant=="yes"),]# removes reads below read count threshold, keeping no variant calls and common variant calls
    n_variants<-length(df$rnames)#total number of reads in each specimen
    total_reads<-sum(df[1])#sum of reads for each specimen
    allele_freq<-df[1]/total_reads# calculate allele freq by line
    names(allele_freq)[1] <- "allele_freq"# rename column for allele_freq
    df<-cbind(df,allele_freq)#add allele_freq column to the working df
    specimen_id<-names(df)[1]# pull the specimen name from the column title
    specimen<-rep_len(specimen_id, n_variants)# make a new column for the specimen id
    df<-cbind(df,specimen)# bind it to the working df
    #fix levels
    index<-as.factor(c(1:length(df$rnames)))#create index string as long as the number of variants, make it as a factor
    index<-factor(index, levels=rev(levels(index)))#reverse the levels
    df<-cbind(df,index)#bind it to the working df
    dir.create(paste0(output_folder,"/",class,"_thresh_",thresh,"_tvaf_",vaf_thresh,"_tpaf_",paf_thresh))
    write.csv(df,file = paste0(output_folder,"/",class,"_thresh_",thresh,"_tvaf_",vaf_thresh,"_tpaf_",paf_thresh,"/",class,"_",specimen_id,".csv"))
  }
  con_vc_all <- as.data.frame(variantCounts(con_crispr_set), stringsAsFactors=FALSE)#big data frame of read counts
  con_vc_list<-list()
  for(i in 1:length(con_mid_names)){con_vc_list[[length(con_vc_list)+1]]<-con_vc_all[i]}#puts individual samples into a list of vectors
  con_vc_prop<-variantCounts(con_crispr_set, result = "proportions")# generates a matrix of variant allele PERCENTS, similar to vc_all
  
  exp_vc_all <- as.data.frame(variantCounts(exp_crispr_set), stringsAsFactors=FALSE)#big data frame of read counts
  exp_vc_list<-list()
  for(i in 1:length(exp_mid_names)){exp_vc_list[[length(exp_vc_list)+1]]<-exp_vc_all[i]}#puts individual samples into a list of vectors
  exp_vc_prop<-variantCounts(exp_crispr_set, result = "proportions")# generates a matrix of variant allele PERCENTS, similar to vc_all
  
  lapply(con_vc_list,generate_v2, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = con_vc_prop,output_folder = output_folder, class = "control")
  lapply(exp_vc_list,generate_v2, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = exp_vc_prop,output_folder = output_folder, class = "experimental")
  
}


####function call####
generalized_func(output_folder = outs,
                 use_UMI = FALSE, 
                 genome_fp = "references/gestalt_pipeline4.fasta", 
                 ref_seq = "references/gestalt2_ref.txt",
                 control_name = "training", 
                 experimental_name = "validation", 
                 bed = "references/gestalt2.bed",
                 tvaf = tvaf,
                 tpaf = tpaf)

#stop cluster
stopCluster(cl)

####bootstrap analysis on control and experimental samples####
#bootstrap function to calculate bootstrap values and generate plots
bootstrap_func<-function(input_folder, seed, output_folder, outfile_cor, outfile_hist, h, w, report){
  #helper function for bootstrap analysis
  bootf<-function(data, indices){
    dt<-data[indices,]
    c(mean(dt[,2]),sd(dt[,3]))
  }
  set.seed(seed)
  
  infiles<-list.files(input_folder, pattern = "\\d.csv", full.names = TRUE)
  
  bootdata_0<-read.csv(infiles[1])
  fi<-1-sum(c(bootdata_0[1,2],bootdata_0[2,2]))/sum(bootdata_0[2])
  vaf<-c(NA,NA,bootdata_0[-c(1:2),2]/sum(bootdata_0[-c(1:2),2]))
  bootdata_1<-bootdata_0[,c(7,2,3)]
  colnames(bootdata_1)<-c("specimen", "counts", "variant")
  bootdata_1$FI<-rep(x = fi,times = nrow(bootdata_1))
  bootdata_1$vaf<-vaf
  bootdata_1$vaf2<-rep(sum(vaf>0.02, na.rm = TRUE), times = nrow(bootdata_1))
  for (i in 2:length(infiles)){
    newdf<-read.csv(infiles[i])
    fi<-1-sum(c(newdf[1,2],newdf[2,2]))/sum(newdf[2])
    vaf<-c(NA,NA,newdf[-c(1:2),2]/sum(newdf[-c(1:2),2]))
    newdf<-newdf[,c(7,2,3)]
    colnames(newdf)<-c("specimen", "counts", "variant")
    newdf$FI<-rep(x = fi,times = nrow(newdf))
    newdf$vaf<-vaf
    newdf$vaf2<-rep(sum(vaf>0.02, na.rm = TRUE), times = nrow(newdf))
    bootdata_1<-rbind(bootdata_1,newdf)
  }
  bootdata_2<-bootdata_1[,c(1,4,6)]#gets rid of the variant-level data for each sample since you don't need it here
  bootdata_3<-unique(bootdata_2)#gets rid of duplicate rows
  #bootstrap analysis starts here
  bootstrap<-boot(bootdata_3,bootf,R=1000)
  boot_df<-as.data.frame(bootstrap$t)
  cor_val<-cor(boot_df$V1,boot_df$V2, method = "pearson")
  cor_return<-cor.test(boot_df$V1,boot_df$V2, method = "pearson")
  boot_cor_plot<-ggplot(boot_df, aes(x = boot_df$V1, y = boot_df$V2))
  boot_cor_plot<-boot_cor_plot+
    geom_density2d(color = "black")+
    theme_cowplot(font_size = 11)+
    labs(x = "Mean Fraction Informative", y = expression(paste("S.D. ",N["Variants"],"\nVAF>0.02",sep = "")))
  boot_cor_plot
  save_plot(plot = boot_cor_plot, filename = paste0(output_folder,"/",outfile_cor), base_height = h, base_width = w)
  #plot distribution
  boot_ci<-boot.ci(bootstrap, index = 1)#makes a list of the bootstrap confidence intervals
  boot_plot<-ggplot(boot_df, aes(x = boot_df$V1))
  boot_plot<-boot_plot+
    geom_histogram(binwidth = 0.01, fill = "white", color = "black")+
    geom_vline(aes(xintercept = mean(boot_df$V1), color = "Estimate"), linetype = "dashed")+#plots bootstrap estimate of mean informative fraction
    geom_vline(aes(xintercept = boot_ci$bca[1,4], color = "Lower\nBound"), linetype = "dashed")+#plots lower bound of 95% ci
    scale_color_manual(name = "", values = c('Lower\nBound' = "#DC0000",Estimate = "#3c5488"))+
    theme_cowplot(font_size = 11)+
    theme(legend.position = c(0.6,0.85))+
    labs(y = "Replicates", x = "Mean Informative Fraction")
  boot_plot
  save_plot(plot = boot_plot, filename = paste0(output_folder,"/",outfile_hist), base_height = h, base_width = w)
  returnlist<-list("bootstrap report:  " = bootstrap, 
                   "correlation of FI and SD over replicates:  " = cor_return, 
                   "95% ci of bootstrap values:  " = boot_ci, 
                   "estimate:  " = mean(boot_df$V1))
  capture.output(print(returnlist), file = paste0(output_folder,"/",report)) 
  return(boot_ci$bca[1,4])
}  

control_boot_lowerbound<-bootstrap_func(input_folder = paste0(outs,"/control_thresh_0_tvaf_",tvaf,"_tpaf_",tpaf), 
               seed = 12345, 
               output_folder = outs,
               outfile_cor = "control_bootstrap_correlation_a5e96.pdf", 
               outfile_hist = "control_bootstrap_estimate_mif_6acfd.pdf", 
               report = "control_bootstrap_report.txt", 
               h = 2.3, w = 2.5)

experimental_boot_lowerbound<-bootstrap_func(input_folder = paste0(outs,"/experimental_thresh_0_tvaf_",tvaf,"_tpaf_",tpaf),
               seed = 54321,
               output_folder = outs,
               outfile_cor = "experimental_bootstrap_correlation_3cfbf.pdf",
               outfile_hist = "experimental_bootstrap_estimate_mif_74d44.pdf",
               report = "experimental_bootstrap_report.txt",
               h = 2.5, w = 3)

####re-run training samples using condition tvaf = 0.1, tpaf = 0.1 and identify informative vs uninformative samples using kmeans####
fi_func<-function(input_folder, output_folder,cutoff, fi_plotname, h, w, informative_vec){
  #load the data
  infiles<-list.files(input_folder, pattern = "\\d.csv", full.names = TRUE)
  fidf_0<-read.csv(infiles[1])
  fi<-1-sum(c(fidf_0[1,2],fidf_0[2,2]))/sum(fidf_0[2])
  fidf_1<-fidf_0[,c(7,2,3)]
  colnames(fidf_1)<-c("specimen", "counts", "variant")
  fidf_1$FI<-rep(x = fi,times = nrow(fidf_0))
  for (i in 2:length(infiles)){
    newdf<-read.csv(infiles[i])
    fi<-1-sum(c(newdf[1,2],newdf[2,2]))/sum(newdf[2])
    newdf<-newdf[,c(7,2,3)]
    colnames(newdf)<-c("specimen", "counts", "variant")
    newdf$FI<-rep(x = fi,times = nrow(newdf))
    fidf_1<-rbind(fidf_1,newdf)
  }
  fidf_2<-fidf_1[,c(1,4)]#gets rid of the variant-level data for each sample since you don't need it here
  fidf_3<-unique(fidf_2)#gets rid of duplicate rows
  fidf_3$rank<-rank(-(fidf_3$FI), ties.method = "first")#ranks samples based on fraction informative with 1 = most informative
  
  # #kmeans
  # km<-kmeans(fidf_3$FI, centers = 2)#generates kmeans of 2 clusters based on fraction informative
  # fidf_3$cluster<-as.factor(km$cluster)#puts cluster vector into the fi dataframe as a new column
  # 
  # top_clust<-fidf_3[which(fidf_3$rank==1),4]#identifies the top cluster since kmeans function assigns cluster 1 at random
  # bottom_clust<-fidf_3[which(fidf_3$rank==length(fidf_3$rank)),4]#identifies bottom cluster
  # fidf_3$cluster<-as.character(fidf_3$cluster)
  # fidf_3$cluster[fidf_3$cluster==top_clust]<-"Informative"#relabel clusters so they are consistent
  # fidf_3$cluster[fidf_3$cluster==bottom_clust]<-"Not Informative"#ditto
  
  #use bootstrap to define clusters
  top_boot_clust<-fidf_3[which(fidf_3$FI>=cutoff),]
  top_boot_clust$boot_clust<-rep("High FI", times = nrow(top_boot_clust))
  bottom_boot_clust<-fidf_3[which(fidf_3$FI<cutoff),]
  bottom_boot_clust$boot_clust<-rep("Low FI", times = nrow(bottom_boot_clust))
  fidf_3.1<-rbind(top_boot_clust,bottom_boot_clust)
 
  #plot the fi distribution with kmeans
  fi_dist<-ggplot(fidf_3.1, aes(x=factor(fidf_3.1$rank), y=fidf_3.1$FI, fill=fidf_3.1$boot_clust))
  fi_dist<-fi_dist+
    geom_bar(stat = "identity", color = "black")+
    scale_fill_manual(values = c("#72ce55","#fde725"),name = "",labels = c(expression(paste("High ",phi),paste("Low ",phi))))+
    scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
    geom_hline(yintercept = cutoff, color = "#DC0000", linetype = "dashed")+
    theme_cowplot(font_size = 11)+
    theme(legend.position = c(0.7,0.95))+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
    labs(x="Sample Index", y=expression(paste(phi)))
  fi_dist
  save_plot(plot = fi_dist, filename = paste0(output_folder,"/",fi_plotname), base_width = w, base_height = h)
  informative_df<-top_boot_clust
  informative_spec_vec<-informative_df$specimen
  return(informative_spec_vec)
}

control_informative<-fi_func(input_folder = paste0(outs,"/control_thresh_0_tvaf_",tvaf,"_tpaf_",tpaf),
                             cutoff = control_boot_lowerbound,
                             output_folder = outs,
                             fi_plotname = "training_fiplot_a132b.pdf",
                             h = 3,w = 4)

experimental_informative<-fi_func(input_folder = paste0(outs,"/experimental_thresh_0_tvaf_",tvaf,"_tpaf_",tpaf),
                             cutoff = experimental_boot_lowerbound,
                             output_folder = outs,
                             fi_plotname = "validation_fiplot_149f4.pdf",
                             h = 2.5,w = 3)

####quantifying alleles from thresholded data in informative samples####
#reassemble the long form data frame of read counts using only informative samples
#set output folder for final plots
output_folder<-outs
control_label<-"Training"
experimental_label<-"Validation"



control_allfiles<-as.data.frame(list.files(paste0(outs,"/control_thresh_0_tvaf_",tvaf,"_tpaf_",tpaf), full.names = TRUE))
control_allfiles$specimen<-as.factor(str_sub(as.character(control_allfiles[,1]),-9,-5))
colnames(control_allfiles)<-c("path","specimen")
control_informative<-as.data.frame(control_informative); colnames(control_informative)<-"specimen"
control_inlist<-left_join(control_informative, control_allfiles, by = "specimen")


experimental_allfiles<-as.data.frame(list.files(paste0(outs,"/experimental_thresh_0_tvaf_",tvaf,"_tpaf_",tpaf), full.names = TRUE))
experimental_allfiles$specimen<-as.factor(str_sub(as.character(experimental_allfiles[,1]),-9,-5))
colnames(experimental_allfiles)<-c("path","specimen")
experimental_informative<-as.data.frame(experimental_informative); colnames(experimental_informative)<-"specimen"
experimental_inlist<-left_join(experimental_informative, experimental_allfiles, by = "specimen")

con_qdf0<-read.csv(as.character(control_inlist$path[1]))
colnames(con_qdf0)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
for (i in 2:nrow(control_inlist)) {
  qdf_new<-read.csv(as.character(control_inlist$path[i]))
  colnames(qdf_new)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
  con_qdf0<-rbind(con_qdf0,qdf_new)
}
con_qdf0$class<-rep("control", times = nrow(con_qdf0))
con_qdf1<-con_qdf0[,-1]

exp_qdf0<-read.csv(as.character(experimental_inlist$path[1]))
colnames(exp_qdf0)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
for (i in 2:nrow(experimental_inlist)) {
  qdf_new<-read.csv(as.character(experimental_inlist$path[i]))
  colnames(qdf_new)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
  exp_qdf0<-rbind(exp_qdf0,qdf_new)
}
exp_qdf0$class<-rep("experimental", times = nrow(exp_qdf0))
exp_qdf1<-exp_qdf0[,-1]

qdf1<-rbind(con_qdf1,exp_qdf1)

spec_df<-unique(qdf1[,c(6,8)])
spec_list<-spec_df$specimen
class_list<-spec_df$class

#generate a vector of shannon entropy values for the informative samples
shannon<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  shannon<-c(shannon,shannon_func(qdf_filt_inf[1]))#calculates off of read count
}
shannon

#generate a vector of simpson diversity indices for the informative samples
simpson<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  simpson<-c(simpson,simpson_func(qdf_filt_inf[1]))#calculates off of read count
}
simpson

#generate a vector of the number of variants over 2%
vaf2<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  qdf_filt_inf$adj_allele_freq<-qdf_filt_inf$read_count/sum(qdf_filt_inf$read_count)#recalculates allele frequencies from read counts without no variant and common variant
  qdf_2percent<-filter(qdf_filt_inf, qdf_filt_inf$adj_allele_freq>=0.02)
  vaf2<-c(vaf2,nrow(qdf_2percent))
}
vaf2

#counts barcodes over 1% raw vaf
vaf1<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  qdf_filt_inf$adj_allele_freq<-qdf_filt_inf$read_count/sum(qdf_filt_inf$read_count)#recalculates allele frequencies from read counts without no variant and common variant
  qdf_1percent<-filter(qdf_filt_inf, qdf_filt_inf$adj_allele_freq>=0.01)
  vaf1<-c(vaf1,nrow(qdf_1percent))
}
vaf1

#compile and summarize data
sumdf0<-as.data.frame(cbind(spec_list,class_list,vaf1,vaf2,shannon,simpson))
colnames(sumdf0)<-c("Sample","Class","variants_over_1%", "variants_over_2%", "Shannon_Entropy", "Simpson_Diversity")
sumdf1<-sumdf0


sumdf_long<-gather(sumdf1, key = "attribute", value = "value_to_plot", -c(Sample,Class))
sumdf_long$value_to_plot<-as.numeric(sumdf_long$value_to_plot)
vaf1_long<-filter(sumdf_long, sumdf_long$attribute=="variants_over_1%")
vaf2_long<-filter(sumdf_long, sumdf_long$attribute=="variants_over_2%")
shannon_long<-filter(sumdf_long, sumdf_long$attribute=="Shannon_Entropy")
simpson_long<-filter(sumdf_long, sumdf_long$attribute=="Simpson_Diversity")

p1<-ggplot(vaf1_long, aes(x = vaf1_long$Class, y = vaf1_long$value_to_plot))
p1<-p1+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=vaf1_long$Class)) + 
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_x_discrete(breaks = c("control","experimental"), labels = c(control_label, experimental_label))+
  labs(y="Clones with VAF>0.01")
p1
save_plot(plot = p1, filename = paste0(output_folder,"/vaf_0.01_60412.pdf"),base_width = 2.5, base_height = 2.3)

p2<-ggplot(vaf2_long, aes(x = vaf2_long$Class, y = vaf2_long$value_to_plot))
p2<-p2+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=vaf2_long$Class)) + 
  #geom_boxplot()+
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  scale_x_discrete(breaks = c("control","experimental"), labels = c(control_label, experimental_label))+
  labs(y="Clones with VAF>0.02")
p2
save_plot(plot = p2, filename = paste0(output_folder,"/vaf_0.02_401af.pdf"),base_width = 3, base_height = 2.5)

p3<-ggplot(shannon_long, aes(x = shannon_long$Class, y = shannon_long$value_to_plot))
p3<-p3+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=shannon_long$Class)) + 
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_x_discrete(breaks = c("control","experimental"), labels = c(control_label, experimental_label))+
  labs(y="Shannon Entropy")
p3
save_plot(plot = p3, filename = paste0(output_folder,"/shannon_ee48c.pdf"),base_width = 2.5, base_height = 2.3)

p4<-ggplot(simpson_long, aes(x = simpson_long$Class, y = simpson_long$value_to_plot))
p4<-p4+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=simpson_long$Class)) + 
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_x_discrete(breaks = c("control","experimental"), labels = c(control_label, experimental_label))+
  labs(y="Simspon Diversity")
p4
save_plot(plot = p4, filename = paste0(output_folder,"/simpson_36475.pdf"),base_width = 2.5, base_height = 2.3)


####stat report####
df_con<-sumdf_long[which(sumdf_long$Class=="control"),]
df_exp<-sumdf_long[which(sumdf_long$Class=="experimental"),]

genstat<-function(x,y,stat,output_folder){
  xvec<-x[which(x$attribute==stat),4]
  yvec<-y[which(y$attribute==stat),4]
  ks<-ks.test(xvec,yvec)
  shapirox<-shapiro.test(xvec)
  shapiroy<-shapiro.test(yvec)
  t_equal<-t.test(xvec,yvec,alternative = "two.sided", var.equal = TRUE)
  t_unequal<-t.test(xvec,yvec,alternative = "two.sided", var.equal = FALSE)
  wilcox<-wilcox.test(xvec,yvec,alternative = "two.sided")
  returnlist<-list(paste0("X is ", unique(x$Class)),xvec,paste0("SD of X is ",sd(xvec)),paste0("SEM of X is ",se(xvec)),
                   paste0("Y is ", unique(y$Class)),yvec,paste0("SD of Y is ",sd(yvec)),paste0("SEM of Y is ",se(yvec)),
                   paste0("Statistic is ",stat),
                   shapirox,shapiroy,ks,t_equal,t_unequal,wilcox)
  sink(paste0(output_folder,"/stat_report_",stat,".txt"))
  return(returnlist)
}

genstat(x = df_con, y = df_exp, stat = "variants_over_1%", output_folder = outs);sink()
genstat(x = df_con, y = df_exp, stat = "variants_over_2%", output_folder = outs);sink()
genstat(x = df_con, y = df_exp, stat = "Shannon_Entropy", output_folder = outs);sink()
genstat(x = df_con, y = df_exp, stat = "Simpson_Diversity", output_folder= outs);sink()

#####generate sharing factor plot for control and experimental sets####
cor_func4<-function(control_or_experimental, informative_df, thresh, tvaf, tpaf, method, tag, output_folder){
  input_directory<-paste0(output_folder,"/",control_or_experimental,"_thresh_",thresh,"_tvaf_",tvaf,"_tpaf_",tpaf)
  #inlist<-list.files(input_directory,full.names = TRUE, pattern = "\\d.csv$")
  inlist<-list.files(paste0(output_folder,"/",control_or_experimental,"_thresh_",thresh,"_tvaf_",tvaf,"_tpaf_",tpaf),full.names = TRUE)
  insample<-str_sub(inlist,-9,-5)
  inframe<-as.data.frame(cbind(inlist,insample))
  inlist1<-left_join(informative_df,inframe,by = c("specimen" = "insample"))
  inlist1$inlist<-as.character(inlist1$inlist)
  df_cor<-data.frame(Variant = character(), Count = integer(), Specimen = character())
  for (i in 1:nrow(inlist1)){
    df_new<-read.csv(file = inlist1$inlist[i], header = TRUE)
    df_new<-df_new[,c(3,2,7)]
    colnames(df_new)<-c("Variant", "Count", "Specimen")
    df_cor<-rbind(df_cor,df_new)
  }
  df_cor_wide0<-reshape(df_cor, idvar = "Variant", timevar = "Specimen", direction = "wide")
  df_cor_wide<-df_cor_wide0
  df_cor_wide[is.na(df_cor_wide)]<-0
  colnames(df_cor_wide)<-substr(colnames(df_cor_wide),7,11)
  rownames(df_cor_wide)<-df_cor_wide[,1]
  df_cor_wide<-df_cor_wide[,-1]
  df_cor_wide1<-df_cor_wide[-c(1,2),]

  stat_matrix<-matrix(vector(),nrow = ncol(df_cor_wide1), ncol = ncol(df_cor_wide1))
  colnames(stat_matrix)<-colnames(df_cor_wide1)
  rownames(stat_matrix)<-colnames(df_cor_wide1)

  for (i in 1:ncol(df_cor_wide1)) {
    for (j in 1:ncol(df_cor_wide1)) {
      a<-pmin(df_cor_wide1[,i],df_cor_wide1[,j])
      b<-pmax(df_cor_wide1[,i],df_cor_wide1[,j])
      c<-(sum(a)*2)/(sum(a)+sum(b))
      stat_matrix[i,j]<-c

    }

  }

  reorder_cormat <- function(x){
    # Use correlation between variables as distance
    dd <- as.dist((1-x)/2)
    hc <- hclust(dd)
    x <-x[hc$order, hc$order]
  }

  get_lower_tri <- function(x){
    x[upper.tri(x)]<- NA
    return(x)
  }

  sm_reordered<-reorder_cormat(stat_matrix)
  sm_lower<-get_lower_tri(sm_reordered)

  msm<-melt(sm_lower)

  msm_plot<-ggplot(data = msm, aes(x = msm$Var1, y = msm$Var2, fill = msm$value))
  msm_plot<-msm_plot+
    geom_tile()+
    #scale_fill_viridis_c(begin = 0, end = 1, guide = "colourbar", aesthetics = "fill", na.value = "white")+
    scale_fill_distiller(palette = "RdYlBu", na.value = "transparent",guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE),limits = c(0,1))+
    #scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
    theme_cowplot(font_size = 11)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))+
    theme(axis.text.y = element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(axis.line.y = element_blank())+
    theme(axis.ticks.y = element_blank())+
    xlab("Sample")+
    theme(axis.title.y = element_blank())+
    theme(plot.title = element_blank())+
    theme(plot.subtitle = element_blank())+
    guides(fill=guide_colorbar(title="Sharing\nFactor", reverse = F))+
    theme(legend.position = c(0,0.7))+
    coord_fixed()
  msm_plot

  #calculate Fraction informative
  FI<-as.data.frame(colSums(df_cor_wide[-c(1,2),])/colSums(df_cor_wide))#ratio of reads excluding no variant and common variant to all reads in each sample.
  colnames(FI)<-"Fraction_Informative"
  FI$sample_x<-rownames(FI)

  h.mean<-function(x){
    harmonic<-1/mean(1/x)
    return(harmonic)}

  avg_share<-round(mean(msm[which(msm$Var1!=msm$Var2),3],na.rm = TRUE),4)#averages pearson correlation excluding identical comparisons

  mean_inf<-round(mean(FI$Fraction_Informative),3)

  save_plot(plot = msm_plot, filename = paste0(output_folder,"/share_factor_",control_or_experimental,"_",tag,".pdf"), base_width = 2.75)
  nubbin<-data.frame(threshold = thresh, tvaf = tvaf, tpaf = tpaf, mean_inf = mean_inf, avg_share = avg_share)
  write.csv(x = nubbin, file = paste0(input_directory,"/nubbin_invormative_only_",method,".csv"))
  #return(inlist1)
}


cor_func4(control_or_experimental = "control", informative_df = control_informative, thresh = 0, tvaf = tvaf, tpaf = tpaf, method = "blaser", tag = "5614d", output_folder = output_folder)
cor_func4(control_or_experimental = "experimental", informative_df = experimental_informative, thresh = 0, tvaf = tvaf, tpaf = tpaf, method = "blaser", tag = "33613", output_folder = output_folder)

####clean up folders####
unlink("temp", recursive = TRUE)
unlink("fastq_control", recursive = TRUE);dir.create("fastq_control")
unlink("fastq_experimental", recursive = TRUE);dir.create("fastq_experimental")