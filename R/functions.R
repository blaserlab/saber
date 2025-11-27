se <- function(x)
  sqrt(var(x, na.rm = TRUE) / length(x))

data_summary <- function(x) {
  m <- median(x)
  ymin <- m - se(x)
  ymax <- m + se(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

data_summary1 <- function(x) {
  m <- median(x)
  ymin <- m - (IQR(x) / 2)
  ymax <- m + (IQR(x) / 2)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

data_summary2 <- function(x) {
  m <- mean(x)
  ymin <- m - se(x)
  ymax <- m + se(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

data_summary3 <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)
  ymax <- m + sd(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

shannon_func <- function(x) {
  prop <- x / sum(x)
  shannon_entropy <- -1 * sum(prop * log2(prop))
  return(shannon_entropy)
}

simpson_func <- function(x) {
  prop <- x / sum(x)
  simpson_diversity <- 1 / sum(prop^2)
  return(simpson_diversity)
}

general_diversity_func <- function(x, q) {
  prop <- x / sum(x)
  general_diversity <- (sum(prop^q))^(1 / (1 - q))
  return(general_diversity)
}


#' @importFrom fs dir_create path_temp path
make_tree <- function(output_folder) {
  temp_folder <- fs::path_temp()
  fs::dir_create(temp_folder)

  fastq_split_R1_folder <- fs::path(temp_folder, "fastq_split_R1")
  fastq_split_R2_folder <- fs::path(temp_folder, "fastq_split_R2")

  fs::dir_create(fastq_split_R1_folder)
  fs::dir_create(fastq_split_R2_folder)

  bam_temp <- fs::path(temp_folder, "bam_temp")

  fs::dir_create(bam_temp)

  fastqc_files_folder <- fs::path(output_folder, "fastqc_files")
  fs::dir_create(fastqc_files_folder)

  bam_output_folder <- fs::path(output_folder, "bam_output")

  fs::dir_create(bam_output_folder)

  PEAR_output_folder <- fs::path(output_folder, "PEAR_output")
  fs::dir_create(PEAR_output_folder)

  fastq_trimmed <- fs::path(temp_folder, "fastq_trimmed")
  fs::dir_create(fastq_trimmed)

  fastq_merged <- fs::path(temp_folder, "fastq_merged")
  fs::dir_create(fastq_merged)

  sam_temp <- fs::path(temp_folder, "sam_temp")
  fs::dir_create(sam_temp)

  sam_fixed <- fs::path(temp_folder, "sam_fixed")
  fs::dir_create(sam_fixed)

  cutadapt1_folder <- fs::path(temp_folder, "cutadapt1")
  cutadapt2_folder <- fs::path(temp_folder, "cutadapt2")
  cutadapt_folder <- fs::path(temp_folder, "cutadapt")
  fs::dir_create(cutadapt1_folder)
  fs::dir_create(cutadapt2_folder)
  fs::dir_create(cutadapt_folder)



}

#' @importFrom fs path_temp dir_ls
#' @importFrom readr read_csv
load_gestalt_fastqs <- function(config,
                                fastq_folder = fs::path_temp("fastq"),
                                fastq_split_R1_folder = fs::path_temp("fastq_split_R1"),
                                fastq_split_R2_folder = fs::path_temp("fastq_split_R2")) {
  for (i in 1:length(config$samples_fp_R1)) {
    cmd <- paste0("cp ", config$samples_fp_R1[i], " ", fastq_split_R1_folder)
    message(cmd, "\n")
    system(cmd)
  }

  for (i in 1:length(config$samples_fp_R2)) {
    cmd <- paste0("cp ", config$samples_fp_R2[i], " ", fastq_split_R2_folder)
    message(cmd, "\n")
    system(cmd)
    return(config)
  }



  # #unzip
  for (i in 1:length(fs::dir_ls(fastq_split_R1_folder))) {
    cmd <- paste0("gzip -d ", fs::dir_ls(fastq_split_R1_folder)[i])
    message(cmd, "\n")
    system(cmd)
  }

  # #unzip
  for (i in 1:length(fs::dir_ls(fastq_split_R2_folder))) {
    cmd <- paste0("gzip -d ", fs::dir_ls(fastq_split_R2_folder)[i])
    message(cmd, "\n")
    system(cmd)
  }




}


#' @importFrom fs path_temp path dir_ls path_file
merge_read_pairs <- function(output_folder,
                             conf,
                             fastq_split_R1_folder = fs::path_temp("fastq_split_R1"),
                             fastq_split_R2_folder = fs::path_temp("fastq_split_R2"),
                             fastq_merged = fs::path_temp("fastq_merged")) {
  PEAR_output_folder <- fs::path(output_folder, "PEAR_output")
  fastq_split_R1_fp <- fs::dir_ls(fastq_split_R1_folder)
  fastq_split_R1_files <- fs::path_file(fastq_split_R1_fp)
  fastq_split_R2_fp <- fs::dir_ls(fastq_split_R2_folder)
  fastq_split_R2_files <- fs::path_file(fastq_split_R2_fp)

  # merge read pairs with PEAR
  mid_names <- conf$sample
  for (i in 1:length(fastq_split_R1_fp)) {
    cmd <-
      paste0(
        "pear -f ",
        fastq_split_R1_fp[i],
        " -r ",
        fastq_split_R2_fp[i],
        " -o ",
        fastq_merged,
        "/",
        mid_names[i],
        " -j 39 > ",
        PEAR_output_folder,
        "/",
        mid_names[i],
        ".PEARreport.txt"
      )
    message(cmd, "\n")
    system(cmd)
  }
  unlink(paste0(fastq_merged, "/*unassembled*"), recursive = TRUE)
  unlink(paste0(fastq_merged, "/*discarded*"), recursive = TRUE)

}


#' @importFrom fs path_temp dir_ls path_file
trim_reads <- function(trimmomatic_path = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
                       conf,
                       fastq_merged = fs::path_temp("fastq_merged"),
                       fastq_trimmed = fs::path_temp("fastq_trimmed")) {
  fastq_merged_fp <-
    fs::dir_ls(fastq_merged)
  fastq_merged_files <-
    fs::path_file(fastq_merged_fp)


  #run trimmomatic on split and ordered files
  mid_names <- conf$sample
  for (i in 1:length(fastq_merged_fp)) {
    cmd <-
      paste0(
        "java -jar ",
        trimmomatic_path,
        " SE ",
        #invoke trimmomatic
        fastq_merged_fp[i],
        " ",
        # merged read input
        fastq_trimmed,
        "/",
        mid_names[i],
        ".trimmed.fastq ",
        # merged trimmed output
        "SLIDINGWINDOW:4:15 MINLEN:100"
      )# trim parameters.
    message(cmd, "\n")
    system(cmd)
  }
}


#' @importFrom fs path_temp dir_ls path_file
run_cutadapt <- function(conf,
                         multiplex,
                         fastq_trimmed = fs::path_temp("fastq_trimmed"),
                         cutadapt1_folder = fs::path_temp("cutadapt1"),
                         cutadapt2_folder = fs::path_temp("cutadapt2"),
                         n_cores = max(1, parallel::detectCores() - 1)) {
  fastq_trimmed_fp <-
    fs::dir_ls(fastq_trimmed)
  fastq_trimmed_files <-
    fs::path_file(fastq_trimmed_fp)
  #run cutadapt to keep only fastqs that have the full flanking primer sequences
  #5 prime adapter = V6/7F:  TCGAGCTCAAGCTTCGG
  mid_names <- conf$sample

  if (multiplex) {
    purrr::walk2(.x = fastq_trimmed_fp, .y = mid_names, .f = \(
      x,
      y,
      rf = system.file("extdata/references/barcodes.fasta", package = "saber")
    ) {
      cmd <- paste0(
        "cutadapt -g file:",
        rf,
        " -o ",
        fs::path_temp("cutadapt"),
        "/{name}.",
        y,
        ".cutadapt.fastq ",
        x
      )
      message(cmd, "\n")
      system(cmd)
    })
    fastq_trimmed_fp <-
      fs::dir_ls(fs::path_temp("cutadapt"))

    # remove unknown sequences
    unknown <- fs::dir_ls(fs::path_temp("cutadapt"))[stringr::str_detect(fs::dir_ls(fs::path_temp("cutadapt")), "unknown")]
    fs::file_delete(unknown)
    mid_names <- stringr::str_remove(fs::path_file(fs::dir_ls(fs::path_temp("cutadapt"))), ".cutadapt.fastq")
  }

  parallel::mclapply(
    X = seq_along(fastq_trimmed_fp),
    FUN = function(i) {
      cmd <-
        paste0(
          "cutadapt -g ^TCGAGCTCAAGCTTCGG --discard-untrimmed -e 0.01 --action=none -o ",
          cutadapt1_folder,
          "/",
          mid_names[i],
          ".cutadapt1.fastq ",
          fastq_trimmed_fp[i]
        )
      message(cmd, "\n")
      system(cmd)
      invisible(NULL)
    },
    mc.cores = n_cores
  )

  fs::file_delete(fs::path_temp("cutadapt1", "NA.cutadapt1.fastq"))

  cutadapt1_fp <-
    fs::dir_ls(cutadapt1_folder)
  cutadapt1_files <-
    fs::path_file(cutadapt1_fp)

  #3 prime adapter = V6/7R:  GACCTCGAGACAAATGGCAG (reverse complement of the primer sequence 5'-3')
  parallel::mclapply(
    X = seq_along(cutadapt1_fp),
    FUN = function(i) {
      cmd <-
        paste0(
          "cutadapt -a GACCTCGAGACAAATGGCAG$ --discard-untrimmed -e 0.01 --action=none -o ",
          cutadapt2_folder,
          "/",
          mid_names[i],
          ".cutadapt2.fastq ",
          cutadapt1_fp[i]
        )
      message(cmd, "\n")
      system(cmd)
      invisible(NULL)
    },
    mc.cores = n_cores
  )

  if (fs::file_exists(fs::path_temp("cutadapt2", "NA.cutadapt2.fastq")))
    fs::file_delete(fs::path_temp("cutadapt2", "NA.cutadapt2.fastq"))


}


#' @importfrom fs path_temp dir_ls path_file path
#' @importFrom parallel mclapply detectCores
run_fastqc <- function(output_folder,
                       multiqc_path = "/workspace/python/anaconda3/bin/multiqc",
                       cutadapt2_folder = fs::path_temp("cutadapt2"),
                       n_cores = max(1, parallel::detectCores() - 1)) {
  cutadapt2_fp <-
    fs::dir_ls(cutadapt2_folder)
  cutadapt2_files <-
    fs::path_file(cutadapt2_fp)
  fastqc_files_folder <- fs::path(output_folder, "fastqc_files")
  parallel::mclapply(
    X = seq_along(cutadapt2_fp),
    FUN = function(i) {
      cmd <-
        paste0("fastqc -o ", fastqc_files_folder, " ", cutadapt2_fp[i])
      message(cmd, "\n")
      system(cmd)
      invisible(NULL)
    },
    mc.cores = n_cores
  )


  # run multiqc to summarize the qc files
  system(paste0(multiqc_path, " -d ", output_folder, " -o ", output_folder))

  invisible(NULL)

}



#' @importFrom fs path_temp path_package dir_ls path_file
#' @importFrom parallel mclapply detectCores
align_reads <- function(cutadapt2_folder = fs::path_temp("cutadapt2"),
                        sam_temp        = fs::path_temp("sam_temp"),
                        output_folder,
                        genome_fp       = fs::path_package("saber", "extdata", "references", "gestalt_pipeline4.fasta"),
                        n_cores         = max(1, parallel::detectCores() - 1)) {
  # list input files
  cutadapt2_fp    <- fs::dir_ls(cutadapt2_folder)
  cutadapt2_files <- fs::path_file(cutadapt2_fp)
  mid_names       <- stringr::str_remove(cutadapt2_files, "\\.cutadapt2\\.fastq$")

  # run needleall in parallel
  parallel::mclapply(
    X = seq_along(cutadapt2_files),
    FUN = function(i) {
      cmd <- paste0(
        "needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence ",
        genome_fp,
        " -bsequence ",
        cutadapt2_fp[i],
        " -outfile ",
        sam_temp,
        "/",
        mid_names[i],
        ".sam",
        " -errfile ",
        output_folder,
        "/con_needleall.error"
      )
      message(cmd, "\n")
      system(cmd)
      invisible(NULL)
    },
    mc.cores = n_cores
  )

  invisible(NULL)
}


#' @importFrom fs path_temp path_package dir_ls path_file
#' @importFrom readr write_tsv
fix_sam_header <- function(sam_temp = fs::path_temp("sam_temp"),
                           sam_fixed = fs::path_temp("sam_fixed"),
                           new_header = fs::path_package(package = "saber", "extdata", "references", "sam_header.csv")) {
  sam_temp_fp <- fs::dir_ls(sam_temp)
  sam_temp_files <-
    fs::path_file(sam_temp_fp)


  # fix sam file header and select only reads matching at position 1
  for (i in 1:length(sam_temp_files)) {
    samdf <-
      read.delim(sam_temp_fp[i], sep = "\t", header = FALSE)#read in sam file
    samdf <- samdf[-(1:2), ]#chop off the old header
    samdf <-
      samdf[which(samdf$V4 == "1"), ]#select only reads mapping to coordinate 1
    sam_header <-
      read.table(
        new_header,
        fill = TRUE,
        header = FALSE,
        sep = ",",
        colClasses = (rep("character", 13))
      )# read in standard sam header
    names(sam_header) <- paste("V", 1:13, sep = "")
    samdf <- rbind(sam_header, samdf)
    readr::write_tsv(
      samdf,
      na = "",
      file = sam_temp_fp[i],
      col_names = FALSE,
      append = FALSE
    )
    cmd <- paste0(
      "awk '{ sub(/[ \\t]+$/, \"\"); print }' ",
      sam_temp_fp[i],
      " > ",
      sam_fixed,
      "/",
      sam_temp_files[i]
    )
    message(cmd, "\n")
    system(cmd)

  }

}


#' @importFrom fs dir_ls path_file
sam_to_bam <- function(sam_fixed = fs::path_temp("sam_fixed"),
                       bam_temp = fs::path_temp("bam_temp"),
                       conf) {
  sam_fixed_fp <- fs::dir_ls(sam_fixed)
  sam_fixed_files <- fs::path_file(sam_fixed_fp)
  mid_names <- stringr::str_remove(sam_fixed_files, ".sam")
  # convert sam to bam
  for (i in 1:length(sam_fixed_files)) {
    cmd <-
      paste0("samtools view -S -b ",
             sam_fixed_fp[i],
             " > ",
             bam_temp,
             "/",
             mid_names[i],
             ".bam")
    message(cmd, "\n")
    system(cmd)
  }

}


#' @importFrom fs dir_ls path_file path
sort_and_index <- function(bam_temp = fs::path_temp("bam_temp"),
                           output_folder) {
  bam_temp_fp <- fs::dir_ls(bam_temp)
  bam_temp_files <- fs::path_file(bam_temp_fp)
  bam_output_folder <- fs::path(output_folder, "bam_output")

  # sort and index bam
  for (i in 1:length(bam_temp_fp)) {
    cmd <-
      paste0("samtools sort ", bam_temp_fp[i], " -o ", bam_temp_fp[i])
    message(cmd, "\n")
    system(cmd)
  }

  for (i in 1:length(bam_temp_fp)) {
    cmd <- paste0("samtools index ", bam_temp_fp[i])
    message(cmd, "\n")
    system(cmd)
  }

  system(paste0("cp -r ", fs::path_temp("bam_temp"), " ", bam_output_folder))


}




#' @importFrom fs dir_ls path
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges resize width
#' @importFrom readr read_file
#' @importFrom Biostrings DNAString
#' @importFrom CrispRVariants readsToTarget
make_crispr_set <- function(bam_temp = fs::path_temp("bam_temp"),
                            conf,
                            output_folder) {
  #build metadata for experiment
  if (fs::file_exists(fs::path_temp("bam_temp", "NA.bam")))
    fs::file_delete(fs::path_temp("bam_temp", "NA.bam"))
  if (fs::file_exists(fs::path_temp("bam_temp", "NA.bam.bai")))
    fs::file_delete(fs::path_temp("bam_temp", "NA.bam.bai"))
  bam_fp <- fs::dir_ls(bam_temp, regexp = ".bai", invert = TRUE)
  bam_files <- fs::path_file(bam_fp)
  mid_names <- stringr::str_remove(bam_files, ".bam")
  group_desig <- rep("group", times = length(mid_names))
  # md <- read.csv("references/blank_metadata.csv", header = TRUE)
  newrow <- data.frame(
    bamfile = bam_fp,
    directory = getwd(),
    Short.name = mid_names,
    Targeting.type = "",
    sgRNA1 = "",
    sgRNA2 = "",
    Group = group_desig
  )
  md <- rbind(md, newrow)


  #create target region
  # gd <- rtracklayer::import("references/gestalt2.bed")
  gdl <-
    GenomicRanges::resize(gd, GenomicRanges::width(gd) + 0, fix = "center") #resize region for analysis
  # reference0 <- readr::read_file("references/gestalt2_ref.txt")
  reference1 <-
    substr(reference0, 1, 310)#this has to be here to remove \n.  Can generalize in future
  reference <- Biostrings::DNAString(reference1)


  # make the crispr set
  crispr_set <- CrispRVariants::readsToTarget(
    bam_fp,
    target = gd,
    reference = reference,
    names = mid_names,
    target.loc = 16,
    #generalize with distance from primer to first sequence of barcode?
    collapse.pairs = FALSE,
    split.snv = TRUE
  )#split.snv=FALSE adds SNVs into no variant count
  save(crispr_set, file = fs::path(output_folder, "crispr_set.rda"))
  return(crispr_set)
}

#' @importFrom CrispRVariants plotVariants
#' @importFrom IRanges IRanges
#' @importFrom grid unit
#' @importFrom fs path
plot_variants <- function(data, conf, output_folder) {
  # group_desig <- conf$group
  group_desig <- rep("group", times = length(data$crispr_runs))
  # plot the variants
  ps <- 37#need to generalize
  pam_seq <- seq(ps, 280, 27)#need to generalize

  while (!is.null(dev.list()))
    dev.off()
  p <- CrispRVariants::plotVariants(
    data,
    col.wdth.ratio = c(1, 1),
    plotAlignments.args = list(
      pam.start = pam_seq,
      #c(37,64), #draws a 3-nt box starting including the position noted
      target.loc = pam_seq - 3,
      #draws a vertical line after the position noted
      guide.loc = IRanges::IRanges(pam_seq -
                                     20, pam_seq + 2),
      #first parameter - beginning of target sequence, second - end of target sequence
      min.count = 200,
      tile.height = 0.9
    ),
    plotFreqHeatmap.args = list(
      min.count = 200,
      plot.text.size = 3,
      x.size = 8,
      group = group_desig,
      legend.text.size = 8,
      legend.key.height = grid::unit(0.5, "lines")
    )
  )
  dev.copy2pdf(
    file = fs::path(output_folder, "crispr_plot.pdf"),
    width = 36,
    height = 36
  )  #for 10 samples use 24 x 24, for 30-40 samples use 48 x 48

  dev.off()


}
