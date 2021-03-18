#!/usr/bin/env Rscript

# NOTE: Script designed to run on neon.hpc.uiowa.edu as a qsub job with 16
#       cores.
library("argparser", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)
library("Rsubread", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE) # For featureCounts
library("DESeq2", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE) # For differential expression analysis
library("dplyr", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)
library("readr", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)

# Expects a data.frame with "old" and "new" sample IDs
clean_sample_names <- function(x, id_map, old="old", new="new") {
  id_map <- id_map %>% rename_("old"=old, "new"=new)
  
  # Update names in x$counts matrix
  counts_old_names <- x$counts %>% colnames()
  
  colnames(x$counts) <- id_map[match(counts_old_names, id_map$old),"new"] %>%
    collect %>% .[["new"]]
  
  # Update names in x$target vector
  x$targets <- id_map[match(x$targets, id_map$old),"new"] %>%
    collect %>% .[["new"]]
  
  # Update names in x$stat data.frame
  stat_old_names <- names(x$stat)[which(names(x$stat) != "Status")]
  
  stat_new_names <- id_map[match(stat_old_names, id_map$old), "new"] %>% 
    collect %>% .[["new"]]
  
  names(x$stat) <- c("Status", stat_new_names)
  
  return(x)
}
# Set up argparser ------------------------

p <- arg_parser("Get counts of reads within genes, exons, or other features using featureCounts")

p <- add_argument(p, "--gtf", 
                  help="GTF file of features", 
                  default = "gencode.vM6.annotation_spikes.gtf")

p <- add_argument(p, "--outdir", 
                  help="number of decimal places", 
                  default = "data/")

p <- add_argument(p, "--temp",
                  help = "path for temporary directory",
                  default = "testing")

p <- add_argument(p, "--threads",
                  help = "number of threads",
                  default = 16)

p <- add_argument(p, "--input",
                  help = "tab-delimited file of bams and sample names",
                  default = "samples_bams_trimmed.txt")

p <- add_argument(p, "--prefix",
                  help = "prefix for sample output",
                  default = "summarized_experiment")

p <- add_argument(p, "--feature",
                  help = "type of feature to count reads",
                  default = "exon")

argv <- parse_args(p)

# Parameters ------------------------
gtf_path <- argv$gtf
output_path <- argv$outdir
temporary_path <- argv$temp
number_threads <- argv$threads
sample_bams_path_file <- argv$input
out_name <- paste(argv$prefix, ".RData", sep="")


# Load data ------------------------
# Load all the sorted bam files
bam_file_df <- read_tsv(sample_bams_path_file) %>%
    mutate(old = make.names(bam_file))
    
bam_file_list <- bam_file_df %>% collect %>% .[["bam_file"]]

working_directory <- getwd()

# Move to temporary directory ------------------------
setwd(temporary_path)

# Note: For strand-specific Illumina TruSeq libraries, set 'strandSpecific=2'.
#       For details, see: https://support.bioconductor.org/p/66733/
summarized_experiment <- featureCounts(files=bam_file_list, 
										 annot.ext=gtf_path, 
										 GTF.featureType="exon", 
										 isPairedEnd=TRUE, 
										 nthreads=number_threads, 
										 isGTFAnnotationFile=TRUE,
										 strandSpecific=0,
										 requireBothEndsMapped=FALSE)

setwd(working_directory) # Return to original directory

# Create the output directory if it does not exist
if(!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

# Clean up the sample IDs
summarized_experiment <- clean_sample_names(summarized_experiment, bam_file_df, 
                                            old="old", new="sample_id")


# Save the read counts ------------------------
save(summarized_experiment, 
     file=file.path(output_path, out_name))

# TODO: Add differntial expression code here.
# DESeqDataSetFromMatrix
