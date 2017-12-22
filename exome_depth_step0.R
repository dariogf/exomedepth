#!/usr/bin/env Rscript

library("optparse")
printf <- function(...) cat(sprintf(...))

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample name", metavar="character"),
  make_option(c("-d", "--bed"), type="character", default=NULL, 
              help="bed_file", metavar="character"),
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="bam_file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sample) || is.null(opt$bed) || is.null(opt$bam) ){
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}


cat(sprintf("Using SAMPLE: %s\n",opt$sample));
cat(sprintf("Using BED: %s\n",opt$bed));


library(GenomicRanges)
library(ExomeDepth)


#############################################################################
########### INPUT ARGUMENTS
############################################################################

wd=getwd();

file_samples='exome_depth.samples'
file_bed=opt$bed 
sample_name=opt$sample
sample_bam=opt$bam


##############################################################################
###########  READING DATA ##############
##############################################################################

########## Load files #################

#Cargamos el nombre de los ficheros que debe leer desde una plantilla
#samples <- read.table(file_samples, header= TRUE, sep= '\t')

#cat(as.vector(samples$Samples))

## load bed #####
#Cargamos el fichero bed 
#mybed.hg19 <- read.table('trusight_one_Covered.bed', header= TRUE, sep= '\t') #al fichero bed le hemos puesto cabecera
mybed.hg19 <- read.table(file_bed, header= FALSE, sep= '\t') #fichero bed sin cabecera


## eliminamos 'chr' de la primera columna del fichero bed
mybed.hg19[,1] <- gsub(as.character(mybed.hg19[,1]),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters

mybed.hg19[,1] <- gsub(as.character(mybed.hg19[,1]),
                                   pattern = 'M',
                                   replacement = 'MT')

## unimos bed y bam y crear objeto ExomeDepht #####
aux_bam <- getBamCounts(bed.frame= mybed.hg19, bam.file= sample_bam)

#bam_count[['bam2']] <- bam2[['bam2']]

setwd(wd);
save(file=paste(basename(sample_bam),'.RData',sep=''),list="aux_bam")
quit()
