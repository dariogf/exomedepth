#!/usr/bin/env Rscript


library("optparse")
printf <- function(...) cat(sprintf(...))

option_list = list(
  make_option(c("-f", "--samples"), type="character", default=NULL, 
              help="samplelist", metavar="character"),
  make_option(c("-d", "--bed"), type="character", default=NULL, 
              help="bed_file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (sample).\n", call.=FALSE)
}

file_samples=opt$samples;
file_bed=opt$bed;

cat(sprintf("Using samplefile: %s\n",file_samples));
cat(sprintf("Using BED: %s\n",file_bed));

library(GenomicRanges)
library(ExomeDepth)

#############################################################################
########### INPUT ARGUMENTS
############################################################################
#file_samples='exome_depth.samples'
#file_bed='bed_file.bed'

##############################################################################
###########  READING DATA ##############
##############################################################################

########## Load files #################

#Cargamos el nombre de los ficheros que debe leer desde una plantilla
samples <- read.table(file_samples, header= TRUE, sep= '\t')

## load bed #####
#Cargamos el fichero bed 
mybed.hg19 <- read.table(file_bed, header= FALSE, sep= '\t') #fichero bed sin cabecera

## eliminamos 'chr' de la primera columna del fichero bed
mybed.hg19[,1] <- gsub(as.character(mybed.hg19[,1]),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters

mybed.hg19[,1] <- gsub(as.character(mybed.hg19[,1]),
                                   pattern = 'M',
                                   replacement = 'MT')

## unimos bed y bam y crear objeto ExomeDepht #####
nfiles=length(samples$Samples)

# load first elem
n=as.vector(samples$Samples[1])
load(file=paste(n,'.RData',sep=''))
my.counts=aux_bam

for (i in 2:nfiles){
  n=as.vector(samples$Samples[i])
  cat(paste(i,") Loading:",n))
  aux_bam=NULL
  load(file=paste(n,'.RData',sep=''))
  
  my.counts[[n]] <- aux_bam[[n]]
}

aux_bam=NULL
#my.counts <- getBamCounts(bed.frame= mybed.hg19, bam.file= as.vector(samples$Samples))

########## Annotation files #################

#### obtenemos un set de anotaciones para ser utilizados posteriormente desde el fichero bed

my.exons.hg19.GRanges <- GRanges(seqnames = mybed.hg19[,1],
                              IRanges(start=mybed.hg19[,2],end=mybed.hg19[,3]),
                              names = mybed.hg19[,4])


########## Analisis CNV #################

### prepare the main matrix of read count data

ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame') #dataframe con la unión d
#ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = 'B.*')]) #dataframe con la unión d
#ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')]) #dataframe con la unión d

save.image(file='base_data.RData')
quit()
