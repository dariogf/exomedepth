#!/usr/bin/env Rscript

library(GenomicRanges)
library(ExomeDepth)

#############################################################################
########### INPUT ARGUMENTS
############################################################################
#args<-commandArgs(TRUE) #Monta un array donde introduce los argumentos de entrada a partir del indice 1
#args = commandArgs(TRUE)

#file_samples = toString(args[1])
#file_bed = toString(args[2])
#cat(file_samples,"\n")
#cat(file_bed)
file_samples='exome_depth.samples'
file_bed='/mnt/home/users/seq_001_genologica/genologica/workflows/templates/trusight_one_Covered.bed'

#cat(typeof(file_samples))
##############################################################################
###########  WORK DIRECTORY ##############
##############################################################################

#t1 <- try(system("pwd", intern = TRUE)) # Se consulta donde se esta ejecutando el script y asi el directorio de trabajo pasa a ser el sitio en el cual se ha invocadp al script
#dirname <- setwd(t1) # Ir al directorio de trabajo

#dirname <- setwd("/Users/rociobm/Desktop/ExomeDetph/bam_file/")







##############################################################################
###########  Load library ##############
##############################################################################
#library(ExomeDepth)



##############################################################################
###########  READING DATA ##############
##############################################################################

########## Load files #################

#Cargamos el nombre de los ficheros que debe leer desde una plantilla
samples <- read.table(file_samples, header= TRUE, sep= '\t')

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

#cat(mybed.hg19[,1])

## unimos bed y bam y crear objeto ExomeDepht #####
my.counts <- getBamCounts(bed.frame= mybed.hg19, bam.file= as.vector(samples$Samples))


########## Annotation files #################

#### obtenemos un set de anotaciones para ser utilizados posteriormente desde el fichero bed

my.exons.hg19.GRanges <- GRanges(seqnames = mybed.hg19[,1],
                              IRanges(start=mybed.hg19[,2],end=mybed.hg19[,3]),
                              names = mybed.hg19[,4])


########## Analisis CNV #################

### prepare the main matrix of read count data

ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame') #dataframe con la unión d
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = 'B.*')]) #dataframe con la unión d
nsamples <- ncol(ExomeCount.mat)
### start looping over each sample
for (i in 1:nsamples) {
  #### Create the aggregate reference set for this sample
  my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                     reference.counts = ExomeCount.mat[,-i],
#                                     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                      bin.length= ExomeCount.dafr$width)
#                                     n.bins.reduced = 10000)
  my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)
  message('Now creating the ExomeDepth object')
  all.exons <- new('ExomeDepth',
                    test = ExomeCount.mat[,i],
                    reference = my.reference.selected,
                    formula = 'cbind(test, reference) ~ 1')
  
  ################ Now call the CNVs
  all.exons <<- CallCNVs(x = all.exons,
                        transition.probability = 10^-4, 
                        chromosome = ExomeCount.dafr$space,
                        start = ExomeCount.dafr$start,
                        end = ExomeCount.dafr$end,
                        name = ExomeCount.dafr$names)
  # transition.probability= Transition probability of the hidden Markov Chain from the normal copy number state to either a deletion or a duplication
  
  ########################### Now annotate the ExomeDepth object

  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = my.exons.hg19.GRanges,
                             min.overlap = 0.0001,
                             column.name = 'exons.hg19')
  output.file <- paste('B_', i, '.csv', sep = '')
  11
  
  write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)
}

#plot (all.exons,
#      sequence = '2',
#      xlim = c(47600592 - 1000, 47672806 + 1000),
#      count.threshold = 20,
#      main = 'EPCAM gene',
#      cex.lab = 0.8,
#      with.gene = TRUE)
## Plotting the gene data

