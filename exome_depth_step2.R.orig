#!/usr/bin/env Rscript

library("optparse")
printf <- function(...) cat(sprintf(...))

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample name", metavar="character"),
  make_option(c("-r", "--rdata"), type="character", default=NULL, 
              help="base_rdata", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sample) || is.null(opt$rdata)){
  print_help(opt_parser)
  stop("All options must be supplied (sample).\n", call.=FALSE)
}


cat(sprintf("Using SAMPLE: %s\n",opt$sample));

library(GenomicRanges)
library(ExomeDepth)

#load('base_data.RData')
load(opt$rdata)

ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')]) #dataframe con la unión d


nsamples <- ncol(ExomeCount.mat)
### start looping over each sample
#for (i in 1:nsamples) {

#i1=match(paste(opt$sample,"_exome_sorted.bam",sep=''),as.matrix(samples))
#i=match(paste(opt$sample,"_exome_sorted.bam",sep=''),names(as.data.frame(ExomeCount.mat[0,])))
i=match(paste(gsub("-",".",opt$sample),"_exome_sorted.bam",sep=''),names(as.data.frame(ExomeCount.mat[0,])))

#print(samples)
print(names(as.data.frame(ExomeCount.mat[0,])))

printf("Indice  in data_table: %d",i);
#if (i1 != i){
#	printf("ERROR sample %s: some other sample failed execution of step1",opt$sample);
#}
  #### Create the aggregate reference set for this sample
  my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                     reference.counts = ExomeCount.mat[,-i],
#                                     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                      bin.length= ExomeCount.dafr$width)
#                                     n.bins.reduced = 10000)
  my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)

  #dump(my.reference.selected);
  printf("Selected_refs: ");
  my.choice$summary.stats["selected"];

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
  #output.file <- paste(opt$sample,'.csv', sep = '')
  output.file <- 'exome_depth.csv'
  write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)
  
#}

#plot (all.exons,
#      sequence = '2',
#      xlim = c(47600592 - 1000, 47672806 + 1000),
#      count.threshold = 20,
#      main = 'EPCAM gene',
#      cex.lab = 0.8,
#      with.gene = TRUE)
## Plotting the gene data

