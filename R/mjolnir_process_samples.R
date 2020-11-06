mjolnir_process_samples <- function(libs,lib,cores){
  libslist <- NULL
  for (prefix in libs) libslist <- paste0(libslist,prefix,".filtered_length_part*.fasta ")
  message("Joining filtered reads into a single file.")
  system(paste0("cat ",libslist," > ",lib,".joined.fasta"),intern=T,wait=T)
  message("Calculating stats per sample.")
  system(paste0("obistat -c sample -a seq_length ",lib,".joined.fasta > sample_stats_",lib,".txt"))
  message("Splitting into individual files for each sample.")
  system(paste0("obisplit -t sample ",lib,".joined.fasta"),intern=T,wait=T)
  sample_db <- read.csv(paste0("sample_stats_",lib,".txt"),sep="\t",head=T)
  sample_list <- sample_db$sample[order(sample_db$sample)]
  message("Grouping unique sequences in every sample")
  library(parallel)
  no_cores <- cores
  clust <- makeCluster(no_cores)
  X <- NULL
  for (i in sample_list) X <- c(X,paste0("obiuniq ",i,".fasta > ",i,".unique.fasta"))
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  message("Converting to vsearch format")
  X <- NULL
  for (i in sample_list) X <- c(X,paste0("i,".unique.fasta"))
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) owi_obisample2vsearch(x))
  message("Removing Chimaeras in every sample")
  X <- NULL
  for (i in sample_list) X <- c(X,paste0("vsearch --uchime_denovo ",i,".unique.vsearch.fasta --sizeout --minh 0.90 --nonchimeras ",i,".nonchimeras.fasta --chimeras ",i,".chimeras.fasta --uchimeout ",i,".uchimeout.log"))
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)
  message("Concatenating non chimaeric sequences")
  system(paste0("cat *.nonchimeras.fasta > ",lib,".no_chimeras.fasta"),intern=T,wait=T)
  message("File ",lib,".no_chimeras.fasta written.")
  system(paste0("sed -i 's/size/ count/g' ",lib,".no_chimeras.fasta"),intern=T,wait=T)
  system(paste0("sed -i 's/;sample=/ sample=/g' ",lib,".no_chimeras.fasta"),intern=T,wait=T)
  message("Finding unique sequences in ",lib,".no_chimeras.fasta file")
  system(paste0("obiuniq -m sample ",lib,".no_chimeras.fasta > ",lib,".unique.fasta"),intern=T,wait=T)
  message("Changing sequence identifiers to a short index")
  system(paste0("obiannotate --seq-rank ",lib,".unique.fasta | obiannotate --set-identifier \'\"\'",lib,"\'_%09d\" % seq_rank\' > ",lib,".new.fasta"),intern=T,wait=T)
  message("Changing format to vsearch, so we can use Swarm.")
  owi_obifasta2vsearch(infile=paste0(lib,".new.fasta"),outfile=paste0(lib,".vsearch.fasta"))
  message("File ",lib,".vsearch.fasta written.")
  message("Obtaining the table file with abundances")
  system(paste0("obitab -o ",lib,".new.fasta >  ",lib,".new.tab"),intern=T,wait=T)
  message("Done.")
}


