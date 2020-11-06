mjolnir_demulti_filter <- function(lib_prefix,cores,Lmin=299,Lmax=320){
  message("Doing paired-end alignment, demultiplexing and length filter.")
  library(parallel)
  no_cores <- cores*length(lib_prefix)
  clust <- makeCluster(no_cores)
  X <- NULL
  for (i in 1:cores) for (j in 1:length(lib_prefix)) {
    X <- c(X,paste0("illuminapairedend -r ",lib_prefix[j],"_R2_part_",sprintf("%02d",i),".fastq ",lib_prefix[j],"_R1_part_",sprintf("%02d",i),".fastq | obigrep -p \'score>40.00\' | ngsfilter -t ngsfilter_",lib_prefix[j],".tsv | obigrep -p \'seq_length>",Lmin,"\' -p \'seq_length<",Lmax,"\' -s \'^[ACGT]+$\' -p \'forward_tag!=None\' -p \'reverse_tag!=None\' --fasta-output > ",lib_prefix[j],".filtered_length_part_",sprintf("%02d",i),".fasta"))
  }
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)
  message("Filtering done.")
}
