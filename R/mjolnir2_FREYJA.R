# FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment  

mjolnir2_FREYJA <- function(lib_prefix,cores,Lmin=299,Lmax=320,obipath=""){
  message("FREYJA will do paired-end alignment, demultiplexing and length filter.")
  suppressPackageStartupMessages(library(parallel))
  no_cores <- cores*length(lib_prefix)
  old_path <- Sys.getenv("PATH")
  clust <- makeCluster(no_cores)
  X <- NULL
  for (i in 1:cores) for (j in 1:length(lib_prefix)) {
    #if (cores<10) {formatted_i <- i} else {formatted_i <- sprintf("%02d",i)}
    formatted_i <- sprintf("%02d",i)
    X <- c(X,paste0("illuminapairedend -r ",lib_prefix[j],"_R2_part_",formatted_i,".fastq ",lib_prefix[j],"_R1_part_",formatted_i,".fastq | obigrep -p \'score>40.00\' | ngsfilter -t ngsfilter_",lib_prefix[j],".tsv | obigrep -p \'seq_length>",Lmin,"\' -p \'seq_length<",Lmax,"\' -s \'^[ACGT]+$\' -p \'forward_tag!=None\' -p \'reverse_tag!=None\' --fasta-output > ",lib_prefix[j],".filtered_length_part_",sprintf("%02d",i),".fasta"))
  }
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)
  message("FREYJA is done.")
}
