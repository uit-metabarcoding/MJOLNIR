# HELA: Hierarchical Elimination of Lurking Artifacts

# This function uses the uchime_denovo algorithm implemented in VSEARCH to remove chimaeric sequences from the dataset.
# HELA works in a sample-by-sample basis. HELA will process all individual fasta files in the current folder matching the pattern XXXX_sample_XXX.fasta.
# This allows for parallel computing, significantly decreasing calculation times.  
# HELA can optionally remove singleton sequences (default: remove_singletons=T), so that the computing time for clustering with ODIN will be considerably reduced. 
# This is possibly a good idea for very large datasets (with > 5 million unique sequences before clustering)
# The final dataset output is in VSEARCH format, so it can be directly fed into SWARM (ODIN).

mjolnir3_HELA <- function(lib, cores, remove_singletons = TRUE, obipath = ""){
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))
  fastq_list <- list.files(pattern="^[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}.fastq$")
  sample_list <- gsub(".fast[aq]","",list.files(pattern="^[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}.fast[aq]$"))
  message("HELA will group unique sequences in every sample")
  suppressPackageStartupMessages(library(parallel))
  no_cores <- cores
  clust <- makeCluster(no_cores)
  X <- NULL
  if (identical(fastq_list,character(0))) {for (i in sample_list) X <- c(X,paste0("obiuniq ",i," > ",i,"_unique.fasta"))} else{
      for (i in 1:length(fastq_list)) X <- c(X,paste0("obiconvert --fasta-output ",fastq_list[i]," | obiuniq > ",sample_list[i],"_unique.fasta"))
  }
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  message("HELA will convert every sample file to vsearch format")
  X <- NULL
  for (i in sample_list) X <- c(X,paste0(i,"_unique.fasta"))
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
  parLapply(clust,X, function(x) owi_obisample2vsearch(x))
  message("HELA will remove chimaeras from each sample")
  X <- NULL
  for (i in sample_list) X <- c(X,paste0("vsearch --uchime_denovo ",i,"_unique_vsearch.fasta --sizeout --minh 0.90 --nonchimeras ",i,"_nonchimeras.fasta --chimeras ",i,"_chimeras.fasta --uchimeout ",i,"_uchimeout.log"))
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)
  message("HELA will put all non-chimaeric sequences together")
  system(paste0("cat *_nonchimeras.fasta > ",lib,"_no_chimeras.fasta"),intern=T,wait=T)
  message("File ",lib,"_no_chimeras.fasta written.")
  system(paste0("sed -i 's/size/ count/g' ",lib,"_no_chimeras.fasta"),intern=T,wait=T)
  system(paste0("sed -i 's/;sample=/ sample=/g' ",lib,"_no_chimeras.fasta"),intern=T,wait=T)
  system(paste0("sed -i 's/;;/;/g' ",lib,"_no_chimeras.fasta"),intern=T,wait=T)
            
  
  # Check if the output of vsearch is in the right format. Otherwise, add ";" at the end of the headers in no_chimeras.fasta file
  check_header <- readLines(paste0(lib,"_no_chimeras.fasta"),1)
  if (substr(check_header,nchar(check_header),nchar(check_header))!=";") {   
    suppressPackageStartupMessages(library(Biostrings))
    file_nochim <- readDNAStringSet(paste0(lib,"_no_chimeras.fasta"))
    file_nochim@ranges@NAMES <- paste0(file_nochim@ranges@NAMES,";")
    writeXStringSet(file_nochim,paste0(lib,"_no_chimeras.fasta"))       
  }
            
  if (!remove_singletons) {system(paste0("obiuniq -m sample ",lib,"_no_chimeras.fasta > ",lib,"_unique.fasta"),intern=T,wait=T)} else {
    system(paste0("obiuniq -m sample ",lib,"_no_chimeras.fasta | obigrep -p 'count>1' > ",lib,"_unique.fasta"),intern=T,wait=T)}

  message("HELA will change sequence identifiers to a short index")
  system(paste0("obiannotate --seq-rank ",lib,"_unique.fasta | obiannotate --set-identifier \'\"\'",lib,"\'_%09d\" % seq_rank\' > ",lib,"_new.fasta"),intern=T,wait=T)
  message("HELA will change the format to vsearch, so ODIN can use it for SWARM.")
  system(paste0("obiannotate --delete-tag=merged_sample --delete-tag=seq_rank ",lib,"_new.fasta | sed 's/ count/;size/g' > ",lib,"_vsearch.fasta"),intern=T,wait=T)
  message("File ",lib,"_vsearch.fasta written.")
  message("HELA is obtaining a table file with abundances of unique sequences in each sample")
  system(paste0("obitab -o ",lib,"_new.fasta >  ",lib,"_new.tab"),intern=T,wait=T)
  message("HELA is done.")
}


