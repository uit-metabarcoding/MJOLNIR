# HELA: Hierarchical Elimination of Lurking Artifacts

# This function uses the uchime_denovo algorithm implemented in VSEARCH to remove chimaeric sequences from the dataset.
# HELA works in a sample-by-sample basis. HELA will process all individual fasta files in the current folder matching the pattern XXXX_sample_XXX.fasta.
# This allows for parallel computing, significantly decreasing calculation times.  
# The final dataset output is in VSEARCH format, so it can be directly fed into SWARM (ODIN).

mjolnir3_HELA <- function(lib,cores,obipath=""){
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))
  sample_list <- gsub(".fasta","",list.files(pattern="^[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}.fasta$"))
  message("HELA will group unique sequences in every sample")
  suppressPackageStartupMessages(library(parallel))
  no_cores <- cores
  clust <- makeCluster(no_cores)
  X <- NULL
  for (i in sample_list) X <- c(X,paste0("obiuniq ",i,".fasta > ",i,".unique.fasta"))
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  message("HELA will convert every sample file to vsearch format")
  X <- NULL
  for (i in sample_list) X <- c(X,paste0(i,".unique.fasta"))
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
  parLapply(clust,X, function(x) owi_obisample2vsearch(x))
  message("HELA will remove chimaeras from each sample")
  X <- NULL
  for (i in sample_list) X <- c(X,paste0("vsearch --uchime_denovo ",i,".unique.vsearch.fasta --sizeout --minh 0.90 --nonchimeras ",i,".nonchimeras.fasta --chimeras ",i,".chimeras.fasta --uchimeout ",i,".uchimeout.log"))
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)
  message("HELA will put all non-chimaeric sequences together")
  system(paste0("cat *.nonchimeras.fasta > ",lib,".no_chimeras.fasta"),intern=T,wait=T)
  message("File ",lib,".no_chimeras.fasta written.")
  system(paste0("sed -i 's/size/ count/g' ",lib,".no_chimeras.fasta"),intern=T,wait=T)
  system(paste0("sed -i 's/;sample=/ sample=/g' ",lib,".no_chimeras.fasta"),intern=T,wait=T)
  
  # Check if the output of vsearch is in the right format. Otherwise, add ";" at the end of the headers in no_chimeras.fasta file
  check_header <- readLines(paste0(lib,".no_chimeras.fasta"),1)
  if (substr(check_header,nchar(check_header),nchar(check_header))!=";") {   
    suppressPackageStartupMessages(library(Biostrings))
    file_nochim <- readDNAStringSet(paste0(lib,".no_chimeras.fasta"))
    file_nochim@ranges@NAMES <- paste0(file_nochim@ranges@NAMES,";")
    writeXStringSet(file_nochim,paste0(lib,".no_chimeras.fasta"))       
  }
            
  system(paste0("obiuniq -m sample ",lib,".no_chimeras.fasta > ",lib,".unique.fasta"),intern=T,wait=T)
  message("HELA will change sequence identifiers to a short index")
  system(paste0("obiannotate --seq-rank ",lib,".unique.fasta | obiannotate --set-identifier \'\"\'",lib,"\'_%09d\" % seq_rank\' > ",lib,".new.fasta"),intern=T,wait=T)
  message("HELA will change the format to vsearch, so ODIN can use it for SWARM.")
  owi_obifasta2vsearch(infile=paste0(lib,".new.fasta"),outfile=paste0(lib,".vsearch.fasta"))
  message("File ",lib,".vsearch.fasta written.")
  message("HELA is obtaining a table file with abundances of unique sequence in each sample")
  system(paste0("obitab -o ",lib,".new.fasta >  ",lib,".new.tab"),intern=T,wait=T)
  message("HELA is done.")
}


