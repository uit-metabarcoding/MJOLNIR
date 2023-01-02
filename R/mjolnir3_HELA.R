# HELA: Hierarchical Elimination of Lurking Artifacts

# This function uses the uchime_denovo algorithm implemented in VSEARCH to remove chimaeric sequences from the dataset.
# HELA works in a sample-by-sample basis. HELA will process all individual fasta files in the current folder matching the pattern XXXX_sample_XXX.fasta.
# This allows for parallel computing, significantly decreasing calculation times.
# HELA can optionally remove singleton sequences (default: remove_singletons=T), so that the computing time for clustering with ODIN will be considerably reduced.
# This is possibly a good idea for very large datasets (with > 5 million unique sequences before clustering)
# The final dataset output is in VSEARCH format, so it can be directly fed into SWARM (ODIN).

mjolnir3_HELA <- function(lib, cores, obipath = ""){

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))
  sample_list <- gsub("_FREYJA_uniq.fasta","",list.files(pattern="^[a-zA-Z0-9]{4}_[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}_FREYJA_uniq.fasta$"))

  message("HELA will remove chimaeras from each sample")
  X <- NULL
  for (i in sample_list) {
    X <- c(X,paste0("vsearch --uchime_denovo ",i,"_FREYJA_uniq.fasta ",
                    "--sizeout --minh 0.90 ",
                    "--nonchimeras ",i,"_HELA_nonchimeras.fasta ",
                    "--chimeras ",i,"_HELA_chimeras.fasta ",
                    "--uchimeout ",i,"_HELA_uchimeout.log"))
  }
  no_cores <- cores
  clust <- makeCluster(no_cores)
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)

  message("HELA is done.")
}


