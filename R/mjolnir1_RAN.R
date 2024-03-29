# RAN: Reads Allotment in N portions
# This function will prepare the FASTQ raw data for parallel processing.
# RAN will be run only if your data consist of multiplexed libraries, which will be split into aliquote parts using obisplit, to be processed by FREYJA.

mjolnir1_RAN <- function(R1_filenames="",cores=1,lib_prefixes="",R1_motif="_R1",R2_motif="_R2",obipath=""){
    message(paste0("RAN will split initial FASTQ files in ",cores," fragments each."))
    filelist <- NULL
    outfilelist <- NULL
    for (file in R1_filenames) filelist <- c(filelist,file,gsub(R1_motif,R2_motif,file))
    for (prefix in lib_prefixes) outfilelist <- c(outfilelist,paste0(prefix,"_R1_part"),paste0(prefix,"_R2_part"))
    suppressPackageStartupMessages((library(parallel)))
    no_cores <- length(filelist)
    old_path <- Sys.getenv("PATH")
    clust <- makeCluster(no_cores)
    X <- NULL
    if (grepl(".gz",filelist[1])) {
      for (i in 1:length(filelist)) {
        X <- c(X,paste0("gzip -dc ",filelist[i]," | obidistribute -n ",cores," -p \'",outfilelist[i],"\'"))
      }
    } else  {
        for (i in 1:length(filelist)) {
          X <- c(X,paste0("cat ",filelist[i]," | obidistribute -n ",cores," -p \'",outfilelist[i],"\'"))
        }
    }
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)
    for (file in outfilelist) for (j in 1:9) if (file.exists(paste0(file,"_",j,".fastq"))) system(paste0("mv ",file,"_",j,".fastq ",file,"_",sprintf("%02d",j),".fastq"))
    message("Splitting done.")
} 
