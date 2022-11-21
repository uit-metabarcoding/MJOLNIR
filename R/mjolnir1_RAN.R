# RAN: Reads Allotment in N portions
# This function will prepare the FASTQ raw data for parallel processing.
# RAN will be run only if your data consist of multiplexed libraries, which will be split into aliquote parts using obisplit, to be processed by FREYJA.

mjolnir1_RAN <- function(R1_filenames="",cores=1,lib_prefixes="",R1_motif="_R1",R2_motif="_R2"){
    message(paste0("RAN will split initial FASTQ files in ",cores," fragments each."))
    filelist <- NULL
    outfilelist <- NULL
    for (file in R1_filenames) filelist <- c(filelist,file,gsub(R1_motif,R2_motif,file)) # here the name of R1 and R2 of all files
    for (prefix in lib_prefixes) outfilelist <- c(outfilelist,paste0(prefix,"_R1_part"),paste0(prefix,"_R2_part")) # ouput file names using lib_prefixes
    suppressPackageStartupMessages((library(parallel)))
    suppressPackageStartupMessages((library(stringr)))
    no_cores <- length(filelist)
    old_path <- Sys.getenv("PATH")
    # unzip files if needed
    X <- NULL

    # check for gz files in filelist and create commands
    thereis_gz <- F
    for (i in 1:length(filelist)) {
      if (grepl(".gz",filelist[i])) {
        thereis_gz <- T
        # create commands to unzip from the system
        X <- c(X,paste0("gzip -dc ",filelist[i]))
        filelist[i] <- str_remove(filelist[i],'.gz')
      }
    }
    # if there is gz run commands
    if (thereis_gz) {
      clust <- makeCluster(no_cores)
      clusterExport(clust, list("X","old_path","obipath"),envir = environment())
      clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
      parLapply(clust,X, function(x) system(x,intern=T,wait=T))
      stopCluster(clust)
    }

    # distribute into <cores> files with new names with prefixes

    # get the number of lines for each file
    X <- NULL
    for (i in 1:length(filelist)) {
          X <- c(X,paste0("wc -l ",filelist[i]," |  cut -f1 -d ' ' "))
        }
    clust <- makeCluster(no_cores)
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
    num_lines <- as.numeric(parLapply(clust,X,
                                      function(x) system(x,intern=T,wait=T)))
    stopCluster(clust)

    # create the command with a vector with lines to create outfiles from filelist
    num_seqs <- num_lines/4
    lines <- list()
    X <- NULL
    for (i in 1:length(filelist)) {
      for (j in 1:cores) {
        lines_vector <- sort(c(seq(1+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4)),
                            seq(2+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4)),
                            seq(3+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4)),
                            seq(4+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4))))
        X <- c(X,
               paste0("perl -e 'while(<>){if(++$l~~[",
                      paste(lines_vector,collapse = ','),
                      "]){print}}' < ", filelist[i]," > ",
                      outfilelist[i],"_",sprintf("%02d",j),".fastq"))
      }
    }

    clust <- makeCluster(no_cores)
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)

    message("Splitting done.")
}
