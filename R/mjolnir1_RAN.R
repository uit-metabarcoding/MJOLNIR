# RAN: Reads Allotment in N portions
# This function will prepare the FASTQ raw data for parallel processing.
# RAN will be run only if your data consist of multiplexed libraries, which will be split into aliquote parts using obisplit, to be processed by FREYJA.

mjolnir1_RAN <- function(R1_filenames="",cores=1,lib_prefixes="",R1_motif="_R1",R2_motif="_R2"){
  setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/adriantich_tests/'); R1_filenames="../example_MJOLNIR_multiplexed_data_metafast/ULO1_R1.fastq.gz"; cores=3;lib_prefixes="ULOY";R1_motif="_R1";R2_motif="_R2"
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
        if (file.exists(filelist[i])) {
          thereis_gz <- T
          # create commands to unzip from the system
          X <- c(X,paste0("gzip -dc ",filelist[i], ' > ',str_remove(filelist[i],'.gz')))
          filelist[i] <- str_remove(filelist[i],'.gz')
        } else {
          print(paste0(filelist[i],' not found or already unziped'))
          filelist[i] <- str_remove(filelist[i],'.gz')
        }

      }
    }
    # if there is gz run commands
    if (thereis_gz) {
      print('unzipping')
      clust <- makeCluster(no_cores)
      clusterExport(clust, list("X","old_path"),envir = environment())
      clusterEvalQ(clust, {Sys.setenv(PATH = old_path)})
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
    clusterExport(clust, list("X","old_path"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = old_path)})
    num_lines <- as.numeric(parLapply(clust,X,
                                      function(x) system(x,intern=T,wait=T)))
    stopCluster(clust)

    # create the command with a vector with lines to create outfiles from filelist
    num_seqs <- num_lines/4
    for (i in 1:length(filelist)) {
      if (i==1) {
        commands <- NULL
        file_commands <- NULL
      }
      for (j in 1:cores) { # cores will be the number of output files per infile
        # get the lines on infile for each outfile
        lines_vector <- sort(c(seq(1+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4)),
                            seq(2+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4)),
                            seq(3+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4)),
                            seq(4+(j-1)*4,num_seqs[i]*4,((cores-1)*4+4))))
        # divide all in sets of 100 lines -> 25 seqs (system command is too long otherwise)
        sqs <- 4*100
        n_sets <- (floor(length(lines_vector)/sqs))
        if (n_sets==0) { # if lines_vector length is lower than 100
          n_sets <- 1
          last_set <- length(lines_vector)
        } else if (n_sets*sqs == length(lines_vector)) { # if lines_vector length is a multiple of 100
          last_set <- sqs
        } else { # if lines_vector length is NOT a multiple of 100
          last_set <- length(lines_vector)-n_sets*sqs
          n_sets <- n_sets+1
        }
        for (sets in 1:4){ #n_sets) {
          if (sets == n_sets) {
            set_endlines <- last_set + (sets-1)*sqs
          } else {
            set_endlines <- sqs + (sets-1)*sqs
          }
          set_beginlines <- 1 + (sets-1)*sqs
          if (sets == 1) {
            X <- c(paste0("perl -e 'while(<>){if(++$l~~[",
                          paste(lines_vector[1:set_endlines],
                                collapse = ','),
                          "]){print}}' < ", filelist[i]," > ",
                          outfilelist[i],"_",sprintf("%02d",j),".fastq"))
          } else {
            X <- c(X,
                   paste0("perl -e 'while(<>){if(++$l~~[",
                          paste(lines_vector[set_beginlines:set_endlines],
                                collapse = ','),
                          "]){print}}' < ", filelist[i]," >> ",
                          outfilelist[i],"_",sprintf("%02d",j),".fastq"))
          }
        }
        commands <- append(commands,list(list(X)))
      }
      file_commands <- append(file_commands,commands)
    }


    clust <- makeCluster(no_cores)
    clusterExport(clust, list("file_commands","old_path"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = old_path)})
    parLapply(clust,file_commands, function(x) lapply(x, function(y) system(y,intern=T,wait=T)))
    stopCluster(clust)

    message("Splitting done.")
}
