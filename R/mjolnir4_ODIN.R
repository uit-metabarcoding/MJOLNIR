# ODIN: OTU Delimitation Inferred by Networks

## ODIN performs MOTU clustering and (optionally) entropy denoising and it is one of the main steps of MJOLNIR.
## Clustering is performed using SWARM, which will produce MOTUs as networks of unique sequences
## After SWARM, ODIN will recalculate the abundances for every MOTU in every sample.
## Then (optionally) ODIN proceeds with within-MOTU denoising, using the DnoisE entropy-ratio algorithm for coding regions to get an ESV table.
## Two obligatory arguments are needed: the name of the library, typically 4 characters, and the number of computing cores.
## Three optional parameters: the clustering distance d (default=13),
## min_reads_MOTU is the minimum number of reads to keep a MOTU in the MOTU table (default=2),
## and the minimum number of reads to keep an ESV in the final ESV file (default=2).
## Two boolean parameters can be selected: run_swarm can be set as FALSE to save time if a SWARM output is already available.
## And generate_ESV=TRUE (default) will use the DnoisE algorithm to produce an ESV table, along with the MOTU table.
## ODIN deprecates the previous owi_recount_swarm script used in old metabarcoding pipelines (e.g. Project Metabarpark 2015).
## By Owen S. Wangensteen



mjolnir4_ODIN <- function(lib,cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,alpha=5,COI=T,
                          algorithm="DnoisE_SWARM",obipath="",python_packages="", swarmpath="", remove_singletons = TRUE){
  # setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/adriantich_tests/');
  # # setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/example_MJOLNIR_demultiplexed_data/');
  # cores=3; d=13; min_reads_MOTU=2; min_reads_ESV=2; obipath="/home/adriantich/obi3-env/bin/"; alpha=5; COI=T
  # algorithm="DnoisE_SWARM";
  # lib='ULOY'
  # # algorithm="SWARM_DnoisE";
  # # algorithm="SWARM";
  # # algorithm="DnoisE";
  # python_packages='/home/adriantich/obi3-env/lib/python3.10/site-packages'
  # remove_singletons = TRUE
  # min_reads_MOTU=5

  # the algorithm has 4 options that are shared between them
  # DnoisE: D
  # DnoisE + SWARM: DS
  # SWARM + DnoisE: SD
  # SWARM: S
  # The outputs have to be in all cases a fasta file divided in many parts and a minimum of one csv file containing the abundances in each sample
  # the resulting csv files will be:
  # D csv of ESVs
  # DS and SD csv of ESV within MOTU and only MOTU representative seqs
  # S csv of MOTU representative seqs

  #####
  # 0: define variables
  #####
  algorithm = tolower(algorithm)
  dnoise_path <- "~/DnoisE/src/"
  swarmpath <- "~/swarm/bin/"
  swarm <- paste0(swarmpath,'swarm')
  dnoise <- paste0("python3 ", dnoise_path, "DnoisE.py") # Change this to where the Dnoise executable is
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(python_packages,old_path, obipath, dnoise_path, sep = ":"))

  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))

  if (!(algorithm=="dnoise_swarm" | algorithm=="dnoise" | algorithm=="swarm_dnoise" | algorithm=="swarm")) {
    message('ERROR: algorithm has to be one of the following:\nDnoisE_SWARM\nSWARM_DnoisE\nSWARM\nDnoisE')
    quit()
  }

  #####
  # 1: D and DS -> denoise the fasta files
  #####

  # inputs
  sample_list <- gsub("_HELA_nonchimeras.fasta","",list.files(pattern="^[a-zA-Z0-9]{4}_[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}_HELA_nonchimeras.fasta$"))
  # ouputs
  # file"_ODIN_Adcorr_denoised_ratio_d.fasta | file"_ODIN_denoised_ratio_d.fasta
  if (algorithm=="dnoise_swarm" | algorithm=="dnoise") {
    if (file.exists("summary_HELA.RData")){
      load("summary_HELA.RData")
      before_1_ODIN <- after_HELA
    } else {
      before_1_ODIN <- mclapply(sample_list,function(file){
        output <- system(paste0("grep '>' ",file,"_HELA_nonchimeras.fasta | wc -l"),intern = T,wait = T)
        value <- as.numeric(output)
        return(data.frame(file=paste0(file,"_HELA_nonchimeras.fasta"),
                          num_seqs=value))
      },mc.cores = cores)
    }
    for (file in sample_list) {
      seq_counts <- before_1_ODIN$num_seqs[before_1_ODIN$file==paste0(file,"_HELA_nonchimeras.fasta")]
      if (seq_counts==1) {
          if (COI) {
            system(paste0("cp ",file,"_HELA_nonchimeras.fasta ",file,"_ODIN_Adcorr_denoised_ratio_d.fasta "),intern = T, wait = T)
          } else  {
            system(paste0("cp ",file,"_HELA_nonchimeras.fasta ",file,"_ODIN_denoised_ratio_d.fasta "),intern = T, wait = T)
          }
        } else {
          if (COI) {
            system(paste0(dnoise," --fasta_input ",file,"_HELA_nonchimeras.fasta ",
                          "--fasta_output ",file,"_ODIN ",
                          "-a ",alpha," -c ",cores," -y -m 313 -r ",min_reads_ESV),intern = T, wait = T)
          } else  {
            system(paste0(dnoise," --fasta_input ",file,"_HELA_nonchimeras.fasta ",
                          "--fasta_output ",file,"_ODIN ",
                          "-a ",alpha," -c ",cores," -r ",min_reads_ESV),intern = T, wait = T)
          }
        }
    }
  }

  #####
  # 2: D,DS,S,S -> cat all fasta and perform obi uniq, annotate
  #####

  # input
  # sample_list + _ODIN_Adcorr_denoised_ratio_d.fasta | _ODIN_denoised_ratio_d.fasta | _HELA_nonchimeras.fasta
  # final outputs
  # DMS as lib'_ODIN/'

  X <- NULL
  for (file in sample_list) {
    if (algorithm=="dnoise_swarm" | algorithm=="dnoise") {
      if (COI) {
        input_file <- paste0(file,'_ODIN_Adcorr_denoised_ratio_d.fasta')
      } else {
        input_file <- paste0(file,'_ODIN_denoised_ratio_d.fasta')
      }
    } else {
      input_file <- paste0(file,'_HELA_nonchimeras.fasta')
    }
    X <- c(X,paste0(
      "sed -i 's/;size/; COUNT/g' ",input_file," ; ",
      "obi import --fasta-input ",input_file, " ", file,"_ODIN/sample ; ",
      "obi annotate -S sample:\"", gsub('^[a-zA-Z0-9]{4}_', '', file), "\" ", file,"_ODIN/sample  ", file,"_ODIN/sample_name "))

  }
  clust <- makeCluster(cores)
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)

  for (i in 1:length(sample_list)) {
    file <- sample_list[i]
    if (i==1) {
      system(paste0("obi cat -c ",
                    file,'_ODIN/sample_name',
                    ' ', lib,"_ODIN/version",i),
             intern = T, wait = T)
    } else if (i==length(sample_list)) {
      system(paste0("obi cat -c ",
                    file,'_ODIN/sample_name',
                    ' -c ', lib,"_ODIN/version",c(i-1),
                    ' ', lib,"_ODIN/samples ; ",
                    "obi rm ", lib,"_ODIN/version",c(i-1)),
             intern = T, wait = T)
    } else {
      system(paste0("obi cat -c ",
                    file,'_ODIN/sample_name',
                    ' -c ', lib,"_ODIN/version",c(i-1),
                    ' ', lib,"_ODIN/version",i," ; ",
                    "obi rm ", lib,"_ODIN/version",c(i-1)),
             intern = T, wait = T)
    }
  }
  system(paste0("obi uniq --merge 'sample' ",
                lib,'_ODIN/samples ',
                lib,"_ODIN/samples_uniq"),
         intern = T, wait = T)
  if (remove_singletons) {
    system(paste0("obi grep -p \"sequence[\'COUNT\'] > 1\" ",
                  lib,"_ODIN/samples_uniq ",
                  lib,"_ODIN/seq_nosing"),
           intern = T, wait = T)
    system(paste0("obi annotate --seq-rank ",
                  lib,"_ODIN/seq_nosing ",
                  lib,"_ODIN/seq_rank"),
           intern = T, wait = T)
  } else {
    system(paste0("obi annotate --seq-rank ",
                  lib,"_ODIN/samples_uniq ",
                  lib,"_ODIN/seq_rank"),
           intern = T, wait = T)

  }
  system(paste0("obi annotate --set-identifier ",
                "\'\"\'",lib,"\'_%09d\" % sequence[\"seq_rank\"]\' ",
                lib,"_ODIN/seq_rank ",
                lib,"_ODIN/seq_id"),
         intern = T, wait = T)

  output <- system(paste0("obi ls ",lib,"_ODIN | grep 'Line count'"),intern = T,wait = T)
  values <- as.numeric(gsub(".*count: ","",output))
  version <- gsub(".*# ","",gsub(": Date.*","",output))
  after_2_ODIN <-data.frame(algorithm=algorithm,
                                         version=version,
                                         num_seqs=values)



  #####
  # 3: D,DS,S,S -> export csv file of sequences (if denoised, they are ESV)
  #####

  # input
  # sample_list + _ODIN_Adcorr_denoised_ratio_d.fasta | _ODIN_denoised_ratio_d.fasta | _HELA_nonchimeras.fasta
  # final outputs
  if (algorithm == 'dnoise'){
    filetab <- paste0(lib, "_ODIN_ESV.csv")
  } else{
    filetab <- paste0(lib, "_ODIN_seqs.csv")
  }

  system(paste0("obi export --tab-output --sep ','  -o ",
                filetab, " ",
                lib, "_ODIN/seq_id"),
         intern = T, wait = T)


  #####
  # 4a: DS,SD,S -> do SWARM and create csv files of MOTUs and ESV or seqs. Also create fasta files from MOTU csv
  # 4ab: SD -> run DnoisE over the csv files
  # 4b: D -> divide fasta file into different parts
  #####

  # input
  # filetab (4a) | lib"_ODIN/seq_id" (4b)
  # final outputs
  outfile <- paste0(lib,"_ODIN")# "_part_",sprintf("%02d",part),'.fasta'
  if (algorithm=="swarm_dnoise" | algorithm=="swarm" | algorithm=="dnoise_swarm"){
    outfile_MOTU <- paste0(outfile,"_counts.tsv")
    outfile_ESV <- paste0(outfile,"_ESV.tsv")
    outfile_DnoisE <- paste0(outfile,"_SWARM_DnoisE") # + _Adcorr_denoised_ratio_d.csv | _Adcorr_denoised_ratio_d.csv
    # intermediate input/ouput
    fileswarm <- paste0(outfile,"_SWARM_output")
  }



  if (algorithm=="swarm_dnoise" | algorithm=="swarm" | algorithm=="dnoise_swarm"){ # 4a
    message("ODIN will cluster sequences into MOTUs with SWARM.")
    system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",lib,"_ODIN/seq_id > ",outfile,".fasta ; ",
                  "sed -i 's/COUNT/size/g' ",outfile,".fasta ; ",
                  "sed -i 's/;//g' ",outfile,".fasta ; ",
                  "sed -E -i 's/(size=[0-9]*).*/\\1;/g' ",outfile,".fasta ; ",
                  "sed -i 's/ /;/g' ",outfile,".fasta "))
    system(paste0(swarm," -d ",d," -z -t ",cores," -o ",outfile,"_SWARM_output -s ",outfile,"_SWARM",d,"nc_stats -w ",outfile,"_SWARM_seeds.fasta ",outfile,".fasta"),intern=T,wait=T)

    get_swarm_size <- function(cadena="="){
      # it gets the positon of the '=' and gets the items after it
      return(as.numeric(gsub(";","",substr(cadena,gregexpr("=",cadena)[[1]][[1]]+1,nchar(cadena)))))
    }

    message("ODIN will recount abundances for every MOTU after Swarm.")

    # Read cluster list database
    message("1. ODIN is reading SWARM results...")
    swarm_db <- readLines(fileswarm)
    total_swarms <- length(swarm_db)
    message("2. ODIN has read ", total_swarms," total MOTUs.")

    message(paste("3. ODIN will now perform a remove of MOTUs with sequences of less than",
                  min_reads_MOTU,"reads."))
    # Calculate reads in each cluster and reduce the dataset if min_reads_MOTU>0
    clusters <- strsplit(swarm_db,"; ")
    if (min_reads_MOTU>0) {
      # first reduction of a proportion of MOTUs
      if (min_reads_MOTU>9) i <- 9 else i <- min_reads_MOTU-1
      clusters <- clusters[!(grepl(paste0('size=[0-',i,']'),clusters) & lengths(clusters)==1)]
      cluster_reads <- NULL
      # second reduction
      cluster_reads <- mclapply(clusters,function(x) sum(as.numeric(lapply(X=(x),FUN=get_swarm_size))), mc.cores = cores)
      # for (i in 1:length(clusters)) cluster_reads[i] <- sum(as.numeric(lapply(X=(clusters[[i]]),FUN=get_swarm_size)))
      clusters <- clusters[cluster_reads>=min_reads_MOTU]
      total_swarms <- length(clusters)
    }

    message("4. ODIN will now calculate the number of reads in every sample for each MOTU.")
    clusters <- mclapply(clusters,function(x){sub(";.*","",x)}, mc.cores = cores)
    names(clusters) <- mclapply(clusters, function(x) x[[1]], mc.cores = cores)

    message("5. ODIN kept only ", total_swarms," MOTUs of size greater than or equal to ",min_reads_MOTU," reads.")
    motu_seqs_names <- stack(clusters) %>% rename(ID = values, MOTU = ind)
    # motu_seqs_names <- unlist(clusters, use.names=F)

    # # Generate a file with the list of ids of non-singleton clusters
    # motulist <- file(paste0(lib,"_non_singleton_motu_list.txt"),"wt")
    # writeLines(id,motulist)
    # message("ODIN has created the file ",paste0(lib,"_non_singleton_motu_list.txt")," with the list of identifiers of non-singleton MOTUs.")

    # Read counts database and keep only the needed clusters
    message("6. ODIN is reading the abundance database. This could take Him a while, since He has just one eye left, after all.")
    db <- read.table(filetab,sep=",",head=T)
    numseqs <- nrow(db)
    db <- merge(motu_seqs_names2,db,by="ID")
    # db <- db[db$ID %in% motu_seqs_names,]
    numseqs_reduced <- nrow(db)
    samples <- sum(grepl('sample',names(db)))
    message("7. ODIN finished reading the Database, which includes ", numseqs," total unique sequences and ",samples," samples.\n",
            "ODIN kept only ", numseqs_reduced," sequences for calculations.")

    message("4. ODIN will now calculate the number of reads in every sample for each MOTU.")
    db.total <- split(db[,grepl('sample',names(db))], db$MOTU)
    db.total <- mclapply(db.total,function(x)as.data.frame(t(as.matrix(c(COUNT=sum(x),colSums(x),CLUST_WEIGHT=dim(x)[1])))), mc.cores = cores)
    db.total <- do.call(rbind,db.total)
    db.total <- cbind(data.frame(ID=rownames(db.total)), db.total)
    db.total <- merge(db.total,db[,grepl('ID|NUC_SEQ',names(db))],by = "ID")

    # order the columns
    col_order <- c('ID', 'COUNT', 'MOTU', names(db)[grepl('sample',names(db))], 'NUC_SEQ' )
    db <- db[,col_order]
    s_opt <- min(grep('sample',names(db)))
    z_opt <- max(grep('sample',names(db)))

    # print datasets
    write.table(db.total,outfile_MOTU,sep="\t",quote=F,row.names=F)
    write.table(db,outfile_ESV,sep="\t",quote=F,row.names=F)

    # divide dataset into different files for THOR
    db.total <- db.total[,c('ID','NUC_SEQ')]
    db.total <- paste(paste0('>',db.total$ID),db.total$NUC_SEQ,sep='\n')
    db.total <- split(db.total, factor(sort(1:length(db.total)%%cores)))
    for (part in 1:length(db.total)) {
      writeLines(paste0(db.total[[part]],collapse = '\n'),paste0(outfile,"_part_",sprintf("%02d",part),'.fasta'))
    }

    message("File ", outfile, " written")

    if (algorithm=='swarm_dnoise') { # 4ab
      message("ODIN will generate now a list of ESVs for every non-singleton MOTU, using DnoisE.")
      if (COI) {
        system(paste0(dnoise," --csv_input ",outfile_ESV," ",
                      "--csv_output ",outfile_DnoisE," ",
                      "-a ",alpha," -c ",cores," -n 'COUNT' -p 1 -q 'NUC_SEQ' ",
                      "-s ", s_opt," -z ", z_opt," ",
                      "-y -w 'MOTU' -r ",min_reads_ESV),intern = T, wait = T)
      } else  {
        system(paste0(dnoise," --csv_input ",outfile_ESV," ",
                      "--csv_output ",outfile_DnoisE," ",
                      "-a ",alpha," -c ",cores," -n 'COUNT' -p 1 -q 'NUC_SEQ' ",
                      "-s ", s_opt," -z ", z_opt," ",
                      "-w 'MOTU' -r ",min_reads_ESV),intern = T, wait = T)
      }
      message("")
    }

  } else { # 4b

    if (cores == 1) {
      system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",lib,"_ODIN/seq_id > ",paste0(outfile,"_part_01.fasta")), intern = T, wait = T)
    } else {
      system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",lib,"_ODIN/seq_id > ",outfile,".fasta "), intern = T, wait = T)
      fasta_file <- readLines(paste0(outfile,'.fasta'))

      seqs <- grep('>',file_read)
      seq_rank <- sort(seqs%%cores)

      for (i in 1:cores) {
        if (i == 1) {
          start_line <- 1
          end_line <- seqs[grep(unique(seq_rank)[i+1],seq_rank)[1]]-1
        } else if (i == cores) {
          start_line <- seqs[grep(unique(seq_rank)[i],seq_rank)[1]]
          end_line <- length(fasta_file)
        } else {
          start_line <- seqs[grep(unique(seq_rank)[i],seq_rank)[1]]
          end_line <- seqs[grep(unique(seq_rank)[i+1],seq_rank)[1]]-1
        }
        writeLines(paste0(fasta_file[start_line:end_line],
                          collapse = '\n'),
                   paste0(outfile,"_part_",sprintf("%02d",i),'.fasta'))
      }
    }
  }
  message("ODIN is done.")
}

