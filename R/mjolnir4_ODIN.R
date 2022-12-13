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



mjolnir4_ODIN <- function(lib,cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,COI=T,
                          algorithm="DnoisE_SWARM",obipath="", swarmpath="", remove_singletons = TRUE){
  setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/adriantich_tests/');
  # setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/example_MJOLNIR_demultiplexed_data/');
  cores=3; d=13; min_reads_MOTU=2; min_reads_ESV=2; obipath="/home/adriantich/obi3-env/bin/"; alpha=5; COI=T
  algorithm="DnoisE_SWARM";
  lib='ULOY'
  # algorithm="SWARM_DnoisE";
  # algorithm="SWARM";
  # algorithm="DnoisE";
  python_packages='/home/adriantich/obi3-env/lib/python3.10/site-packages'
  remove_singletons = TRUE
  min_reads_MOTU=5

  algorithm = tolower(algorithm)
  dnoise_path <- "~/DnoisE/src/"
  swarmpath <- "~/swarm/bin/"
  swarm <- paste0(swarmpath,'swarm')
  dnoise <- paste0("python3 ", dnoise_path, "DnoisE.py") # Change this to where the Dnoise executable is
  old_path <- Sys.getenv("PATH")
  filetab=paste0(lib, "_ODIN_seqs.csv")
  Sys.setenv(PATH = paste(python_packages,old_path, obipath, dnoise_path, sep = ":"))

  if (!(algorithm=="dnoise_swarm" | algorithm=="dnoise" | algorithm=="swarm_dnoise" | algorithm=="swarm")) {
    message('ERROR: algorithm has to be one of the following:\nDnoisE_SWARM\nSWARM_DnoisE\nSWARM\nDnoisE')
    quit()
  }

  sample_list <- gsub("_HELA_nonchimeras.fasta","",list.files(pattern="^[a-zA-Z0-9]{4}_[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}_HELA_nonchimeras.fasta$"))

  if (algorithm=="dnoise_swarm" | algorithm=="dnoise") {
    for (file in sample_list) {
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
    system(paste0("obi grep-p \"sequence[\'score\'] > 1\" ",
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


  system(paste0("obi export --tab-output --sep ','  -o ",
                filetab, " ",
                lib, "_ODIN/seq_id"),
         intern = T, wait = T)

  outfile <- paste0(lib,"_ODIN")

  if (algorithm=="swarm_dnoise" | algorithm=="swarm" | algorithm=="dnoise_swarm"){
    message("ODIN will cluster sequences into MOTUs with SWARM.")
    system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",lib,"_ODIN/seq_id > ",outfile,".fasta ; ",
                  "sed -i 's/COUNT/size/g' ",outfile,".fasta ; ",
                  "sed -i 's/;//g' ",outfile,".fasta ; ",
                  "sed -E -i 's/(size=[0-9]*).*/\\1;/g' ",lib,".fasta ; ",
                  "sed -i 's/ /;/g' ",outfile,".fasta "))
    system(paste0(swarm," -d ",d," -z -t ",cores," -o ",outfile,"_SWARM_output -s ",outfile,"_SWARM",d,"nc_stats -w ",outfile,"_SWARM_seeds.fasta ",outfile,".fasta"),intern=T,wait=T)
    message("ODIN will recount abundances for every MOTU after Swarm.")

    fileswarm <- paste0(outfile,"_SWARM_output")
    outfile_MOTU <- paste0(outfile,"_counts.tsv")
    outfile_ESV <- paste0(outfile,"_ESV.tsv")
    outfile_DnoisE <- paste0(outfile,"_SWARM_DnoisE")
    get_swarm_size <- function(cadena="="){
      # it gets the positon of the '=' and gets the items after it
      return(as.numeric(substr(cadena,gregexpr("=",cadena)[[1]][[1]]+1,nchar(cadena))))
    }

    # Read cluster list database
    message("ODIN is reading SWARM results...")
    swarm_db <- readLines(fileswarm)
    # remove the last ';' of each cluster
    swarm_db <- gsub(';$','',swarm_db)
    total_swarms <- length(swarm_db)
    message("ODIN has read ", total_swarms," total MOTUs.")

    # Calculate reads in each cluster and reduce the dataset if min_reads_MOTU>0
    if (min_reads_MOTU>0) {
      message("ODIN will now calculate the number of reads in every sample for each MOTU.")
      clusters <- strsplit(swarm_db,"; ")
      # first reduction of a proportion of MOTUs
      if (min_reads_MOTU>9) i <- 9 else i <- min_reads_MOTU-1
      clusters <- clusters[!(grepl(paste0('size=[0-',i,']'),clusters) & lengths(clusters)==1)]
      cluster_reads  <- NULL
      # second reduction
      for (i in 1:length(clusters)) cluster_reads[i] <- sum(as.numeric(lapply(X=(clusters[[i]]),FUN=get_swarm_size)))
      clusters <- clusters[cluster_reads>=min_reads_MOTU]
      total_swarms <- length(clusters)
    }

    ID <- NULL
    for (i in 1:total_swarms) for (j in 1:length(clusters[[i]])) {
      clusters[[i]][[j]] <- sub(";.*","",clusters[[i]][[j]])
      ID[i] <- clusters[[i]][1]
    }

    names(clusters) <- ID

    message("ODIN kept only ", total_swarms," MOTUs of size greater than or equal to ",min_reads_MOTU," reads.")
    motu_names <- unlist(clusters, use.names=F)

    # # Generate a file with the list of ids of non-singleton clusters
    # motulist <- file(paste0(lib,"_non_singleton_motu_list.txt"),"wt")
    # writeLines(id,motulist)
    # message("ODIN has created the file ",paste0(lib,"_non_singleton_motu_list.txt")," with the list of identifiers of non-singleton MOTUs.")

    # Read counts database and keep only the needed clusters
    message("ODIN is reading the abundance database. This could take Him a while, since He has just one eye left, after all.")
    db <- read.table(filetab,sep=",",head=T)
    numseqs <- nrow(db)
    db <- db[db$ID %in% motu_names,]
    numseqs_reduced <- nrow(db)
    samples <- sum(grepl('sample',names(db)))
    message("ODIN finished reading the Database, which includes ", numseqs," total unique sequences and ",samples," samples.")
    message("ODIN kept only ", numseqs_reduced," sequences for calculations.")
    # create a database with the seeds of each cluster
    db.total <- merge(data.frame(ID),db,by="ID")

    # add MOTU name to each ESV and add all reads of the same MOTU in the db.total
    db$MOTU <- 'no_motu'
    for (motu in motu_names) {
      motu_seqs <- unlist(clusters[names(clusters)==motu])
      db$MOTU[db$ID %in% motu_seqs] <- motu
      db.total[db.total$ID==motu,grepl('sample',names(db))] <- colSums(db[db$MOTU==motu,grepl('sample',names(db))])
    }
    db.total$COUNT <- rowSums(db.total[,grepl('sample',names(db.total))])

    # remove all ESV not assigned to a retained MOTU
    db <- db[!db$MOTU=='no_motu',]

    # order the columns
    col_order <- c('ID', 'COUNT', 'MOTU', names(db)[grepl('sample',names(db))], 'NUC_SEQ' )
    db <- db[,col_order]
    s_opt <- min(grep('sample',names(db)))
    z_opt <- max(grep('sample',names(db)))

    col_order_total <- c('ID', 'COUNT', names(db)[grepl('sample',names(db))], 'NUC_SEQ' )
    db.total <- db.total[,col_order_total]

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

    if (algorithm=='swarm_dnoise') {
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

  } else {

  }






  message("ODIN will now remove MOTUs with total abundance less than ",min_reads_MOTU," from the fasta output file, to decrease THOR's workload.")
  system(paste0("sed -i 's/;size/ size/g' ",lib,"_SWARM_seeds.fasta"),intern=T,wait=T)
  system(paste0("obigrep -p 'size>",(min_reads_MOTU-1),"' ",lib,"_SWARM_seeds.fasta > ",lib,"_seeds_abundant.fasta"),intern=T,wait=T)
  message("ODIN is done.")
}

