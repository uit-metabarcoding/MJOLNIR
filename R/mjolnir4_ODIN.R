# ODIN: OTU Delimitation Inferred by Networks

## ODIN performs MOTU clustering and (optionally) entropy denoising and it is one of the main steps of MJOLNIR.
## Clustering is performed using SWARM, which will produce MOTUs as networks of unique sequences 
## After SWARM, ODIN will recalculate the abundances for every MOTU in every sample.
## Then (optionally) ODIN proceeds with within-MOTU denoising, using the DnoisE entropy-ratio algorithm for coding regions to get an ESV table. 
## Two obligatory arguments are needed: the name of the library, typically 4 characters, and the number of computing cores.
## Three optional parameters: the clustering distance d (default=13), the minimum number of reads to keep a MOTU in the MOTU table (default=2),
## and the minimum number of reads to keep an ESV in the final ESV file (default=2).
## Two boolean parameters can be selected: run_swarm can be set as FALSE to save time if a SWARM output is already available.
## And generate_ESV=TRUE (default) will use the DnoisE algorithm to produce an ESV table, along with the MOTU table.
## ODIN deprecates the previous owi_recount_swarm script used in old metabarcoding pipelines (e.g. Project Metabarpark 2015).
## By Owen S. Wangensteen

mjolnir4_ODIN <- function(lib,cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,run_swarm=TRUE,generate_ESV=TRUE){
  # Minimum MOTU abundance to be kept in the output file can be selected. It is 2 by default.
  dnoise_path <- "~/DnoisE/src/DnoisE.py"
  if (run_swarm){
    message("ODIN will cluster sequences into MOTUs with SWARM.")
    system(paste0("swarm -d ",d," -z -t ",cores," -o ",lib,".SWARM_output -s ",lib,".SWARM",d,"nc_stats -w ",lib,".SWARM_seeds.fasta ",lib,".vsearch.fasta"),intern=T,wait=T)
    message("ODIN will recount abundances for every MOTU after Swarm.")
  }
  fileswarm=paste0(lib,".SWARM_output")
  filetab=paste0(lib,".new.tab")
  outfile <-paste(fileswarm,".counts.csv",sep="")
  outfile_ESV <-paste(fileswarm,".ESV.csv",sep="")
  
  get_swarm_size <- function(cadena="="){
    return(as.numeric(substr(cadena,gregexpr("=",cadena)[[1]][[1]]+1,nchar(cadena))))
  }

  # Read cluster list database
  message("ODIN is reading SWARM results...")
  swarm_db <- readLines(fileswarm)
  total_swarms <- length(swarm_db)
  message("ODIN has read ", total_swarms," total MOTUs.")

  # Calculate reads in each cluster
  message("ODIN will now calculate the number of reads in every sample for each MOTU.")
  clusters <- strsplit(swarm_db,"; ")
  for (i in 1:total_swarms) for (j in 1:length(clusters[[i]])) if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))==";"){
    clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-1)
  }
  cluster_reads  <- NULL
  for (i in 1:total_swarms) cluster_reads[i] <- sum(as.numeric(lapply(X=(clusters[[i]]),FUN=get_swarm_size)))
  swarm_db_reduced <- swarm_db[cluster_reads>=min_reads_MOTU]
  clusters <- strsplit(swarm_db_reduced,"; ")
  total_swarms_reduced <- length(swarm_db_reduced)

  id <- NULL
  for (i in 1:total_swarms_reduced) for (j in 1:length(clusters[[i]])) {
    clusters[[i]][[j]] <- sub(";.*","",clusters[[i]][[j]])
    id[i] <- clusters[[i]][1]
  }

  names(clusters) <- id

  message("ODIN kept only ", total_swarms_reduced," MOTUs of size greater than or equal to ",min_reads_MOTU," reads.")
  necesarios <- unlist(clusters, use.names=F)
  
  # Generate a file with the list of ids of non-singleton clusters
  motulist <- file(paste0(lib,"_non_singleton_motu_list.txt"),"wt")
  writeLines(id,motulist)
  message("ODIN has created the file ",paste0(lib,"_non_singleton_motu_list.txt")," with the list of identifiers of non-singleton MOTUs.")
  
  # Read counts database and keep only the needed clusters
  message("ODIN is reading the abundance database. This could take Him a while, since He has just one eye left, after all.")
  db <- read.table(filetab,sep="\t",head=T)
  numseqs <- nrow(db)
  db$id <- gsub(";","",db$id)
  db <- db[db$id %in% necesarios,]
  numseqs_reduced <- nrow(db)
  samples <- length(names(db)[substr(names(db),1,6)=="sample"])
  message("ODIN finished reading the Database, which includes ", numseqs," total unique sequences and ",samples," samples.")
  message("ODIN kept only ", numseqs_reduced," sequences for calculations.")

  db.total <- merge(data.frame(id),db,by="id") # This will keep just the heads
  id <- db.total$id
  numclust <- nrow(db.total)

  if (generate_ESV) dir.create("MOTU_tsv", showWarnings = FALSE)
  for (fila in 1:numclust){
    head <- id[fila]
    tails <- unlist(clusters[names(clusters)==head])
    db.reduced <- db[db$id %in% tails,]
    if (generate_ESV) write.table(db.reduced,paste0("MOTU_tsv/",head),sep="\t",quote=F,row.names=F)
    suma <- colSums(db.reduced[,substr(names(db.total),1,6)=="sample"])
    db.total[fila,substr(names(db.total),1,6)=="sample"] <- suma
    db.total$cluster_weight[fila] <- nrow(db.reduced)
    message("Cluster ", fila, " / ",numclust, " ready, including ", db.total$cluster_weight[fila]," sequences.","\r",appendLF = FALSE)
  }
  db.total$total_reads <- rowSums(db.total[,substr(names(db.total),1,6)=="sample"])
  names(db.total[substr(names(db.total),1,6)=="sample"]) <- substr(names(db.total[substr(names(db.total),1,6)=="sample"]),8,nchar(names(db.total[substr(names(db.total),1,6)=="sample"])))
  write.table(db.total[,c(1:(ncol(db.total)-3),(ncol(db.total)-1):ncol(db.total),(ncol(db.total)-2))],outfile,sep=";",quote=F,row.names=F)
  message("File ", outfile, " written")
  
  if (generate_ESV) {
    message("ODIN will generate now a list of ESVs for every non-singleton MOTU, using DnoisE.")
    sample_cols <- (1:ncol(db.total))[substr(names(db.total),1,6)=="sample"]
    start_samp <- sample_cols[1]
    end_samp <- sample_cols[length(sample_cols)]
    suppressPackageStartupMessages(library(parallel))
    clust <- makeCluster(cores)
    X <- NULL
    for (motu in id) X <- c(X,paste0("python3 ",dnoise_path," -i MOTU_tsv/",motu," -o ",motu," -f F -F F -s ",start_samp," -z ",end_samp," -n 'count' -p 1 -y T"))
    clusterExport(clust, "X",envir = environment())
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)
    message("ODIN will now merge all ESVs into a final ESV table.")
    ESV_tables <- NULL
    for (i in 1:length(id)) {
      ESV_tables[[i]] <- read.table(paste0(id[i],"_Adcorr_denoised_ratio_d.csv"),sep=",",head=T)
      ESV_tables[[i]]$definition <- id[i]
    }
    ESV_table <- rbind(ESV_tables)
    write.table(ESV_table,outfile_ESV,sep=";",quote=F,row.names=F)
    message("File ", outfile_ESV, " written with ",nrow(ESV_table)," ESVs in ",length(id), " MOTUs")
  }
    
  message("ODIN is now removing singleton MOTUs from the fasta output file, to make THOR's work an easier task.")
  system(paste0("sed -i 's/;size/ size/g' ",lib,".SWARM_seeds.fasta"),intern=T,wait=T)
  system(paste0("obigrep -p 'size>1' ",lib,".SWARM_seeds.fasta > ",lib,".seeds_nonsingleton.fasta"),intern=T,wait=T)
  message("ODIN is done.")
}
