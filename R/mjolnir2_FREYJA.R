# FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment  

# FREYJA will use OBITools commands to merge paired-end reads, demultiplex libraries into samples (if needed),trim primer sequences, filter by length.
# In case the data are already demultiplexed and consist of individual fastq files for each sample, use the option demultiplexed=TRUE.
# When demultiplexed=TRUE, FREYJA will read the names of each individual R1 fastq files from a column in the LIBX_metadata.tsv file, called fastq_name_R1
# In the metadata table, each sample in the original_samples column must have a matching fastq_name_R1 and a matching mjolnir_agnomen (LIBX_sample_XXX).  
# When demultiplexed=TRUE, you must also specify the R1_motif and R2_motif strings in the options input to FREYJA. 
# When demultiplexed=TRUE, you must also specify the primer_F and primer_R sequences in the options input to FREYJA. COI Leray-XT primers are specified by default.
# Otherwise, when demultiplexed=FALSE, the primers information must be already written in the LIBX_ngsfilter.tsv files.

mjolnir2_FREYJA <- function(lib_prefix="",cores=1,Lmin=299,Lmax=320,lib="", fasta_output=T,
                            demultiplexed=F,primer_F="GGWACWRGWTGRACWNTNTAYCCYCC",primer_R="TANACYTCNGGRTGNCCRAARAAYCA",R1_motif="_R1",R2_motif="_R2",obipath=""){
  message("FREYJA will do paired-end alignment, demultiplexing and length filter.")  
  suppressPackageStartupMessages(library(parallel))
  no_cores <- cores*length(lib_prefix)
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":")) 
  clust <- makeCluster(no_cores)
  X <- NULL
  libslist <- NULL
  ext_file <- ifelse(fasta_output,".fasta",".fastq")
  fasta_modifier <- ifelse(fasta_output," --fasta-output","")
  if (!demultiplexed){
    for (i in 1:cores) for (j in 1:length(lib_prefix)) {
      formatted_i <- sprintf("%02d",i)
      X <- c(X,paste0("illuminapairedend -r ",lib_prefix[j],"_R2_part_",formatted_i,".fastq ",lib_prefix[j],"_R1_part_",formatted_i,
                                        ".fastq | obigrep -p \'score>40.00\' | ngsfilter -t ngsfilter_",lib_prefix[j],".tsv | obigrep -p \'seq_length>",
                                        Lmin,"\' -p \'seq_length<",Lmax,"\' -s \'^[ACGT]+$\' -p \'forward_tag!=None\' -p \'reverse_tag!=None\'",fasta_modifier," > ",
                                        lib_prefix[j],"_filtered_length_part_",sprintf("%02d",i),ext_file)) }
    for (prefix in lib_prefix) libslist <- paste0(libslist,prefix,"_filtered_length_part*",ext_file," ")
  } else {
      metadata <- read.table(paste0(lib,"_metadata.tsv"),sep="\t",header=T)
      fastqR1_list <- metadata$fastq_name_R1
      agnomens <-  metadata$mjolnir_agnomens
      # Create ngsfilter files
      for (ag in agnomens) writeLines(paste(lib,ag,":",primer_F,primer_R,sep="\t"),paste0("ngsfilter_",ag,".tsv"))
      # Create obitool commands
      for (i in 1:length(agnomens)) {
         X <- c(X,paste0("illuminapairedend -r ",gsub(R1_motif,R2_motif,fastqR1_list[i]), " ",fastqR1_list[i]," | obigrep -p \'score>40.00\' | ngsfilter -t ngsfilter_",agnomens[i],".tsv | obigrep -p \'seq_length>",Lmin,"\' -p \'seq_length<",Lmax,"\' -s \'^[ACGT]+$\' -p \'forward_tag!=None\' -p \'reverse_tag!=None\' --fasta-output > ",agnomens[i],".fasta"))
         libslist <- paste0(libslist,agnomens[i],".filtered_length_part_",sprintf("%02d",i),".fasta ")
      }
  }
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)
    # If not demultiplexed, then join all parts into a joined file and then split it into samples
    if (!demultiplexed){
      message("FREYJA is joining filtered reads into a single file.")
      system(paste0("cat ",libslist," > ",lib,"_joined",ext_file),intern=T,wait=T)
      message("File ",lib,".joined",ext_file," written.")
      message("FREYJA will create individual files for each sample.")
      system(paste0("obisplit -t sample ",lib,"_joined",ext_file),intern=T,wait=T)
    }
    message("FREYJA is done.")
}
