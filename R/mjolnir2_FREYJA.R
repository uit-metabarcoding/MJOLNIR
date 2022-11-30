# FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment

# FREYJA will use OBITools commands to merge paired-end reads, demultiplex libraries into samples (if needed),trim primer sequences, filter by length.
# In case the data are already demultiplexed and consist of individual fastq files for each sample, use the option demultiplexed=TRUE.
# When demultiplexed=TRUE, FREYJA will read the names of each individual R1 fastq files from a column in the LIBX_metadata.tsv file, called fastq_name_R1
# In the metadata table, each sample in the original_samples column must have a matching fastq_name_R1 and a matching mjolnir_agnomen (LIBX_sample_XXX).
# When demultiplexed=TRUE, you must also specify the R1_motif and R2_motif strings in the options input to FREYJA.
# When demultiplexed=TRUE, you must also specify the primer_F and primer_R sequences in the options input to FREYJA. COI Leray-XT primers are specified by default.
# Otherwise, when demultiplexed=FALSE, the primers information must be already written in the LIBX_ngsfilter.tsv files.

mjolnir2_FREYJA <- function(lib_prefix="",cores=1,Lmin=299,Lmax=320,lib="", fasta_output=F,fastq_output=F,score_obialign=40,
                            demultiplexed=F,primer_F="GGWACWRGWTGRACWNTNTAYCCYCC",primer_R="TANACYTCNGGRTGNCCRAARAAYCA",R1_motif="_R1",R2_motif="_R2",obipath=""){

  message("FREYJA will do paired-end alignment, demultiplexing and length filter.")
  suppressPackageStartupMessages(library(parallel))
  no_cores <- cores*length(lib_prefix)
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))
  clust <- makeCluster(no_cores)
  X <- NULL
  libslist <- NULL
  if (!demultiplexed){
    for (i in 1:cores) for (j in 1:length(lib_prefix)) {
      formatted_i <- sprintf("%02d",i)
      X <- c(X,paste0(
        "obi import --fastq-input ",lib_prefix[j],"_R2_part_",formatted_i,".fastq ", lib_prefix[j],"_",formatted_i,"_FREYJA/reads2 ; ",
        "obi import --fastq-input ",lib_prefix[j],"_R1_part_",formatted_i,".fastq ", lib_prefix[j],"_",formatted_i,"_FREYJA/reads1 ; ",
        "obi import --ngsfilter-input ngsfilter_",lib_prefix[j],".tsv ", lib_prefix[j],"_",formatted_i,"_FREYJA/ngsfile ; ",
        "obi alignpairedend -R ", lib_prefix[j],"_",formatted_i,"_FREYJA/reads2 ", lib_prefix[j],"_",formatted_i,"_FREYJA/reads1 ", lib_prefix[j],"_",formatted_i,"_FREYJA/aligned_seqs ; ",
        "obi grep -p \"sequence[\'score\'] > ",score_obialign,"\" ", lib_prefix[j],"_",formatted_i,"_FREYJA/aligned_seqs ", lib_prefix[j],"_",formatted_i,"_FREYJA/good_seqs ; ",
        "obi ngsfilter -t ", lib_prefix[j],"_",formatted_i,"_FREYJA/ngsfile -u ", lib_prefix[j],"_",formatted_i,"_FREYJA/unidentified_seqs ", lib_prefix[j],"_",formatted_i,"_FREYJA/good_seqs ",lib_prefix[j],"_",formatted_i,"_FREYJA/identified_seqs ; ",
        "obi grep -p \"len(sequence)>",Lmin," and len(sequence)<",Lmax," and sequence[\'forward_tag\']!=None and sequence[\'reverse_tag\']!=None\" -S \"^[ACGT]+$\" ",lib_prefix[j],"_",formatted_i,"_FREYJA/identified_seqs ",lib_prefix[j],"_",formatted_i,"_FREYJA/filtered_seqs"))
      libslist <- paste0(libslist,paste0("-c ",lib_prefix[j],"_",formatted_i,"_FREYJA/filtered_seqs "))
    }
  } else {
      metadata <- read.table(paste0(lib,"_metadata.tsv"),sep="\t",header=T)
      fastqR1_list <- metadata$fastq_name_R1
      agnomens <-  metadata$mjolnir_agnomens
      # Create ngsfilter files
      for (ag in agnomens) writeLines(paste(lib,ag,"None:None",primer_F,primer_R,sep="\t"),paste0("ngsfilter_",ag,".tsv"))
      # Create obitool commands
      for (i in 1:length(agnomens)) {
        X <- c(X,paste0(
          "obi import --fastq-input ",gsub(R1_motif,R2_motif,fastqR1_list[i]), " ", lib,"_",agnomens[i],"_FREYJA/reads2 ; ",
          "obi import --fastq-input ",fastqR1_list[i], " ", lib,"_",agnomens[i],"_FREYJA/reads1 ; ",
          "obi import --ngsfilter-input ngsfilter_",agnomens[i],".tsv ", lib,"_",agnomens[i],"_FREYJA/ngsfile ; ",
          "obi alignpairedend -R ", lib,"_",agnomens[i],"_FREYJA/reads2 ", lib,"_",agnomens[i],"_FREYJA/reads1 ", lib,"_",agnomens[i],"_FREYJA/aligned_seqs ; ",
          "obi grep -p \"sequence[\'score\'] > ",score_obialign,"\" ", lib,"_",agnomens[i],"_FREYJA/aligned_seqs ", lib,"_",agnomens[i],"_FREYJA/good_seqs ; ",
          "obi ngsfilter --no-tags -t ", lib,"_",agnomens[i],"_FREYJA/ngsfile -u ", lib,"_",agnomens[i],"_FREYJA/unidentified_seqs ", lib,"_",agnomens[i],"_FREYJA/good_seqs ",lib,"_",agnomens[i],"_FREYJA/identified_seqs_nosample ; ",
          "obi annotate -S \"sample:'", agnomens[i], "'\" ", lib,"_",agnomens[i],"_FREYJA/identified_seqs_nosample ",lib,"_",agnomens[i],"_FREYJA/identified_seqs ; ",
          "obi grep -p \"len(sequence)>",Lmin," and len(sequence)<",Lmax,"\" -S \"^[ACGT]+$\" ",lib,"_",agnomens[i],"_FREYJA/identified_seqs ",lib,"_",agnomens[i],"_FREYJA/filtered_seqs ; "))
      }
  }
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)
    # If not demultiplexed, then join all parts into a joined file and then split it into samples
    if (!demultiplexed){
      message("FREYJA is joining filtered reads into a single file.")
      # this step is necessary to join all sets of sequences to split them into different files
      system(paste0("obi cat ",libslist," ",lib,"_FREYJA/concatenated"),intern=T,wait=T)
      message("Sequence sets concatenated")
      message("FREYJA will create individual files for each sample in the dms directory, this could take a while..")
      # here the concatenated file is split into different samples
      system(paste0("obi split -t \"sample\" -p ",lib,"_ ",lib,"_FREYJA/concatenated"),intern=T,wait=T)
      # in order to make the next steps parallelizable it is nesary to export each file into a new dms
      files <- system(paste0("obi ls ",lib,"_FREYJA | cut -f1 -d ':' | cut -f4 -d ' '"),intern=T,wait=T)
      files <- files[grep(lib,files)]
      for (file in files) {
        system(paste0("obi cat -c ",
                      lib,"_FREYJA/",file, " ",
                      file,"_FREYJA/filtered_seqs "),intern=T,wait=T)
      }
    }
    # obi uniq vas performed in HELA in previous versions but now is computed here
    files <- list.dirs(recursive = F)
    files <- files[grepl('sample',files)&grepl('FREYJA.obidms',files)]
    files <- gsub('./','',gsub('.obidms','',files))
    X <- NULL
    for (file in files) {
      X <- c(X,paste0("obi uniq ",file,"/filtered_seqs ",file,"/uniq ; ",
                      "obi export --fasta-output --only-keys \"sample\" --only-keys \"COUNT\" ",file,"/uniq > ",file,"_uniq.fasta ; ",
                      "sed -i 's/COUNT/size/g' ",file,"_uniq.fasta ; ",
                      "sed -i 's/;//g' ",file,"_uniq.fasta ; ",
                      "sed -i 's/ /;/g' ",file,"_uniq.fasta ; ",
                      "sed -i -E 's/(.*);.*$/\\1;/' ",file,"_uniq.fasta "))
    }
    clust <- makeCluster(no_cores)
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)

    # if required, export to fasta or fastq. These options is to return the files without joining amplicons with obi uniq
    if (fasta_output | fastq_output) {
      files <- list.dirs(recursive = F)
      files <- files[grepl('sample',files)&grepl('FREYJA.obidms',files)]
      files <- gsub('./','',gsub('.obidms','',files))
      X <- NULL
      if (fasta_output) {
        print('Exporting fasta files, this could take a while.')
        for (file in files) {
          X <- c(X,paste0("obi export --fasta-output ",
                          file,"/filtered_seqs > ",
                          file,".fasta "))
        }
      }
      if (fastq_output) {
        print('Exporting fastq files, this could take a while')
        for (file in files) {
          X <- c(X,paste0("obi export --fastq-output ",
                          file,"/filtered_seqs > ",
                          file,".fastq "))
        }
      }
      clust <- makeCluster(no_cores)
      clusterExport(clust, list("X","old_path","obipath"),envir = environment())
      clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
      parLapply(clust,X, function(x) system(x,intern=T,wait=T))
      stopCluster(clust)

    }
    message("FREYJA is done.")
}
