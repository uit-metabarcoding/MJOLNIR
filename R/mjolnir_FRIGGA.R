# FRIGGA: Final Recount and Integration of Generated Genealogies and Abundances

## Script for combining an abundance CSV file from ODIN with a taxonomy-annotated CSV file from THOR.
## The script will have the name of the library (typically 4 characters, e.g. LIBR) as input.
## The taxonomy file must be called LIBR.ecotag.fasta.annotated.csv
## The abundance file must be called LIBR.SWARM_output.counts.csv
## The separator characters for both the taxonomy and abundance CSV files can be also used as input.
## By default, they are both are semicolons ";;".
## The abundances file must include sample columns names starting with "sample."
## The output file is a CSV file called LIBR.All_MOTUs.csv
## FRIGGA deprecates the owi_combine function from previous pipelines (Project Metabarpark, 2016).
## By Owen S. Wangensteen 

mjolnir_FRIGGA <- function(lib=NULL,sept=";;"){
  # sept = separator characters used in taxonomy-annotated file and abundances file, respectively (default: ';;' )
  message("Producing the combined file.")

  infile=paste0(lib,".ecotag.fasta.annotated.csv")
  abundances=paste0(lib,".SWARM_output.counts.csv")
  outfile=paste0(lib,".All_MOTUs.csv"))

  message("Reading ecotag-annotated database from THOR...")
  ecotag_db <- read.table(infile,sep=substr(sept,1,1),head=T,stringsAsFactors=F)
  message("Ecotag database read including ", nrow(ecotag_db)," total MOTUs.")
  # Delete "None" from the taxonomy database
  ecotag_db[ecotag_db=="None"] <- ""

  message("Reading abundances database from ODIN...")
  abun_db <- read.table(abundances,sep=substr(sept,2,2),head=T,stringsAsFactors=F)
  n_samples <- ncol(abun_db[,substr(names(abun_db),1,7)=="sample."])
  message("Abundances database read including ", nrow(abun_db)," total MOTUs and ",n_samples," samples.")

  # Merge databases
  db <- merge(ecotag_db,abun_db,by="id")
  db$sequence <- db$sequence.y
  db <- db[substr(names(db),1,9)!="sequence."]
  names(db)[substr(names(db),1,7)=="sample."] <- substr(names(db)[substr(names(db),1,7)=="sample."],8,nchar(names(db)[substr(names(db),1,7)=="sample."]))

  # Delete unnecessary columns
  db <- db[,!(names(db) %in% c("definition","ali_length","avg_quality","count","direction","experiment","forward_match","forward_primer","forward_score",
                             "forward_tag","goodali","head_quality","mid_quality","mode","position","reverse_match","reverse_primer","reverse_score","reverse_tag","score","score_norm",
                             "seq_a_deletion","seq_a_insertion","seq_a_mismatch","seq_a_single","seq_ab_match","seq_b_deletion","seq_b_insertion","seq_b_mismatch","seq_b_single",
                             "seq_length_ori","seq_rank","status","tail_quality"))]

  # Reorder total_counts and cluster weight columns
  db <- db[,c(1:4,ncol(db)-1,ncol(db)-2,5:ncol(db)-3,ncol(db))]

  write.table(db,outfile,sep=";",quote=F,row.names=F)
  message("File ", outfile, " written, including ",nrow(db)," MOTUs with ",sum(db$total_reads)," total reads in ",n_samples," samples.")
  message("(",nrow(db[db$total_reads>1,])," non-singletons MOTUs).")
}

