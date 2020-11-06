## Script for combining an abundance file with an ecotag annotated file
## The script will read three arguments from the command line: two input file names and the output file name.
## The abundances file must include sample columns names beginning with "sample."
## If the output file name is empty, it will just add ".merged.csv" at the end of the name of the annotated input file.
## By Owen S. Wangensteen - Project Metabarpark  2016

owi_combine <- function(infile,abundances,outfile=NULL,sept=";;") {
  # sept = separator characters used in annotated file and abundances file, respectively (default: ';;' )
  if (is.null(outfile)) outfile <- paste0(substr(infile,1,nchar(infile)-3),"combined.csv")

  # Read cluster list database
  message("Reading ecotag database...")
  ecotag_db <- read.table(infile,sep=substr(sept,1,1),head=T,stringsAsFactors=F)
  message("Ecotag database read including ", nrow(ecotag_db)," total MOTUs.")
  # Delete "None" from the taxonomy database
  ecotag_db[ecotag_db=="None"] <- ""

  message("Reading abundance database...")
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
