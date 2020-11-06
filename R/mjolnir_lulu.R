mjolnir_lulu <- function(lib,min_id=.84){
  message("Producing a pairwise match list for LULU.")
  system(paste0("vsearch --usearch_global ",lib,".seeds_nonsingleton.fasta --db ",lib,".seeds_nonsingleton.fasta --self --id ",min_id," --iddef 1 --userout ",lib,".match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"),intern=T,wait=T)
  message("Removing pseudogenes with LULU.")
  system(paste0("owi_LULU ",lib),intern=T,wait=T)
  message("Done!")
}

