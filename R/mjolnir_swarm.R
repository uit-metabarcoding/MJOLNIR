mjolnir_swarm <- function(lib,cores,d=13){
  message("Clustering with Swarm.")
  system(paste0("swarm -d ",d," -z -t ",cores," -o ",lib,".SWARM_output -s ",lib,".SWARM",d,"nc_stats -w ",lib,".SWARM_seeds.fasta ",lib,".vsearch.fasta"),intern=T,wait=T)
  message("Recounting after Swarm.")
  owi_recount_swarm(fileswarm=paste0(lib,".SWARM_output"),filetab=paste0(lib,".new.tab"))
  message("Removing singletons.")
  system(paste0("sed -i 's/;size/ size/g' ",lib,".SWARM_seeds.fasta"),intern=T,wait=T)
  system(paste0("obigrep -p 'size>1' ",lib,".SWARM_seeds.fasta > ",lib,".seeds_nonsingleton.fasta"),intern=T,wait=T)
  message("Done.")
}
