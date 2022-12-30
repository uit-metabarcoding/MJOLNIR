# THOR: Taxonomy with Higher-than-Order Ranks
# This is a wrapper of ecotag
# After assignment with ecotag, higher taxa at ranks higher than order are added from custom CSV files.
# THOR replaces owi_add_taxonomy

mjolnir5_THOR <- function(lib,cores,tax_dir,tax_dms=NULL,ref_db=NULL,taxo_db=NULL,obipath="",run_ecotag=T){
  setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/adriantich_tests/');
  # setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/example_MJOLNIR_demultiplexed_data/');
  cores=3; obipath="/home/adriantich/obi3-env/bin/";
  lib='ULOY'
  tax_dir='taxo_dir'
  tax_db='taxdump.tar.gz'
  ref_db='DUFA_proves_small.fasta'
  # tax_dms=NULL
  tax_dms=paste0(lib, "_THOR_taxo")


  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))

  if (run_ecotag) {
    message("THOR will assign the taxonomy to the order level with ecotag.")

    if (is.null(tax_dms)) {
      system(paste0("obi import ",tax_dir,"/",ref_db," ",lib, "_THOR_taxo/ref_seqs ; ",
                    "obi import --taxdump ",tax_dir,"/",tax_db," ",lib, "_THOR_taxo/taxonomy/my_tax ; ",
                    "obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy ",lib, "_THOR_taxo/taxonomy/my_tax ",lib, "_THOR_taxo/ref_seqs ",lib, "_THOR_taxo/ref_seqs_clean ; ",
                    "obi uniq --taxonomy ",lib, "_THOR_taxo/taxonomy/my_tax ",lib, "_THOR_taxo/ref_seqs_clean ",lib, "_THOR_taxo/ref_seqs_uniq ; ",
                    "obi build_ref_db -t 0.95 --taxonomy ",lib, "_THOR_taxo/taxonomy/my_tax ",lib, "_THOR_taxo/ref_seqs_uniq ",lib, "_THOR_taxo/ref_seqs_db "),
             intern = T, wait = T)
      tax_dms <- paste0(lib, "_THOR_taxo")
    }

    # it is necessary to run ecotag within a new directory for each part.
    # this is because dms can not be calles from two processes at the same time
    # and can not change the name of the dms so make a copy in each directory and
    # run there the ecotag
    X <- NULL
    for (i in 1:cores) {
      system(paste0("mkdir ",lib, "_THOR_",sprintf("%02d",i)," ; cp -r ",tax_dms,".obidms ",lib, "_THOR_",sprintf("%02d",i),"/. ; "),intern = T, wait = T)
      X <- c(X,paste0("cd ",lib, "_THOR_",sprintf("%02d",i)," ; ",
                      "obi import --fasta-input ../",lib,"_ODIN_part_",sprintf("%02d",i),".fasta ",tax_dms,"/seqs ; ",
                      # "obi import --fasta-input ",lib,"_ODIN_part_",sprintf("%02d",i),".fasta ",lib, "_THOR_",sprintf("%02d",i),"/seqs ; ",
                      "obi ecotag --taxonomy ",tax_dms,"/taxonomy/my_tax -R ",tax_dms,"/ref_seqs_db ",tax_dms,"/seqs ",tax_dms,"/assigned_seqs ; ",
                      "obi export --fasta-output ",tax_dms,"/assigned_seqs >../",lib,"_THOR_part_",sprintf("%02d",i),".fasta"))
      # "obi ecotag --taxonomy ",lib, "_THOR_",sprintf("%02d",i),"/taxonomy/my_tax -R ",lib, "_THOR_",sprintf("%02d",i),"/ref_seqs_db ",lib, "_THOR_",sprintf("%02d",i),"/seqs ",lib, "_THOR_",sprintf("%02d",i),"/assigned_seqs"))
    }

    suppressPackageStartupMessages(library(parallel))
    no_cores <- cores
    clust <- makeCluster(no_cores)
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)
  }

  message("THOR will add higher taxonomic ranks now.")
  filefasta <-paste0(lib,"_THOR.fasta")
  system(paste0("cat ",lib,"_THOR_part_??.fasta > ",filefasta),intern=T,wait=T)
  outfile <- paste0(lib,"_THOR_annotated.tsv")
  # Here old owi_add_taxonomy starts
  suppressPackageStartupMessages(library("Biostrings"))
  length_id <- 14 # This is the total length of the MOTU IDs in filefasta. It can be changed if needed.
  # Load external information on higher taxa:
  taxo_names <- read.table(paste0(tax_dir,"/","order.complete.csv"),sep=",",head=T)
  class_to_sk <- unique(taxo_names[,2:5])
  phylum_to_sk <- unique(class_to_sk[,2:4])
  kingdom_to_sk <-unique(phylum_to_sk[,2:3])
  showlines <- 1

  family_to_order <- read.table(paste0(tax_dir,"/","family_to_order.csv"),sep=",",head=T)
  genus_to_family <- read.table(paste0(tax_dir,"/","genus_to_family.csv"),sep=",",head=T)

  get_rank <- function(cadena){
    if (gregexpr("rank=",cadena)[[1]]>0)
    {cadena <- substr(cadena,gregexpr("rank=",cadena)[[1]]+5,nchar(cadena))
    cadena <- substr(cadena,1,regexpr(";",cadena)[[1]]-1)
    }  else cadena <- "None"
    return(cadena)
  }

  # The following are taxonomic ranks that can be found in the "scientific_name" field and should be treated as exceptions, by function fix_exceptions.
  exceptions <- c("Bacillariophycidae","Bilateria","Collembola","Chelicerata","Crustacea","Deuterostomia","Ecdysozoa",
                  "Eleutherozoa","Echinozoa","Asterozoa","Endopterygota","Eumalacostraca","Eumetazoa","Eurotiomycetidae","Euteleostomi",
                  "Euthyneura","Hydroidolina","hypocreomyceta","Hypocreomycetidae","Hypsogastropoda","leotiomyceta","Leotiomycetidae","Lophotrochozoa",
                  "Myxogastria","Neocopepoda","Podoplea","Neoptera","Octocorallia","Hexacorallia","Scleractinia","Opisthokonta","Panarthropoda",
                  "Pancrustacea","Peracarida","Protostomia","PX clade","Rhodophyta","Sarcoptiformes","Scolecida","Aciculata","Palpata",
                  "sordariomyceta","Sordariomycetidae","Stramenopiles","Stylommatophora","Trombidiformes","Acariformes","Acari",
                  "Pterygota","Paleoptera","Nudipleura","Panpulmonata","Heterobranchia","Caenogastropoda","Eucarida","Hoplocarida",
                  "Streptophytina","Entomobryomorpha","Poduromorpha","Symphypleona","Fungi incertae sedis","Haptophyceae","Jakobida","Alveolata",
                  "Palaeonemertea","Euopisthobranchia","Euteleosteomorpha","Littorinimorpha","Littorinoidea","Fragilariophycidae","Mandibulata",
                  "Actinopterygii","Batoidea","Boreoeutheria","Cephalaspidea","Clitellata","Clupeocephala","Decapodiformes","Echinacea","Euechinoidea",
                  "Embryophyta","Eupercaria","Percomorphaceae","Ophiuridea","Patellogastropoda","Piroplasmida","Pteriomorphia","Synurophyceae",
                  "Vetigastropoda","CirripFedia","Dikarya","Euacanthomorphacea","Galeoidea","Gnathostomata","Intramacronucleata","Neogastropoda",
                  "Phascolosomatidea","Trachylinae","Hexapoda","Oligochaeta","Sacoglossa","Thoracica","Buccinoidea","Digenea",
                  "Agaricomycetes incertae sedis","Agaricomycetidae","Agaricomycotina","Amniota","Amoebozoa","Anystina","Apusozoa","Armophorea",
                  "Bryophyta","Centroheliozoa","Cercozoa","Enoplia","asterids","campanulids","lamiids","rosids","Magnoliophyta","Mesangiospermae",
                  "Pentapetalae","Pezizomycotina","Pinidae","Poduroidea","Prostigmata","Rhabditophora","Scuticociliatia","Silicofilosea","Sipuncula",
                  "Spermatophyta","Tubulinea","Ustilaginomycotina","Xanthophyceae","Petrosaviidae","Peritrichia","Moniliformopses",
                  "Lecanoromycetidae","Labyrinthulomycetes","Hydracarina","Gregarinasina","Dactylopodida","Cymbellales","eudicotyledons",
                  "Pucciniomycotina","Taphrinomycotina","Dorylaimia","Teleostei","Stichotrichia","Peritrichia","Oligotrichia",
                  "Choreotrichia","Ciliophora","Cryptophyta","Trichostomatia","Seriata","Tunicata","Tubulinida","Paraneoptera","saccharomyceta",
                  "Euphyllophyta","Coscinodiscophycidae","Corallinophycidae","Biddulphiophycidae","Bdelloidea","Acanthomorphata","Leotiomycetes incertae sedis",
                  "Ctenosquamata","Euteleosteomorpha","Heteroscleromorpha","Hexanauplia","Imparidentia","Neoloricata","Ochrophyta","Prymnesiophyceae",
                  "Rhodymeniophycidae","Sedentaria","Sar","Nemaliophycidae","Heteroconchia","Euheterodonta","Palaeoheterodonta",
                  "Picobiliphyte sp. MS584-11","beta proteobacterium CB","Chrysophyceae sp.","Bacteria","uncultured bacterium","uncultured actinobacterium"
  )

  fix_exceptions <- function(scientific_name){
    if (scientific_name %in% c("asterids","campanulids","lamiids","rosids","Pentapetalae","Mesangiospermae","eudicotyledons","Magnoliophyta")) {
      matrix.data["class_name",2] <- "Magnoliopsida"
      matrix.data["phylum_name",2] <- "Tracheophyta"
      matrix.data["kingdom_name",2] <- "Viridiplantae"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
    }
    if (scientific_name %in% c("Pinidae")) {
      matrix.data["class_name",2] <- "Pinopsida"
      matrix.data["phylum_name",2] <- "Tracheophyta"
      matrix.data["kingdom_name",2] <- "Viridiplantae"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
    }
    if (scientific_name %in% c("Moniliformopses","Euphyllophyta")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "Tracheophyta"
      matrix.data["kingdom_name",2] <- "Viridiplantae"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
    }
    if (scientific_name %in% c("Bacteria","uncultured bacterium","uncultured actinobacterium")) {
      matrix.data["kingdom_name",2] <- "Bacteria"
      matrix.data["superkingdom_name",2] <- "Prokaryota"
    }
    if (scientific_name %in% c("beta proteobacterium CB")) {
      matrix.data["class_name",2] <- "Betaproteobacteria"
      matrix.data["phylum_name",2] <- "Proteobacteria"
      matrix.data["kingdom_name",2] <- "Bacteria"
      matrix.data["superkingdom_name",2] <- "Prokaryota"
    }
    if (scientific_name %in% c("Stichotrichia","Peritrichia","Oligotrichia","Choreotrichia","Ciliophora","Trichostomatia")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "Ciliophora"
      matrix.data["kingdom_name",2] <- "Alveolata"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }

    if (scientific_name %in% c("Petrosaviidae")) {
      matrix.data["class_name",2] <- "Liliopsida"
      matrix.data["phylum_name",2] <- "Tracheophyta"
      matrix.data["kingdom_name",2] <- "Viridiplantae"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
    }
    if (scientific_name %in% c("Centroheliozoa")) {
      matrix.data["class_name",2] <- "Centrohelea"
      matrix.data["phylum_name",2] <- "Heliozoa"
      matrix.data["kingdom_name",2] <- "Hacrobia"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Sar")) {
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Cryptophyta")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "Cryptophyta"
      matrix.data["kingdom_name",2] <- "Hacrobia"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Prymnesiophyceae")) {
      matrix.data["class_name",2] <- "Prymnesiophyceae"
      matrix.data["phylum_name",2] <- "Haptophyta"
      matrix.data["kingdom_name",2] <- "Hacrobia"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Picobiliphyte sp. MS584-11")) {
      matrix.data["class_name",2] <- "Picomonadea"
      matrix.data["phylum_name",2] <- "Picozoa"
      matrix.data["kingdom_name",2] <- "Hacrobia"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Chrysophyceae sp.")) {
      matrix.data["class_name",2] <- "Chrysophyceae"
      matrix.data["phylum_name",2] <- "Ochrophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Ochrophyta")) {
      matrix.data["phylum_name",2] <- "Ochrophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "phylum"
    }
    if (scientific_name %in% c("Cymbellales")) {
      matrix.data["class_name",2] <- "Bacillariophyceae"
      matrix.data["phylum_name",2] <- "Bacillariophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Coscinodiscophycidae")) {
      matrix.data["class_name",2] <- "Coscinodiscophyceae"
      matrix.data["phylum_name",2] <- "Bacillariophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }

    if (scientific_name %in% c("Biddulphiophycidae")) {
      matrix.data["class_name",2] <- "Mediophyceae"
      matrix.data["phylum_name",2] <- "Bacillariophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }

    if (scientific_name %in% c("Enoplia","Dorylaimia")) {
      matrix.data["class_name",2] <- "Enoplea"
      matrix.data["phylum_name",2] <- "Nematoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Rhabditophora","Seriata")) {
      matrix.data["class_name",2] <- "Rhabditophora"
      matrix.data["phylum_name",2] <- "Platyhelminthes"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Gregarinasina")) {
      matrix.data["class_name",2] <- "Conoidasida"
      matrix.data["phylum_name",2] <- "Apicomplexa"
      matrix.data["kingdom_name",2] <- "Alveolata"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name %in% c("Tubulinea","Tubulinida")) {
      matrix.data["class_name",2] <- "Tubulinea"
      matrix.data["phylum_name",2] <- "Tubulinea"
      matrix.data["kingdom_name",2] <- "Lobosa"
      matrix.data["superkingdom_name",2] <- "Amoebozoa"
    }
    if (scientific_name %in% c("Dactylopodida")) {
      matrix.data["class_name",2] <- "Flabellinia"
      matrix.data["phylum_name",2] <- "Discosea"
      matrix.data["kingdom_name",2] <- "Lobosa"
      matrix.data["superkingdom_name",2] <- "Amoebozoa"
    }
    if (scientific_name %in% c( "Amoebozoa","Apusozoa")) {
      matrix.data["superkingdom_name",2] <- scientific_name
    }
    if (scientific_name %in% c("Pezizomycotina","Taphrinomycotina")) {
      matrix.data["phylum_name",2] <- "Ascomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("saccharomyceta")) {
      matrix.data["class_name",2] <- "Saccharomycetes"
      matrix.data["phylum_name",2] <- "Ascomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Lecanoromycetidae")) {
      matrix.data["class_name",2] <- "Lecanoromycetes"
      matrix.data["phylum_name",2] <- "Ascomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name %in% c("Agaricomycotina","Pucciniomycotina","Ustilaginomycotina")) {
      matrix.data["phylum_name",2] <- "Basidiomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name %in% c("Agaricomycetes incertae sedis","Agaricomycetidae")) {
      matrix.data["class_name",2] <- "Agaricomycetes"
      matrix.data["phylum_name",2] <- "Basidiomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Bryophyta")) {
      matrix.data["phylum_name",2] <- "Bryophyta"
      matrix.data["kingdom_name",2] <- "Viridiplantae"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
    }

    if (scientific_name %in% c("Armophorea")) {
      matrix.data["class_name",2] <- "Armophorea"
      matrix.data["phylum_name",2] <- "Ciliophora"
      matrix.data["kingdom_name",2] <- "Alveolata"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name == "Cercozoa") {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "Cercozoa"
      matrix.data["kingdom_name",2] <- "Rhizaria"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name == "Silicofilosea") {
      matrix.data["class_name",2] <- "Imbricatea"
      matrix.data["phylum_name",2] <- "Cercozoa"
      matrix.data["kingdom_name",2] <- "Rhizaria"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }

    if (scientific_name == "Jakobida") {
      matrix.data["order_name",2] <- "Jakobida"
      matrix.data["class_name",2] <- "Jakobea"
      matrix.data["phylum_name",2] <- "Loukozoa"
      matrix.data["kingdom_name",2] <- "Loukozoa"
      matrix.data["superkingdom_name",2] <- "Excavata"
    }
    if (scientific_name == "Piroplasmida") {
      matrix.data["order_name",2] <- "Piroplasmorida"
      matrix.data["class_name",2] <- "Aconoidasida"
      matrix.data["phylum_name",2] <- "Apicomplexa"
      matrix.data["kingdom_name",2] <- "Alveolata"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name == "Intramacronucleata") {
      matrix.data["class_name",2] <- "Intramacronucleata"
      matrix.data["phylum_name",2] <- "Ciliophora"
      matrix.data["kingdom_name",2] <- "Alveolata"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name == "Scuticociliatia") {
      matrix.data["class_name",2] <- "Oligohymenophorea"
      matrix.data["phylum_name",2] <- "Ciliophora"
      matrix.data["kingdom_name",2] <- "Alveolata"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
    }
    if (scientific_name == "Heteroscleromorpha") {
      matrix.data["class_name",2] <- "Demospongiae"
      matrix.data["phylum_name",2] <- "Porifera"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name == "Palaeonemertea") {
      matrix.data["order_name",2] <- "Palaeonemertea"
      matrix.data["class_name",2] <- "Anopla"
      matrix.data["phylum_name",2] <- "Nemertea"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "order"
    }
    if (scientific_name %in%  c("Phascolosomatidea")) {
      matrix.data["order_name",2] <- "Phascolosomatiformes"
      matrix.data["class_name",2] <- "Sipuncula"
      matrix.data["phylum_name",2] <- "Annelida"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Sipuncula")) {
      matrix.data["class_name",2] <- "Sipuncula"
      matrix.data["phylum_name",2] <- "Annelida"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Bdelloidea")) {
      matrix.data["class_name",2] <- "Eurotatoria"
      matrix.data["phylum_name",2] <- "Rotifera"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name == "Digenea") {
      matrix.data["class_name",2] <- "Trematoda"
      matrix.data["phylum_name",2] <- "Platyhelminthes"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }


    if (scientific_name == "Alveolata") {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "None"
      matrix.data["kingdom_name",2] <- "Alveolata"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "kingdom"
    }
    if (scientific_name == "Haptophyceae") {
      matrix.data["class_name",2] <- "Haptophyceae"
      matrix.data["phylum_name",2] <- "Haptophyta"
      matrix.data["kingdom_name",2] <- "Hacrobia"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "class"
    }
    if (scientific_name == "Synurophyceae") {
      matrix.data["class_name",2] <- "Chrysophyceae"
      matrix.data["phylum_name",2] <- "Ochrophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "class"
    }
    if (scientific_name == "Labyrinthulomycetes") {
      matrix.data["class_name",2] <- "Labyrinthulea"
      matrix.data["phylum_name",2] <- "Bigyra"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "class"
    }

    if (scientific_name == "Xanthophyceae") {
      matrix.data["class_name",2] <- "Xanthophyceae"
      matrix.data["phylum_name",2] <- "Ochrophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "class"
    }

    if (scientific_name %in% c("Fungi incertae sedis","Dikarya")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "None"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "kingdom"
    }
    if (scientific_name == "Opisthokonta") {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "None"
      matrix.data["kingdom_name",2] <- "None"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "superkingdom"
    }
    if (scientific_name %in% c("Streptophytina","Embryophyta")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "Streptophyta"
      matrix.data["kingdom_name",2] <- "Viridiplantae"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
    }
    if (scientific_name %in% c("Spermatophyta")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "None"
      matrix.data["kingdom_name",2] <- "Viridiplantae"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
    }
    if (scientific_name %in% c("Bilateria","Ecdysozoa","Eumetazoa","Lophotrochozoa","Panarthropoda","Protostomia","Deuterostomia")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "None"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Scolecida","Aciculata","Palpata","Sedentaria")) {
      matrix.data["class_name",2] <- "Polychaeta"
      matrix.data["phylum_name",2] <- "Annelida"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Clitellata")) {
      matrix.data["class_name",2] <- "Clitellata"
      matrix.data["phylum_name",2] <- "Annelida"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Oligochaeta")) {
      matrix.data["class_name",2] <- "Oligochaeta"
      matrix.data["phylum_name",2] <- "Annelida"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name %in% c("Eleutherozoa","Echinozoa","Asterozoa")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "Echinodermata"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Euechinoidea","Echinacea","Gnathostomata")) {
      matrix.data["class_name",2] <- "Echinoidea"
      matrix.data["phylum_name",2] <- "Echinodermata"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Ophiuridea")) {
      matrix.data["class_name",2] <- "Ophiuroidea"
      matrix.data["phylum_name",2] <- "Echinodermata"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name == "Rhodophyta") {
      matrix.data["phylum_name",2] <- "Rhodophyta"
      matrix.data["kingdom_name",2] <- "Rhodophyta"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
      matrix.data["rank",2] <- "phylum"
    }
    if (scientific_name %in% c("Corallinophycidae","Rhodymeniophycidae","Nemaliophycidae")) {
      matrix.data["class_name",2] <- "Florideophyceae"
      matrix.data["phylum_name",2] <- "Rhodophyta"
      matrix.data["kingdom_name",2] <- "Rhodophyta"
      matrix.data["superkingdom_name",2] <- "Archaeplastida"
      matrix.data["rank",2] <- "phylum"
    }
    if (scientific_name == "PX clade") {
      matrix.data["phylum_name",2] <- "Ochrophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "phylum"
    }
    if (scientific_name == "Stramenopiles") {
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "kingdom"
    }
    if (scientific_name %in% c("Euteleostomi","Teleostei","Euteleosteomorpha","Actinopterygii","Clupeocephala","Eupercaria","Percomorphaceae",
                               "Acanthomorphata","Euacanthomorphacea","Ctenosquamata","Euteleosteomorpha")) {
      matrix.data["class_name",2] <- "Actinopterygii"
      matrix.data["phylum_name",2] <- "Chordata"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Batoidea","Galeoidea")) {
      matrix.data["class_name",2] <- "Chondrichthyes"
      matrix.data["phylum_name",2] <- "Chordata"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Boreoeutheria")) {
      matrix.data["class_name",2] <- "Mammalia"
      matrix.data["phylum_name",2] <- "Chordata"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Amniota","Tunicata")) {
      matrix.data["phylum_name",2] <- "Chordata"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name %in% c("Chelicerata","Crustacea","Pancrustacea","Mandibulata")) {
      matrix.data["class_name",2] <- "None"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Eumalacostraca","Peracarida","Eucarida","Hoplocarida")) {
      matrix.data["class_name",2] <- "Malacostraca"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Neocopepoda","Podoplea","Cirripedia","Thoracica","Hexanauplia")) {
      matrix.data["class_name",2] <- "Maxillopoda"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Collembola","Entomobryomorpha","Poduromorpha","Symphypleona","Poduroidea")) {
      matrix.data["order_name",2] <- "Collembola"
      matrix.data["class_name",2] <- "Hexapoda"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "order"
    }
    if (scientific_name %in% c("Hexapoda")) {
      matrix.data["class_name",2] <- "Hexapoda"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "class"
    }
    if (scientific_name %in% c("Acariformes","Acari")) {
      matrix.data["class_name",2] <- "Arachnida"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "order"
    }
    if (scientific_name=="Sarcoptiformes") {
      matrix.data["order_name",2] <- "Sarcoptiformes"
      matrix.data["class_name",2] <- "Arachnida"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "order"
    }
    if (scientific_name %in% c("Trombidiformes","Anystina","Hydracarina","Prostigmata")) {
      matrix.data["order_name",2] <- "Trombidiformes"
      matrix.data["class_name",2] <- "Arachnida"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "order"
    }
    if (scientific_name %in% c("Endopterygota","Neoptera","Pterygota","Paleoptera","Paraneoptera")) {
      matrix.data["class_name",2] <- "Insecta"
      matrix.data["phylum_name",2] <- "Arthropoda"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Stylommatophora","Euthyneura","Nudipleura","Panpulmonata","Heterobranchia",
                               "Euopisthobranchia","Hypsogastropoda","Caenogastropoda")) {
      matrix.data["class_name",2] <- "Gastropoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Littorinimorpha","Littorinoidea")) {
      matrix.data["order_name",2] <- "Littorinimorpha"
      matrix.data["class_name",2] <- "Gastropoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Sacoglossa")) {
      matrix.data["order_name",2] <- "Sacoglossa"
      matrix.data["class_name",2] <- "Gastropoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
      matrix.data["rank",2] <- "order"
    }
    if (scientific_name %in%  c("Patellogastropoda")) {
      matrix.data["order_name",2] <- "Patellogastropoda"
      matrix.data["class_name",2] <- "Gastropoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Vetigastropoda")) {
      matrix.data["order_name",2] <- "Vetigastropoda"
      matrix.data["class_name",2] <- "Gastropoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Neogastropoda","Buccinoidea")) {
      matrix.data["order_name",2] <- "Neogastropoda"
      matrix.data["class_name",2] <- "Gastropoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Neoloricata")) {
      matrix.data["class_name",2] <- "Polyplacophora"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name %in%  c("Cephalaspidea")) {
      matrix.data["order_name",2] <- "Cephalaspidea"
      matrix.data["class_name",2] <- "Gastropoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in%  c("Pteriomorphia","Imparidentia","Heteroconchia","Euheterodonta","Palaeoheterodonta")) {
      matrix.data["class_name",2] <- "Bivalvia"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name %in% c("Decapodiformes")) {
      matrix.data["class_name",2] <- "Cephalopoda"
      matrix.data["phylum_name",2] <- "Mollusca"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name == "Myxogastria") {
      matrix.data["class_name",2] <- "Myxomycetes"
      matrix.data["phylum_name",2] <- "Mycetozoa"
      matrix.data["kingdom_name",2] <- "Conosa"
      matrix.data["superkingdom_name",2] <- "Amoebozoa"
      matrix.data["rank",2] <- "class"
      matrix.data["scientific_name",2] <- "Myxomycetes"
    }

    if (scientific_name == "Fragilariophycidae") {
      matrix.data["class_name",2] <- "Fragilariophyceae"
      matrix.data["phylum_name",2] <- "Bacillariophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "class"
      matrix.data["scientific_name",2] <- "Bacillariophyceae"
    }
    if (scientific_name == "Bacillariophycidae") {
      matrix.data["class_name",2] <- "Bacillariophyceae"
      matrix.data["phylum_name",2] <- "Bacillariophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "class"
      matrix.data["scientific_name",2] <- "Bacillariophyceae"
    }
    if (scientific_name %in% c("Octocorallia","Hexacorallia","Scleractinia")) {
      matrix.data["class_name",2] <- "Anthozoa"
      matrix.data["phylum_name",2] <- "Cnidaria"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Hydroidolina")) {
      matrix.data["class_name",2] <- "Hydrozoa"
      matrix.data["phylum_name",2] <- "Cnidaria"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("Trachylinae")) {
      matrix.data["order_name",2] <- "Trachylina"
      matrix.data["class_name",2] <- "Hydrozoa"
      matrix.data["phylum_name",2] <- "Cnidaria"
      matrix.data["kingdom_name",2] <- "Metazoa"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }

    if (scientific_name %in% c("sordariomyceta","Sordariomycetidae")) {
      matrix.data["class_name",2] <- "Sordariomycetes"
      matrix.data["phylum_name",2] <- "Ascomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("leotiomyceta","Leotiomycetidae","Leotiomycetes incertae sedis")) {
      matrix.data["class_name",2] <- "Leotiomycetes"
      matrix.data["phylum_name",2] <- "Ascomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("eurotiomyceta","Eurotiomycetidae")) {
      matrix.data["class_name",2] <- "Eurotiomycetes"
      matrix.data["phylum_name",2] <- "Ascomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    if (scientific_name %in% c("hypocreomyceta","Hypocreomycetidae")) {
      matrix.data["class_name",2] <- "Sordariomycetes"
      matrix.data["phylum_name",2] <- "Ascomycota"
      matrix.data["kingdom_name",2] <- "Fungi"
      matrix.data["superkingdom_name",2] <- "Opisthokonta"
    }
    return(matrix.data)
  }

  #Read  fasta file
  message("Reading ecotagged fasta file")
  db <- readDNAStringSet(filefasta)
  info <- db@ranges@NAMES
  #Initialize
  inicio <- T
  message("Read ",length(db)," records")
  for (fila in 1:length(db)) {
    infofasta <- info[fila]
    if (length(gregexpr("ORDER_NAME=",infofasta)[[1]])==2) infofasta <- sub(pattern="; ORDER_NAME=",replacement="",x=infofasta)
    if (length(gregexpr("ORDER=",infofasta)[[1]])==2) infofasta <- sub(pattern="ORDER=",replacement="",x=infofasta)
    id <- strsplit(infofasta,split=c(" "))[[1]][1]
    infofasta2 <- substr(infofasta,nchar(id)+2,nchar(infofasta))
    infofasta2 <- (strsplit(sub(pattern=";'",replacement="'",x=infofasta2),";"))
    buenos <- unlist(infofasta2)
    buenos <- buenos[grepl("=",buenos)]
    matrix.data <- data.frame(matrix(ncol=2,nrow=length(buenos)))

    for (i in 1:length(buenos)) matrix.data[i,] <- strsplit(buenos[i],"=")[[1]]
    rownames(matrix.data) <- gsub("^ ", "", matrix.data$X1)
    if (inicio) {
      db.out <- data.frame(names=c("id",matrix.data$X1,"class_name","phylum_name","kingdom_name",
                                   "superkingdom_name","sequence"))
      rownames(db.out) <- c("id",rownames(matrix.data),"class_name","phylum_name","kingdom_name",
                            "superkingdom_name","sequence")
      inicio<-F
    }
    if (sum(grepl("ramosa",matrix.data$X2))!=0){
      submatrix.data <- matrix.data$X2[grep(":",matrix.data$X2)]
      dobles <- grep(":",matrix.data$X2)
      for (j in 1:length(submatrix.data)) submatrix.data[j]  <-
        substr(strsplit(submatrix.data[j],":")[[1]][2],2,nchar(strsplit(submatrix.data[j],":")[[1]][2])-1)
      matrix.data$X2[dobles] <- submatrix.data
    }
    matrix.data["id",2] <- substr(id,1,length_id)
    # matrix.data["seq_length",2] <- db@ranges@width[fila]
    matrix.data["sequence",2] <- as.character(db[fila])

    if (matrix.data["order_name",2]=="None" & matrix.data["family_name",2]=="None" & matrix.data["genus_name",2]!="None"){
      if (matrix.data["genus_name",2] %in% genus_to_family$genus_name){
        matrix.data["family_name",2] <- as.character(genus_to_family$family_name[genus_to_family$genus_name==matrix.data["genus_name",2]])
        matrix.data["order_name",2] <- as.character(family_to_order$order_name[family_to_order$family_name==matrix.data["family_name",2]])} else{
          matrix.data["order_name",2] <- "Correct_manually"
          matrix.data["family_name",2] <- "Correct_manually"
        }
    }


    if (matrix.data["order_name",2]=="None" & matrix.data["family_name",2]!="None"){
      if (matrix.data["family_name",2] %in% family_to_order$family_name){
        matrix.data["order_name",2] <- as.character(family_to_order$order_name[family_to_order$family_name==matrix.data["family_name",2]])} else{
          matrix.data["order_name",2] <- "Correct_manually"
        }
    }
    if (length(as.character(taxo_names[taxo_names$order_name==matrix.data["order_name",2],2]))!=0){
      matrix.data["class_name",2] <- as.character(taxo_names[taxo_names$order_name==matrix.data["order_name",2],2])
      matrix.data["phylum_name",2] <- as.character(taxo_names[taxo_names$order_name==matrix.data["order_name",2],3])
      matrix.data["kingdom_name",2] <- as.character(taxo_names[taxo_names$order_name==matrix.data["order_name",2],4])
      matrix.data["superkingdom_name",2] <- as.character(taxo_names[taxo_names$order_name==matrix.data["order_name",2],5])
    } else {
      matrix.data["class_name",2] <- "Correct_manually"
      matrix.data["phylum_name",2] <- "Correct_manually"
      matrix.data["kingdom_name",2] <- "Correct_manually"
      matrix.data["superkingdom_name",2] <- "Correct_manually"
    }
    if (matrix.data["rank",2]=="phylum" & (length(as.character(phylum_to_sk[phylum_to_sk$phylum_name==matrix.data["scientific_name",2],2]))!=0) ){
      matrix.data["phylum_name",2] <- matrix.data["scientific_name",2]
      matrix.data["kingdom_name",2] <- as.character(phylum_to_sk[phylum_to_sk$phylum_name==matrix.data["phylum_name",2],2])
      matrix.data["superkingdom_name",2] <- as.character(phylum_to_sk[phylum_to_sk$phylum_name==matrix.data["phylum_name",2],3])

    }
    if (matrix.data["rank",2]=="class" & (length(as.character(class_to_sk[class_to_sk$class_name==matrix.data["scientific_name",2],2]))!=0) ){
      matrix.data["class_name",2] <- matrix.data["scientific_name",2]
      matrix.data["phylum_name",2] <- as.character(class_to_sk[class_to_sk$class_name==matrix.data["class_name",2],2])
      matrix.data["kingdom_name",2] <-as.character(class_to_sk[class_to_sk$class_name==matrix.data["class_name",2],3])
      matrix.data["superkingdom_name",2] <-as.character(class_to_sk[class_to_sk$class_name==matrix.data["class_name",2],4])
    }

    if (matrix.data["rank",2]=="kingdom" & (length(as.character(kingdom_to_sk[kingdom_to_sk$kingdom_name==matrix.data["scientific_name",2],2]))!=0)){
      matrix.data["kingdom_name",2] <- matrix.data["scientific_name",2]
      matrix.data["superkingdom_name",2] <-as.character(kingdom_to_sk[kingdom_to_sk$kingdom_name==matrix.data["kingdom_name",2],2])
    }
    if (matrix.data["rank",2]=="phylum" & matrix.data["scientific_name",2]=="Phaeophyceae") {
      matrix.data["class_name",2] <- "Phaeophyceae"
      matrix.data["phylum_name",2] <- "Ochrophyta"
      matrix.data["kingdom_name",2] <- "Stramenopiles"
      matrix.data["superkingdom_name",2] <- "Chromalveolata"
      matrix.data["rank",2] <- "class"
    }
    if (matrix.data["scientific_name",2] %in% exceptions) matrix.data <- fix_exceptions(matrix.data["scientific_name",2])


    db.out <- cbind(db.out,row=matrix.data[match(rownames(db.out),rownames(matrix.data)),2])

    if (fila %% showlines == 0) message(fila,"/",length(db)," lines written.","\r",appendLF = FALSE)

  }
  db.outt <- data.frame(t(db.out[-1]))

  # Substitute commas by "|" in species_list and remove special characters
  db.outt$species_list <- gsub(", ","|",db.outt$species_list)
  db.outt$species_list <- gsub("\\[","",db.outt$species_list)
  db.outt$species_list <- gsub("\\]","",db.outt$species_list)
  db.outt$species_list <- gsub("\'","",db.outt$species_list)
  db.outt$species_list <- gsub("#","",db.outt$species_list)
  db.outt$best_match <- gsub("\'","",db.outt$best_match)
  db.outt$best_match <- gsub("#","",db.outt$best_match)

  db.outt$species_name <- gsub("#","",db.outt$species_name)
  db.outt$genus_name <- gsub("#","",db.outt$genus_name)
  db.outt$scientific_name <- gsub("#","",db.outt$scientific_name)

  db.outt$species_name <- gsub("\'","",db.outt$species_name)
  db.outt$genus_name <- gsub("\'","",db.outt$genus_name)
  db.outt$scientific_name <- gsub("\'","",db.outt$scientific_name)



  db.outt$merged_sample <- NULL
  reordered_db.outt <- db.outt[,c("id","rank","scientific_name","best_identity",
                                  "superkingdom_name","kingdom_name","phylum_name","class_name",
                                  "order_name","family_name","genus_name","species_name",
                                  "best_match","species_list","taxid","sequence")]

  write.table(reordered_db.outt,outfile,row.names=F,col.names=T,sep="\t",quote = FALSE)
  message("THOR is done. He wrote the output file ",outfile," with ",length(db), " assigned sequences.")
}

