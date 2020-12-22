# MJOLNIR
![MJOLNIR](https://github.com/uit-metabarcoding/MJOLNIR/blob/main/mjolnir_blue_mid.png)

<b>MJOLNIR: Metabarcoding Joining Obitools &amp; Linkage Networks In R</b>

MJOLNIR is an R package to run modular metabarcoding pipelines from the R environment. MJOLNIR runs on Linux and Mac systems. It is not clear to me if it may run in Windows 10 using the Ubuntu Linux subsystem. In any case, the extensive use of package parallel and several dependencies that are designed primarily for Linux systems (see below) makes the success of Windows installations highly improbable. Users are welcome to try to install and run MJOLNIR on Windows Linux Subsystem, but I would not recommended that.

MJOLNIR depends on the following dependencies, which must be installed in the system and properly working:

- OBITools (Boyer et al. 2016):
  Original information about OBITools here: https://git.metabarcoding.org/obitools/obitools/wikis/home
  Help on installing OBITools: http://rleca.pbworks.com/w/file/fetch/124098201/tuto_obitools_install_W10OS.html
  If this does not work, try: "sudo apt-get install obitools".
  Note that OBITools currently runs on Python 2.7. It is not working in Python 3. So Python 2.7 is required for the instalation.

- VSEARCH (Rognes et al. 2016): 
  Help on installing VSEARCH: https://github.com/torognes/vsearch
  
- SWARM (Mahé et al. 2015):
  Help on installing SWARM: https://github.com/torognes/swarm
  
- LULU (Frøslev et al. 2017):
  Help on installing LULU:
  https://github.com/tobiasgf/lulu

Installing MJOLNIR:

1. If the devtools library is properly installed in the system: MJOLNIR can be installed directly from the R console using:
library(devtools)
install_github("uit-metabarcoding/MJOLNIR")

2. If the devtools library is not installed, then MJOLNIR can be downloaded as a package from this link: https://github.com/uit-metabarcoding/MJOLNIR/archive/main.zip
Then the file must be unzipped and MJOLNIR can be installed offline from the R console using:
install.packages("MJOLNIR-main", repos=NULL)

MJOLNIR is currently optimized by default to process COI metabarcoding data (the Leray Fragment). For other metabarcoding markers, some settings must be changed as follows:

The following settings are recommended for COI Leray/Leray-XT primers (Leray et al. 2013; Wangensteen et al. 2018): 
- In mjolnir2_FREYJA: Lmin=299,Lmax=320 
- In mjolnir4_ODIN: d=13

The following settings are recommended for 12S MiFish primers (Miya et al. 2015): 
- In mjolnir2_FREYJA: Lmin=140,Lmax=190 
- In mjolnir4_ODIN: d=1

The following settings are recommended for 18S All_shorts primers (Guardiola et al. 2015): 
- In mjolnir2_FREYJA: Lmin=75,Lmax=180 
- In mjolnir4_ODIN: d=1

The following settings are recommended for 16S Bacterial F515/R806 primers (Caporaso et al. 2011): 
- In mjolnir2_FREYJA: Lmin=215,Lmax=299 
- In mjolnir4_ODIN: d=1

The MJOLNIR Pipeline

0. Input data

MJOLNIR is optimized to process paired-end FASTQ files from Illumina sequencing of multiplexed libraries prepared using the METAFAST procedure. This procedure adds sample-tags on both ends of the amplicons, at 5' from the metabarcoding primers. Several samples (usually 96 or more) can be multiplexed into a single metafast library. Each pair of fastq files belonging to a library usually contain multiple samples, identified by unique combination of forward:reverse sample-tags. MJOLNIR can process several such libraries simultaneously, spanning hundreds or thousands of samples, which will be joined together into a final combined dataset before the clustering step.

1. RAN: Reads Allotment in N portions 

MJOLNIR starts with a call to mjolnir1_RAN() which will split initial FASTQ files in aliquot parts for parallel processing. This will distribute the initial paired-end fastq files in fragments of equal number of reads, to be processed in parallel. 

2. FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment 

The function mjolnir2_FREYJA() calls OBITools for the following three steps: (1) paired-end alignment, (2) demultiplexing and (3) read-length filtering, wrapped together in a single function. Paired-end alignment is done using illuminapairedend. Demultiplexing and primer-removal are done using ngsfilter. Removal of low-quality reads is based on the individual read quality results from both steps. A length filter using obigrep follows, based on the length of the metabarcoding fragment (excluding the primers). All reads having bases different from [ACGT] are also removed. The ngsfilter function needs an input table (ngsfilter-table), containing information about the internal sample names (agnomens), the sample-tags used to identify them at both ends, and the sequences of the metabarcoding primers used. ngsfilter-tables must be provided for each library.

3. HELA: Hierarchical Elimination of Lurking Artifacts

MJOLNIR uses the mjolnir3_HELA() function to call the uchime_denovo algorithm implemented in VSEARCH, to remove chimaeric sequences from individual sample files, in a sample-by-sample basis. This procedure is much faster than removing the chimaeras from the whole dataset together. After removing the chimaeras from every sample, HELA joins all the samples again, and de-replicates the sequences. Then creates new sequence identifiers for all different sequence variants found, and generates two outputs: (1) a fasta file with total abundances (ending in .vsearch.fasta) ready to be fed into Swarm for the clustering step, and (2) a table including the abundance information of all sequence variants in every sample. 

4. ODIN: OTU Delimitation Inferred by Networks

The function mjolnir4_ODIN() uses the swarm algorithm to delimit MOTUs, based on linkage-networks created by step-by-step agregation. This clustering algorithm is not based on an arbitrary, constant, absolute identity threshold. Conversely, swarm is based on an iterative aggregation of sequences that differ less than a given distance d. This strategy results into linkage-networks of different sizes, which can have different effective values for the within-MOTU identity threshold, depending on the complexity of the natural variability of the sequences present in the sample. This procedure is very convenient in the case of hypervariable metabarcoding markers, such as COI, which usually feature extremely high levels of natural diversity, which is added to the random sequencing errors. Moreover, the resulting networks are very convenient for further processing of metaphylogeography datasets (Turon et al. 2019). mjolnir4_ODIN() takes the .vsearch.fasta file produced by HELA, as an input. Then it used the table with the read abundances by sample to compute the combined abundance of the resulting MOTUs in every individual sample. The output is a CSV file with contains a row for each MOTU, with its abundances in every sample. Also, a fasta file including the representative sequences of each MOTU (usually only non-singleton MOTUs are included) is generated, which will be used for taxonomic assignment in the next step.

5. THOR: Taxonomy with Higher-than-Order Ranks

Taxonomic assignment is performed with the mjolnir5_THOR() function, which is a wrapper of ecotag (Boyer et al. 2016). This step depends on the availability of a taxonomic database in ecoPCR format (from which the phylogenetic tree information is retrieved) and a reference database file including sequences for the exclusively metabarcoded fragment, conveniently identified by a taxonomic identifier (taxid), in fasta format. These can be downloaded from our DUFA repository: https://github.com/uit-metabarcoding/DUFA. mjolnir5_THOR then uses a custom procedure to assign higher taxonomic ranks (higher than order). This taxonomic ranks are stored in a csv file called "order_complete-csv", which is customizable to meet the particular preferences of the user. This is specially useful for assigning the higher ranks of microeukaryotes, where there is not a universally-accepted consensus. After taxonomy assignment, THOR joins the taxonomic information with the abundances obtained by ODIN, and creates an output CSV file ending in .All_MOTUs.csv 

6. LOKI: LULU Overseeing with Kinship Identification 

The next step of MJOLNIR is the removal of pseudogenes using the mjolnir6_LOKI() function. This is a wrapper of the LULU algorithm. The information about the removed putative pseudogene MOTUs is also stored as an additional output file, together with the information of their possible mother sequences. This file can be checked to assess the taxonomic coherence of the results. A .curated.csv file is created.

7. RAGNAROC: Replace AGnomens with Names And Recover Original Codification

mjolnir7_RAGNAROC() is the last step of the MJOLNIR pipeline. It requires a table of original sample names with the internal agnomes used along the rest of the pipeline. RAGNAROC also reorders the sample columns by alphabetical order and produces the .final_dataset.csv file.
