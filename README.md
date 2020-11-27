# MJOLNIR
MJOLNIR Metabarcoding Joining Obitools &amp; Linkage Networks In R

MJOLNIR is an R package to run modular metabarcoding pipelines from the R environment. MJOLNIR runs on Linux and Mac systems. It is not clear to me if it may run in Windows 10 using the Ubuntu Linux subsystem. In any case, the extensive use of package parallel and several dependencies that are designed primarily for Linux systems (see below) makes the success of Windows installations highly improbable. Users are welcome to try to install and run MJOLNIR on Windows Linux Subsystem, but I do not recommended that.

MJOLNIR depends on the following dependencies, which must be installed and properly working:

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

MJOLNIR is currently optimized by default to process COI metabarcoding data (the Leray Fragment). For other metabarcoding markers, some settings must be changed as follows:

The following settings are recommended for COI Leray/Leray-XT primers (Leray et al. 2013; Wangensteen et al. 2018): 
- In mjolnir_demulti_filter: Lmin=299,Lmax=320 
- In mjolnir_swarm: d=13

The following settings are recommended for 12S MiFish primers (Miya et al. 2015): 
- In mjolnir_demulti_filter: Lmin=140,Lmax=190 
- In mjolnir_swarm: d=1

The following settings are recommended for 18S All_shorts primers (Guardiola et al. 2015): 
- In mjolnir_demulti_filter: Lmin=75,Lmax=180 
- In mjolnir_swarm: d=1

The following settings are recommended for 16S Bacterial F515/R806 primers (Caporaso et al. 2011): 
- In mjolnir_demulti_filter: Lmin=215,Lmax=299 
- In mjolnir_swarm: d=1

The MJOLNIR Pipeline

0. Input data

MJOLNIR is optimized to process paired-end FASTQ files from Illumina sequencing of multiplexed libraries prepared using the METAFAST procedure. This procedure adds sample-tags on both ends of the amplicons, at 5' from the metabarcoding primers. Several samples (usually 96 or more) are multiplexed into a single metafast library. Each pair of fastq files belongs to a library (usually with several samples, identified by unique combination of forward:reverse sample-tags). MJOLNIR can process several such libraries simultaneously, spanning hundreds or thousands of samples, which will be joined together into a final combined dataset before the clustering step.

1. Splitting of FASTQ files for parallel processing

MJOLNIR starts with a call to mjolnir_distribute(). This will distribute the initial paired-end fastq files in fragments of equal number of reads, to be processed in parallel. 

2. Paired-end alignment, demultiplexing and quality control

MJOLNIR uses OBITools for these three steps, wrapped together in a single function called mjolnir_demulti_filter(). Paired-end alignment is done using illuminapairedend. Demultiplexing and primer-removal are done using ngsfilter and removal of low-quality reads is based on the individual read quality results from both steps. A length filter using obigrep follows, based on the length of the metabarcoding fragment (excluding the primers). All reads having bases different from [ACGT] are also removed.  ngsfilter needs an input table containing information on the sample-tags at both ends, the names of the samples and the sequences of the metabarcoding primers used. Ngsfilter tables must be provided for each library (fastq file).

3. Sample processing and chimaera removal

MJOLNIR uses the uchime_denovo algorithm implemented in VSEARCH to remove chimaeric sequences from individual sample files, in a sample-by-sample basis. This procedure is much faster than removing the chimaeras from the whole dataset, and it is performed using the mjolnir_process_samples() function.

4. Clustering into linkage-network MOTUs

MJOLNIR uses the SWARM 2.0 algorithm to delimit MOTUs, wrapped into the mjolnir_swarm() function. This clustering algorithm is not based on an arbitrary, constant, absolute identity threshold. Conversely, Swarm is based on an iterative step-by-step aggregation of sequences that differ less than a given Hamming distance d. This strategy results into linkage-networks of different sizes, which can have different effective values for the within-MOTU identity threshold, depending on the complexity of the natural variability of the sequences present in the sample. This procedure is very convenient in the case of hypervariable metabarcoding markers, such as COI, which usually feature extremely high levels of natural diversity, which is added to the random sequencing errors. Moreover, the resulting networks are very convenient for further processing of metaphylogeography datasets (Turon et al. 2019). All samples are combined together into a single dataset to be clustered by mjolnir_swarm(), and then the read abundances of the resulting MOTUs in each individual sample are computed.

5. Taxonomic assignment

Taxonomic assignment is performed with the mjolnir_ecotag() function, which is a wrapper of ecotag (Boyer et al. 2016). This step depends on the availability of a taxonomic database in ecoPCR format (from which the phylogenetic tree information is retrieved) and a reference database file including sequences for the exclusively metabarcoded fragment, conveniently identified by a taxonomic identifier (taxid), in fasta format. These can be downloaded from our DUFA repository: https://github.com/uit-metabarcoding/DUFA.

6. Removal of pseudogenes

The last step of MJOLNIR is the removal of pseudogenes using the mjolnir_lulu() function. This is a wrapper of the LULU algorithm. The information about the removed putative pseudogene MOTUs is provided in an output file, together with the information of their possible mother sequences. This file can be checked to assess the taxonomic coherence of the results.
