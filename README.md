# MJOLNIR
MJOLNIR Metabarcoding Joining Obitools &amp; Linkage Networks In R

MJOLNIR is an R package to run modular metabarcoding pipelines from the R environment. MJOLNIR runs on Linux and Mac systems. It is not clear to me if it can work in Windows 10 using the Ubuntu linux subsystem. In any case, the extensive use of package parallel and several dependencies that are designed primarily for Linux systems (see below) makes the success of Windows installations highly improbable. Users are welcome to try to install and run MJOLNIR on Windows systems, but it is not recommended at all.

MJOLNIR depends on the following dependences, which must be installed and properly working:

- OBITools (Boyer et al. 2016): 
  Help on installing OBITools: http://rleca.pbworks.com/w/file/fetch/124098201/tuto_obitools_install_W10OS.html
  If this does not work, try: "sudo apt-get install obitools"

- VSEARCH (Rognes et al. 2016): 
  Help on installing VSEARCH: https://github.com/torognes/vsearch
  
- SWARM 2.0 (Mahé et al. 2015):
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

MOTU delimitation is performed using SWARM 2.0. This iterative algorithm will result in a linkage network for every MOTU, clustering all sequences with a Hamming distance less than d. This clustering is very convenient for further processing of metaphylogeography datasets (Turon et al. 2019). 

Taxonomic assignment is performed with ecotag (Boyer et al. 2016) and it depends on the availability of a taxonomic database in ecoPCR format and a reference database file in fasta format. These can be downloaded from our DUFA repository: https://github.com/uit-metabarcoding/DUFA.
