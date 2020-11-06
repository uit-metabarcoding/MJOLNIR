# MJOLNIR
MJOLNIR Metabarcoding Joining Obitools &amp; Linkage Networks In R

MJOLNIR is an R package to do metabarcoding pipelines from R. MJOLNIR depends on the following dependences, which must be installed and properly working:

- OBITools: 
  Help on installing OBITools: http://rleca.pbworks.com/w/file/fetch/124098201/tuto_obitools_install_W10OS.html
  If this does not work, try: "sudo apt-get install obitools"

- VSEARCH: 
  Help on installing VSEARCH: https://github.com/torognes/vsearch
  
- SWARM 2.0:
  Help on installing SWARM: https://github.com/torognes/swarm
  
-LULU:
  Help on installing LULU:
  https://github.com/tobiasgf/lulu

MJOLNIR is currently optimized to process COI metabarcoding data (the Leray Fragment). Swarm is performed with d=13.

Taxonomic assignment is performed with ecotag and it depends on the availability of a taxonomic database in ecoPCR format and a reference database file in fasta format.
These can be downloaded from our DUFA repository: https://github.com/uit-metabarcoding/DUFA.
