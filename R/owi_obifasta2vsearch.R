## Script for formatting an obitools fasta file into a vsearch fasta file
## The script will read the input file name, the output file name, and the name of the attribute containing abundances
## If the output file argument is empty, it will just add ".vsearch.fasta" at the end of the name of the input file.
## By Owen S. Wangensteen - Project Metabarpark  2016

owi_obifasta2vsearch <- function(infile,outfile=NULL,attrib="count"){
  # if outfile is not provided, calculate it from the infile
  if (is.null(outfile)) outfile <- paste(substr(infile,1,nchar(infile)-6),"_vsearch.fasta",sep="")

    suppressPackageStartupMessages(library(Biostrings))

  # Function to get the abundances for each sequence in the fasta
  abundance <- function(texto) {
    comienzo <- as.vector((regexpr(paste(" ",attrib,"=",sep=""),texto)))+nchar(attrib)+2
    medio <- substr(texto, comienzo,nchar(texto))
    fin <- as.vector(regexpr(";",medio))
    return(substring(medio,1,fin-1))
  }

    # main loop
  db <- readDNAStringSet(infile)
  db@ranges@NAMES <- paste(substr(db@ranges@NAMES,1,14),";size=",abundance(db@ranges@NAMES),";" ,sep="")
  writeXStringSet(db,outfile)
  message("File ",outfile," written")
}
