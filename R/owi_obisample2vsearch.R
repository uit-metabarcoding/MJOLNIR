## Script for formatting a obitools fasta file into a vsearch fasta file
## The script will read the input file name, the output file name, and the name of the attribute containing abundances
## The name of the sample will be kept
## If the output file argument is empty, it will just add ".vsearch.fasta" at the end of the name of the input file.
## By Owen S. Wangensteen - Project Metabarpark  2016

owi_obisample2vsearch <- function(infile,outfile=NULL,attrib=count){

    # if outfile is not provided, calculate it from the infile
    if (is.null(outfile)) outfile <- paste(substr(infile,1,nchar(infile)-6),".vsearch.fasta",sep="")

    suppressPackageStartupMessages(library(Biostrings))

        # Function to get abundances for each obifasta header
    abundance <- function(texto){
      comienzo <- as.vector((regexpr(paste(" ",attrib,"=",sep=""),texto)))+nchar(attrib)+2
      medio <- substr(texto, comienzo,nchar(texto))
      fin <- as.vector(regexpr(";",medio))
      return(substring(medio,1,fin-1))
    }

    # Function to get sample names for each obifasta header
    sample_name <- function(texto){
      comienzo <- as.vector((regexpr(" sample=",texto)))+nchar(attrib)+2
      medio <- substr(texto, comienzo,nchar(texto))
      fin <- as.vector(regexpr(";",medio))
      return(substring(medio,1,fin-1))
    }

    # main loop
    db <- readDNAStringSet(infile)
    db@ranges@NAMES <- paste0(as.character(lapply(strsplit(db@ranges@NAMES," "),function(l) l[[1]])),";sample",sample_name(db@ranges@NAMES),";size=",abundance(db@ranges@NAMES),";")
    writeXStringSet(db,outfile)
    message("File ",outfile," written")
}
