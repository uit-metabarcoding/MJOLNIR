# RAGNAROC: Replace AGnomens with Names And Recover Original Codification
# Sort_MOTUs options available: "id" (default), "taxonomy", "abundance"
# To remove bacterial MOTUs use remove_bacteria = T
# To remove contaminant MOTUs use remove_contamination = T. 
# When remove_contamination = T a list of contaminants must be provided in contamination_file (a text file with one line for each contaminant scientific_name)


mjolnir8_RAGNAROC <- function(lib,metadata_table="",output_file="",min_reads=2,min_relative=1/50000,sort_MOTUs="id",remove_bacteria=F,remove_contamination=F,contamination_file="contaminants.txt"){
    message("RAGNAROC is coming. Original sample names will be recovered.")
    if (output_file == "") output_file <- paste0(lib,"_final_dataset.tsv")
    #Load the dataset
        # If LULU file exists, load it, otherwise load the All_MOTUs file
        file_to_load <- ifelse(file.exists(paste0(lib,".Curated_LULU.tsv")),paste0(lib,".Curated_LULU.tsv"),paste0(lib,".All_MOTUs.tsv"))
        db <- read.table(file_to_load,sep="\t",head=T,stringsAsFactors = F)
    # Select sample abundance columns
        sample_cols <- (1:ncol(db))[tolower(substr(names(db),6,11))=="sample"]
        no_sample_cols <- (1:ncol(db))[tolower(substr(names(db),6,11))!="sample"]
        sample_names <- names(db[sample_cols])
    # Load the metadata_table
        if (metadata_table=="") metadata_table <- paste0(lib,"_metadata.tsv")
        sample_db <- read.table(metadata_table,sep="\t",head=T,stringsAsFactors = F)
        new_sample_names <- rep("",length(sample_names))
    # Change agnomens by original names
        agnomens_true <- sample_names %in% sample_db$mjolnir_agnomens
        for (i in 1:length(sample_names)) if (agnomens_true[i]) new_sample_names[i] <- sample_db$original_samples[sample_db$mjolnir_agnomens==sample_names[i]]
        empty_samples <- new_sample_names==""
        if (sum(empty_samples)>0) for (i in 1:sum(empty_samples)) new_sample_names[empty_samples][i] <- paste0("EMPTY",i)
        for (i in 1:length(sample_cols)) names(db)[names(db)==sample_names[i]] <- new_sample_names[i]
    # Apply relative abundance filter per sample
        message("RAGNAROC is applying a relative abundance filter. MOTUs with less than ",min_relative," relative abundance will be removed from each sample.")
        for (col in sample_cols) {
            tot_reads <- sum(db[,col])
            min_reads_samp <- ceiling(min_relative*tot_reads)
            db[,col][db[,col] < min_reads_samp] <- 0
        }
    # Reorder sample columns
        db_sample <- db[,sample_cols]
        db_sample_ordered <- db_sample[,sort(names(db_sample))]
        new_total_reads <- rowSums(db_sample_ordered[,substr(names(db_sample_ordered),1,5)!="EMPTY"])
        db_new <- cbind(db[,no_sample_cols[-length(no_sample_cols)]],db_sample_ordered,sequence=db[,length(db)])
    # Remove empty and duplicated columns
        db_new <- db_new[,substr(names(db_new),1,5)!="EMPTY" & substr(names(db_new),(nchar(names(db_new))-1),nchar(names(db_new)))!=".1"]
    # Replace total_reads
        db_new$total_reads <- new_total_reads
    # Remove MOTUs with new total_reads less than or equal than min_reads
        message("RAGNAROC is removing MOTUs with less than ",min_reads," total reads.")
        db_new <- db_new[db_new$total_reads >= min_reads,]
    # Remove bacteria
        if (remove_bacteria) {
            message("RAGNAROC is removing bacterial MOTUs now.")
            db_new <- db_new[(db_new$superkingdom_name != "Prokaryota" & db_new$scientific_name != "root"),]
            }
    # Remove contamination
        if (remove_contamination){
        message("RAGNAROC is removing contaminant MOTUs now.")
        contamination <- readLines(contamination_file)
              db_new <- db_new[!((db_new$scientific_name %in% contamination) |
                      (db_new$phylum_name %in% contamination) |
                      (db_new$class_name %in% contamination) |
                      (db_new$order_name %in% contamination) |
                      (db_new$family_name %in% contamination) |
                      (db_new$genus_name %in% contamination)) ,]       
        }
    # Order by taxonomy
        if (sort_MOTUs == "taxonomy"){
        message("RAGNAROC is ordering MOTUs by taxonomy.")
        db_sort <- db_new[,substr(names(db_new),nchar(names(db_new))-4,nchar(names(db_new)))=="_name"]
        db_sort[db_sort==""] <- "ZZZZ"
        db_sort$abundance <- sum(db_new$total_reads) - db_new$total_reads
        sort_vector <- paste(db_sort$superkingdom_name,db_sort$kingdom_name,db_sort$phylum_name,db_sort$class_name,db_sort$order_name,
                     db_sort$family_name,db_sort$genus_name,db_sort$species_name,db_sort$scientific_name,db_sort$abundance)
        db_new <- db_new[order(sort_vector),]
    }
    # Order by abundance
        if (sort_MOTUs == "abundance") {
            message("RAGNAROC is ordering MOTUs by abundance.")
            db_new <- db_new[order(db_new$total_reads,decreasing = T),]
    }
    # Write final table
        if (output_file=="") output_file <- paste0(lib,".final_dataset.tsv") 
        write.table(db_new,output_file,row.names = F,sep="\t",quote = F)
        message("After RAGNAROC, MJOLNIR is done. File: ",output_file, " written with ",nrow(db_new), " MOTUs and ",sum(db_new$total_reads)," total reads.")
}
