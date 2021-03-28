# RAGNAROC: Replace AGnomens with Names And Recover Original Codification


mjolnir8_RAGNAROC <- function(lib,sample_table,output_file){
    message("RAGNAROC is coming. Original sample names will be recovered.")
    #Load the dataset
        # If LULU file exists, load it, otherwise load the All_MOTUs file
        file_to_load <- ifelse(file.exists(paste0(lib,".Curated_LULU.csv")),paste0(lib,".Curated_LULU.csv"),paste0(lib,".All_MOTUs.csv"))
        db <- read.csv(file_to_load,sep=";",head=T,stringsAsFactors = F)
    # Select sample abundance columns
        sample_cols <- (1:ncol(db))[tolower(substr(names(db),6,11))=="sample"]
        no_sample_cols <- (1:ncol(db))[tolower(substr(names(db),6,11))!="sample"]
        sample_names <- names(db[sample_cols])
    # Load the sample_table
        sample_db <- read.csv(sample_table,sep=";",head=T,stringsAsFactors = F)
        new_sample_names <- rep("",length(sample_names))
    # Change agnomens by original names
        agnomens_true <- sample_names %in% sample_db$mjolnir_agnomens
        for (i in 1:length(sample_names)) if (agnomens_true[i]) new_sample_names[i] <- sample_db$original_samples[sample_db$mjolnir_agnomens==sample_names[i]]
        empty_samples <- new_sample_names==""
        if (sum(empty_samples)>0) for (i in 1:sum(empty_samples)) new_sample_names[empty_samples][i] <- paste0("EMPTY",i)
        for (i in 1:length(sample_cols)) names(db)[names(db)==sample_names[i]] <- new_sample_names[i]
    # Reorder sample columns
        db_sample_ordered <- db[,sort(names(db[,sample_cols]))]
        new_total_reads <- rowSums(db_sample_ordered[,substr(names(db_sample_ordered),1,5)!="EMPTY"])
        db_new <- cbind(db[,no_sample_cols[-length(no_sample_cols)]],db_sample_ordered,sequence=db[,length(db)])
    # Remove empty and duplicated columns
        db_new <- db_new[,substr(names(db_new),1,5)!="EMPTY" & substr(names(db_new),(nchar(names(db_new))-1),nchar(names(db_new)))!=".1"]
    # Replace total_reads
        db_new$total_reads <- new_total_reads
    # Remove MOTUs with new total-reads==0
        db_new <- db_new[db_new$total_reads>0,]
    # Write final table
        if (output_file=="") output_file <- paste0(lib,".final_dataset.csv") 
        write.table(db_new,output_file,row.names = F,sep=";",quote = F)
        message("After RAGNAROC, MJOLNIR is done. File: ",output_file, " written with ",nrow(db_new), " MOTUs and ",sum(db_new$total_reads)," total reads.")
}
