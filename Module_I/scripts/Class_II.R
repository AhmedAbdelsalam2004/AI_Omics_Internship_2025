classify <- function(logFC, padj){
  
  if(logFC > 1 & padj < 0.05){
    return("Upregulated")
  }
  
  else if(logFC < -1 & padj < 0.05){
    return("Downregulated")
  }
  
  else {
    return("Not_Significant")
  }
  
}


Process <- function(file){
  
  data <- read.csv(file) # Load CSV file
  
  data <- na.omit(data) # Remove missing values
  
  data$classification <- mapply(classify, data$logFC, data$padj) # Classification 
  
  output <- paste0("clean_data/", sub(".csv", "_processed.csv", basename(file))) # Save processed data to clean_data
  write.csv(data, output, row.names = FALSE)
  
  
  print(table(data$classification)) # summary
  
}

Process("raw_data/DEGs_Data_1.csv")
Process("raw_data/DEGs_Data_2.csv")


save.image(file = "AhmedAbdelsalam_Class_2_Assignment.RData") # Save the entire R workspace


