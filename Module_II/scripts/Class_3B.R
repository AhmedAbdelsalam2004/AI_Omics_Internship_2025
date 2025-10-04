# Microarray Preprocessing and Quality Control Pipeline
# GSE79973 Analysis - AI & Omics Internship 2025 - Module II
# Updated version with corrected outlier detection

# ======================
# INITIALIZATION SECTION
# ======================

# Package Management Function
setup_environment <- function() {
  # Bioconductor packages
  bioc_packages <- c("GEOquery", "affy", "arrayQualityMetrics", "genefilter", "AnnotationDbi")
  
  # CRAN packages  
  cran_packages <- c("dplyr", "matrixStats", "ggplot2", "reshape2")
  
  # Install and load Bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing Bioconductor package: ", pkg)
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Install and load CRAN packages
  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing CRAN package: ", pkg)
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  options(stringsAsFactors = FALSE)
  message("All packages loaded successfully")
}

# Initialize environment
setup_environment()

# ======================
# CONFIGURATION SECTION
# ======================

# Analysis parameters
CONFIG <- list(
  gse_accession = "GSE79973",
  results_dir = "Analysis_Results",
  raw_data_dir = "Raw_Data_Files",
  intensity_threshold = 3.5,
  variance_cutoff = 0.5,
  outlier_sd_threshold = 2
)

# Create directory structure
create_directories <- function() {
  dirs_to_create <- c(CONFIG$results_dir, CONFIG$raw_data_dir)
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      message("Created directory: ", dir)
    }
  }
}

create_directories()

# ======================
# DATA LOADING SECTION
# ======================

load_geo_data <- function() {
  message("Downloading GEO dataset: ", CONFIG$gse_accession)
  
  geo_data <- getGEO(CONFIG$gse_accession, GSEMatrix = TRUE)
  
  if (length(geo_data) == 0) {
    stop("Failed to download GEO data for: ", CONFIG$gse_accession)
  }
  
  # Use first dataset if multiple are available
  expression_set <- geo_data[[1]]
  if (length(geo_data) > 1) {
    message("Multiple datasets found, using the first one")
  }
  
  return(expression_set)
}

# Load the main dataset
expression_set <- load_geo_data()

# Extract data components
expression_matrix <- exprs(expression_set)
feature_info <- fData(expression_set)
sample_info <- pData(expression_set)

# Check data quality
message("Dataset dimensions: ", nrow(expression_matrix), " probes × ", 
        ncol(expression_matrix), " samples")
message("Missing values in sample source: ", 
        sum(is.na(sample_info$source_name_ch1)))

# ======================
# GET EXACT COUNTS
# ======================

# Total samples
total_samples <- ncol(expression_matrix)
message("=== ANSWER: Total samples: ", total_samples)

# Disease vs Normal counts
if (!is.null(sample_info$source_name_ch1)) {
  sample_types <- table(trimws(as.character(sample_info$source_name_ch1)))
  disease_count <- sample_types["gastric adenocarcinoma"]
  normal_count <- sample_types["gastric mucosa"]
  
  message("=== ANSWER: Disease samples: ", disease_count)
  message("=== ANSWER: Normal samples: ", normal_count)
} else {
  message("=== ANSWER: Sample type information not available in source_name_ch1")
}

# Probes before filtering
probes_before <- nrow(expression_matrix)
message("=== ANSWER: Probes before filtering: ", probes_before)

# ======================
# RAW DATA PROCESSING
# ======================

process_raw_data <- function() {
  # Look for CEL files
  cel_files <- list.files(
    path = CONFIG$raw_data_dir,
    pattern = "\\.CEL(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(cel_files) > 0) {
    message("Found ", length(cel_files), " CEL files")
    cel_directory <- dirname(cel_files[1])
    
    # Read Affymetrix data
    raw_affy_data <- ReadAffy(celfile.path = cel_directory)
    message("Raw data loaded: ", ncol(raw_affy_data), " arrays")
    
    # Quality control on raw data
    perform_quality_control(raw_affy_data, "Raw_Data_QC")
    
    # Normalize using RMA
    message("Performing RMA normalization...")
    normalized_set <- rma(raw_affy_data)
    
    return(normalized_set)
  } else {
    message("No CEL files found - using preprocessed expression set")
    return(expression_set)
  }
}

# Process raw data if available
normalized_expression_set <- process_raw_data()
normalized_matrix <- exprs(normalized_expression_set)

# Save normalized data
write.csv(normalized_matrix, 
          file.path(CONFIG$results_dir, "normalized_expression_data.csv"))
message("Normalized data saved")

# ======================
# QUALITY CONTROL
# ======================

perform_quality_control <- function(expression_data, qc_name) {
  qc_directory <- file.path(CONFIG$results_dir, qc_name)
  dir.create(qc_directory, recursive = TRUE, showWarnings = FALSE)
  
  message("Running quality control: ", qc_name)
  arrayQualityMetrics(
    expressionset = expression_data,
    outdir = qc_directory,
    force = TRUE
  )
}

# QC on normalized data
perform_quality_control(normalized_expression_set, "Normalized_Data_QC")

# ======================
# DATA FILTERING
# ======================

filter_low_intensity <- function(expression_data, threshold) {
  probe_medians <- rowMedians(expression_data, na.rm = TRUE)
  
  # Create intensity histogram
  intensity_plot <- ggplot(data.frame(MedianIntensity = probe_medians), 
                           aes(x = MedianIntensity)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    geom_vline(xintercept = threshold, color = "red", linewidth = 1) +
    labs(title = "Probe Intensity Distribution",
         x = "Median Intensity",
         y = "Frequency") +
    theme_minimal()
  
  ggsave(file.path(CONFIG$results_dir, "intensity_distribution.png"), 
         intensity_plot, width = 8, height = 6)
  
  # Apply filter
  keep_probes <- probe_medians > threshold
  filtered_data <- expression_data[keep_probes, ]
  
  message("Intensity filtering: ", sum(keep_probes), "/", nrow(expression_data), 
          " probes retained")
  
  return(filtered_data)
}

filter_by_variance <- function(expression_data, cutoff = 0.5) {
  if (requireNamespace("genefilter", quietly = TRUE)) {
    variance_filtered <- genefilter::varFilter(
      expression_data, 
      var.cutoff = cutoff, 
      filterByQuantile = TRUE
    )
    message("Variance filtering: ", nrow(variance_filtered), " probes retained")
    return(variance_filtered)
  } else {
    message("genefilter package not available - skipping variance filtering")
    return(expression_data)
  }
}

# Apply filtering pipeline
intensity_filtered <- filter_low_intensity(normalized_matrix, 
                                           CONFIG$intensity_threshold)

# ANSWER 1: Transcripts after intensity filtering
message("=== ANSWER: Transcripts remaining after intensity filtering: ", nrow(intensity_filtered))

final_filtered_data <- filter_by_variance(intensity_filtered, 
                                          CONFIG$variance_cutoff)

# Save filtered data
write.csv(final_filtered_data, 
          file.path(CONFIG$results_dir, "filtered_expression_matrix.csv"))

# ======================
# SAMPLE ANNOTATION
# ======================

setup_sample_groups <- function(sample_metadata) {
  if (!is.null(sample_metadata$source_name_ch1)) {
    sample_sources <- trimws(as.character(sample_metadata$source_name_ch1))
    
    # Create factor for analysis groups
    analysis_groups <- factor(
      sample_sources,
      levels = c("gastric mucosa", "gastric adenocarcinoma"),
      labels = c("Normal", "Cancer")
    )
    
    message("Sample groups created: ", 
            paste(levels(analysis_groups), collapse = " vs "))
    
    return(analysis_groups)
  } else {
    warning("Sample source information not available")
    return(NULL)
  }
}

sample_groups <- setup_sample_groups(sample_info)

# ======================
# OUTLIER DETECTION
# ======================

detect_sample_outliers <- function(expression_data) {
  sample_distances <- as.matrix(dist(t(expression_data)))
  mean_distances <- rowMeans(sample_distances)
  
  outlier_threshold <- mean(mean_distances) + 
    CONFIG$outlier_sd_threshold * sd(mean_distances)
  
  outlier_samples <- names(mean_distances)[mean_distances > outlier_threshold]
  
  return(outlier_samples)
}

outlier_samples <- detect_sample_outliers(final_filtered_data)

# ANSWER 2: Outlier detection results
if (length(outlier_samples) > 0) {
  message("=== ANSWER: Outlier arrays detected after normalization: Yes")
  message("=== ANSWER: Number of outliers flagged: ", length(outlier_samples))
  message("=== ANSWER: Outlier samples: ", paste(outlier_samples, collapse = ", "))
} else {
  message("=== ANSWER: Outlier arrays detected after normalization: No")
  message("=== ANSWER: Number of outliers flagged: 0")
}

# ======================
# PROBE ANNOTATION
# ======================

annotate_probes <- function(expression_set, filtered_data) {
  platform_id <- annotation(expression_set)
  annotation_package <- paste0(platform_id, ".db")
  
  message("Platform detected: ", platform_id)
  
  if (requireNamespace(annotation_package, quietly = TRUE)) {
    library(annotation_package, character.only = TRUE)
    
    probe_ids <- rownames(filtered_data)
    gene_annotations <- tryCatch({
      select(get(annotation_package), 
             keys = probe_ids,
             columns = c("SYMBOL", "GENENAME", "ENTREZID"),
             keytype = "PROBEID")
    }, error = function(e) {
      message("Annotation failed: ", e$message)
      return(NULL)
    })
    
    if (!is.null(gene_annotations)) {
      message("Probe annotation completed: ", nrow(gene_annotations), " mappings")
      return(gene_annotations)
    }
  } else {
    message("Annotation package not available: ", annotation_package)
  }
  
  return(NULL)
}

probe_annotations <- annotate_probes(normalized_expression_set, final_filtered_data)

# ======================
# RESULTS SUMMARY
# ======================

generate_summary_report <- function() {
  summary_file <- file.path(CONFIG$results_dir, "analysis_summary.txt")
  
  summary_content <- c(
    "MICROARRAY ANALYSIS SUMMARY REPORT",
    "===================================",
    paste("GEO Accession:", CONFIG$gse_accession),
    paste("Analysis Date:", Sys.Date()),
    "",
    "DATA DIMENSIONS:",
    paste("Original probes × samples:", nrow(expression_matrix), "×", ncol(expression_matrix)),
    paste("Normalized data dimensions:", nrow(normalized_matrix), "×", ncol(normalized_matrix)),
    paste("Final filtered probes:", nrow(final_filtered_data)),
    "",
    "QUALITY METRICS:",
    paste("Missing sample annotations:", sum(is.na(sample_info$source_name_ch1))),
    paste("Outlier samples detected:", length(outlier_samples)),
    "",
    "FILTERING PARAMETERS:",
    paste("Intensity threshold:", CONFIG$intensity_threshold),
    paste("Variance cutoff:", CONFIG$variance_cutoff),
    "",
    "SAMPLE INFORMATION:",
    paste("Total samples:", ncol(final_filtered_data))
  )
  
  if (!is.null(sample_groups)) {
    group_summary <- table(sample_groups)
    for (group in names(group_summary)) {
      summary_content <- c(summary_content, 
                           paste("  ", group, ":", group_summary[group]))
    }
  }
  
  writeLines(summary_content, summary_file)
  message("Analysis summary saved to: ", summary_file)
}

generate_summary_report()

# ======================
# FINAL MESSAGE
# ======================

message("\n", strrep("=", 50))
message("ANALYSIS COMPLETED SUCCESSFULLY")
message("Results saved in: ", CONFIG$results_dir)
message("Files generated:")
message("  - normalized_expression_data.csv")
message("  - filtered_expression_matrix.csv") 
message("  - intensity_distribution.png")
message("  - analysis_summary.txt")
message("  - Quality control reports (HTML)")
message(strrep("=", 50))