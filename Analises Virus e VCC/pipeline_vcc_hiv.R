# ==============================================================================
# PIPELINE: 3D Spatial Anchoring and Quantitative Phenotyping of VCCs (HIV-1)
# 
# Authors and Affiliations:
#   MSc. Felippe Egon Castrignano Camargo (felippe.camargo@sou.ufmt.br)
#   MSc. Kamilla Batista da Silva Souza (kamillasouza@usp.br)
#   Profa. Dra. Bruna Cunha de Alencar Bargieri (bruna.alencar@usp.br)
#   Laboratório de Biologia Celular do Sistema Imune
#   Instituto de Ciências Biomédicas IV (ICB-IV), Universidade de São Paulo (USP)
#
# Facility Collaborators:
#   Dr. Mário Costa Cruz (confocal@usp.br) - CEFAP, Core CONFOCAL (USP)
#   Dr. Alexandre Freitas da Silva (alexandredasilva@usp.br) - CEFAP, Core GENIAL (USP)
#
# Scientific Hypothesis:
#   Beta-2 integrin (CD18) sustains a biomechanotransducer axis required to maintain
#   the structural integrity and compartmentalization of HIV-1 VCCs in myeloid cells.
#
# Biological Model:
#   U937 human monocytic cell line stably expressing CD169 (Siglec-1)
#   WT vs. CRISPR-NT1 vs. CRISPR-CD18 (EX3 & EX4) | Extensible to Dendritic Cells (DCs)
# ==============================================================================

# Ensures reproducibility of stochastic jitter rendering
set.seed(42)

# Check and install required packages
required_packages <- c("tidyverse", "rstudioapi", "ggpubr", "rstatix", "scales")
missing_packages  <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing_packages) > 0) {
  message(paste("Installing required packages:", paste(missing_packages, collapse = ", ")))
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstudioapi)
  library(ggpubr)
  library(rstatix)
  library(scales)
})

# ==============================================================================
# GLOBAL BIOPHYSICAL & BIOLOGICAL PARAMETERS
# ==============================================================================
# 1. Membrane contact / internalization threshold (µm):
#    Values <= 0 indicate membrane-docked or internalized virus; <= 0.5 captures glycocalyx/filopodia.
MEMBRANE_CONTACT_THRESHOLD <- 0.50

# 2. VCC Volume threshold (µm³):
#    Empirically calibrated via blind ground-truth visual inspection (32 VCCs vs. 692 Spreads; 100% sensitivity & specificity).
VCC_THRESHOLD <- 0.50

# 3. Literature comparative references:
GARCIA_PARAJO_THRESHOLD <- 2.00

# ==============================================================================
# 1. DIRECTORY MAPPING & SELECTION
# ==============================================================================
message(">>> [1/7] Selecting master directory...")

master_dir <- if (rstudioapi::isAvailable()) {
  rstudioapi::selectDirectory(caption = "Select MASTER FOLDER containing Imaris exports")
} else {
  getwd()
}

if (is.null(master_dir) || is.na(master_dir) || master_dir == "") {
  master_dir <- getwd()
}

message(paste("Master directory:", master_dir))

all_files   <- list.files(master_dir, recursive = TRUE, full.names = TRUE)
all_folders <- unique(dirname(all_files))

# Filter folders following the 'CampoXX' naming pattern
valid_folders <- all_folders[str_detect(basename(all_folders), "(?i)Campo\\d+")]
field_ids     <- str_extract(basename(valid_folders), "(?i)Campo\\d+") %>% unique()

if (length(field_ids) == 0) {
  stop("ERROR: No folder matching 'CampoXX' pattern was found. Please check folder names.")
}

message(paste("Identified fields (", length(field_ids), "):", paste(field_ids, collapse = ", ")))

# ==============================================================================
# 2. RESILIENT IMARIS IMPORT FUNCTION
# ==============================================================================
# Imports strictly the 3 essential metrics: Position, Volume, and Shortest Distance.
# Automatically discards whole-image aggregate histograms (_Img=1.csv).
read_imaris_metric <- function(folder, metric) {
  csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  file_names <- basename(csv_files)
  
  # Exclude summary/histogram tables where ID is blank
  object_files <- file_names[!str_detect(file_names, "(?i)_Img=\\d+\\.csv$")]
  
  if (metric == "Position") {
    target_idx <- which(str_detect(object_files, "(?i)(^|.*_)Position\\.csv$"))
  } else if (metric == "Shortest_Distance_to_Surfaces") {
    target_idx <- which(str_detect(object_files, "(?i)(^|.*_)Shortest_Distance_to_Surfaces.*\\.csv$"))
  } else if (metric == "Volume") {
    target_idx <- which(str_detect(object_files, "(?i)(^|.*_)Volume\\.csv$"))
  } else {
    target_idx <- integer(0)
  }
  
  # Fallback if filtered list did not match
  if (length(target_idx) == 0) {
    if (metric == "Shortest_Distance_to_Surfaces") {
      target_idx <- which(str_detect(file_names, "(?i)(^|.*_)Shortest_Distance_to_Surfaces.*\\.csv$"))
    } else {
      target_idx <- which(str_detect(file_names, paste0("(?i)(^|.*_)", metric, ".*\\.csv$")))
    }
  }
  
  if (length(target_idx) > 0) {
    chosen_file <- if (length(object_files) >= target_idx[1]) {
      file.path(folder, object_files[target_idx[1]])
    } else {
      csv_files[target_idx[1]]
    }
    
    df <- tryCatch({
      read.csv(chosen_file, skip = 3, stringsAsFactors = FALSE) %>% 
        select(-matches("^X(\\.\\d+)?$")) %>%
        filter(!is.na(ID) & ID != "") %>%
        mutate(ID = as.integer(ID))
    }, error = function(e) NULL)
    
    return(df)
  } else {
    return(NULL)
  }
}

# ==============================================================================
# 3. PAIRED IMPORT (CELLS & VIRUSES)
# ==============================================================================
message("\n>>> [2/7] Importing metric tables (Position, Volume, Shortest Distance)...")

cell_list  <- list()
virus_list <- list()

for (id in field_ids) {
  field_folders <- valid_folders[str_detect(basename(valid_folders), regex(id, ignore_case = TRUE))]
  
  cell_folder  <- field_folders[str_detect(basename(field_folders), "(?i)c[eé]lulas")]
  virus_folder <- field_folders[str_detect(basename(field_folders), "(?i)v[ií]rus")]
  
  # --- Cells (Requires Position & Volume) ---
  if (length(cell_folder) > 0) {
    pos <- read_imaris_metric(cell_folder[1], "Position")
    vol <- read_imaris_metric(cell_folder[1], "Volume")
    
    if (!is.null(pos) && !is.null(vol)) {
      df_cel <- pos %>% left_join(vol %>% select(ID, Volume), by = "ID")
      
      cell_list[[id]] <- df_cel %>% 
        mutate(
          Field_ID      = id, 
          Original_Name = basename(cell_folder[1]),
          Clean_Name    = str_remove(Original_Name, "(?i)(_Statistics)+$")
        )
      message(paste("  ->", id, "| Cells loaded:", nrow(df_cel)))
    } else {
      warning(paste("  ->", id, "| Failed to read Position or Volume in:", basename(cell_folder[1])))
    }
  }
  
  # --- Virus (Requires Position, Volume & Shortest Distance) ---
  if (length(virus_folder) > 0) {
    pos <- read_imaris_metric(virus_folder[1], "Position")
    vol <- read_imaris_metric(virus_folder[1], "Volume")
    dst <- read_imaris_metric(virus_folder[1], "Shortest_Distance_to_Surfaces")
    
    if (!is.null(pos) && !is.null(vol) && !is.null(dst)) {
      col_dst <- names(dst)[str_detect(names(dst), "(?i)Distance|Shortest")][1]
      if (is.na(col_dst)) col_dst <- names(dst)[1]
      
      dst_clean <- dst %>% select(ID, Shortest_Distance = all_of(col_dst))
      
      df_vir <- pos %>% 
        left_join(vol %>% select(ID, Volume), by = "ID") %>%
        left_join(dst_clean, by = "ID")
      
      virus_list[[id]] <- df_vir %>% 
        mutate(
          Field_ID      = id, 
          Original_Name = basename(virus_folder[1]),
          Clean_Name    = str_remove(Original_Name, "(?i)(_Statistics)+$")
        )
      message(paste("  ->", id, "| Viral particles loaded:", nrow(df_vir)))
    } else {
      warning(paste("  ->", id, "| Failed to read viral metrics in:", basename(virus_folder[1])))
    }
  }
}

global_cells_table <- bind_rows(cell_list)
global_virus_table <- bind_rows(virus_list)

if (nrow(global_cells_table) == 0) {
  stop("ERROR: No cells were imported. Please check Imaris folder names and files.")
}

message(paste0("\nTotal consolidated: ", nrow(global_cells_table), " cells and ", nrow(global_virus_table), " viral particles."))

# ==============================================================================
# 4. TAXONOMIC DECOMPOSITION (Condition_Lineage_Field_Target)
# ==============================================================================
message("\n>>> [3/7] Parsing metadata (Condition_Lineage_Field_Target)...")

global_cells_table <- global_cells_table %>%
  separate(Clean_Name, into = c("Condition", "Lineage", "Field", "Target"), sep = "_", remove = FALSE, extra = "merge")

global_virus_table <- global_virus_table %>%
  separate(Clean_Name, into = c("Condition", "Lineage", "Field", "Target"), sep = "_", remove = FALSE, extra = "merge")

# ==============================================================================
# 5. 3D EUCLIDEAN SPATIAL ANCHORING & CONTACT FILTER
# ==============================================================================
message("\n>>> [4/7] Performing 3D Euclidean spatial anchoring and membrane contact filtering...")

coords_cells <- global_cells_table %>%
  select(Field, Condition, Lineage, Cell_ID = ID, 
         X_Cell = Position.X, Y_Cell = Position.Y, Z_Cell = Position.Z, Cell_Volume = Volume)

coords_virus <- global_virus_table %>%
  select(Field, Condition, Lineage, Virus_ID = ID, 
         X_Vir = Position.X, Y_Vir = Position.Y, Z_Vir = Position.Z, 
         Viral_Volume = Volume, Dist_Surface = Shortest_Distance)

# Step 1: 3D Euclidean distance to the nearest cell centroid in the same field/condition
anchored_viruses <- coords_virus %>%
  left_join(coords_cells, by = c("Field", "Condition", "Lineage"), relationship = "many-to-many") %>%
  mutate(Euclidean_Distance = sqrt((X_Vir - X_Cell)^2 + (Y_Vir - Y_Cell)^2 + (Z_Vir - Z_Cell)^2)) %>%
  group_by(Field, Condition, Lineage, Virus_ID) %>%
  slice_min(Euclidean_Distance, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 2: Biophysical membrane contact filter (Shortest Distance <= 0.50 µm)
associated_viruses <- anchored_viruses %>%
  filter(Dist_Surface <= MEMBRANE_CONTACT_THRESHOLD)

message(paste("Validated cell-associated viruses (Dist <=", MEMBRANE_CONTACT_THRESHOLD, "µm):", 
              nrow(associated_viruses), "of", nrow(anchored_viruses)))

# ==============================================================================
# 6. CELLULAR PHENOTYPIC CLASSIFICATION (VCC+ vs. Spread vs. Uninfected)
# ==============================================================================
message("\n>>> [5/7] Computing phenotypic classification per individual cell...")

safe_max <- function(x) {
  val <- x[!is.na(x)]
  if (length(val) == 0) return(0)
  max(val)
}

safe_sum <- function(x) {
  val <- x[!is.na(x)]
  if (length(val) == 0) return(0)
  sum(val)
}

cell_virus_metrics <- associated_viruses %>%
  group_by(Condition, Lineage, Field, Cell_ID) %>%
  summarise(
    Associated_Particles_Count   = n(),
    Total_Viral_Volume           = safe_sum(Viral_Volume),
    Largest_Compartment_Volume   = safe_max(Viral_Volume),
    Mean_Membrane_Distance       = mean(Dist_Surface, na.rm = TRUE),
    .groups = "drop"
  )

# Retain all cells (preserving uninfected cell denominator)
cell_phenotypic_profile <- coords_cells %>%
  left_join(cell_virus_metrics, by = c("Condition", "Lineage", "Field", "Cell_ID")) %>%
  mutate(
    Associated_Particles_Count = replace_na(Associated_Particles_Count, 0),
    Total_Viral_Volume         = replace_na(Total_Viral_Volume, 0),
    Largest_Compartment_Volume = replace_na(Largest_Compartment_Volume, 0),
    
    # Phenotype definition:
    Phenotype = case_when(
      Associated_Particles_Count == 0          ~ "Uninfected",
      Largest_Compartment_Volume >= VCC_THRESHOLD ~ "VCC+",
      TRUE                                     ~ "Spread"
    ),
    Phenotype = factor(Phenotype, levels = c("Uninfected", "Spread", "VCC+"))
  )

# Standardized condition order
condition_order <- c("WT", "CRISPRV2-NT1", "CRISPRV2-CD18-EX3", "CRISPRV2-CD18-EX4")
active_conditions <- intersect(condition_order, unique(cell_phenotypic_profile$Condition))
if (length(active_conditions) > 0) {
  cell_phenotypic_profile$Condition <- factor(cell_phenotypic_profile$Condition, levels = active_conditions)
  associated_viruses$Condition      <- factor(associated_viruses$Condition, levels = active_conditions)
}

message("\nOverall Phenotypic Distribution:")
print(table(cell_phenotypic_profile$Condition, cell_phenotypic_profile$Phenotype))

# ==============================================================================
# 7. STATISTICAL ANALYSIS & PUBLICATION GRAPHICS (IN ENGLISH)
# ==============================================================================
message("\n>>> [6/7] Generating statistical tables and publication figures...")

phenotype_palette <- c("Uninfected" = "#95a5a6", "Spread" = "#3498db", "VCC+" = "#e74c3c")
condition_palette <- c("WT" = "#2c3e50", "CRISPRV2-NT1" = "#7f8c8d", 
                       "CRISPRV2-CD18-EX3" = "#e67e22", "CRISPRV2-CD18-EX4" = "#d35400")

# --- Descriptive Statistical Summary of Volumes ---
statistical_summary <- associated_viruses %>%
  group_by(Condition) %>%
  summarise(
    N_Particles   = n(),
    Mean_Volume   = round(mean(Viral_Volume, na.rm = TRUE), 4),
    SD_Volume     = round(sd(Viral_Volume, na.rm = TRUE), 4),
    Median_Volume = round(median(Viral_Volume, na.rm = TRUE), 4),
    P75_Volume    = round(quantile(Viral_Volume, 0.75, na.rm = TRUE), 4),
    P90_Volume    = round(quantile(Viral_Volume, 0.90, na.rm = TRUE), 4),
    P95_Volume    = round(quantile(Viral_Volume, 0.95, na.rm = TRUE), 4),
    Max_Volume    = round(max(Viral_Volume, na.rm = TRUE), 4),
    .groups = "drop"
  )

# --- Dynamic Bimodal Cutoff Estimation ---
log_vols_temp  <- log10(associated_viruses$Viral_Volume[associated_viruses$Viral_Volume > 0 & !is.na(associated_viruses$Viral_Volume)])
dens_temp      <- density(log_vols_temp, n = 512, adjust = 1.2)
valleys_temp   <- which(diff(sign(diff(dens_temp$y))) == 2) + 1
real_valleys   <- 10^(dens_temp$x[valleys_temp])
real_valleys   <- real_valleys[real_valleys >= 0.2 & real_valleys <= 3.0]
bimodal_cutoff <- if (length(real_valleys) > 0) round(real_valleys[1], 3) else round(quantile(10^log_vols_temp, 0.90), 3)

# --- FIGURE 1: Phenotypic Proportions across Conditions ---
df_prop <- cell_phenotypic_profile %>%
  count(Condition, Phenotype) %>%
  group_by(Condition) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

fig_01_phenotypes <- ggplot(df_prop, aes(x = Condition, y = Proportion, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.65, color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phenotype_palette) +
  theme_classic(base_size = 13) +
  labs(
    title = "Phenotypic Distribution of HIV-1 Infection in U937-CD169",
    subtitle = paste0("VCC+ (Compartment >= ", VCC_THRESHOLD, " µm³) vs. Surface Spread vs. Uninfected"),
    x = "Condition",
    y = "Relative Frequency (%)",
    fill = "Phenotype"
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"))

# --- FIGURE 2: Compartment Volumes with Means & 3 Reference Thresholds ---
fig_02_volume <- ggplot(associated_viruses, aes(x = Condition, y = Viral_Volume, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.18, alpha = 0.30, size = 1.2, color = "#2c3e50") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "#c0392b", linewidth = 1.0, linetype = "solid") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.8, fill = "#c0392b", color = "black") +
  geom_hline(yintercept = VCC_THRESHOLD, color = "#27ae60", linetype = "solid", linewidth = 0.9) +
  geom_hline(yintercept = bimodal_cutoff, color = "#e74c3c", linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = GARCIA_PARAJO_THRESHOLD, color = "#2980b9", linetype = "dotdash", linewidth = 0.9) +
  annotate("text", x = 0.55, y = VCC_THRESHOLD * 1.15, 
           label = paste0("Visual: ", VCC_THRESHOLD, " µm³"), color = "#27ae60", fontface = "bold", size = 3.2, hjust = 0) +
  annotate("text", x = 0.55, y = bimodal_cutoff * 1.15, 
           label = paste0("Bimodal: ", bimodal_cutoff, " µm³"), color = "#e74c3c", fontface = "bold", size = 3.2, hjust = 0) +
  annotate("text", x = 0.55, y = GARCIA_PARAJO_THRESHOLD * 1.15, 
           label = paste0("Garcia-Parajo: ", GARCIA_PARAJO_THRESHOLD, " µm³"), color = "#2980b9", fontface = "bold", size = 3.2, hjust = 0) +
  scale_y_log10(labels = scales::label_log()) +
  scale_fill_manual(values = condition_palette) +
  theme_classic(base_size = 13) +
  labs(
    title = "Viral Compartment Volume and VCC Thresholds",
    subtitle = "Green: Visual Ground Truth (0.50 µm³) | Red: Bimodal | Blue: Garcia-Parajo et al. (2.00 µm³)",
    x = "Condition",
    y = expression("Viral Volume (" * mu * "m"^3 * " - Log10 Scale)"),
    fill = "Condition"
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"))

# --- FIGURE 3: Total Viral Load per Infected Cell ---
infected_cells <- cell_phenotypic_profile %>% filter(Phenotype %in% c("VCC+", "Spread"))

fig_03_viral_load <- ggplot(infected_cells, aes(x = Condition, y = Total_Viral_Volume, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "#c0392b", linewidth = 1.0) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.8, fill = "#c0392b", color = "black") +
  scale_y_log10(labels = scales::label_log()) +
  scale_fill_manual(values = condition_palette) +
  theme_classic(base_size = 13) +
  labs(
    title = "Total Cell-Associated Viral Load",
    subtitle = "Infected cells only | Red diamonds: Condition mean",
    x = "Condition",
    y = expression("Total Viral Volume per Cell (" * mu * "m"^3 * " - Log10)"),
    fill = "Condition"
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"))

# --- FIGURE 4: Cumulative Distribution Function (ECDF) ---
fig_04_ecdf <- ggplot(associated_viruses, aes(x = Viral_Volume, color = Condition)) +
  stat_ecdf(geom = "step", linewidth = 1.1) +
  geom_vline(xintercept = VCC_THRESHOLD, color = "#27ae60", linetype = "solid", linewidth = 0.9) +
  geom_vline(xintercept = bimodal_cutoff, color = "#e74c3c", linetype = "dashed", linewidth = 0.9) +
  geom_vline(xintercept = GARCIA_PARAJO_THRESHOLD, color = "#2980b9", linetype = "dotdash", linewidth = 0.9) +
  annotate("text", x = VCC_THRESHOLD * 1.1, y = 0.25, label = paste0("Visual: ", VCC_THRESHOLD, " µm³"), 
           color = "#27ae60", fontface = "bold", size = 3.2, hjust = 0) +
  annotate("text", x = bimodal_cutoff * 1.1, y = 0.50, label = paste0("Bimodal: ", bimodal_cutoff, " µm³"), 
           color = "#e74c3c", fontface = "bold", size = 3.2, hjust = 0) +
  annotate("text", x = GARCIA_PARAJO_THRESHOLD * 1.1, y = 0.75, label = paste0("Garcia-Parajo: ", GARCIA_PARAJO_THRESHOLD, " µm³"), 
           color = "#2980b9", fontface = "bold", size = 3.2, hjust = 0) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = condition_palette) +
  theme_classic(base_size = 13) +
  labs(
    title = "Empirical Cumulative Distribution Function (ECDF) of Viral Volumes",
    subtitle = "Green: Visual Ground Truth (0.50 µm³) | Red: Bimodal | Blue: Garcia-Parajo et al. (2.00 µm³)",
    x = expression("Viral Volume (" * mu * "m"^3 * " - Log10 Scale)"),
    y = "Cumulative Percentage (%)",
    color = "Condition"
  )

# --- Statistical Hypothesis Testing ---
kruskal_test_res <- tryCatch({
  associated_viruses %>% kruskal_test(Viral_Volume ~ Condition)
}, error = function(e) NULL)

pairwise_test_res <- tryCatch({
  associated_viruses %>% wilcox_test(Viral_Volume ~ Condition, p.adjust.method = "fdr")
}, error = function(e) NULL)

# ==============================================================================
# 8. EXPORTING RESULTS & FIGURES (UNIVERSAL ENGLISH NAMING)
# ==============================================================================
message("\n>>> [7/7] Saving output tables and figures in English...")

output_dir <- file.path(master_dir, "Results_VCC_Analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(cell_phenotypic_profile, file.path(output_dir, "Table_Cellular_Phenotypic_Profile_HIV.csv"), row.names = FALSE)
write.csv(associated_viruses, file.path(output_dir, "Table_Associated_Anchored_Viruses.csv"), row.names = FALSE)
write.csv(statistical_summary, file.path(output_dir, "Table_Statistical_Summary_Volumes.csv"), row.names = FALSE)

if (!is.null(pairwise_test_res)) {
  write.csv(pairwise_test_res, file.path(output_dir, "Table_Pairwise_Statistics_Volume.csv"), row.names = FALSE)
}

ggsave(file.path(output_dir, "Figure_01_Phenotypic_Proportions.png"), fig_01_phenotypes, width = 7.5, height = 5.0, dpi = 300)
ggsave(file.path(output_dir, "Figure_02_Compartment_Volume_Trend.png"), fig_02_volume, width = 8.0, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "Figure_03_Cellular_Viral_Load.png"), fig_03_viral_load, width = 8.0, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "Figure_04_ECDF_Volume_Distribution.png"), fig_04_ecdf, width = 8.5, height = 5.5, dpi = 300)

message(paste0("\n=== PROCESSING COMPLETED SUCCESSFULLY! ==="))
message(paste0("Files and figures saved to: ", output_dir))
