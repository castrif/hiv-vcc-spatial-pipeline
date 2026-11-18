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
# 1. DIRECTORY MAPPING & MULTI-EXPERIMENT SELECTION
# ==============================================================================
message(">>> [1/7] Mapping master directory and experiment folders...")

master_dir <- if (rstudioapi::isAvailable()) {
  tryCatch(rstudioapi::selectDirectory(caption = "Select MASTER FOLDER containing Imaris exports or experiments"),
           error = function(e) NULL)
} else {
  NULL
}

if (is.null(master_dir) || is.na(master_dir) || master_dir == "") {
  # Smart fallback: if working inside 'Analises Virus e VCC', check parent dir first
  if (basename(getwd()) == "Analises Virus e VCC" && 
      length(list.files(dirname(getwd()), pattern = "(?i)Campo[0-9]+", recursive = TRUE)) > 0) {
    master_dir <- dirname(getwd())
  } else {
    master_dir <- getwd()
  }
}

message(paste("Master directory:", master_dir))

all_files   <- list.files(master_dir, recursive = TRUE, full.names = TRUE)
all_folders <- unique(dirname(all_files))

# Exclude internal git/rproj/temporary folders
valid_folders <- all_folders[
  str_detect(basename(all_folders), "(?i)Campo[0-9]+") &
  !str_detect(all_folders, "\\.git|\\.Rproj\\.user")
]

if (length(valid_folders) == 0) {
  stop("ERROR: No folder matching 'CampoXX' pattern was found. Please check folder names.")
}

# Distinguish if folders are inside sub-experiments (e.g. Exp1, Exp2, Exp3) or in root
norm_master <- normalizePath(master_dir, winslash = "/", mustWork = FALSE)

folder_metadata <- tibble(folder_path = valid_folders) %>%
  mutate(
    norm_path  = normalizePath(folder_path, winslash = "/", mustWork = FALSE),
    rel_path   = sub(norm_master, "", norm_path, fixed = TRUE) %>% str_remove("^/"),
    parent_dir = dirname(rel_path),
    # If the folder is in the root, assign Exp01; otherwise, use the subfolder name
    Experiment = if_else(parent_dir == "." | parent_dir == "", "Exp01", basename(parent_dir)),
    Folder_Name = basename(folder_path),
    Clean_Name  = str_remove(Folder_Name, "(?i)(_Statistics)+$")
  ) %>%
  # Exclude GitHub backup/example folder if actual Exp folders are present
  filter(!(Experiment == "EXEMPLO_GITHUB" & any(str_detect(Experiment, "(?i)^Exp[0-9]+")))) %>%
  separate(Clean_Name, into = c("Condition", "Lineage", "Field", "Target"), sep = "_", remove = FALSE, extra = "merge") %>%
  mutate(
    Target_Type = case_when(
      str_detect(Target, "(?i)c[eé]lula") ~ "Cells",
      str_detect(Target, "(?i)v[ií]rus")   ~ "Virus",
      TRUE                                ~ "Other"
    ),
    # Unique composite key: Experiment + Condition + Lineage + Field
    Group_ID = paste(Experiment, Condition, Lineage, Field, sep = "___")
  ) %>%
  filter(Target_Type %in% c("Cells", "Virus"))

experiments_found <- unique(folder_metadata$Experiment)
group_ids_found   <- unique(folder_metadata$Group_ID)

message(paste0("Identified experiments (", length(experiments_found), "): ", paste(experiments_found, collapse = ", ")))
message(paste0("Identified unique acquisition fields: ", length(group_ids_found)))

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
# 3. PAIRED IMPORT BY GROUP_ID (EXPERIMENT + CONDITION + LINEAGE + FIELD)
# ==============================================================================
message("\n>>> [2/7] Importing metric tables (Position, Volume, Shortest Distance)...")

cell_list  <- list()
virus_list <- list()

for (gid in group_ids_found) {
  field_rows <- folder_metadata %>% filter(Group_ID == gid)
  
  cell_row  <- field_rows %>% filter(Target_Type == "Cells") %>% slice(1)
  virus_row <- field_rows %>% filter(Target_Type == "Virus") %>% slice(1)
  
  # --- Cells (Requires Position & Volume) ---
  if (nrow(cell_row) > 0) {
    pos <- read_imaris_metric(cell_row$folder_path, "Position")
    vol <- read_imaris_metric(cell_row$folder_path, "Volume")
    
    if (!is.null(pos) && !is.null(vol) && nrow(pos) > 0 && nrow(vol) > 0) {
      df_cel <- pos %>% left_join(vol %>% select(ID, Volume), by = "ID")
      
      cell_list[[gid]] <- df_cel %>% 
        mutate(
          Group_ID      = gid,
          Experiment    = cell_row$Experiment,
          Condition     = cell_row$Condition,
          Lineage       = cell_row$Lineage,
          Field         = cell_row$Field,
          Original_Name = cell_row$Folder_Name
        )
      message(paste("  ->", cell_row$Experiment, "|", cell_row$Field, "| Cells loaded:", nrow(df_cel)))
    } else {
      warning(paste("  -> Skipping invalid/incomplete cell folder:", cell_row$Folder_Name))
    }
  }
  
  # --- Virus (Requires Position, Volume & Shortest Distance) ---
  if (nrow(virus_row) > 0) {
    pos <- read_imaris_metric(virus_row$folder_path, "Position")
    vol <- read_imaris_metric(virus_row$folder_path, "Volume")
    dst <- read_imaris_metric(virus_row$folder_path, "Shortest_Distance_to_Surfaces")
    
    if (!is.null(pos) && !is.null(vol) && !is.null(dst) && nrow(pos) > 0 && nrow(vol) > 0) {
      col_dst <- names(dst)[str_detect(names(dst), "(?i)Distance|Shortest")][1]
      if (is.na(col_dst)) col_dst <- names(dst)[1]
      
      dst_clean <- dst %>% select(ID, Shortest_Distance = all_of(col_dst))
      
      df_vir <- pos %>% 
        left_join(vol %>% select(ID, Volume), by = "ID") %>%
        left_join(dst_clean, by = "ID")
      
      virus_list[[gid]] <- df_vir %>% 
        mutate(
          Group_ID      = gid,
          Experiment    = virus_row$Experiment,
          Condition     = virus_row$Condition,
          Lineage       = virus_row$Lineage,
          Field         = virus_row$Field,
          Original_Name = virus_row$Folder_Name
        )
      message(paste("  ->", virus_row$Experiment, "|", virus_row$Field, "| Viral particles loaded:", nrow(df_vir)))
    } else {
      warning(paste("  -> Skipping invalid/incomplete virus folder:", virus_row$Folder_Name))
    }
  }
}

global_cells_table <- bind_rows(cell_list)
global_virus_table <- bind_rows(virus_list)

if (nrow(global_cells_table) == 0) {
  stop("ERROR: No cells were imported. Please check Imaris folder names and files.")
}

message(paste0("\nTotal consolidated: ", nrow(global_cells_table), " cells and ", nrow(global_virus_table), " viral particles across all experiments."))

# ==============================================================================
# 4. 3D EUCLIDEAN SPATIAL ANCHORING & CONTACT FILTER
# ==============================================================================
message("\n>>> [3/7] Performing 3D Euclidean spatial anchoring and membrane contact filtering...")

coords_cells <- global_cells_table %>%
  select(Group_ID, Experiment, Condition, Lineage, Field, Cell_ID = ID, 
         X_Cell = Position.X, Y_Cell = Position.Y, Z_Cell = Position.Z, Cell_Volume = Volume)

if (nrow(global_virus_table) > 0) {
  coords_virus <- global_virus_table %>%
    select(Group_ID, Experiment, Condition, Lineage, Field, Virus_ID = ID, 
           X_Vir = Position.X, Y_Vir = Position.Y, Z_Vir = Position.Z, 
           Viral_Volume = Volume, Dist_Surface = Shortest_Distance)
  
  # Step 1: 3D Euclidean distance strictly within the same acquisition field (Group_ID)
  anchored_viruses <- coords_virus %>%
    inner_join(coords_cells, by = c("Group_ID", "Experiment", "Condition", "Lineage", "Field"), relationship = "many-to-many") %>%
    mutate(Euclidean_Distance = sqrt((X_Vir - X_Cell)^2 + (Y_Vir - Y_Cell)^2 + (Z_Vir - Z_Cell)^2)) %>%
    group_by(Group_ID, Virus_ID) %>%
    slice_min(Euclidean_Distance, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Step 2: Biophysical membrane contact filter (Shortest Distance <= 0.50 µm)
  associated_viruses <- anchored_viruses %>%
    filter(Dist_Surface <= MEMBRANE_CONTACT_THRESHOLD)
  
  message(paste("Validated cell-associated viruses (Dist <=", MEMBRANE_CONTACT_THRESHOLD, "µm):", 
                nrow(associated_viruses), "of", nrow(anchored_viruses)))
} else {
  associated_viruses <- tibble()
}

# ==============================================================================
# 5. CELLULAR PHENOTYPIC CLASSIFICATION (VCC+ vs. Spread vs. Uninfected)
# ==============================================================================
message("\n>>> [4/7] Computing phenotypic classification per individual cell...")

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

if (nrow(associated_viruses) > 0) {
  cell_virus_metrics <- associated_viruses %>%
    group_by(Group_ID, Experiment, Condition, Lineage, Field, Cell_ID) %>%
    summarise(
      Associated_Particles_Count = n(),
      Total_Viral_Volume         = safe_sum(Viral_Volume),
      Largest_Compartment_Volume = safe_max(Viral_Volume),
      Mean_Membrane_Distance     = mean(Dist_Surface, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  cell_virus_metrics <- tibble(
    Group_ID = character(), Experiment = character(), Condition = character(),
    Lineage = character(), Field = character(), Cell_ID = integer(),
    Associated_Particles_Count = integer(), Total_Viral_Volume = numeric(),
    Largest_Compartment_Volume = numeric(), Mean_Membrane_Distance = numeric()
  )
}

# Retain all cells (preserving uninfected cell denominator)
cell_phenotypic_profile <- coords_cells %>%
  left_join(cell_virus_metrics, by = c("Group_ID", "Experiment", "Condition", "Lineage", "Field", "Cell_ID")) %>%
  mutate(
    Associated_Particles_Count = replace_na(Associated_Particles_Count, 0),
    Total_Viral_Volume         = replace_na(Total_Viral_Volume, 0),
    Largest_Compartment_Volume = replace_na(Largest_Compartment_Volume, 0),
    
    # Phenotype definition:
    Phenotype = case_when(
      Associated_Particles_Count == 0             ~ "Uninfected",
      Largest_Compartment_Volume >= VCC_THRESHOLD ~ "VCC+",
      TRUE                                        ~ "Spread"
    ),
    Phenotype = factor(Phenotype, levels = c("Uninfected", "Spread", "VCC+"))
  )

# Dynamic and resilient condition ordering (supports WT, KOCD18, CRISPR, and any user condition)
found_conditions <- unique(as.character(cell_phenotypic_profile$Condition))
preferred_order  <- c("WT", "CRISPRV2-NT1", "NT1", "CRISPR-NT1", 
                      "CRISPRV2-CD18-EX3", "CRISPRV2-CD18-EX4", 
                      "KOCD18", "KO-CD18", "CD18KO", "CD18-KO", "KO")
ordered_levels   <- c(intersect(preferred_order, found_conditions), 
                      setdiff(found_conditions, preferred_order))

cell_phenotypic_profile$Condition <- factor(cell_phenotypic_profile$Condition, levels = ordered_levels)
if (nrow(associated_viruses) > 0) {
  associated_viruses$Condition <- factor(associated_viruses$Condition, levels = ordered_levels)
}

message("\nOverall Phenotypic Distribution:")
print(table(cell_phenotypic_profile$Condition, cell_phenotypic_profile$Phenotype))

# ==============================================================================
# 6. STATISTICAL ANALYSIS & PUBLICATION GRAPHICS (IN ENGLISH)
# ==============================================================================
message("\n>>> [5/7] Generating statistical tables and publication figures...")

phenotype_palette <- c("Uninfected" = "#95a5a6", "Spread" = "#3498db", "VCC+" = "#e74c3c")

# Dynamic condition color palette
base_colors <- c(
  "WT" = "#2c3e50", "CRISPRV2-NT1" = "#7f8c8d", "NT1" = "#7f8c8d", "CRISPR-NT1" = "#7f8c8d",
  "CRISPRV2-CD18-EX3" = "#e67e22", "CRISPRV2-CD18-EX4" = "#d35400",
  "KOCD18" = "#e67e22", "KO-CD18" = "#e67e22", "CD18KO" = "#d35400", "KO" = "#e67e22"
)
extra_palette <- c("#16a085", "#8e44ad", "#2980b9", "#f39c12", "#27ae60", "#c0392b")
condition_palette <- base_colors[ordered_levels]
unassigned <- which(is.na(condition_palette))
if (length(unassigned) > 0) {
  condition_palette[unassigned] <- extra_palette[seq_along(unassigned)]
}
names(condition_palette) <- ordered_levels

# --- Descriptive Statistical Summary of Volumes (Overall & By Experiment) ---
if (nrow(associated_viruses) > 0) {
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
  
  statistical_summary_exp <- associated_viruses %>%
    group_by(Experiment, Condition) %>%
    summarise(
      N_Particles   = n(),
      Mean_Volume   = round(mean(Viral_Volume, na.rm = TRUE), 4),
      Median_Volume = round(median(Viral_Volume, na.rm = TRUE), 4),
      P75_Volume    = round(quantile(Viral_Volume, 0.75, na.rm = TRUE), 4),
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
} else {
  statistical_summary <- tibble()
  statistical_summary_exp <- tibble()
  bimodal_cutoff <- 0.50
}

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
    subtitle = paste0("VCC+ (Compartment >= ", VCC_THRESHOLD, " µm³) vs. Surface Spread vs. Uninfected | Pooled Experiments (N=", length(experiments_found), ")"),
    x = "Condition",
    y = "Relative Frequency (%)",
    fill = "Phenotype"
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"))

# FIGURE 1B: Phenotypic Proportions by Experiment (Biological Replicate Quality Control)
df_prop_exp <- cell_phenotypic_profile %>%
  count(Experiment, Condition, Phenotype) %>%
  group_by(Experiment, Condition) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

fig_01b_exp <- ggplot(df_prop_exp, aes(x = Condition, y = Proportion, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.65, color = "black", linewidth = 0.3) +
  facet_wrap(~ Experiment) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phenotype_palette) +
  theme_classic(base_size = 12) +
  labs(
    title = "Biological Replicate Reproducibility: Phenotype Proportions by Experiment",
    x = "Condition",
    y = "Relative Frequency (%)",
    fill = "Phenotype"
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
        strip.background = element_rect(fill = "#ecf0f1", color = NA),
        strip.text = element_text(face = "bold"))

# --- FIGURE 2: Compartment Volumes with Means & 3 Reference Thresholds ---
if (nrow(associated_viruses) > 0) {
  fig_02_volume <- ggplot(associated_viruses, aes(x = Condition, y = Viral_Volume, fill = Condition)) +
    geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(aes(shape = Experiment), width = 0.18, alpha = 0.35, size = 1.3, color = "#2c3e50") +
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
      subtitle = paste0("Green: Ground Truth (0.50 µm³) | Red: Bimodal | Blue: Lit. (2.00 µm³) | Shapes denote Experiments (N=", length(experiments_found), ")"),
      x = "Condition",
      y = expression("Viral Volume (" * mu * "m"^3 * " - Log10 Scale)"),
      fill = "Condition",
      shape = "Experiment"
    ) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"))
  
  # --- FIGURE 3: Total Viral Load per Infected Cell ---
  infected_cells <- cell_phenotypic_profile %>% filter(Phenotype %in% c("VCC+", "Spread"))
  
  fig_03_viral_load <- ggplot(infected_cells, aes(x = Condition, y = Total_Viral_Volume, fill = Condition)) +
    geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(aes(shape = Experiment), width = 0.15, alpha = 0.4, size = 1.5) +
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
      fill = "Condition",
      shape = "Experiment"
    ) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"))
  
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
  pairwise_test_res <- tryCatch({
    if (length(unique(associated_viruses$Condition)) > 1) {
      associated_viruses %>% wilcox_test(Viral_Volume ~ Condition, p.adjust.method = "fdr")
    } else NULL
  }, error = function(e) NULL)
}

# --- FIGURE 5A & 5B: Cellular Phenotypes Separated (NI vs. Spread vs. VCC+) ---
df_cells_plot <- cell_phenotypic_profile %>%
  mutate(
    Phenotype_Display = case_when(
      Phenotype == "Uninfected" ~ "NI (Uninfected)",
      Phenotype == "Spread"     ~ "Spread (< 0.50 µm³)",
      Phenotype == "VCC+"       ~ "VCC+ (≥ 0.50 µm³)",
      TRUE                      ~ as.character(Phenotype)
    ),
    Phenotype_Display = factor(Phenotype_Display, levels = c("NI (Uninfected)", "Spread (< 0.50 µm³)", "VCC+ (≥ 0.50 µm³)"))
  )

palette_display <- c(
  "NI (Uninfected)"     = "#95a5a6",
  "Spread (< 0.50 µm³)" = "#3498db",
  "VCC+ (≥ 0.50 µm³)"   = "#e74c3c"
)

fig_05a_phenotype_load <- ggplot(df_cells_plot, aes(x = Phenotype_Display, y = Total_Viral_Volume, fill = Phenotype_Display)) +
  geom_boxplot(width = 0.45, alpha = 0.65, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(aes(shape = Experiment), width = 0.18, alpha = 0.60, size = 2.4, color = "#2c3e50") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4.0, fill = "yellow", color = "black") +
  stat_compare_means(comparisons = list(c("Spread (< 0.50 µm³)", "VCC+ (≥ 0.50 µm³)")), 
                     method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = palette_display) +
  theme_classic(base_size = 13) +
  labs(
    title = "Cellular Viral Load by Separated Infection Phenotype",
    subtitle = "Distinct categories: NI vs. Spread vs. VCC+ | Yellow diamond: Mean | Shapes denote Experiment",
    x = "Cellular Phenotype",
    y = expression("Total Viral Volume per Cell (" * mu * "m"^3 * ")"),
    fill = "Phenotype",
    shape = "Experiment"
  ) +
  theme(legend.position = "right", axis.text.x = element_text(face = "bold"))

fig_05b_phenotype_largest <- ggplot(df_cells_plot %>% filter(Phenotype != "Uninfected"), 
                                    aes(x = Phenotype_Display, y = Largest_Compartment_Volume, fill = Phenotype_Display)) +
  geom_boxplot(width = 0.40, alpha = 0.65, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(aes(shape = Experiment), width = 0.18, alpha = 0.60, size = 2.4, color = "#2c3e50") +
  geom_hline(yintercept = VCC_THRESHOLD, color = "#27ae60", linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 1.5, y = VCC_THRESHOLD * 1.15, label = paste0("VCC Cutoff: ", VCC_THRESHOLD, " µm³"), 
           color = "#27ae60", fontface = "bold", size = 3.8) +
  scale_y_log10(labels = scales::label_log()) +
  scale_fill_manual(values = palette_display) +
  theme_classic(base_size = 13) +
  labs(
    title = "Largest Viral Compartment Volume per Cell",
    subtitle = paste0("Dashed green line: Empirical Ground Truth cutoff (", VCC_THRESHOLD, " µm³)"),
    x = "Cellular Phenotype",
    y = expression("Largest Compartment Volume (" * mu * "m"^3 * " - Log10)"),
    fill = "Phenotype",
    shape = "Experiment"
  ) +
  theme(legend.position = "right", axis.text.x = element_text(face = "bold"))

# ==============================================================================
# 7. EXPORTING RESULTS & FIGURES (UNIVERSAL ENGLISH NAMING)
# ==============================================================================
message("\n>>> [6/7] Saving output tables and figures in English...")

output_dir <- file.path(master_dir, "Results_VCC_Analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(cell_phenotypic_profile, file.path(output_dir, "Table_Cellular_Phenotypic_Profile_HIV.csv"), row.names = FALSE)
if (nrow(associated_viruses) > 0) {
  write.csv(associated_viruses, file.path(output_dir, "Table_Associated_Anchored_Viruses.csv"), row.names = FALSE)
  write.csv(statistical_summary, file.path(output_dir, "Table_Statistical_Summary_Volumes.csv"), row.names = FALSE)
  write.csv(statistical_summary_exp, file.path(output_dir, "Table_Statistical_Summary_By_Experiment.csv"), row.names = FALSE)
  
  if (!is.null(pairwise_test_res)) {
    write.csv(pairwise_test_res, file.path(output_dir, "Table_Pairwise_Statistics_Volume.csv"), row.names = FALSE)
  }
  
  ggsave(file.path(output_dir, "Figure_02_Compartment_Volume_Trend.png"), fig_02_volume, width = 8.5, height = 5.5, dpi = 300)
  ggsave(file.path(output_dir, "Figure_03_Cellular_Viral_Load.png"), fig_03_viral_load, width = 8.5, height = 5.5, dpi = 300)
  ggsave(file.path(output_dir, "Figure_04_ECDF_Volume_Distribution.png"), fig_04_ecdf, width = 8.5, height = 5.5, dpi = 300)
}

ggsave(file.path(output_dir, "Figure_01_Phenotypic_Proportions.png"), fig_01_phenotypes, width = 7.5, height = 5.0, dpi = 300)
ggsave(file.path(output_dir, "Figure_01B_Phenotypic_Proportions_By_Experiment.png"), fig_01b_exp, width = 9.0, height = 5.0, dpi = 300)
ggsave(file.path(output_dir, "Figure_05A_Cellular_Viral_Load_By_Phenotype.png"), fig_05a_phenotype_load, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "Figure_05B_Largest_Compartment_By_Phenotype.png"), fig_05b_phenotype_largest, width = 8.5, height = 5.5, dpi = 300)

message(paste0("\n=== PROCESSING COMPLETED SUCCESSFULLY! ==="))
message(paste0("Files and figures saved to: ", output_dir))
