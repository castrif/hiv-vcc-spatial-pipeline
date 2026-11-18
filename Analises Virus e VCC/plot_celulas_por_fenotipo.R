# ==============================================================================
# SCRIPT: Visualização Celular por Fenótipo Separado (NI vs. Spread vs. VCC+)
# Boxplot com Dotplot / Jitter por Célula Individual
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(scales)
})

# 1. Localizar tabela de perfil fenotípico celular
caminho_csv <- file.path(dirname(getwd()), "Results_VCC_Analysis", "Table_Cellular_Phenotypic_Profile_HIV.csv")
if (!file.exists(caminho_csv)) {
  caminho_csv <- file.path(getwd(), "Results_VCC_Analysis", "Table_Cellular_Phenotypic_Profile_HIV.csv")
}

if (!file.exists(caminho_csv)) {
  stop("Arquivo 'Table_Cellular_Phenotypic_Profile_HIV.csv' não encontrado. Execute o 'pipeline_vcc_hiv.R' primeiro.")
}

message(paste("Carregando dados de:", caminho_csv))
df_cells <- read.csv(caminho_csv, stringsAsFactors = FALSE)

# 2. Formatar rótulos de fenótipo e ordem
df_cells <- df_cells %>%
  mutate(
    Phenotype_Clean = case_when(
      Phenotype == "Uninfected" ~ "NI (Não-Infectada)",
      Phenotype == "Spread"     ~ "Spread (< 0.50 µm³)",
      Phenotype == "VCC+"       ~ "VCC+ (≥ 0.50 µm³)",
      TRUE                      ~ Phenotype
    ),
    Phenotype_Clean = factor(Phenotype_Clean, levels = c("NI (Não-Infectada)", "Spread (< 0.50 µm³)", "VCC+ (≥ 0.50 µm³)"))
  )

phenotype_palette <- c(
  "NI (Não-Infectada)"  = "#95a5a6",
  "Spread (< 0.50 µm³)" = "#3498db",
  "VCC+ (≥ 0.50 µm³)"   = "#e74c3c"
)

pasta_saida <- dirname(caminho_csv)

# ------------------------------------------------------------------------------
# FIGURA A: Carga Viral Celular Total (Total Viral Volume) por Fenótipo Separado
# ------------------------------------------------------------------------------
fig_carga <- ggplot(df_cells, aes(x = Phenotype_Clean, y = Total_Viral_Volume, fill = Phenotype_Clean)) +
  geom_boxplot(width = 0.45, alpha = 0.65, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(aes(shape = Experiment), width = 0.18, alpha = 0.60, size = 2.4, color = "#2c3e50") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4.0, fill = "yellow", color = "black") +
  stat_compare_means(comparisons = list(c("Spread (< 0.50 µm³)", "VCC+ (≥ 0.50 µm³)")), 
                     method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = phenotype_palette) +
  theme_classic(base_size = 13) +
  labs(
    title = "Carga Viral Total por Célula Individual (Por Fenótipo)",
    subtitle = "Barras separadas: NI vs. Spread vs. VCC+ | Diamante amarelo: Média | Formas: Réplicas/Experimentos",
    x = "Fenótipo Celular",
    y = expression("Volume Viral Total por Célula (" * mu * "m"^3 * ")"),
    fill = "Fenótipo",
    shape = "Experimento"
  ) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

# Se houver mais de uma condição (ex: WT vs KO), adiciona visualização comparativa
if (length(unique(df_cells$Condition)) > 1) {
  fig_carga_cond <- ggplot(df_cells, aes(x = Phenotype_Clean, y = Total_Viral_Volume, fill = Condition)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               alpha = 0.6, size = 2.0, aes(shape = Experiment)) +
    theme_classic(base_size = 13) +
    labs(
      title = "Carga Viral por Célula: Comparativo entre Condições Genéticas",
      x = "Fenótipo Celular",
      y = expression("Volume Viral Total (" * mu * "m"^3 * ")"),
      fill = "Condição Genética",
      shape = "Experimento"
    ) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(pasta_saida, "Figure_05_Carga_Viral_Por_Fenotipo_e_Condicao.png"), 
         fig_carga_cond, width = 9.5, height = 5.5, dpi = 300)
}

# ------------------------------------------------------------------------------
# FIGURA B: Volume do Maior Compartimento (Largest Compartment Volume)
# ------------------------------------------------------------------------------
fig_maior_comp <- ggplot(df_cells %>% filter(Phenotype != "Uninfected"), 
                         aes(x = Phenotype_Clean, y = Largest_Compartment_Volume, fill = Phenotype_Clean)) +
  geom_boxplot(width = 0.40, alpha = 0.65, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(aes(shape = Experiment), width = 0.18, alpha = 0.60, size = 2.4, color = "#2c3e50") +
  geom_hline(yintercept = 0.50, color = "#27ae60", linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 1.5, y = 0.58, label = "Corte VCC (0.50 µm³)", color = "#27ae60", fontface = "bold", size = 3.8) +
  scale_y_log10(labels = scales::label_log()) +
  scale_fill_manual(values = phenotype_palette) +
  theme_classic(base_size = 13) +
  labs(
    title = "Volume do Maior Compartimento Viral por Célula",
    subtitle = "Linha tracejada verde: Limiar de corte empírico (0.50 µm³)",
    x = "Fenótipo Celular",
    y = expression("Volume do Maior Compartimento (" * mu * "m"^3 * " - Log10)"),
    fill = "Fenótipo",
    shape = "Experimento"
  ) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

# ------------------------------------------------------------------------------
# FIGURA C: Contagem de Partículas Virais Associadas por Célula
# ------------------------------------------------------------------------------
fig_particulas <- ggplot(df_cells, aes(x = Phenotype_Clean, y = Associated_Particles_Count, fill = Phenotype_Clean)) +
  geom_boxplot(width = 0.45, alpha = 0.65, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(aes(shape = Experiment), width = 0.18, alpha = 0.60, size = 2.4, color = "#2c3e50") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4.0, fill = "yellow", color = "black") +
  scale_fill_manual(values = phenotype_palette) +
  theme_classic(base_size = 13) +
  labs(
    title = "Número de Partículas Virais Ancoradas por Célula",
    subtitle = "Contagem total por célula | Diamante amarelo: Média",
    x = "Fenótipo Celular",
    y = "Número de Partículas Virais por Célula",
    fill = "Fenótipo",
    shape = "Experimento"
  ) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

# Salvar gráficos
ggsave(file.path(pasta_saida, "Figure_05A_Carga_Viral_Celular_Por_Fenotipo.png"), fig_carga, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(pasta_saida, "Figure_05B_Volume_Maior_Compartimento_Por_Fenotipo.png"), fig_maior_comp, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(pasta_saida, "Figure_05C_Contagem_Particulas_Por_Fenotipo.png"), fig_particulas, width = 8.5, height = 5.5, dpi = 300)

message("\n=== GRÁFICOS GERADOS COM SUCESSO! ===")
message(paste("Salvos em:", pasta_saida))
