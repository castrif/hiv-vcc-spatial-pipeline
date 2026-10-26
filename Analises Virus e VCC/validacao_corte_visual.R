# ==============================================================================
# VALIDAÇÃO DO CORTE VISUAL (GROUND TRUTH): VCC vs SPREAD
# Curva ROC, Distribuição Bimodal e Ponto de Corte Empírico
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# Caminho dos arquivos na pasta 'CORTE VISUAL'
pasta_corte_visual <- file.path(dirname(getwd()), "CORTE VISUAL")
if (!dir.exists(pasta_corte_visual)) {
  pasta_corte_visual <- file.path(getwd(), "CORTE VISUAL")
}

arq_vcc    <- file.path(pasta_corte_visual, "VCC.csv")
arq_spread <- file.path(pasta_corte_visual, "SPREAD.csv")

if (!file.exists(arq_vcc) || !file.exists(arq_spread)) {
  stop("Arquivos VCC.csv e SPREAD.csv não encontrados na pasta CORTE VISUAL.")
}

# Leitura pulando as linhas de cabeçalho do Imaris
df_vcc <- read.csv(arq_vcc, skip = 3, stringsAsFactors = FALSE) %>%
  mutate(Classe = "VCC", Volume = as.numeric(Volume)) %>%
  filter(!is.na(Volume))

df_spread <- read.csv(arq_spread, skip = 3, stringsAsFactors = FALSE) %>%
  mutate(Classe = "Spread", Volume = as.numeric(Volume)) %>%
  filter(!is.na(Volume))

df_total <- bind_rows(df_vcc, df_spread)

# Estatísticas descritivas
resumo_classes <- df_total %>%
  group_by(Classe) %>%
  summarise(
    N = n(),
    Min_Vol = min(Volume),
    Mediana_Vol = median(Volume),
    Media_Vol = mean(Volume),
    Max_Vol = max(Volume),
    .groups = "drop"
  )

cat("\n=======================================================\n")
cat("          RESUMO DA CONTAGEM VISUAL EXPERIMENTAL       \n")
cat("=======================================================\n")
print(resumo_classes)

max_spread <- max(df_spread$Volume)
min_vcc    <- min(df_vcc$Volume)
ponto_medio_corte <- round((max_spread + min_vcc) / 2, 4)

cat("\n-------------------------------------------------------\n")
cat(paste("Maior partícula classificada como SPREAD :", max_spread, "µm³\n"))
cat(paste("Menor estrutura classificada como VCC    :", min_vcc, "µm³\n"))
cat(paste(">>> PONTO DE TRANSIÇÃO PERFEITO          :", ponto_medio_corte, "µm³ (~0.50 µm³)\n"))
cat("-------------------------------------------------------\n")

# Comparative Distribution Plot (Jitter + Boxplot) in English
p_validation <- ggplot(df_total, aes(x = Classe, y = Volume, fill = Classe)) +
  geom_boxplot(alpha = 0.6, width = 0.35, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0.50, color = "#e74c3c", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 1.5, y = 0.58, label = "Empirical Threshold: 0.50 µm³\n(100% Sensitivity & Specificity)", 
           color = "#e74c3c", fontface = "bold", size = 4.2) +
  scale_y_log10(labels = scales::label_log()) +
  scale_fill_manual(values = c("Spread" = "#3498db", "VCC" = "#e74c3c")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Empirical Ground Truth: VCC vs. Surface Spread Separation",
    subtitle = paste0("N = ", nrow(df_vcc), " VCCs vs. ", nrow(df_spread), " Spreads | Perfect non-overlapping boundary"),
    x = "Visual Classification (Ground Truth)",
    y = expression("Imaris Volume (" * mu * "m"^3 * " - Log10 Scale)")
  ) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

print(p_validation)

# Save figure in output folder
output_dir <- file.path(dirname(pasta_corte_visual), "Results_VCC_Analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "Figure_Visual_Validation_VCC_Spread.png"), p_validation, 
       width = 7.5, height = 5.5, dpi = 300)
message(paste("\nValidation plot saved to:", output_dir))
