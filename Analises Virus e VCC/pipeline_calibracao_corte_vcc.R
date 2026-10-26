# ==============================================================================
# PIPELINE: Calibração Estatística e Biofísica do Ponto de Corte para VCC
# Projeto: Quantificação e Fenotipagem de VCCs (HIV-1) em Macrófagos U937
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstudioapi)
  library(scales)
})

message(">>> [1/4] Localizando dados para calibração...")

# Tenta carregar a tabela já processada pelo pipeline principal
pasta_padrao <- getwd()
caminho_tabela_existente <- file.path(pasta_padrao, "Resultados_Analise_VCC", "Tabela_Virus_Associados_Ancorados.csv")

if (!file.exists(caminho_tabela_existente)) {
  # Se não estiver no diretório local, busca pelo seletor de pasta
  if (rstudioapi::isAvailable()) {
    pasta_selecionada <- rstudioapi::selectDirectory(caption = "Selecione a PASTA MÃE ou a pasta 'Resultados_Analise_VCC'")
    if (dir.exists(file.path(pasta_selecionada, "Resultados_Analise_VCC"))) {
      caminho_tabela_existente <- file.path(pasta_selecionada, "Resultados_Analise_VCC", "Tabela_Virus_Associados_Ancorados.csv")
    } else if (file.exists(file.path(pasta_selecionada, "Tabela_Virus_Associados_Ancorados.csv"))) {
      caminho_tabela_existente <- file.path(pasta_selecionada, "Tabela_Virus_Associados_Ancorados.csv")
    }
  }
}

if (file.exists(caminho_tabela_existente)) {
  message(paste("Carregando dados processados de:", caminho_tabela_existente))
  df_virus <- read.csv(caminho_tabela_existente, stringsAsFactors = FALSE)
} else if (exists("virus_associados")) {
  message("Utilizando objeto 'virus_associados' presente na memória do RStudio.")
  df_virus <- virus_associados
} else {
  stop("ERRO: Execute primeiro o 'pipeline_vcc_hiv.R' para gerar a tabela de partículas ancoradas ou selecione a pasta correta.")
}

# ==============================================================================
# 2. ANÁLISE DESCRITIVA POPULACIONAL E PERCENTIS
# ==============================================================================
message("\n>>> [2/4] Calculando distribuição e percentis dos volumes virais...")

vols <- df_virus$Vol_Vir[df_virus$Vol_Vir > 0 & !is.na(df_virus$Vol_Vir)]

tabela_percentis <- data.frame(
  Metrica = c("Mínimo", "Mediana (P50)", "P75", "P85", "P90", "P95", "P97.5", "P99", "Máximo"),
  Volume_um3 = c(
    round(min(vols), 4),
    round(quantile(vols, 0.50), 4),
    round(quantile(vols, 0.75), 4),
    round(quantile(vols, 0.85), 4),
    round(quantile(vols, 0.90), 4),
    round(quantile(vols, 0.95), 4),
    round(quantile(vols, 0.975), 4),
    round(quantile(vols, 0.99), 4),
    round(max(vols), 4)
  ),
  Diametro_Equivalente_um = round(2 * ((3 * c(
    min(vols),
    quantile(vols, 0.50),
    quantile(vols, 0.75),
    quantile(vols, 0.85),
    quantile(vols, 0.90),
    quantile(vols, 0.95),
    quantile(vols, 0.975),
    quantile(vols, 0.99),
    max(vols)
  )) / (4 * pi))^(1/3), 3)
)

print(tabela_percentis)

# ==============================================================================
# 3. DETECÇÃO MATEMÁTICA DO PONTO DE INFLEXÃO (VALE BIMODAL DE KERNEL)
# ==============================================================================
message("\n>>> [3/4] Estimando curva contínua de densidade (Kernel Density Estimation)...")

log_vols <- log10(vols)
dens <- density(log_vols, n = 1024, adjust = 1.2)

df_dens <- data.frame(
  log_vol     = dens$x,
  densidade   = dens$y,
  Volume_Real = 10^(dens$x)
)

# Identifica mínimos locais (vales) e máximos locais (picos)
diff_sinal <- diff(sign(diff(dens$y)))
vales_indices <- which(diff_sinal == 2) + 1
picos_indices <- which(diff_sinal == -2) + 1

# Filtra vales que estejam na faixa plausível de transição de VCC (> 0.2 µm³)
vales_candidatos <- df_dens$Volume_Real[vales_indices]
vales_filtrados  <- vales_candidatos[vales_candidatos >= 0.2 & vales_candidatos <= 3.0]

if (length(vales_filtrados) > 0) {
  corte_bimodal <- round(vales_filtrados[1], 3)
  tipo_corte <- "Vale Bimodal (Mínimo Local de Densidade)"
} else {
  # Se a distribuição for contínua decrescente, adota o Percentil 90
  corte_bimodal <- round(quantile(vols, 0.90), 3)
  tipo_corte <- "Percentil 90 (Filtro Top 10% de Maiores Compartimentos)"
}
corte_visual_usuario <- 0.500
corte_garcia_parajo  <- 2.000

cat("\n=======================================================\n")
cat("          DIAGNÓSTICO E CORTES COMPARATIVOS            \n")
cat("=======================================================\n")
cat(paste("1. Corte Visual Empírico (Ground Truth) :", corte_visual_usuario, "µm³\n"))
cat(paste("2. Corte Estatístico Bimodal            :", corte_bimodal, "µm³ (Baseado em:", tipo_corte, ")\n"))
cat(paste("3. Literatura (Garcia-Parajo et al.)    :", corte_garcia_parajo, "µm³\n"))
cat("=======================================================\n\n")

# ==============================================================================
# 4. PLOTAGEM DO GRÁFICO DE CALIBRAÇÃO COM OS 3 CORTES
# ==============================================================================
message(">>> [4/4] Generating comparative calibration plot...")

df_lines <- data.frame(
  Cutoff = c(
    paste0("Empirical Visual Ground Truth (", corte_visual_usuario, " µm³)"),
    paste0("Bimodal Kernel Cutoff (", corte_bimodal, " µm³)"),
    paste0("Literature Reference / Garcia-Parajo et al. (", corte_garcia_parajo, " µm³)")
  ),
  Value = c(corte_visual_usuario, corte_bimodal, corte_garcia_parajo),
  Color = c("#27ae60", "#e74c3c", "#2980b9"),
  Type  = c("solid", "dashed", "dotdash")
)

p_calibration <- ggplot(df_dens, aes(x = Volume_Real, y = densidade)) +
  geom_area(fill = "#3498db", alpha = 0.2) +
  geom_line(color = "#2c3e50", linewidth = 1.1) +
  geom_vline(data = df_lines, aes(xintercept = Value, color = Cutoff, linetype = Cutoff), 
             linewidth = 0.9, show.legend = TRUE) +
  scale_color_manual(
    name = "Reference Thresholds:",
    values = setNames(df_lines$Color, df_lines$Cutoff)
  ) +
  scale_linetype_manual(
    name = "Reference Thresholds:",
    values = setNames(df_lines$Type, df_lines$Cutoff)
  ) +
  scale_x_log10(
    labels = scales::label_log(),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Viral Volume Density Distribution and VCC Calibration",
    subtitle = "Single virions / spread (left peak) vs. sequestered VCCs (right tail)",
    x = expression("Viral Volume (" * mu * "m"^3 * " - Log10 Scale)"),
    y = "Kernel Population Density"
  ) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 11, color = "#444444"),
    legend.position    = "bottom",
    legend.direction   = "vertical",
    legend.box         = "vertical",
    legend.title       = element_text(face = "bold", size = 11),
    legend.text        = element_text(size = 10.5),
    legend.background  = element_rect(fill = "#f8f9fa", color = "#dcdde1", linewidth = 0.4),
    plot.margin        = margin(t = 15, r = 25, b = 15, l = 15)
  )

# Display plot in RStudio
print(p_calibration)

# Save copy in output directory
output_dir <- dirname(caminho_tabela_existente)
if (dir.exists(output_dir)) {
  write.csv(tabela_percentis, file.path(output_dir, "Table_Calibration_Percentiles_Volumes.csv"), row.names = FALSE)
  ggsave(file.path(output_dir, "Figure_Calibration_VCC_Threshold.png"), p_calibration, 
         width = 8.5, height = 6.0, dpi = 300)
  message(paste("\nCalibration plot and percentiles table saved to:", output_dir))
}
