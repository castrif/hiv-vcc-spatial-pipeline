# ==============================================================================
# PIPELINE: Quantificação Espacial 3D e Fenotipagem de VCCs (HIV-1)
#
# Autoria e Desenvolvimento:
#   MSc. Felippe Egon Castrignano Camargo (felippe.camargo@sou.ufmt.br)
#   MSc. Kamilla Batista da Silva Souza (kamillasouza@usp.br)
#   Profa. Dra. Bruna Cunha de Alencar Bargieri (bruna.alencar@usp.br)
#   Laboratório de Biologia Celular do Sistema Imune
#   Instituto de Ciências Biomédicas IV (ICB-IV), Universidade de São Paulo (USP)
#
# Colaboradores (Facilities):
#   Dr. Mário Costa Cruz (confocal@usp.br) - CEFAP, Core CONFOCAL (USP)
#   Dr. Alexandre Freitas da Silva (alexandredasilva@usp.br) - CEFAP, Core GENIAL (USP)
#
# Hipótese Biológica:
#   O CD18 (integrina beta-2) sustenta um eixo biomecanotransdutor essencial para
#   a manutenção da integridade de compartimentos virais (VCCs) em células mieloides.
#
# Modelo Celular:
#   Linhagem monocítica U937 expressando estavelmente CD169 (Siglec-1)
#   WT vs. CRISPR-NT1 vs. CRISPR-CD18 (EX3 e EX4) | Extensível a Células Dendríticas
# ==============================================================================

# 🔬 3D Spatial Anchoring and Quantitative Phenotyping of HIV-1 Capture and VCC Formation in a Stable U937-CD169 Monocytic Model

Pipeline bioinformático automatizado em R para importação, higienização, ancoragem espacial tridimensional (3D) e classificação fenotípica de compartimentos contendo vírus maduros (**VCCs**) e dispersão viral superficial (**spread**) de HIV-1 a partir de dados multidimensionais exportados do Imaris.

---

## 🏛️ Contexto Científico e Filiação Institucional

O desenvolvimento deste pipeline e a caracterização do modelo celular foram realizados no **Laboratório de Biologia Celular do Sistema Imune**, sediado no **Instituto de Ciências Biomédicas IV (ICB-IV)** da **Universidade de São Paulo (USP)**.

### 💡 Racional e Hipótese Científica
> A integrina β₂ (CD18) atua sustentando um **eixo biomecanotransdutor** essencial para a manutenção da integridade arquitetural, estabilidade e tráfego de compartimentos contendo vírus (**VCCs**) em células mieloides infectadas pelo HIV-1.

### 🧪 Modelo Biológico
* **Modelo Celular Estável:** Células de linhagem monocítica humana (**U937**) manipuladas para expressar de forma estável o receptor **CD169 (Siglec-1)** — principal ligante de captura do HIV-1 —, estabelecendo um modelo padronizado e reprodutível para investigar a captura viral.
* **Grupos Genéticos:**
  * **WT** (*Wild-Type* U937-CD169)
  * **CRISPR-NT1** (Controle transduzido não-alvo)
  * **CRISPR-CD18-EX3** e **EX4** (Nocautes para integrina β₂/CD18)
* **Escopo Estendido:** O pipeline também é aplicado para investigar e comparar o comportamento de invaginação e formação de VCCs em **Células Dendríticas (DCs)**.

---

## 👥 Equipe e Colaboradores

### Laboratório de Biologia Celular do Sistema Imune (ICB-IV / USP):
* **MSc. Felippe Egon Castrignano Camargo**
  *Desenvolvedor do pipeline e experimentador principal*
  📧 felippe.camargo@sou.ufmt.br
* **MSc. Kamilla Batista da Silva Souza**
  📧 kamillasouza@usp.br
* **Profa. Dra. Bruna Cunha de Alencar Bargieri**
  *Coordenação e Orientação do Projeto*
  📧 bruna.alencar@usp.br

### Centro de Facilidades de Apoio à Pesquisa (CEFAP / USP):
* **Dr. Mário Costa Cruz**
  *Core CONFOCAL*
  📧 confocal@usp.br
* **Dr. Alexandre Freitas da Silva**
  *Core GENIAL*
  📧 alexandredasilva@usp.br

---

## 📌 Algoritmo de Ancoragem e Classificação

### 1. Reconstrução da Relação Hospedeiro-Patógeno (Ancoragem Euclidiana 3D)
Como o Imaris exporta as superfícies virais sem o identificador da célula hospedeira (ID_Cel), o pipeline calcula a **Distância Euclidiana 3D** entre o centro de massa de cada vírus e todas as células presentes no mesmo campo confocal:

  d = √[(X_vir - X_cel)² + (Y_vir - Y_cel)² + (Z_vir - Z_cel)²]

A partícula é associada ao centroide celular mais próximo (slice_min).

### 2. Filtro Biofísico de Envelope Celular
Para evitar que vírus livres no meio extracelular sejam computados, aplica-se o filtro:

  Shortest Distance to Surfaces ≤ 0.50 μm

* **Distância negativa:** Vírus internalizado no envelope/citoplasma.
* **Distância zero:** Vírus na membrana plasmática.
* **Distância até +0.50 μm:** Contato com glicocálice, invaginações e filopódios.
* **Distância > 0.50 μm:** Descartado como partícula livre no meio.

### 3. Padrão de Corte de VCC
* **Limiar Empírico Validado (*Ground Truth*):** **0.50 μm³**
  Calibrado por inspeção visual cega (32 VCCs vs. 692 partículas de spread), resultando em **100% de Sensibilidade** e **100% de Especificidade** (AUC = 1.000).
* **Diâmetro Equivalente:** ≈ 1.0 μm (empacotamento de 15 a 30 vírions maduros em um único microdomínio confocal).
* **Triangulação Visual:** Os gráficos exibem simultaneamente o corte visual (0.50 μm³), o corte bimodal estatístico e o corte conservador da literatura (**Garcia-Parajo et al.: 2.00 μm³**).

---

## 📂 Nomenclatura Taxonômica das Pastas

Para viabilizar o processamento automatizado em lote via `tidyr::separate()`, as pastas das superfícies do Imaris devem seguir o formato de 4 blocos separados por underscores (`_`):

  Condicao_Linhagem_Campo_Alvo

### Exemplos:

* WT_U937CD169_Campo01_Celulas e WT_U937CD169_Campo01_Virus
* CRISPRV2-NT1_U937CD169_Campo01_Celulas e CRISPRV2-NT1_U937CD169_Campo01_Virus
* CRISPRV2-CD18-EX3_U937CD169_Campo01_Celulas e CRISPRV2-CD18-EX3_U937CD169_Campo01_Virus
* CRISPRV2-CD18-EX4_U937CD169_Campo01_Celulas e CRISPRV2-CD18-EX4_U937CD169_Campo01_Virus

---

## 📊 Métricas do Imaris Empregadas

O pipeline requer estritamente 3 métricas do Imaris:

1. **Position.csv:** Coordenadas (X, Y, Z) para o cálculo da distância euclidiana.
2. **Shortest_Distance_to_Surfaces...csv:** Distância real à membrana celular (filtro ≤ 0.50 μm).
3. **Volume.csv:** Volume 3D reconstruído para classificar VCC (≥ 0.50 μm³) vs Spread (< 0.50 μm³) e computar a carga viral celular total.

---

## 🚀 Como Executar

### Pré-requisitos

* R ≥ 4.0
* Pacotes: tidyverse, rstudioapi, ggpubr, rstatix, scales (verificados e instalados automaticamente).

### Execução no RStudio

1. Abra o arquivo de projeto: `Analises Virus e VCC.Rproj`.
2. Abra e execute o script: `pipeline_vcc_hiv.R`.
3. Selecione a pasta mãe quando a janela do explorador abrir.

---

## 📈 Resultados Gerados (Resultados_Analise_VCC/)

* **Tabela_Perfil_Celular_HIV.csv:** Perfil fenotípico por célula individual (VCC+, Spread ou Não Infectada).
* **Tabela_Virus_Associados_Ancorados.csv:** Partículas virais validadas com identificador de célula hospedeira (ID_Cel).
* **Tabela_Sumario_Estatistico_Volumes.csv:** Estatística descritiva e percentis populacionais (P75, P90, P95).
* **Estatistica_Pairwise_Volume.csv:** Testes de Kruskal-Wallis e Wilcoxon com ajuste de FDR.
* **Figuras em 300 DPI (.png):**
  * **Plot_01_Proporcao_Fenotipos.png:** Proporções fenotípicas populacionais.
  * **Plot_02_Volume_VCC_Tendencia.png:** Violin plots com médias conectadas por linha de tendência e os 3 marcos de corte.
  * **Plot_03_Carga_Celular.png:** Carga viral total por célula infectada.
  * **Plot_04_Curva_Tendencia_ECDF.png:** Curva acumulada de distribuição de volumes (ECDF).
