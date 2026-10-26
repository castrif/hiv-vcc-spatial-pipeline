# 🔬 3D Spatial Anchoring and Quantitative Phenotyping of HIV-1 Capture and VCC Formation in a Stable U937-CD169 Monocytic Model

Pipeline bioinformático automatizado em R para importação, higienização, ancoragem espacial tridimensional (3D) e classificação fenotípica de compartimentos contendo vírus maduros (**VCCs**) e dispersão viral superficial (**spread**) de HIV-1 a partir de superfícies reconstruídas no Imaris.

---

## 🏛️ Contexto Científico e Filiação Institucional

O desenvolvimento deste pipeline e a caracterização do modelo celular foram realizados no **Laboratório de Biologia Celular do Sistema Imune**, sediado no **Instituto de Ciências Biomédicas IV (ICB-IV)** da **Universidade de São Paulo (USP)**.

### 💡 Racional e Hipótese Científica
> A integrina $\beta_2$ (CD18) atua sustentando um **eixo biomecanotransdutor** essencial para a manutenção da integridade arquitetural, estabilidade e tráfego de compartimentos contendo vírus (**VCCs**) em células mieloides infectadas pelo HIV-1.

### 🧪 Modelo Biológico
* **Modelo Celular Estável:** Células de linhagem monocítica humana (**U937**) manipuladas para expressar de forma estável o receptor **CD169 (Siglec-1)** — principal ligante de captura do HIV-1 —, estabelecendo um modelo padronizado e reprodutível para investigar o sequestro viral.
* **Grupos Genéticos:**
  * **WT** (*Wild-Type* U937-CD169)
  * **CRISPR-NT1** (Controle transduzido não-alvo)
  * **CRISPR-CD18-EX3** e **EX4** (Nocautes para integrina $\beta_2$/CD18)
* **Escopo Estendido:** O pipeline também é aplicado para investigar e comparar o comportamento de invaginação e formação de VCCs em **Células Dendríticas (DCs)**.

---

## 👥 Equipe e Colaboradores

### Laboratório de Biologia Celular do Sistema Imune (ICB-IV / USP):
* **MSc. Felippe Egon Castrignano Camargo**  
  *Desenvolvedor do pipeline e experimentador principal*  
  📧 [felippe.camargo@sou.ufmt.br](mailto:felippe.camargo@sou.ufmt.br)
* **MSc. Kamilla Batista da Silva Souza**  
  📧 [kamillasouza@usp.br](mailto:kamillasouza@usp.br)
* **Profa. Dra. Bruna Cunha de Alencar Bargieri**  
  *Coordenação e Orientação do Projeto*  
  📧 [bruna.alencar@usp.br](mailto:bruna.alencar@usp.br)

### Centro de Facilidades de Apoio à Pesquisa (CEFAP / USP):
* **Dr. Mário Costa Cruz**  
  *CEFAP, Core CONFOCAL*  
  📧 [confocal@usp.br](mailto:confocal@usp.br)
* **Dr. Alexandre Freitas da Silva**  
  *CEFAP, Core GENIAL*  
  📧 [alexandredasilva@usp.br](mailto:alexandredasilva@usp.br)

---

## 📌 Algoritmo de Ancoragem e Classificação

### 1. Reconstrução da Relação Hospedeiro-Patógeno (Ancoragem Euclidiana 3D)
Como o Imaris exporta as superfícies virais sem o identificador da célula hospedeira (`ID_Cel`), o pipeline calcula a **Distância Euclidiana 3D** entre o centro de massa de cada vírus e todas as células presentes no mesmo campo confocal:
$$d = \sqrt{(X_{vir} - X_{cel})^2 + (Y_{vir} - Y_{cel})^2 + (Z_{vir} - Z_{cel})^2}$$
A partícula é associada ao centroide celular mais próximo (`slice_min`).

### 2. Filtro Biofísico de Envelope Celular
Para evitar que vírus livres no meio extracelular sejam computados, aplica-se o filtro:
$$\text{Shortest Distance to Surfaces} \le 0.50\,\mu m$$
* **Distância negativa:** Vírus internalizado no envelope/citoplasma.
* **Distância zero:** Vírus na membrana plasmática.
* **Distância até $+0.50\,\mu m$:** Contato com glicocálice, invaginações e filopódios.
* **Distância $> 0.50\,\mu m$:** Descartado como partícula livre no meio.

### 3. Padrão de Corte de VCC
* **Limiar Empírico Validado (*Ground Truth*):** **$0.50\,\mu m^3$**  
  Calibrado por inspeção visual cega (32 VCCs vs. 692 partículas de spread), resultando em **100% de Sensibilidade** e **100% de Especificidade** ($AUC = 1.000$).
* **Diâmetro Equivalente:** $\approx 1.0\,\mu m$ (empacotamento de 15 a 30 vírions maduros em um único microdomínio confocal).
* **Triangulação Visual:** Os gráficos exibem simultaneamente o corte visual ($0.50\,\mu m^3$), o corte bimodal estatístico e o corte conservador da literatura (**Garcia-Parajo et al.: $2.00\,\mu m^3$**).

---

## 📂 Nomenclatura Taxonômica das Pastas

Para viabilizar o processamento automatizado em lote via `tidyr::separate()`, as pastas das superfícies do Imaris devem seguir o formato de 4 blocos separados por underscores (`_`):

```text
Condicao_Linhagem_Campo_Alvo
```

### Exemplos:
* `WT_U937CD169_Campo01_Celulas` e `WT_U937CD169_Campo01_Virus`
* `CRISPRV2-NT1_U937CD169_Campo01_Celulas` e `CRISPRV2-NT1_U937CD169_Campo01_Virus`
* `CRISPRV2-CD18-EX3_U937CD169_Campo01_Celulas` e `CRISPRV2-CD18-EX3_U937CD169_Campo01_Virus`
* `CRISPRV2-CD18-EX4_U937CD169_Campo01_Celulas` e `CRISPRV2-CD18-EX4_U937CD169_Campo01_Virus`

> **Nota:** O script trata automaticamente acentuações (`Células`/`Celulas`, `Vírus`/`Virus`) e remove sufixos injetados pelo Imaris (ex: `_Statistics`).

---

## 📊 Métricas do Imaris Empregadas

O pipeline requer estritamente 3 métricas do Imaris:
1. **`Position.csv`**: Coordenadas $(X, Y, Z)$ para o cálculo da distância euclidiana.
2. **`Shortest_Distance_to_Surfaces...csv`**: Distância real à membrana celular (filtro $\le 0.50\,\mu m$).
3. **`Volume.csv`**: Volume 3D reconstruído para classificar VCC ($\ge 0.50\,\mu m^3$) vs Spread ($< 0.50\,\mu m^3$) e computar a carga viral celular total.

---

## 🚀 Como Executar

### Pré-requisitos
* R $\ge 4.0$
* Pacotes: `tidyverse`, `rstudioapi`, `ggpubr`, `rstatix`, `scales` (verificados e instalados automaticamente).

### Execução no RStudio
1. Abra o arquivo de projeto: `Analises Virus e VCC.Rproj`.
2. Abra e execute o script: `pipeline_vcc_hiv.R`.
3. Selecione a pasta mãe quando a janela do explorador abrir.

---

## 📈 Saídas Geradas em Inglês (`Results_VCC_Analysis/`)

Ao final da execução, o pipeline salva automaticamente:

### 1. Tabelas de Dados (`.csv`):
* `Table_Cellular_Phenotypic_Profile_HIV.csv`: Perfil fenotípico por célula individual (**VCC+**, **Spread** ou **Uninfected**), contendo contagem de partículas e volume total de carga viral.
* `Table_Associated_Anchored_Viruses.csv`: Partículas virais validadas com identificador de célula hospedeira (`Cell_ID`).
* `Table_Statistical_Summary_Volumes.csv`: Estatística descritiva e percentis populacionais (P75, P90, P95).
* `Table_Pairwise_Statistics_Volume.csv`: Testes de Kruskal-Wallis e Wilcoxon com ajuste de FDR.

### 2. Figuras em Resolução de Publicação (300 DPI, `.png`):
* `Figure_01_Phenotypic_Proportions.png`: Proporções fenotípicas populacionais.
* `Figure_02_Compartment_Volume_Trend.png`: Violin plots com médias conectadas por linha de tendência e os 3 marcos de corte.
* `Figure_03_Cellular_Viral_Load.png`: Carga viral total por célula infectada.
* `Figure_04_ECDF_Volume_Distribution.png`: Curva acumulada de distribuição de volumes (ECDF).

---

## 🔬 Scripts Complementares

* `validacao_corte_visual.R`: Executa a comparação direta entre o *Ground Truth* visual de VCC vs. Spread e gera a figura `Figure_Visual_Validation_VCC_Spread.png`.
* `pipeline_calibracao_corte_vcc.R`: Realiza a modelagem de densidade bimodal de Kernel e gera a figura `Figure_Calibration_VCC_Threshold.png` e a tabela `Table_Calibration_Percentiles_Volumes.csv`.

---

## 📜 Licença e Reprodutibilidade
O código possui semente estocástica fixa (`set.seed(42)`) para garantir reprodutibilidade exata de jittering em qualquer sistema operacional (macOS, Windows ou Linux).
