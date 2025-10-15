# Temperature-mediated-supply-drives-organism-size
Temperature-dependent supply and demand drive shrinking body sizes across levels of ecological organisation

# Simulations and Analyses of Monocultures and Polycultures

This repository contains the **R code and simulated data** used in the analyses for the manuscript:  

Temperature-dependent supply and demand drive shrinking body sizes across levels of ecological organisation

The project models **monocultures and polycultures** to analyze community dynamics, diversity, and functional responses.

---

## Repository Contents

### 📊 Simulated Data
- **BEEP_Final-communityniche.csv** – Final simulated community niches.  
- **Beep_Monosimu2024.csv** – Simulation output for monocultures.  
- **BEEP_Poly2024-2.csv** – Simulation output for polycultures.  

These CSV files can be used directly for plotting and visualization with the scripts below.  

---

### 💻 Main Analysis Scripts
- **BEEP_monocultures20251001.R** – Main script for monoculture simulations.  
- **BEEP_polyculture20251001.R** – Main script for polyculture simulations.  

Both scripts rely on **BeepSupplyDemandFunctions.R**, which defines the model and other core simulation functions.  

---

### 📈 Plotting Scripts
- **Monocultures_plot.R** – Generates figures for monoculture simulations.  
  - Includes additional analysis of hump- and U-shaped responses in the parameter space.  
  - Uses helper function from `localextrema2.R`.  

- **Polycultures_plot.R** – Generates figures for polyculture simulations.  

---

### 🛠 Supporting Scripts
- **BeepSupplyDemandFunctions.R** – Functions implementing the supply-demand model.  
- **localextrema2.R** – Function to detect local minima and maxima in response curves.  
- **MatrixSpecAnalysis.R** – Preliminary script used to generate 1,000,000 interaction matrices and classify communities by specialization/generalization (see Methods section of manuscript).  

---

## 🔧 Usage Instructions

1. Install **R**.  
2. Install required R packages:  
   ```r
   install.packages(c("ggplot2", "dplyr", "tidyr", ...))  # complete list here
   ```
3. Run monoculture simulations:  
   ```r
   source("BEEP_monocultures20251001.R")
   ```
4. Run polyculture simulations:  
   ```r
   source("BEEP_polyculture20251001.R")
   ```
5. To reproduce manuscript figures directly from the provided CSVs:  
   ```r
   source("Monocultures_plot.R")
   source("Polycultures_plot.R")
   ```

---

## ⚠️ Notes
- **File naming**: Some systems (Linux/Mac) are case-sensitive — ensure you respect file names exactly (e.g., `BEEP_` vs. `Beep_`).  
- **Reproducibility**: Consider setting a random seed (e.g., `set.seed(123)`) for exact reproducibility of simulations.  
- **Large files**: All large-scale results are provided as CSVs; re-running simulations may take substantial compute time.  

---

## 📖 Citation
If you use this code or data, please cite:  

> [Authors] ([Year]). *Title of manuscript.* Journal. DOI  

A Zenodo DOI for this repository will be provided once archived.  

---

## 👥 Contact
For questions about the code or data, please contact:  
**Dr. Mikael Pontarp** – Department of Biology, Lund University – mikael.pontarp@biol.lu.se  

---
