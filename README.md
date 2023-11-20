# Artificial Intelligence for Unified Analysis of Historical and Landscape Influences on Genetic Diversity

## Overview

This repository contains scripts and resources related to the paper titled "Artificial Intelligence Enables Unified Analysis of Historical and Landscape Influences on Genetic Diversity." The goal of this project is to provide a unified platform for analyzing the influences of historical and landscape factors on genetic diversity, leveraging artificial intelligence techniques.


## Getting Started

To initiate the analysis follow these steps:

1. **Locality and Demographic Information Extraction:**
   - Run `Species_information_L_troglodytes.R` to extract information about localities and demographic history for the species Leptodactylus troglodytes.
   - Run `Species_information_R_granulosa.R` to extract information about localities and demographic history for the species Rhinella granulosa.

2. **Genomic Image Creation:**
   - Execute `Observed_image_L_troglodytes.R` to generate images from genomic information for the species Leptodactylus troglodytes.
   - Execute `Observed_image_R_granulosa.R` to generate images from genomic information for the species Rhinella granulosa.

3. **Simulated Dataset Generation and HPC Submission Scripts:**
   - Utilize `Main_script_script_L_troglodytes.R` to generate scripts for running simulated datasets and creating PBS scripts for HPC submission for the species Leptodactylus troglodytes.
   - Utilize `Main_script_script_R_granulosa.R` for the same purpose but for the species Rhinella granulosa.

4. **Convolutional Neural Network (CNN) Model Execution:**
   - Run `CNN_script_L_troglodytes.py` to execute the CNN model for the species Leptodactylus troglodytes.
   - Run `CNN_script_R_granulosa.py` to execute the CNN model for the species Rhinella granulosa.

Make sure to review and customize the input parameters within each script according to your specific requirements. Additionally, ensure that all dependencies are installed before running the scripts by using the provided requirements files or installing them manually.
