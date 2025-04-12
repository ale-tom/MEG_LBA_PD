[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ale-tom/MEG_LBA_PD/HEAD)
# Parkinson’s disease impairs cortical sensori-motor decision-making cascades  

This repository contains code and data to reproduce the figures in the article [“**Parkinson’s disease impairs cortical sensori-motor decision-making cascades**”](https://academic.oup.com/braincomms/article/6/2/fcae065/7628357?login=false) from Tomassini et al., published on Brain Communications (2023) where we combine Magnetoencephalography (MEG), Electroencephalography (EEG) and computational modelling to show that Parkinson's disease disrupts the beta-frequency activity mediating the accumulation of evidence for decision-making, leading to inefficient processing and transfer of perceptual information to the frontal cortex and atypical decision-making. These findings highlight the connection between changes in neural dynamics and visuomotor function in neurodegenerative disease.
The included code and data are to allow the reproduction of the manuscript's results and figures. For the sake of reproducibility, a renv.lock file including metadata and dependencies has been included.

![Graphical_abstract](https://github.com/ale-tom/MEG_LBA_PD/assets/30290119/ae6bf36f-96c1-4458-b189-ca03d90b3bba)


## Repository Structure

Below is the current layout of the repository:

```plaintext
MEG_LBA_PD/
├── Data/                  # Neurophysiological and behavioral data
├── R_analyses/            # R scripts and analysis notebooks for reproducing figures
├── .gitignore             
├── LICENSE                # MIT License file for this repository
├── renv.lock              # Dependency lock file for reproducible R environments
└── runtime.txt            # R runtime environment specification
```

## Prerequisites and Installation
We recommend using [renv](https://rstudio.github.io/renv/) to manage dependencies.
  1. Install R from [CRAN](https://cran.r-project.org/).
  2. Install `renv` and restore the environment:
     ```R
     install.packages("renv")
     renv::restore()
     ```
## Usage
### Reproducing Analyses:
Navigate to the R_analyses/ directory.
Open the main analysis script or R Markdown file and follow the instructions to run analyses and generate figures.
### Using Binder:
Click the Binder badge at the top to launch an interactive session in your web browser and run the analysis without local setup.
