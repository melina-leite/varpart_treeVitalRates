# Major axes of variation in tree demography across global forests

 Analysis on variance partitioning on tree vital rates. 
 
 **Mansucript in preparation.**

This repository uses Git LFS (Large Files Storage) to deal with the large files (>100Mb) of the example data analysis. If you wish to run the example workflow make sure you have (git lfs installed)[https://git-lfs.github.com/] in your computer before cloning or downloading the repository.

# Data availability and example of workflow

The forest census data are available from the ForestGEO network. For some of the sites, the data is publicly available at https://forestgeo.si.edu/explore-data. Restrictions apply, however, to the availability of the data from other sites, which were used under license for the current study, and so are not publicly available. Data are however available from the authors upon reasonable request and with permission of the principal investigators of the ForestGEO sites.

To illustrate the analysis workflow for data cleaning, preparation and modeling we used the Barro Colorado Islando (BCI) forest plot, which data is avaliable at the Dryad repository:

> Condit, Richard et al. (2019), Complete data from the Barro Colorado 50-ha plot: 423617 trees, 35 years, Dryad, Dataset, https://doi.org/10.15146/5xcp-0d46.

The example workflow is separate in the `workflow_example` directory.

What we present here are the results of the models for growth mortality and recruitment applied for all the 21 forest plots and the analysis using thse results.

The `data` directory contain information for the forest plots, as the table of classification of the species in rare and common, mean density, richness, etc.

# Models outputs

The models outputs are in the `models_outpus` directory separated by the models with (`models_time`) and without (`models_notime`) temporal organizing principles (OPs). In each, here are the models applied to the entire dataset (`all`), excluding rare species (`exclude`) or regrouping rare species (`regroup`). Inside each type of analysis there are the results for each vital rate: growth (`grow`), mortality (`mort`) and recruitment (`rec`). 


Models results are organized in 3 subfolders:

- `summary`: the summary ouptus `summary`function applied to the model objet

- `table`: the .Rdata table for the summary output with the estimates, and information of the model. This is the data that will be used for the subsequente analysis and figures.

- `traceplots`: pdf files for the posterior distribution of each parameter and traceplots showing the mixture of the 3 chains.

Given the very large sizes, models' objects were not saved and stored in the repository.


# Analysis

The scripts for the subsequent analysis using the models outputs are numbered from 0 to 6 and are available as html files with the code and results.The original .Rmd files are in the `scripts` directory.


# `renv` environment

This repository uses `renv` package to improve reproducibility. It contains the `renv` folder with the R library for all the packages used at the current version of publication.


