# NDE Notebooks
All the notebook is in github https://github.com/uc-cdis/ndh-demo/tree/master/demo

There are five notebooks (including python and r notebook), each developed at one phase of NIAID project.

## DAIDS:
DAIDS_notebook has python notebook, the analysis is consist of three parts:
1. Compare the CD4 positive cells and the viral load between the last HIV seronegative timepoint and first seronegative timepoint. Perform box plot and calculate the p value.

2. Compare the CD4 positive cells and the viral load among the HAART treatement negative timepoint and the HAART treatment for one, two and three years. Perform box plot and calculate the p value.

3. Compare the survival curve between HAART negative group and HAART positive group. Plot Kaplan-Meier curve and calculate p value using logrank test.

This notebook is tested and ready for demo.

## DIMD:
DMID_notebook has python and R notebook, this notebook is developed at the second phase of NIAID project. The purpose is to reproduce the result published in “pathogenic influenza viruses and coronaviruses utilize similar and contrasting approaches to control interferon-stimulated gene responses”.

```"nde_dmid_pynb.ipynb"``` query structure data from niaid data ecosystem flu enviroment, and perform data analysis on data from FLU-LHV project. The analysis shows influenza virus titer and influenza virus RNA copy number in different time points post infection in cell line/mouse model, and weight loss post influenza virus infection in mouse model. The aim is to show the successful establishment of the infection model in cell line/mouse and compares the virus virulence in mouse model.

```"nde_dmid_rnb.ipynb"``` query structure data and download object data from niaid data ecosystem flu enviroment, and perform data analysis on data from FLU-LHV project. The notebook analyzed RNA microarray datasets using limma package to compare gene expression for different virus infection model and plot heatmap to show the activation/suppression of type I inferon signal pathway upon virus infection. In addition, the notebook compared type I inferon signal pathway genes at protein level by G test and T test using processed protein expression data from Mass Spec datasets upon virus infection.

This notebook is tested and ready for demo.


## DAIT:

DAIT notebook has python and r notebook, this notebook is developed at the third phase of NIAID project. The purpose is to show the interoperability of analyzing data using data in multiple projects. 

```"nde_dait_pynb.ipynb"``` query structure data from niaid data ecosystem aids enviroment and perform data analysis on data from HIV-CHARLIE and DAIT-immune_controls project. The analysis shows the total cholesterol and high density lipid cholesteroal level across age in female and male groups or in different races. The analysis shows the advantage of performing analysis across projects. Combining data from different projects enriched datasets so that the measured variable is normally distributed to represent the attributes from population.

```"nde_dait_rnb.ipynb"``` query structure data and download object data from niaid data ecosystem microbiome enviroment and perform data analysis on data from DAIT-microbiome project. The analysis used phyloseq package to calculate alpha and beta diversity from the samples in different organ or samples from pregnant woman at different trimester to reveal the complexity of microbacterial in different organ or at different stages of pregnancy.

This notebook is tested and ready for demo.


## TB:
TB notebook has python and r notebook, this notebook is developed at the fifth phase of NIAID project. The python notebook shows the analysis of drug resistance of tuberculosis bacterial by using Ariba and Mykrobe software. The r notebook shows the accuracy of analysis result derived from Ariba and Mykrobe software.

```"nde_tb_pynb.ipynb"``` query structure data and download object file from niaid data ecosystem TB enviroment, and perform data analysis on data from TB-PATRIC project. The notebook run Ariba and Mykrobe software on tuberculosis raw sequencing files and extract the drug resistance phenotype predicted by the softwares. The prediction result will be submitted to TB-PATRIC project in sequencing result node.

```"nde_tb_rnb.ipynb"``` visualized the specificity and sensitivity of predicting the resistance to four first line drugs by Ariba and Mykrobe softwares. It also compares the prediction accuracy between Ariba and Mykrobe softwares.

The r notebook is tested and ready for demo. The python notebook is not ready for demo, the “Run Mykrobe for drug resistance prediction” and “Submission of Ariba and Mykrobe to Sheepdog” need to be fixed. Currently, this work is blocked by the software issue in NDE, the Mykrobe software is not properly installed in NDE. 



## Inflammatory:

Inflammatory r script run in R studio per customer’s request. This r script is developed at the sixth phase of NIAID project. The script compared demographic attribute and inflammatory bio markers across three cohorts (LTNP, PTC and EC). In order to run the R studio, the user should create a directory "/home/rstudio/pd/nb_output/immune". Upload json files download from NDH HIV Classifier App for Post Treatment Controller (PTC), Elite Controller (EC) or Long Term Non Progressor (LTNP) under "/home/rstudio/pd/nb_output/immune". Type setwd("/home/rstudio/inflammatory/") in console to set the working enviroment and run R code "immune_marker.r".

"immune_marker.r" queries the structured data from the AIDS environment in the NIAID Data Ecosystem and performs the actual data analysis using data from HIV-CHARLIE project. The R code compares the demographic attributes such as age, race, HIV status, drug use, etc. across 3 cohorts, Post Treatment Controllers (PTC), Elite Controllers (EC), and Long Term Non-Progressers (LTNP). The analysis also compares important inflammatory marker levels such as TNF alpha, IL2, IL4, IFN gamma etc. across the PTC, EC, and LTNP cohorts to reveal inflammation in response to HIV virus infection

This notebook is tested and ready for demo.
