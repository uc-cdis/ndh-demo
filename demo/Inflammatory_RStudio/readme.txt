Steps to run mounted R code in NIAID Data ecosystem RStudio.

1. Login in https://niaiddata.org/ using you user id and password.

2. Navigate to the enviroment that you will run the notebook (eg: https://aids.niaiddata.org/) and login.

3. Navigate to Profile and click Create API key. Then click Download json to save the credential.json file locally.

4. Navigate to Workspace, select the mounted script you would like to run.

5. Upload your credential.json file under /home/rstudio/pd.

6. Double click to open the RStudio. Create directory "/home/rstudio/pd/nb_output/immune". Upload json files download from NDH HIV Classifier App for Post Treatment Controller (PTC), Elite Controller (EC) or Long Term Non Progressor (LTNP) under "/home/rstudio/pd/nb_output/immune". Type setwd("/home/rstudio/inflammatory/") in console to set the working enviroment and run R code "immune_marker.r".


R Studio (InflammatoryR) overview:

Under "/home/rstudio/inflammatory" folder, there are 3 files: "immune_marker.r", "Gen3AuthHelper.R" and "Gen3Submission.R". "Gen3AuthHelper.R" and "Gen3Submission.R" are the R packages developed by Gen3 team to facilitate users interacting with Gen3 API endpoints. "immune_marker.r" will use the two packages to interact with Fence and Peregrine API endpoits in NIAID Data Ecosystem.

"immune_marker.r" queries the structured data from the AIDS environment in the NIAID Data Ecosystem and performs the actual data analysis using data from HIV-CHARLIE project. The R code compares the demographic attributes such as age, race, HIV status, drug use, etc. across 3 cohorts, Post Treatment Controllers (PTC), Elite Controllers (EC), and Long Term Non-Progressers (LTNP). The analysis also compares important inflammatory marker levels such as TNF alpha, IL2, IL4, IFN gamma etc. across the PTC, EC, and LTNP cohorts to reveal inflammation in response to HIV virus infection
