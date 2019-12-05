Steps to run mounted R code in NIAID Data ecosystem RStudio.

1. Login in https://niaiddata.org/ using you user id and password.

2. Navigate to the enviroment that you will run the notebook (eg: https://aids.niaiddata.org/) and login.

3. Navigate to Profile and click Create API key. Then click Download json to save the credential.json file locally.

4. Navigate to Workspace, select the mounted script you would like to run.

5. Upload your credential.json file under /home/rstudio/pd.

6. Double click to open the RStudio. Create directory "/home/rstudio/pd/nb_output/immune". Upload json files download from NDH HIV Classifier App for Post Treatment Controller (PTC), Elite Controller (EC) or Long Term Non Progressor (LTNP) under "/home/rstudio/pd/nb_output/immune". Type setwd("/home/rstudio/inflammatory/") in console to set the working enviroment and run R code "immune_marker.R".


R Studio (InflammatoryR) overview:

Under "/home/rstudio/inflammatory" folder, there are "immune_marker.r", "Gen3AuthHelper.R" and "Gen3Submission.R". "Gen3AuthHelper.R" and "Gen3Submission.R" are the R packages developed by Gen3 team to facilitate users interacting with Gen3 API endpoints. "immune_marker.R" will use the two packages to interact with Fence and Peregrine API endpoits in NIAID data ecosystem.

"immune_marker.r" query structure data from niaid data ecosystem aids enviroment, and perform data analysis on data from HIV-CHARLIE project. The R code compares the demographic attributes such as age, race, hcv statas, drug use, etc. among PTC, EC and LTNP groups. The analysis also compares important inflammatory markers level such as TNF alpha, IL2, IL4, IFN gamma etc. amonng PTC, EC and LTNP groups to reveal the inflammation in response to HIV virus infection.
