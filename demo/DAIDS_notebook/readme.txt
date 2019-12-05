Steps to run mounted notebook in NIAID Data Ecosystem.

1. Login in https://niaiddata.org/ using you user id and password.

2. Navigate to the enviroment that you will run the notebook (eg: https://aids.niaiddata.org/) and login.

3. Navigate to Profile and click Create API key. Then click Download json to save the credential.json file locally.

4. Navigate to Workspace, select the mounted notebook you would like to run.

5. Upload your credential.json file under /home/jovyan/pd.

6. Double click to open the notebook. From the navigation bar on the top left, click Run and then select Run All Cells to execute the notebook.


Jupyter Notebook-DAIDS(Lab Edition) overview:

Under daids-notebook directory, there are "nde_aids_pynb.ipynb" python notebook and its dependent library "ndh_aids_function.py"

"nde_aids_pynb.ipynb" query structure data from niaid data ecosystem aids enviroment, and perform data analysis on data from CHARLIE project.

The analysis is consist of three parts:
1. Compare the CD4 positive cells and the viral load between the last HIV seronegative timepoint and first seronegative timepoint. Perform box plot and calculate the p value.

2. Compare the CD4 positive cells and the viral load among the HAART treatement negative timepoint and the HAART treatment for one, two and three years. Perform box plot and calculate the p value.

3. Compare the survival curve between HAART negative group and HAART positive group. Plot Kaplan-Meier curve and calculate p value using logrank test.
