Steps to run mounted notebook in NIAID Data Ecosystem.

1. Login in https://niaiddata.org/ using you user id and password.

2. Navigate to the enviroment that you will run the notebook (eg: https://aids.niaiddata.org/) and login.

3. Navigate to Profile and click Create API key. Then click Download json to save the credential.json file locally.

4. Navigate to Workspace, select the mounted notebook you would like to run.

5. Upload your credential.json file under /home/jovyan/pd.

6. Double click to open the notebook. From the navigation bar on the top left, click Run and then select Run All Cells to execute the notebook.

Jupyter Notebook-DAIT(Lab Edition) overview:

Under DAIT-notebook folder, there are "nde_dait_pynb.ipynb" python notebook and its dependent library "nde_dait_function.py", and "nde_dait_rnb.ipynb" R notebook and its dependent library "nde_dait_function.r".

"nde_dait_pynb.ipynb" query structure data from niaid data ecosystem aids enviroment and perform data analysis on data from HIV-CHARLIE and DAIT-immune_controls project. The analysis shows the total cholesterol and high density lipid cholesteroal level across age in female and male groups or in different races. The analysis shows the advantage of performing analysis across projects. Combining data from different projects enriched datasets so that the measured variable is normally distributed to represent the attributes from population.

"nde_dait_rnb.ipynb" query structure data and download object data from niaid data ecosystem microbiome enviroment and perform data analysis on data from DAIT-microbiome project. The analysis used phyloseq package to calculate alpha and beta diversity from the samples in different organ or samples from pregnant woman at different trimester to reveal the complexity of microbacterial in different organ or at different stages of pregnancy.
