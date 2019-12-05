Steps to run mounted notebook in NIAID Data Ecosystem.

1. Login in https://niaiddata.org/ using you user id and password.

2. Navigate to the enviroment that you will run the notebook (eg: https://aids.niaiddata.org/) and login.

3. Navigate to Profile and click Create API key. Then click Download json to save the credential.json file locally.

4. Navigate to Workspace, select the mounted notebook you would like to run.

5. Upload your credential.json file under /home/jovyan/pd.

6. Double click to open the notebook. From the navigation bar on the top left, click Run and then select Run All Cells to execute the notebook.


Jupyter Notebook- DMID (Lab Edition) overview:

Under DMID-notebook folder, there are "nde_dmid_pynb.ipynb" python notebook and its dependent library "nde_dmid_function.py", and "nde_dmid_rnb.ipynb" R notebook and its dependent library "nde_dmid_function.r" with a list of inteferon pathway genes "CalU3_TypeI_ISG.txt"

"nde_dmid_pynb.ipynb" query structure data from niaid data ecosystem flu enviroment, and perform data analysis on data from FLU-LHV project. The analysis shows influenza virus titer and influenza virus RNA copy number in different time points post infection in cell line/mouse model, and weight loss post influenza virus infection in mouse model. The aim is to show the successful establishment of the infection model in cell line/mouse and compares the virus virulence in mouse model.

"nde_dmid_rnb.ipynb" query structure data and download object data from niaid data ecosystem flu enviroment, and perform data analysis on data from FLU-LHV project. The notebook analyzed RNA microarray datasets using limma package to compare gene expression for different virus infection model and plot heatmap to show the activation/suppression of type I inferon signal pathway upon virus infection. In addition, the notebook compared type I inferon signal pathway genes at protein level by G test and T test using processed protein expression data from Mass Spec datasets upon virus infection.
