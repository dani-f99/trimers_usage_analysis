-----------------------------
------Daniel-Fridman---------
--- System-Immunology-Lab --- 
------Haifa-University ------
-----------2025--------------
-----------------------------
-----------------------------
--- General-Information -----
This program aim is to create, preprocess, analyze and visualize trimers usage across the BCR heavy chain region
of multiple immune repotiores. To do so it utilize the LPA algorithm (see bottom LPA section).
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Required-Python-Modules -
1. pandas
2. numpy 
3. scipy
4. sklearn
5. seaborn
6. natsort
7. tqdm
8. pymysql
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Guide -------------------
See 'trimers_usage_analysis.ipynb' example notebook
1. Trimer usage input file to be used with the LPA algorithm
   - scripts.TrimerSQL -> contating scripts which can be used to import tables from MySQL server. 
                          used to import the labs ImmuneDB tables from our MySQL serever. Please
						  complete the input example list to use (in the example noteb   
   - If the tables were already exported we skip the MySQL import and go ahead and create the trimer 
     input file via the scripts.TrimerCreate.Trimer function. see documentation and example in the 
	 notebook.
	 
2. Trimer usage LPA results visualization
   - Use the file that was created in the previous step to with the scripts.TrimerLPA.TrimerLPA class:
     a. perform PCA on the KLDe distances with the class method pca applied on the TrimerLPA object.
	 b. Visualize the PC clustring via the code in cell block 7.
	
3. Trimer usage analysis
   - In the notebook we take into consideration the results from our prevoius analysis (diversity analysis)
   - We use the scripts.TrimerAnalyze.TrimerAnalysis class to both study the diversity and presence of trimers
     across the different immune repotiore.
   - We initiate the medians_info class and use it's visualization methods plot_median_hist and plot_median_presence
     to see the distribution of trimers between antigen specific sub-repotiores.

4. Genetic motif 
   - In the last step of this pipeline we are looking for trimers that answer two conditions:
     a. They are in the extreamly high or low median distance from the domain
	 b. They are common, meaning their diveristy score is at least 1.
   - The notebook has preprocessing step followed by trimers visualization. The code can be modified as needed.
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Subfolders --------------
1. scripts -> Containing all the relevent code for this program.
2. trimers_data -> Datasets containg the data for the LPA analysis.
3. trimers_figures -> Results Figures.
4. trimers_medians -> Medians tables and figure of the treims usage KLDe distances.
5. trimers_medians_analysis -> Tables of the trimes KLDe analysis tables output.
6. trimers_processed_tables -> Processed sequences tables. 
7. trimers_raw_tables -> Unprocessed sequences tables from the ImmuneBD MySQL server.
8. trimers_results -> Diversity and presance of trimers tables (results).
-----------------------------
-----------------------------
-----------------------------
-----------------------------
---LPA-----------------------
Publication: https://link.springer.com/article/10.1007/s11257-021-09295-7
GitHub: https://github.com/ScanLab-ossi/LPA
-----------------------------

