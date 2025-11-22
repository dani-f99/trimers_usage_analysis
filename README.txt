-----------------------------
------Daniel-Fridman---------
--- System-Immunology-Lab --- 
------Haifa-University ------
-----------2025--------------
-----------------------------
-----------------------------
--- General-Information -----
This program aim is to analyze and visualize the results of trimers usage in the BCR heavy variable region of multiple immune repotiore
which derived from the ImmuneDB database. 
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Required-Python-Modules -
1.  pandas
2.  numpy
3.  datetime
4.  functools
5.  pathlib 
6.  typing
7.  copy
8.  importlib
8.  itertools
9.  scipy
10. sklearn
11. math
12. seaborn
13. natsort
14. os
15. tqdm
16. pymysql
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Pipeline --------------
follow the steps in trimers_tutorial.ipynb to perform LPA and get pca plot:
1. Configure the MySQL connector (block no.2 - trimers_tutorial.ipynb) or import existing tables
   the required tables (ImmuneBD format) are: sequences, sequence_collapse and sample_metadata.
2. Create trimer class object via the scripts.TrimerLPA.TrimerLPA class method.
3. perform PCA on the signatures matrix with the class method trimerLPA.pca.
4. after performing the PCA, can plot heatmap with the trimerLPA.plot_heatmap method.

follow the steps in the steps in trimers_study_v3.pynb to get information about the trimers properties:
1. Ceate LPA object and perform LPA
2. Visualize the KLDe distances from the domain, analyze the medians via the Trimers_Analyze.medians_info
   and medians_hist class methods.
 
follow the steps in the steps in trimers_motif.pynb to get information about specific found genetic motif:
1. define the motif in the "motif_list" var and follow the steps.
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Subfolders --------------
1. scripts -> Contaning all the relevent code for this program.
2. trimers_data -> Datasets containg the data for the LPA analysis.
3. trimers_figures -> Results Figures.
4. trimers_medians -> Medians tables and figure of the treims usage KLDe distances.
5. trimers_medians_analysis -> Tables of the trimes KLDe analysis tables output.
6. trimers_processed_tables -> Processed sequences tables. 
7. trimers_raw_tables -> Unprocessed sequences tables from the ImmuneBD MySQL server.
8. trimers_results -> Diversity and presance of trimers tables (results).
-----------------------------


