Requirements
------

	The analysis performed in this study can be grouped into four distinct steps: data downloading, data preprocessing, marker regression analysis, allele transition bias analysis, machine learning analysis. The analysis requires R, and python3.

The following python3 packages should be installed:
	* numpy,
	* pandas,
	* bokeh,
	* matplotlib,
	* scipy,
	* sklearn.


Analysis workflow
-----

	1. To download data ./download_data.sh shell script should be done. All the required data will be downloaded into the ./data directory.

	2. Before any analysis can be done the jupyter notebook ./analysis/src/result_generation/do_preprocessing.ipynb should be executed to generate any auxiliary tables/files needed for the analysis. These files will be saved in the ./data/product directory. The subsequent analysis steps can be done in any order.

	3. To perform the marker regression analysis the jupyter notebook MarkerRegressionResultGeneration.ipynb should be executed. The results are ROC curve graphs, and these shall be saved into ./results/tf_ROC.png (de novo TF prediction ROC curve), and ./results/SNP/\*.png (de novo target prediction ROC curves for each SNP).

	4. To perform the allele bias transition analysis the jupyter notebook AlleleTransitionBias.ipynb should be run. The results are ROC curve graphs and these shall be saved into  ./results/BIAS/SNP/\*.png (de novo target prediction ROC curves for each SNP).

	5. To perform the machine learning analysis the jupyter notebook RandomForestForTFandTargets.ipynb should be run. The results are the ./results/rf_importances_snps_joined_and_sorted.csv table, which contains two columns - feature name (SNP id or ORF id) and its importance for the random forest regressor; and the ./results/rf_training_scores.csv

