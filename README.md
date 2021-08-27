# Subcellular-3D-spatially-resolved-transcriptional-density-and-function-of-single-stem-cells

SOFTWARE PACKAGES IN JUPYTER NOTEBOOK, PYTHON .PY, AND MATLAB .M FILES

FOR THE MANUSCRIPT TITLED:

Subcellular 3D spatially resolved transcriptional density and function of single stem cells
Zhou Fang1,#,  Ryan Abramowitz1,2,#, Sungwoong Kim2, Thomas Hu1,3, Mayar Allam1, Levi Wood1,2,3, Ankur Singh1,2, and Ahmet F. Coskun1,4,5*
1 Wallace H. Coulter Department of Biomedical Engineering, Georgia Institute of Technology and Emory University, Atlanta, GA, USA
2 Woodruff School of Mechanical Engineering, Georgia Institute of Technology, Atlanta, GA, USA
3 School of Electrical and Computer Engineering, Georgia Institute of Technology, Atlanta, GA, USA
4 Interdisciplinary Bioengineering Graduate Program, Georgia Institute of Technology, Atlanta, GA, USA
5 Parker H. Petit Institute for Bioengineering and Bioscience, Georgia Institute of Technology, Atlanta, GA 30332
*Corresponding author
#Co-first authors 
Ahmet F. Coskun, Ph.D. (ahmet.coskun@bme.gatech.edu)


The "Single_image_clustering_package" includes code and example dataset to conduct RNA density analysis from projected images.
The code extracts pixels within a pre-defined cell mask, and conducts k-means clustering analysis. k values were determined by finding the "elbow point" of sum of within-cluster error vs k-values plot. The noise within the final clustering results were then masked out via the cell mask. All plots are then saved in the same directory as the input data.

The "Agglomerative_clustering_package" includes code that conducts RNA density analysis through agglomerative clustering algorithm.
The code uses "2D_7gene_bm_msc" and "2D_7gene_uc_msc" dataset, and generates RNA density clustering via agglomerative clustering. The heatmap of cluster centers and visualization image are saved.

The "3D_spatial_clustering_package" includes code and example dataset to conduct RNA density analysis in original stacked images.
The code conducts k-means clustering analysis. k values were determined by finding the "elbow point" of sum of within-cluster error vs k-values plot. The clustering heatmap and clustering results were saved as individual image stack for each cluster.

The "Counting_and_visualization_package" includes code that counts dots in RNA images, and generate visualization in figures.
Python code MakeTXTfiles.py counts dots in images, and generate .txt files containing information regarding the dot counts. "RNAcountsBM.txt" and "RNAcountsUC.txt" files are example output of "MakeTXTfiles.py". Other scripts takes the dot counts .txt file, and generate the visualization presented in each figure.

The "Luminex_data_package" includes code that read raw luminex assay cytokine measurements, and generate heatmap visualization of secreted cytokine concentration.

The "Differentiation_prediction_package" generates differentiation time prediction.
The algorithm reads bulk RNAseq data from MSC to chondrocytes differentiation assay, and generate predictions based on simulated RNA expression level and raw RNA image.
