# DeconICA

### This is a package for TranscRiptome deconvolution using unsupervised blind source deconvolution methods

Functions : 

* **run_fastica**  : is a wrapper ove fastica that runs by default 100 components  : status = works 
* * when there is less than 120 columns in the data, number of components will be selected by Keiser rule
  * PCA for n = ncol is performed every time so that dimensions match
  * if with.names option selected row names (genes) are exported in $names if duplicated genes with higher variance are kept
  * save sample names
* **correlate_metagenes** : function that  correlates components with metagenes (or ranked list)
* * gives out correlation matrix between metagenes and components
  * there is a possibility to delete data of certain threshold from ICA matrix
* **assign_metagenes** : selects the reciprocal matches between metagenes and ICs
* * examples missing
  * tests missing
* **identify_immune_ic** defines components possibly of cell types 
* * examples missing
  * tests missing
* **get_enrichment** : fisher test with immgene db : status = working progress
* * examples missing
  * tests missing
  * still results are quite fishy - hard to adjust threshold for number of genes, poor overlap with signatures
* **cell_voting_immgene** : among n top results of gene enrichment, the code counts percentage of given cell type
* * examples missing
  * tests missing
* **helpers** : 
  * .orient_components 
  * .center_rowmeans
  * .rowVars
  * .cumVar
  * .remove_duplicates
  * .intersect.genes
  * .corr_matrix
  * .verify.n


* make_coeff_matrix : status = none
  * transform ranks to counts (more or less) 
  * how many genes select for the final matrix
* compute_FEV : status = experimental
  * for all assigned components


* import_gmt status  done
* import_rnk

PLOTS

lolypop plot for signal plot

tables

This would be the very basic first set of features with some plotting options added 

The next step would be to

* estimate immune cell /tumor absolute fraction
* apply a solver (many options to test, first the one used by competitors)  and find immune cell proportions and confront it with the competition
* test on methylome