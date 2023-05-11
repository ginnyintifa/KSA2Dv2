KSA2Dv2

## 1 Introduction 

This package aims to help users identify differentially expressed proteins and phosphorylation sites with significantly increased/reduced levels. 
The main algorithm for 1-dimensional discovery is based on
>"Efron, Bradley, et al. "Empirical Bayes analysis of a microarray experiment." Journal of the American statistical association 96.456 (2001): 1151-1160.
Efron, B., Tibshirani, R., Storey, J. "

We extended the algorithm to 2-dimensional setting to discover joint variation of kinase and substrates (KSA2D) . Later we adapted the tool to analyse unpaired group comparisons in 2D setting (KSA2Dv2).

## 2 Installation and preparation 
KSA2Dv2 can be downloaded and installed in R. Installation requires devtools as a prerequisite:

```{r}
install.packages("devtools")
```
Next, install KSA2D by:

```{r}
library("devtools")
devtools::install_github("ginnyintifa/KSA2Dv2")
library(KSA2Dv2)
```

Please load these packages 

```{r}
library(data.table)
library(dplyr)
library(magrittr)
library(limma)
library(stringr)
library(splines)
library(MASS)
library(scatterplot3d)
library(limma)
library(ggplot2)
library(ggrepel)
library(gridExtra)
```

Suppose we wish to identify differentially expressed proteins and phosphorylation sites with different abundance levels between 2 groups. Firstly we need to prepare 3 input files 


* protein data: Data frame of proteome data, each row is a protein, columns are samples.  
* phosphorylation data: Data frame of site leve phosphoproteome data, each row is a protein, columns are samples
* Files for annotating proteins and mapping kinase-substrate relationships can be downloaded from PhosphoSitePlus or OmniPath

Please notice that there are specific format requirements for both protData and psiteData. 

protData.tsv:

```
Gene_name   Sample1 Sample2 Sample3 ...
A0AVF1  0.8 0.4 -0.2 ...
AASS    1   0.7 -0.5 ...
AATF    -0.5    0.8 NA ...
```

psiteData.tsv

```
Gene_site   SequenceWindow	 Sample1 Sample2 Sample3 ...
A0AUZ9_S714 QEPVHLDsPAIKHQF -0.004  0.3 -0.2
A0AVT1_S743 VGSLTPPsSPKTQRA 0.3 0.1 -0.01
A0FGR8_S676 ASLNKSKsATTTPSG NA  0.03    0.1
```
Please note that these are the default data format of tidied up FragPipe TMT-I results using the [nesviLab_scripts](https://github.com/ginnyintifa/nesviLab_scripts)


## 3 Functions 

### 3.1 kinase substrate information mapping
```ks_map``` mapps the user's proteome and phosphoproteome to the known kinase substrate mapping information 

```{r}
ks_list_sub = ks_map(
  protData_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ncc_ratio_gene_MD_tidy_mapFiile.tsv",
  psiteData_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ncc_normalized_phospho_single-site_MD_tidy_mapFile.tsv",
  ksNetwork_filename = "/Users/ginny/My Drive/nonCCRCC20211020/merge_ks_omni.tsv",
  ks_outputName = "/Users/ginny/My Drive/KSA2D_v2/test_output/KSA2D_ks_map.tsv")
```

Return value of this function is a list ```ks_list_sub``` and can be used in the following main function:

### 3.2 Main effect model 

```balance_sample_purity_bandwidth``` performs the empirical bayes modelling using the resultant lists from the ```ks_map``` function. At the same time, if purity adjustment is wanted, an annotation file should be prepared. In the annotation file, the columna name for the samples IDs should be "caseID" and the column name for purity should be "purity". Samples should be mappable to the column names of the proteome and phospho site data files. 


```{r}
balance_sample_purity_bandwidth(d1_data = ks_list_sub[[1]],
d2_data = ks_list_sub[[2]],
s1_col_name = col_ncc_wgii_high,
s2_col_name = col_ncc_wgii_low,
balance_flag = F,
nna_cutoff = 6,    
working_dir = "/Users/ginny/My Drive/KSA2D_v2/test_output/ncc_105_1/",
compare_name = "ncc_wgii_105H_unfiltered",
permute_time =200,
bandwidth_factor = 1.5,                                   
adjust_purity_flag = T,
annot_file = "/Users/ginny/My Drive/nonCCRCC20211020/clinical/col_annot_meta_20220526.tsv")
```


here ```col_ncc_wgii_high``` and ```col_ncc_wgii_low``` are the sample IDs of the two comparison groups of interest. 

Important input parameters:

* ```nna_cutoff```: the minimun number of non-missing values acceptable for each feature, should be at least 6. 
* ```permute_time```: the number of permutations needed to generate the null distribution, recommend 200.
* ```compare_name```: a tag for the output files.
* ```bandwidth_factor```: a numerical value for smoothing the empirical distribution, the larger the bandwidth, the smoother the distribution gets. Default at 1. 
* ```adjust_purity_flag```: a boolean value indicate if purity needs to be adjusted when comparing the two groups. 

Return value of this function are files deposited to the ```working_dir```. If the sample sizes of the two comparison groups are very different from each other, the comparison will be done in slices. For each slice, you will find the empirical distribution of the overall signal, the fitted null distribution and the calculated differential distribution. A text file (.tsv) detailing the results will be generated as well. 

### 3.3 Prepare for cytoscape results visualization

```ksa2d_result_cytoscape``` post-precess the results and extract kinase-substrate mappings with desired cutoffs and prepare edge and node files for cytoscapte visualization. 

```{r}

ksa2d_result_cytoscape(working_dir = "/Users/ginny/My Drive/KSA2D_v2/test_output/ncc_105_1/",
                       output_dir = "/Users/ginny/My Drive/KSA2D_v2/test_output/ncc_105_1/test1/",
                       freq_cutoff = NULL,
                       fdr_cutoff = 0.1,
                       kinase_abs_fc_cutoff = 0.05,
                       subsite_abs_fc_cutoff = 0.5,
                       ksNetwork_filename = "/Users/ginny/My Drive/nonCCRCC20211020/merge_ks_omni.tsv")

```

Return values of this function are four tsv files. 
* ```merge_all.tsv```: merged kinase-substrate relationships from all the slices. 
* ```merge_sig.tsv```: the significant subset of the above file. 
* ```edge.tsv```: edges extracted from the merge table for cytoscape illustration. 
* ```node.tsv```: nodes extracted from the merge table for cytoscape illustration. 

## 4 Visualization 
Recommended settings in Cytoscape:

* size: degree
* label Font size: degree
* shape: kinase_bl (diamond or ellipse)
* fill color: kinase_fc
* border width: site_num
* border paint: substrate_av_fc
* edge transparentcy: edge_sig_bl
* edge stroke color: edge_dir 
* edge width: edge_num




