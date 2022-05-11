### Microarray data analysis

Hepatocarcinogenesis (HCC) is a multi-step process. The aim of this study (Torrecilla et al.,
2017) was to characterize early stages of HCC by defining genes altered in early lesions. For
this reason tissue samples of early stages of HCC represented by dysplastic nodules and small
lesions were analyzed by targeted-sequencing and a single nucleotide polymorphism array.
Genes altered in early lesions were defined as trunk. Their distribution was explored in different
regions of large tumors and in different nodules of the same patient. ​ TERT , ​ ​ TP53 and ​ CTNNB1
mutations were the most frequent trunk events and were ubiquitous across different regions of
the same tumor and between primary and metastatic nodules (Torrecilla et al., 2017).

![Hepatocarcinogenesis_progression:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/Hepatocarcinogenesis_progression.png)

Used dataset
SubSeries: ​ GSE102451​ of the Superseries: GSE98620
The dataset consists of tissue samples from 30 different patients and from each patient
Non-tumor cirrhotic liver (shown in green) and HCC tissue (shown in orange) or Pre-malignant
dysplastic nodule (shown in yellow) (in the case of two from the 29 patients) from cirrhotic liver
have been isolated . One patient was excluded from the analysis because the dataset had only
one sample from him: Non-tumor cirrhotic liver.


Boxplot​ of the expression dataset:

![Boxplot_expression_dataset:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/Boxplot%E2%80%8B_expression_dataset.png)


The arrayQualityMetrics() function was used to visualise initial statistics of the dataset. The
density plot​ shows the distribution of expression values in the expression dataset:

![density_plot:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/density_plot.png)


The dataset downloaded with the function getGEO() was already normalised so the data were
log2 transformed. No gene duplicates and empty gene names were found.
Histogram​ and ​ boxplot ​ of the expression data​ ​ after log2 transformation:

![Histogramlog2transformation:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/Histogram%E2%80%8B_log2transformation.png)

![Boxplotlog2transformation:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/boxplotlog2_transformation.png)

Differential Expression analysis
Since we have samples with dependencies (tissue samples from the same patient) the limma
package was used for the DE Analysis because it allows to make comparisons within and
between subjects. This requires that the person/patient is treated as a random effect. In this
case, the differential expression analysis was performed on the Tumor/ Pre-malignant dysplastic
VS Cirrhotic tissues. The comparison is done within each subject because each person returns
both tissues.
The table returning the top genes was filtered for: p-value < 0.05 and logFC > 1 or logFC < -1
and 208 genes passed the filtering - statistically significant (shown in blue).
Volcano plot​ of the Differential Expression analysis:

![volcanoplot:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/Volcano_plot.png)



Heatmap​ of the 208 statistically significant genes (top genes):

![heatmap:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/Heatmap.png)

QSage plot​ :

![QSage_plot:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/QSage_plot.png)


A simple ​ network was created using GeneMANIA in order to visualise the 208 genes. The
following genes were not clustered with the others: Rho GTPase activating protein 11B (on the
wright), IRX3, ALG1L, ETV3L, KRTCAP3 (bottom).

![network:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/network.png)

The 208 top genes were compared with the top genes outputs from the other members of the
HCC group and a venn diagram was created using VENNY online tool. 28 genes were found
mutually in all 3 dataset. The 28 common genes were visualised in a network created using
geneMANIA. Using DAVID online tool no clusters of the genes could be identified.
Network​ of the 28 common genes:

![network_28:](https://github.com/ourtheol/microarray_analysis_study/blob/main/Pictures/network_28.png)




