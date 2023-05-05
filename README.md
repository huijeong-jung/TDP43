# Analysis of single nuclei RNAseq data to identify TDP-43 and its role in Cryptic Exon inclusion 

All the results from this project have been displayed on the shiny server https://huijeong-jung.shinyapps.io/app-2/

## Files 
- Differential_Expression.Rmd: This file was used to conduct edgeR differential expression analysis 
between the raw counts for genes with corresponding cryptic exons of interest. The raw counts were TMM
normalized then compared between the three different patient diagnoses (control, dementia, and AD) as 
well as between the three Braak stage groups (Braak stage 0-1, 2-4, 5-6). The multiple test corrected
p-values from these differential expression analyses are displayed on top of the figures c and d of the shiny. 

- app.R: This file contains all the code that was utilized to read in, categorize, 
and visualize normalized cryptic exon read counts on an R shiny. The R packages ggplot2 and Shiny were used to visualize the results from cryptic exon 
quantification and the summarized results are 
deployed on the shinyapps server https://huijeong-jung.shinyapps.io/app-2/


# TODO : INCLUDE FILE FOR JUST CRYPTIC EXON 
