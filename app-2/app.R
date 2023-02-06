library('edgeR')
library('dplyr')
library(shiny)
library(ggplot2)
library(DT)
library(ggbeeswarm)
library(ggpubr)
library(gridExtra)
library(tidyverse)

## APP-2

crExon_geneID <- read.csv('CrypticExon_geneID.csv')
cellTypes = c('Astrocytes', 'Endothelial', 'GABAergic', 'Glutamatergic', 'Microglia', 'OPCs', 'Oligodendrocytes', 'Pericytes')

# metadata contains information regarding patients that provided the samples
metadata <- read.table('metadata.txt', sep='\t', header=TRUE)

# deep white matter patient ID 
dwm <- metadata[metadata[,'tissue']=='Deep white matter', 'sample_id']
# temporal cortex patient ID 
tc <- metadata[metadata[,'tissue']=='Temporal cortex', 'sample_id']

# Braak stage dependent patient ID 
group1 <- metadata[metadata$braak_stage == 0 | metadata$braak_stage==1, 'sample_id']
group2 <- metadata[metadata$braak_stage == 2 | metadata$braak_stage == 3 | metadata$braak_stage == 4, 
                   'sample_id']
group3 <- metadata[metadata$braak_stage == 5 | metadata$braak_stage == 6, 'sample_id']

# Diagnosis dependent patient ID 
AD <- metadata[metadata$diagnosis=='AD', 'sample_id']
Dementia <- metadata[metadata$diagnosis=='Dementia', 'sample_id']
Control <- metadata[metadata$diagnosis=='Ctrl', 'sample_id']

# tc_df and dwm_df contains cryptic exon counts separated according to brain region
tc_df <- read.csv('tc_df.csv', row.names= 1)
dwm_df <- read.csv('dwm_df.csv', row.names = 1)
# contains the regular gene counts corresponding to cryptic exon counts
# separated according to brain region 
regGenes_tc <- read.csv('regGenes_tc.csv', row.names = 1)
regGenes_dwm <- read.csv('regGenes_dwm.csv', row.names = 1)
# tardbp and grn expression separated according to brain region
tardbp_tc <- read.csv('tardbp_tc.csv', row.names=1)
tardbp_dwm <- read.csv('tardbp_dwm.csv', row.names=1)
grn_tc <- read.csv('grn_tc.csv', row.names=1)
grn_dwm <- read.csv('grn_dwm.csv', row.names= 1)

# putting all the read counts into a list format for easy plotting for temporal cortex data
reads_tc_list <- list()
geneNames_tc <- c()
for (i in 1:nrow(tc_df)) {
  cr <- tc_df[i,]
  reg <- regGenes_tc[i, ]
  # loop through all the columns, plot case vs control. braak stage. 
  reads_tc <- matrix(ncol=7)
  colnames(reads_tc) <- c('pid', 'Group', 'Proportion', 'Diagnosis', 'cellType', 'countE', 'countG')
  crEx <- rownames(tc_df)[i]
  geneID <- substr(crEx, 1, 15)
  geneName <- crExon_geneID[crExon_geneID$ensembleID== geneID, ]$Genes[1]
  for (j in 1:ncol(tc_df)) {
    countG <- as.numeric(reg[j])
    countE <- as.numeric(cr[j])
    colname <- colnames(tc_df)[j]
    pid <- substr(colname, 1, 5)
    # categorizing according to braak stage
    if (pid %in% group1) {
      groupt <- '0-1'
    }
    else if (pid %in% group2) {
      groupt <- '2-4'
    }
    else {
      groupt <- '5-6'
    }
    # categorizing according to diagnosis 
    if (pid %in% AD) {
      diagnosis <- 'AD'
    }
    else if (pid %in% Dementia) {
      diagnosis <- 'Dementia'
    }
    else {
      diagnosis <- 'Control'
    }
    idx <- grep(substr(colname, 7, 10), cellTypes)
    celltype <- cellTypes[idx]
    proportion <- countE / (countG + countE)
    insert <- c(pid, groupt, proportion, diagnosis, celltype, countE, countG)
    reads_tc <- rbind(reads_tc, insert)
  }
  
  reads_tc <- data.frame(reads_tc[-c(1),])
  reads_tc$Proportion <- as.numeric(reads_tc$Proportion)
  reads_tc$countE <- as.numeric(reads_tc$countE)
  reads_tc$countG <- as.numeric(reads_tc$countG)
  
  reads_tc <- cbind(reads_tc, t(tardbp_tc))
  reads_tc <- cbind(reads_tc, t(grn_tc))
  reads_tc <- cbind(reads_tc, rep(geneName, nrow(reads_tc)))
  
  # only extracting reads that have a finite proportion
  idx_fin <- is.finite(reads_tc[,'Proportion'])
  reads_tc <- data.frame(reads_tc[idx_fin,])
  
  colnames(reads_tc) <- c('pid', 'Group', 'Proportion', 'Diagnosis', 'cellType', 'countE', 'countG',
                          'TARDBP', 'GRN', 'gene')
  
  myComparisons <- list(c('0-1', '2-4'), c('2-4', '5-6'), c('0-1', '5-6'))
  
  reads_tc_list[[geneName]] <- reads_tc
  add <- c(geneName, nrow(reads_tc))
  geneNames_tc <- rbind(geneNames_tc, add)
  
}
colnames(geneNames_tc) <- c('Gene Name', 'Occurrences')
rownames(geneNames_tc) <- NULL

# putting all the read counts into a list format for easy plotting for deep white matter data
reads_dwm_list <- list()
geneNames_dwm <- c()
for (i in 1:nrow(dwm_df)) {
  cr <- dwm_df[i,]
  reg <- regGenes_dwm[i, ]
  # loop through all the columns, plot case vs control. braak stage. 
  reads_dwm <- matrix(ncol=7)
  colnames(reads_dwm) <- c('pid', 'Group', 'Proportion', 'Diagnosis', 'cellType', 'countE', 'countG')
  
  crEx <- rownames(dwm_df)[i]
  geneID <- substr(crEx, 1, 15)
  geneName <- crExon_geneID[crExon_geneID$ensembleID== geneID, ]$Genes[1]
  for (j in 1:ncol(dwm_df)) {
    countG <- as.numeric(reg[j])
    countE <- as.numeric(cr[j])
    colname <- colnames(dwm_df)[j]
    pid <- substr(colname, 1, 5)
    if (pid %in% group1) {
      groupt <- '0-1'
    }
    else if (pid %in% group2) {
      groupt <- '2-4'
    }
    else {
      groupt <- '5-6'
    }
    
    if (pid %in% AD) {
      diagnosis <- 'AD'
    }
    else if (pid %in% Dementia) {
      diagnosis <- 'Dementia'
    }
    else {
      diagnosis <- 'Control'
    }
    idx <- grep(substr(colname, 7, 10), cellTypes)
    celltype <- cellTypes[idx]
    proportion <- countE / (countG + countE)
    insert <- c(pid, groupt, proportion, diagnosis, celltype, countE, countG)
    reads_dwm <- rbind(reads_dwm, insert)
  }
  
  reads_dwm <- data.frame(reads_dwm[-c(1),])
  reads_dwm$Proportion <- as.numeric(reads_dwm$Proportion)
  reads_dwm$countE <- as.numeric(reads_dwm$countE)
  reads_dwm$countG <- as.numeric(reads_dwm$countG)
  
  reads_dwm <- cbind(reads_dwm, t(tardbp_dwm))
  reads_dwm <- cbind(reads_dwm, t(grn_dwm))
  reads_dwm <- cbind(reads_dwm, rep(geneName, nrow(reads_dwm)))
  
  reads_dwm <- reads_dwm[is.finite(reads_dwm[,'Proportion']),]
  
  idx_fin <- is.finite(reads_dwm[,'Proportion'])
  reads_dwm <- data.frame(reads_dwm[idx_fin,])
  
  colnames(reads_dwm) <- c('pid', 'Group', 'Proportion', 'Diagnosis', 'cellType', 'countE', 'countG',
                           'TARDBP', 'GRN', 'gene')
  
  myComparisons <- list(c('0-1', '2-4'), c('2-4', '5-6'), c('0-1', '5-6'))
  
  reads_dwm_list[[geneName]] <- reads_dwm
  add <- c(geneName, nrow(reads_dwm))
  geneNames_dwm <- rbind(geneNames_dwm, add)
  
}
colnames(geneNames_dwm) <- c('Gene Name', 'Occurrences')
rownames(geneNames_dwm) <- NULL




# each list contains the plots displayed in the rshiny 
# dwm1: cryptic exon proportion comparison between control and case
# dwm2: cryptic exon proportion comparison between braak stages
# dwm3: cryptic exon proportion comparison between TARDBP expression 
# dwm4: cryptic exon proportion comparison between GRN expression 
# dwm5: regular gene expression comparison between control and case
# dwm6: regular gene expression comparison between braak stages
# pval_dwm: used for storing p-value obtained from cryptic exon proportion comparison between control and case (plot1)
# geneNames_dwm2: the table that will be used as the reactive table that shows each cryptic exon's p-value as well as the relevant plots
dwm1 <- list()
dwm2 <- list()
dwm3 <- list()
dwm4 <- list()
dwm5 <- list()
dwm6 <- list()
pval_dwm <- vector()
geneNames_dwm2 <- c()

# loading p.adj values from edgeR differential expression analysis 
diffexp_dwm_diag <- read_csv("DE_dwm_diagnosis.csv")
diffexp_dwm_braak <- read_csv("DE_dwm_braak.csv")
diffexp_tc_diag <- read_csv("DE_tc_diagnosis.csv")
diffexp_tc_braak <- read_csv("DE_tc_braak.csv")


#have another loop for generating the plots
for (i in 1:length(reads_dwm_list)) {
  data <- reads_dwm_list[[i]]
  g <- geneNames_dwm[i, 1]
  
  data2 <- data %>% group_by(cellType) %>% top_n(1, Proportion) %>%
    slice_max(Proportion, with_ties= FALSE) %>%
    select(cellType, "y.position"= Proportion) %>% as.data.frame() 
  data3 <- data %>% group_by(cellType) %>% top_n(1, countG) %>%
    slice_max(countG, with_ties= FALSE) %>%
    select(cellType, "y.position"= countG) %>% as.data.frame() 
  
  # load the corrected p values
  # corr1 and 2 is between the proportions 
  corr1 <- read.csv(paste0("dwm_fig1", g, ".csv")) %>%
    as.data.frame() %>% merge(data2, by= "cellType", all= TRUE) %>% 
    mutate(p.adj = signif(p.adj, digits= 2))
  corr2 <- read.csv(paste0("dwm_fig2", g, ".csv")) %>%
    as.data.frame() %>% merge(data2, by= "cellType", all= TRUE) %>% 
    mutate(p.adj= signif(p.adj, digits=2))
  
  # corr3 and 4 is between the raw gene counts 
  corr3 <- diffexp_dwm_diag %>% filter(symbol == g) %>% 
    mutate(p.adj = signif(FDR, digits= 2)) %>% 
    merge(data3, by= "cellType", all= TRUE) 
  
  corr4 <- diffexp_dwm_braak %>% filter(symbol == g) %>% 
    mutate(p.adj = signif(FDR, digits= 2)) %>% 
    merge(data3, by= "cellType", all= TRUE)
  
  pval1 <- corr1[is.finite(corr1$p.adj),'p.adj'] %>% min() %>% as.character()
  
  # cells can only be specific for the first two plots
  celldata <- data[data$Proportion != 0, ]$cellType
  cells <- sort(unique(celldata))
  
  # plots proportion between diagnoses 
  p_dwm <- ggplot(data, aes(x=factor(Diagnosis, levels= c("Control", "Dementia", "AD")), y=Proportion, fill= cellType)) +
    geom_quasirandom(size=0.1) +
    geom_boxplot(alpha=0.1, col= 'black', width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab('Proportion of Counts') +
    xlab('Patient Diagnosis') + stat_pvalue_manual(corr1, label= 'p.adj') +
    ggtitle(paste0("Gene ", g, ' White Matter Proportion vs Diagnosis')) +
    theme(plot.title=element_text(size=10)) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # plots proportion between braak stages
  p_dwm_2 <- ggplot(data, aes(x=factor(Group, levels= c("0-1", "2-4", "5-6")), y= Proportion, fill= cellType)) +
    geom_quasirandom(size=0.1) + geom_boxplot(alpha=0.1, col="black", width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab("Proportion of Counts") + xlab("Patient Braak Stage") +
    stat_pvalue_manual(corr2, label= 'p.adj', step.group.by= 'cellType', step.increase= 0.15) +
    ggtitle(paste0('Gene ', g, ' White Matter Proportion vs Braak Stage')) +
    theme(plot.title= element_text(size=10)) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # plots gene counts between diagnoses
  p_dwm_diagnosis <- ggplot(data, aes(x=factor(Diagnosis, levels= c("Control", "Dementia", "AD")), y=countG, fill= cellType)) +
    geom_quasirandom(size=0.1) +
    geom_boxplot(alpha=0.1, col= 'black', width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab('Counts Per Million') +
    xlab('Patient Diagnosis') + 
    stat_pvalue_manual(corr3, label= 'p.adj', step.group.by = "cellType", step.increase= 0.15) +
    ggtitle(paste0("Gene ", g, ' White Matter Proportion vs Diagnosis')) +
    theme(plot.title=element_text(size=10)) 
  
  # plots gene counts between braak stages 
  p_dwm_stage <- ggplot(data, aes(x=factor(Group, levels= c("0-1", "2-4", "5-6")), y= countG, fill= cellType)) +
    geom_quasirandom(size=0.1) + geom_boxplot(alpha=0.1, col="black", width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab("Counts Per Million") + xlab("Patient Braak Stage") +
    stat_pvalue_manual(corr4, label= 'p.adj', step.group.by= 'cellType', step.increase= 0.15) +
    ggtitle(paste0('Gene ', g, ' White Matter Proportion vs Braak Stage')) +
    theme(plot.title= element_text(size=10)) 
  
  p_dwm_tardbp <- ggplot(data, aes(x= TARDBP, y= Proportion)) + geom_point(aes(colour=Diagnosis)) +
    facet_wrap(~cellType, nrow=2) + geom_smooth(color = "blue", method= 'lm') +
    theme_classic() + ylab('Proportion of Counts') + xlab('TMM Normalized TARDBP count (log10)') + scale_x_log10() +
    ggtitle(paste0('Gene ', g, ' White Matter Proportion vs GRN')) +
    theme(plot.title =element_text(size=10), legend.position= 'bottom') + stat_cor(cor.coef.name = "R")
  
  p_dwm_grn <- ggplot(data, aes(x= GRN, y=Proportion)) + geom_point(aes(colour=Diagnosis)) +
    facet_wrap(~cellType, nrow=2) + geom_smooth(color= "blue", method= 'lm') +
    theme_classic() + ylab('Proportion of Counts') + xlab('TMM Normalized GRN count (log10)') + scale_x_log10() +
    ggtitle(paste0('Gene ', g, ' White Matter Proportion vs GRN')) +
    theme(plot.title =element_text(size=10), legend.position= 'bottom') + stat_cor(cor.coef.name = "R")
  
  # requires cell data to be updated bc some read counts exist for raw counts but not proportion
  celldata <- data[data$countG != 0, ]$cellType
  cells <- sort(unique(celldata))
  max_countG <- max(data$countG[is.finite(data$countG)])
  
  
  
  dwm1[[g]] <- p_dwm
  dwm2[[g]] <- p_dwm_2
  dwm3[[g]] <- p_dwm_tardbp
  dwm4[[g]] <- p_dwm_grn
  dwm5[[g]] <- p_dwm_diagnosis 
  dwm6[[g]] <- p_dwm_stage
  
  temp <- c(g, pval1)
  geneNames_dwm2 <- rbind(geneNames_dwm2, temp)
}
sortidx <- sort(as.numeric(geneNames_dwm2[,2]), index.return= TRUE)$ix

# sortidx <- sort(pval_dwm, index.return=TRUE)$ix
dwm_final <- cbind(dwm1, dwm2, dwm3, dwm4, dwm5, dwm6)
dwm_final <- dwm_final[sortidx, ]

geneNames_dwm2 <- geneNames_dwm2[sortidx, ]

colnames(geneNames_dwm2) <- c('Gene Name', 'P Value')
rownames(geneNames_dwm2) <- NULL

# temporal cortex 
# each list contains the plots displayed in the rshiny 
# tc1: cryptic exon proportion comparison between control and case
# tc2: cryptic exon proportion comparison between braak stages
# tc3: cryptic exon proportion comparison between TARDBP expression 
# tc4: cryptic exon proportion comparison between GRN expression 
# tc5: regular gene expression comparison between control and case
# tc6: regular gene expression comparison between braak stages
# pval_tc: used for storing p-value obtained from cryptic exon proportion comparison between control and case (plot1)
# geneNames_tc2: the table that will be used as the reactive table that shows each cryptic exon's p-value as well as the relevant plots
tc1 <- list()
tc2 <- list()
tc3 <- list()
tc4 <- list()
tc5 <- list()
tc6 <- list()
pval_tc <- vector()
geneNames_tc2 <- c()

for (j in 1:length(reads_tc_list)) { 
  data <- reads_tc_list[[j]]
  g <- geneNames_tc[j, 1]
  
  # Load the y.position from the original dataset, merge to the 
  # dataset containing the corrected p values and set as y.position 
  # data2 is for proportion
  # data3 is for gene counts
  data2 <- data %>% group_by(cellType) %>% top_n(1, Proportion) %>%
    slice_max(Proportion, with_ties= FALSE) %>%
    select(cellType, "y.position"= Proportion) %>% as.data.frame()
  data3 <- data %>% group_by(cellType) %>% top_n(1, countG) %>%
    slice_max(countG, with_ties= FALSE) %>%
    select(cellType, "y.position"= countG) %>% as.data.frame()
  
  # load the corrected p values
  corr1 <- read.csv(paste0("tc_fig1", g, ".csv")) %>%
    as.data.frame() %>% merge(data2, by= "cellType", all= TRUE) %>% 
    mutate(p.adj = signif(p.adj, digits=2))
  corr2 <- read.csv(paste0("tc_fig2", g, ".csv")) %>% as.data.frame() %>%
    merge(data2, by= "cellType", all= TRUE) %>% 
    mutate(p.adj = signif(p.adj, digits= 2))
  # corr3 and 4 is between the raw gene counts 
  corr3 <- diffexp_tc_diag %>% filter(symbol == g) %>% 
    mutate(p.adj = signif(FDR, digits= 2)) %>% 
    merge(data3, by= "cellType", all= TRUE) 
  
  corr4 <- diffexp_tc_braak %>% filter(symbol == g) %>% 
    mutate(p.adj = signif(FDR, digits= 2)) %>% 
    merge(data3, by= "cellType", all= TRUE)


  pval2 <- corr1[is.finite(corr1$p.adj),'p.adj'] %>% min()
  
  celldata <- data[data$Proportion != 0, ]$cellType
  cells <- sort(unique(celldata))
  
  p_tc <- ggplot(data, aes(x=factor(Diagnosis, levels= c("Control", "Case")), y=Proportion, fill= cellType)) +
    geom_quasirandom(size=0.1) +
    geom_boxplot(alpha=0.1, col= 'black', width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab('Proportion of Counts') +
    xlab('Patient Diagnosis') + stat_pvalue_manual(corr1, label= 'p.adj') +
    ggtitle(paste0('Gene ', g, ' Grey Matter Proportion vs Diagnosis')) +
    theme(plot.title =element_text(size=10)) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 
  
  p_tc_2 <- ggplot(data, aes(x=factor(Group, levels=c("0-1", "2-4", "5-6")), y=Proportion,fill= cellType)) +
    geom_quasirandom(size=0.1) + geom_boxplot(alpha=0.1, col="black", width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab('Proportion of Counts') + xlab("Patient Braak Stage") +
    stat_pvalue_manual(corr2, label= 'p.adj', step.group.by= 'cellType', step.increase= 0.15) +
    ggtitle(paste0('Gene ', g, ' Grey Matter Proportion vs Braak Stage')) +
    theme(plot.title =element_text(size=10)) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  p_tc_diagnosis <- ggplot(data, aes(x=factor(Diagnosis, levels= c("Control", "Dementia", "AD")), y=countG, fill= cellType)) +
    geom_quasirandom(size=0.1) +
    geom_boxplot(alpha=0.1, col= 'black', width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab('Counts Per Million') + xlab('Patient Diagnosis') + 
    stat_pvalue_manual(corr3, label= 'p.adj', step.group.by= 'cellType', step.increase= 0.15) +
    ggtitle(paste0('Gene ', g, ' Grey Matter Proportion vs Diagnosis')) +
    theme(plot.title =element_text(size=10)) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  p_tc_stage <- ggplot(data, aes(x=factor(Group, levels=c("0-1", "2-4", "5-6")), y=countG, fill= cellType)) +
    geom_quasirandom(size=0.1) + geom_boxplot(alpha=0.1, col="black", width=0.25) +
    facet_wrap(~cellType, nrow=2) + theme_classic() +
    ylab('Counts Per Million') + xlab("Patient Braak Stage") +
    stat_pvalue_manual(corr4, label= 'p.adj', step.group.by= 'cellType', step.increase= 0.15) +
    ggtitle(paste0('Gene ', g, ' Grey Matter Proportion vs Braak Stage')) +
    theme(plot.title =element_text(size=10)) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  p_tc_tardbp <- ggplot(data, aes(x= TARDBP, y= Proportion)) + geom_point(aes(colour=Diagnosis)) + 
    facet_wrap(~cellType, nrow=2) + geom_smooth(color = "blue", method= 'lm') +
    theme_classic() + ylab('Proportion of Counts') + xlab('TMM Normalized TARDBP count (log10)') + 
    scale_x_log10() +
    ggtitle(paste0('Gene ', g, ' Grey Matter Proportion vs TARDBP')) +
    theme(plot.title =element_text(size=10), legend.position= 'bottom') + stat_cor(cor.coef.name = "R")
  
  p_tc_grn <- ggplot(data, aes(x= GRN, y=Proportion)) + geom_point(aes(colour=Diagnosis)) +
    facet_wrap(~cellType, nrow=2) + geom_smooth(color= "blue", method= 'lm') +
    theme_classic() + ylab('Proportion of Counts') + xlab('TMM Normalized GRN count (log10)') + 
    scale_x_log10() +
    ggtitle(paste0('Gene ', g, ' Grey Matter Proportion vs GRN')) +
    theme(plot.title =element_text(size=10), legend.position= 'bottom') + stat_cor(cor.coef.name = "R")
  
  # requires cell data to be updated bc some read counts exist for raw counts but not proportion
  celldata <- data[data$countG != 0, ]$cellType
  cells <- sort(unique(celldata))
  max_countG <- max(data$countG[is.finite(data$countG)])
  

  
  tc1[[g]] <- p_tc
  tc2[[g]] <- p_tc_2
  tc3[[g]] <- p_tc_tardbp
  tc4[[g]] <- p_tc_grn
  tc5[[g]] <- p_tc_diagnosis
  tc6[[g]] <- p_tc_stage
  
  temp <- c(g, pval2)
  geneNames_tc2 <- rbind(geneNames_tc2, temp)
  
}
sortidx2 <- sort(as.numeric(geneNames_tc2[,2]), index.return= TRUE)$ix

# sortidx2 <- sort(pval_tc, index.return=TRUE)$ix
tc_final <- cbind(tc1, tc2, tc3, tc4, tc5, tc6)
# rownames(tc_final) <- geneNames_tc[,1]
tc_final <- tc_final[sortidx2, ]

geneNames_tc2 <- geneNames_tc2[sortidx2, ]

colnames(geneNames_tc2) <- c('Gene Name', 'P Value')
rownames(geneNames_tc2) <- NULL

# sorting the reads to match sorted index
reads_dwm_list <- reads_dwm_list[sortidx]
reads_tc_list <- reads_tc_list[sortidx2]

# r shiny code
ui <- fluidPage(
  # App title ----
  titlePanel("Cryptic Exon Expression in Single Cells"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      radioButtons("which_matter", "Brain Region",  
                   choices= c("Grey Matter" = "part1", 
                              "White Matter" = "part2")), 
      DT::dataTableOutput("mytable"), width= 3
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      uiOutput("myUIOutput"), 
      width= 9
    )
    
  )
  
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  # a reactive expression uses a widget input and returns a value 
  crExon_list_tc <- reactive({
    return(geneNames_tc2)
  })
  
  crExon_list_dwm <- reactive({
    return(geneNames_dwm2)
  })
  
  output$mytable <- renderDataTable({
    region <- input$which_matter
    if (region=='part2') {
      crExon_list_dwm()
    }
    else if (region== 'part1') {
      crExon_list_tc()
    }
    
  }, selection= list(mode='single', selected=1L))
  
  output$p1 <- renderPlot({
    # use plots in the tc_plots list?
    region <- input$which_matter
    idx <- input$mytable_rows_selected
    # geneName <- geneNames_dwm[idx, 1]
    # data <- reads_dwm_list[[idx]]
    if (region=='part2') {
      # use the variable gene to select for gene of interest
      # get the index of the selected gene, obtain the reads from reads_tc_list,
      # plot with ggplot
      p_dwm <- dwm_final[idx, 1][[1]]
      p_dwm_2 <- dwm_final[idx, 2][[1]]
      p_dwm_tardbp <- dwm_final[idx, 3][[1]]
      p_dwm_grn <- dwm_final[idx, 4][[1]]
      p_dwm_diagnosis <- dwm_final[idx, 5][[1]]
      p_dwm_stage <- dwm_final[idx, 6][[1]]
      
      grid.arrange(p_dwm, p_dwm_2, p_dwm_diagnosis, p_dwm_stage, p_dwm_tardbp, p_dwm_grn, ncol=2)
    }
    else if (region=='part1') {
      p_tc <- tc_final[idx, 1][[1]]
      p_tc_2 <- tc_final[idx, 2][[1]]
      p_tc_tardbp <- tc_final[idx, 3][[1]]
      p_tc_grn <- tc_final[idx, 4][[1]]
      p_tc_diagnosis <- tc_final[idx, 5][[1]]
      p_tc_stage <- tc_final[idx, 6][[1]]
      
      grid.arrange(p_tc, p_tc_2, p_tc_diagnosis, p_tc_stage, p_tc_tardbp, p_tc_grn, ncol=2)
    }
  }, height= 1500, width= 1100)
  
  output$myUIOutput <- renderUI({
    region <- input$which_matter
    if (is.null(input$mytable_rows_selected)==TRUE) {
      p("No genes selected")
    }
    else {
      idx <- input$mytable_rows_selected
      if (region=='part1') {
        data <- reads_dwm_list[[idx]]
      }
      else {
        data <- reads_tc_list[[idx]]
      }
      
      if (nrow(data) == 0) {
        p("No expression")
      }
      else {
        plotOutput("p1")
      }
      
    }
  })
  
  
}

shinyApp(ui= ui, server=server)

