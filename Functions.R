
# mytheme -----------------------------------------------------------------

mytheme <- theme_classic() +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),  # 标题居中
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    axis.title = element_text(size = 17),# 添加黑色方框边界
    legend.position = "right"
  )

# tcga_download///Obtaining the data of RNA-sequencing expression profile.----
tcga_download <- function(name){
  x <- str_c("TCGA-",name)
  query_expr <- TCGAbiolinks::GDCquery(project = x,
                                       data.category = "Transcriptome Profiling",
                                       data.type = "Gene Expression Quantification",
                                       workflow.type = "STAR - Counts"
  )
  TCGAbiolinks::GDCdownload(query_expr)
  expr <- TCGAbiolinks::GDCprepare(query_expr)
  se <- expr
  rowdata <- SummarizedExperiment::rowData(se)
  se_mrna <- se[rowdata$gene_type == "protein_coding",]
  se_lnc <- se[rowdata$gene_type == "lncRNA"]
  
  mrna_counts <- SummarizedExperiment::assay(se_mrna,"unstranded")
  
  lnc_counts <- SummarizedExperiment::assay(se_lnc,"unstranded")
  mrna <- SummarizedExperiment::rowData(se_mrna)$gene_name
  lnc <- SummarizedExperiment::rowData(se_lnc)$gene_name
  symbol_mrna <- data.frame(gene=mrna) %>%
    cbind(.,mrna_counts) %>%
    rownames_to_column(var = "ensembl") %>%
    select(1:2)
  symbol_lnc <- data.frame(gene=lnc) %>%
    cbind(.,lnc_counts) %>%
    rownames_to_column(var = "ensembl") %>%
    select(1:2)
  message("矩阵下载成功！(nkyqgong)")
  expr_data <- list(mrna=mrna_counts,lnc=lnc_counts,symbol_lnc=symbol_lnc,symbol_mrna=symbol_mrna)
  return(expr_data)
}
# tcga_clinical///Obtaining the clinical data matched with the data of expression profile.----
tcga_clinical <- function(x){
  y <- str_c("TCGA-",x)
  clin.coad <- TCGAbiolinks::GDCquery_clinic(y, "clinical")
  clinical <- clin.coad
  message("临床数据下载成功！(nkyqgong)")
  return(clinical)
}
# tcga_snp///Obtaining the data of somatic mutation.----
tcga_snp <- function(x){
  y <- str_c("TCGA-",x)
  query_snp <- TCGAbiolinks::GDCquery(
    project = y,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  TCGAbiolinks::GDCdownload(query_snp)
  snp <- TCGAbiolinks::GDCprepare(query_snp)
  message("体细胞突变数据下载成功！(nkyqgong)")
  return(snp)
}
# tcga_cbind///The fuction was created for converting names in ensembl format to symbol format.----
tcga_cbind <- function(counts){
  mrna <- counts$mrna
  lnc <- counts$lnc
  mrna <- cbind(mrna,counts$symbol_mrna) %>% 
    select(-ensembl) %>% 
    select(gene,everything())
  lnc <- cbind(lnc,counts$symbol_lnc) %>% 
    select(-ensembl) %>% 
    select(gene,everything())
  result <- list(mrna=mrna,lnc=lnc)
  return(result)
}

# tcga_clinical_normalization ---------------------------------------------
# The function was written for the normalization of clinical data
tcga_clinical_normalization <- function(clinical){
  clinical <- clinical %>% 
    dplyr::select(project,submitter_id,days_to_birth,days_to_death,
                  days_to_last_follow_up,days_to_diagnosis,
                  year_of_birth,year_of_death,
                  gender,vital_status,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_n,ajcc_pathologic_m,
                  age_at_diagnosis,
                  age_at_index)
  clinical$days_to_death[is.na(clinical$days_to_death)] <- 0   #缺失值标记为0
  clinical$days_to_last_follow_up[is.na(clinical$days_to_last_follow_up)] <- 0
  clinical <- clinical %>% 
    mutate(os.time=days_to_death+days_to_last_follow_up,.before = days_to_birth) %>% 
    mutate(os=ifelse(vital_status=='Alive',0,1),.before = os.time) %>% 
    drop_na(ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_n,ajcc_pathologic_m)
  return(clinical)
}

# tcga_normalization ------------------------------------------------------
# The function can help you quickly for obtaining the expr standarded and grouping the samples of patients.
tcga_normalization <- function(exprname){
  expr <- exprname
  expr <- expr %>% 
    mutate(number=rowSums(.==0))%>%
    select(number,everything())%>%
    filter(number<=((ncol(.)-2)/2))%>%
    select(-number)%>%
    aggregate(.~gene,mean)
  rownames(expr) <- NULL
  expr <- expr %>% 
    column_to_rownames(var = "gene")
  tumor_normal <- expr %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "sample") %>% 
    filter(str_sub(sample,14,15)=="01"|str_sub(sample,14,15)=="11") %>%
    mutate(group=ifelse(str_sub(sample,14,15)=="01","Tumor","Normal")) %>% 
    select(group,everything()) %>% 
    column_to_rownames(var = "sample") %>% 
    t() %>% 
    as.data.frame()
  tumor <- tumor_normal %>% 
    t() %>% 
    as.data.frame() %>% 
    filter(group=="Tumor") %>% 
    select(-group) %>% 
    t() %>% 
    as.data.frame()
  normal <- tumor_normal %>% 
    t() %>% 
    as.data.frame() %>% 
    filter(group=="Normal") %>% 
    select(-group) %>% 
    t() %>% 
    as.data.frame()
  normalization <- list(tumor_normal=tumor_normal,tumor=tumor,normal=normal)
  return(normalization)
}

# gtex_normalization ------------------------------------------------------

gtex_normalization <- function(type="Kidney"){
  gtex_clinical <- read_tsv("GTEX_phenotype") 
  gtex_id <- read.csv("probeMap_gencode.v23.annotation.gene.csv")
  gtex <- fread("gtex_gene_expected_count")
  gtex_clinical <- gtex_clinical %>% 
    filter(`_primary_site`==type)
  sample <- gtex_clinical %>% 
    select("Sample") %>% 
    as.data.frame() %>% 
    dplyr::rename(sample=Sample)
  gtex <- gtex %>% 
    rename(id=sample) %>% 
    column_to_rownames(var = "id") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "sample") %>% 
    inner_join(.,sample,by="sample") %>% 
    column_to_rownames(var = "sample") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "id") %>% 
    inner_join(.,gtex_id,by="id") %>% 
    select(-id) %>% 
    select(gene,everything()) %>%
    aggregate(.~gene,mean) %>% 
    column_to_rownames(var = "gene") %>% 
    rownames_to_column(var = "id") %>% 
    mutate(id = if_else(str_detect(id, "\\."), str_extract(id, "^[^.]+"), id)) %>% 
    aggregate(.~id,mean) %>% 
    column_to_rownames(var = "id")
  gtex <- 2^gtex-1
  gtex <- round(gtex,3)
  result <- list(gtex_clinical=gtex_clinical,gtex=gtex)
}


#  tcga_gtex_cbind---------------------------------------
# Merging the data of TCGA and GTEX
tcga_gtex_cbind <- function(tcga,gtex){
  tcga <- tcga$tumor_normal %>% 
    rownames_to_column(var = "id")
  gtex <- gtex$gtex  
  gtex <- gtex %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(group="Normal") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "id")
  tcga_gtex <- gtex %>% 
    inner_join(.,tcga,by="id") %>% 
    column_to_rownames(var = "id")
  return(tcga_gtex)
}


# expr_numeric ------------------------------------------------------------

expr_numeric <- function(expr){
  rownames <- rownames(expr)
  expr <- apply(expr,2,as.numeric)
  row.names(expr) <- rownames
  expr <- as.data.frame(expr)
  return(expr)
}


