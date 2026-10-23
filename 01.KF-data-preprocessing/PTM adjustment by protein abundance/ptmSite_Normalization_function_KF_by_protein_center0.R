# site_df = site_df
# prot_df = prot_df
# site.id = 'id'
# protein.id = 'protein'

ptmSite_Normalization.by.protein <- function(site_df = NULL,
                                             prot_df = NULL,
                                             site.id = 'id',
                                             protein.id = 'protein'
                                             # ,
                                             # sample_cols = NULL
                                             # ,
                                             # output_filename
)
{
  library(dplyr)
  
  print("apply PTM site-protein normalization")
  # prot_df = proteome_mat
  # site_df = ptm_site_mat
  
  ### identify site name /  protein name / annotation columns / sample columns 
  
  site_df$proteinID = site_df[[protein.id]]
  site_df$Index = site_df[[site.id]]
  prot_df$Index = prot_df[[protein.id]]
  
  id.anno.site = names(which(sapply(site_df,typeof)=='character'))
  id.anno.prot = names(which(sapply(prot_df,typeof)=='character'))
  
  sample_cols = intersect(names(which(sapply(site_df,typeof)=='double')),
                          names(which(sapply(prot_df,typeof)=='double')))
  
  ### calculate protein / site means and center protein / sites to mean 0 vectors
  prot_rmeans = rowMeans(prot_df[,sample_cols],na.rm = T)
  site_rmeans = rowMeans(site_df[,sample_cols],na.rm = T)
    
  names(prot_rmeans) = prot_df[[protein.id]]
  names(site_rmeans) = site_df[[site.id]]
    
  prot_center = data.frame(prot_df[,id.anno.prot],
                      prot_df[,sample_cols]-prot_rmeans,check.names = F)
  site_center = data.frame(site_df[,id.anno.site],
                      site_df[,sample_cols]-site_rmeans,check.names = F)
  
  # site_df$proteinID <- gsub("_[^_]*","",site_df$Index)
  
  ### filter protein / site by shared proteins 
  
  proteins <- intersect(prot_center$Index, site_center$proteinID)
  site_flt <- site_center[site_center$proteinID %in% proteins, ]
  prot_flt <- prot_center[prot_center$Index %in% proteins, ]
  print(paste0(nrow(site_flt), " out of ",nrow(site_center), " sites will be used for normalization"))
  
  rm(site_df)
  rm(prot_df)
  
  ### trans form protein / site to longlist and merge them together
  
  site_melt <- reshape2::melt(site_flt,
                              # id.vars = c("Index", "Gene_site", "SequenceWindow", "proteinID"),
                              id.vars = id.anno.site,
                              measure.vars = sample_cols,
                              value.name = "PTM_site",
                              variable.name = "sample")
  
  prot_melt <- reshape2::melt(prot_flt,
                              # id.vars = c("Index", "Gene"),
                              id.vars = id.anno.prot,
                              measure.vars = sample_cols,
                              value.name = "prot",
                              variable.name = "sample")
  
  rm(site_flt)
  rm(prot_flt)
  
  
  site_prot <- merge(site_melt, 
                     prot_melt[,c("Index","sample","prot")], 
                     by.x = c("proteinID","sample"), 
                     by.y = c("Index","sample"), 
                     all.x = T)
  
  site_prot <- na.omit(site_prot)
  
  rm(site_melt)
  rm(prot_melt)
  
  
  ### fit protein specific linear models and get residuals for each regression
  
  # p = lm(PTM_site ~ prot, data = site_prot)
  # res = p$residuals
  
  res = protein_lm_res(site_prot)
  
  site_prot = site_prot %>%
    dplyr::mutate(residual = res)
  
  # rm(p)
  rm(res)
  
  result_cor = cor(site_prot$residual, site_prot$prot)
  cat("correlation with prot after normalization: ",  result_cor, "\n")

  ### reformat residuals to matrix and add back corresponding site mean values
    
  site_norm <- reshape2::dcast(site_prot, 
                               formula = as.formula(paste(paste(setdiff(colnames(site_prot),c('PTM_site','prot','sample','residual')),collapse = '+'),'~sample')),   
                               # formula = Index + Gene_site + SequenceWindow ~ sample, 
                               value.var = "residual")
  
  site_norm_abundance = site_norm
  site_norm_abundance[,sample_cols] = site_norm_abundance[,sample_cols]+site_rmeans[match(site_norm$Index,names(site_rmeans))]
    
  rm(site_prot)
  rm(result_cor)
  
  # write.table(site_norm, file = output_filename, 
  #             quote = F, row.names = F, sep = "\t")
  return(site_norm_abundance)
  
}



# function to calculate protein-ptm linear regression residuals stratified by proteins

protein_lm_res = function(site_prot){
  protein.list = unique(site_prot$protein)
  
  residual = rep(NA,nrow(site_prot))  
  for(pp in protein.list)
  {
    ind.pp = (site_prot$protein == pp)
    residual[ind.pp] = (lm(site_prot$PTM_site[ind.pp]~site_prot$prot[ind.pp]))$residuals
  }
  return(residual)
}



