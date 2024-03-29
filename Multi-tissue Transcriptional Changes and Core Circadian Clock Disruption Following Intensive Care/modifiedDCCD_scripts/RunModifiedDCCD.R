#This script takes genes expression matrices and runs a modified version of the CalcDeltaCCD function from the
#deltaccd library. In the modified function, data from different centers/genders are permuted separately
#in order to account for batch effects. The functions are defined and then called at the bottom of the script.

library(rstudioapi)
library(tidyverse)
library(doParallel)
library(doRNG)
library(deltaccd)
library(openxlsx)
library(edgeR)
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path) #setwd to path of this file
source("./deltaccd_utils.R")
source("./custom_CalcDeltaCCD.R")
registerDoParallel(cores = 4)
in_directory = paste0(path, "/TMMNormedCountsEdgeR/")
corr_plots = paste0(path, "/CorrelationPlots/")
setwd(in_directory)

gene_list = c("PER3" , "CIART" ,"NPAS2", "PER2",  "NR1D2", "CLOCK" ,"ARNTL", "CRY2" , "CRY1",  "PER1" , "NR1D1" ,"DBP","TEF")
gtex_files=list.files(pattern = "*.xlsx")

try(if(!assertthat::not_empty(gtex_files)) stop(paste("No .xlsx files found in", getwd())))

refCor = getRefCor('human', 'pan', FALSE)


######### helper functions #########
get_tissue_name = function(str){
  tissue_name = str_extract(str, pattern = "(.*)(?=_TMM)")
  return(tissue_name)
}

age_to_num = function(age_str){
  str_list = str_split(age_str, pattern = "-", n =2)
  num = mean(sapply(str_list, as.numeric))
  return(num)
}
char_to_binary = function(vect, control = "A", treatment = "V"){
  vect = matrix(vect)
  vect[grepl(control,vect )] = 0
  vect[grepl(treatment, vect)] = 1
  vect = sapply(vect, as.numeric)
  return(vect)
}

######### heavy lifting functions ########

run_perms = function(file, min_in_group = 50, custom_func = T, nperm = 5000, correct = 'both', combat = F, plot = F){
  #function takes the tmm filename and runs the modified DCCD procedure
  data = read.xlsx(file, sheet = "TMM")
  rownames(data) = data[,1]
  coldata = read.xlsx(file, sheet = "subjectInfo")
  data = data[,colnames(data) %in% coldata[,1]] #only keep columns that are V, or A dthhrdy
  tiss_name = get_tissue_name(file)  #get tissue name from filename
  print(tiss_name)
  tiss_name = str_replace_all(tiss_name, fixed(" "), "") #replace any spaces in name
  start = Sys.time()
  
  sex = as.character(coldata$SEX)
  dthhrdy = as.character(coldata$dthhrdy)
  sex[sex == '1'] = "M"
  sex[sex == '2'] = "F"
  dthhrdy[dthhrdy == '1' | dthhrdy == '2'] = 'A'
  dthhrdy[dthhrdy == '0'] = 'V'

  dthhrdy_vec = paste0(coldata$center, dthhrdy, sex)
  dthhrdy_bin = char_to_binary(dthhrdy)
  center = char_to_binary(coldata$center, control = "B1", treatment = "C1")  
  age = coldata$AGE
  if(all(is.na(coldata))){return(NA)}
  if(any(table(coldata$dthhrdy) < min_in_group)){ 
    return(NA)
  }
  
  if(combat){
    print("using Combat corrected TMMS")
    combat_tmms = read.xlsx(file, sheet = "combat_corrected")
    rownames(combat_tmms) = combat_tmms[,1]
    data = combat_tmms[, -1]
  }
  
  emat = data[which(rownames(data) %in% gene_list),]

  
  
  if(custom_func){
    print("custom dccd")
    p_test = customCalcDeltaCCD(refCor, emat, nPerm = nperm , dthhrdy_vec, "A", correct = correct)
  }else{
    p_test = calcDeltaCCD(refCor, emat, nPerm = nperm, binarize_labels_cust(dthhrdy_vec, 3, 3), "A" )
  }
  
  if(plot){
    setwd(corr_plots)
    plotHeatmap(rownames(refCor), emat, groupVec =get_sub_str(dthhrdy_vec, start = 3, end = 3))
    ggsave(paste0(tiss_name,"CorrelationPlot.pdf"))
    setwd(in_directory)
  }
  
  
  p_test$group1_count = sum(!as.logical(char_to_binary(dthhrdy_vec)))
  p_test$group2_count = sum(as.logical(char_to_binary(dthhrdy_vec)))
  return(p_test)
}


get_ccd = function(file, min_in_group = 50, nPerm = 1000 ){
  #This function calculates the basic ccd
  data = read.xlsx(file, sheet = "TMM")
  rownames(data) = data[,1]
  coldata = read.xlsx(file, sheet = "subjectInfo")

  data = data[,colnames(data) %in% coldata[,1]] #only keep columns that are V, or A dthhrdy
  
  tiss_name = get_tissue_name(file)  #get tissue name from filename
  print(tiss_name)
  tiss_name = str_replace_all(tiss_name, fixed(" "), "") #replace any spaces in name
  dthhrdy = as.character(coldata$dthhrdy)
  dthhrdy[dthhrdy == '1' | dthhrdy == '2'] = 'A'
  dthhrdy[dthhrdy == '0'] = 'V'
  
  dthhrdy_bin = char_to_binary(dthhrdy)

  if(any(table(coldata$dthhrdy) < min_in_group)){ 
    return(NA)
  }
  emat = data[which(rownames(data) %in% gene_list),]
  
  acute_ccd= calcCCD(refCor, emat[,!as.logical(dthhrdy_bin)], nPerm = 1000)  #acute ccd
  vent_ccd= calcCCD(refCor, emat[,as.logical(dthhrdy_bin)], nPerm = 1000)  #acute ccd
  
  names= c("A", "V", "A_pval", "V_pval")
  out = c(acute_ccd$CCD,vent_ccd$CCD,  acute_ccd$Pvalue, vent_ccd$Pvalue)
  names(out) = names
  out = as.data.frame(out)
  return(out)
  
}

###### RUN the functions #######

#uncomment for correcting for gender only:
# gender_shuffle = sapply(gtex_files, run_perms, nperm = 5000, correct = 'gender')

#uncomment for correcting for center only:
#center_shuffle = sapply(gtex_files, run_perms, nperm = 5000, correct = 'center')

#correcting for both gender and center
GC_shuffle = sapply(gtex_files, run_perms, nperm = 1000, correct = 'both')

#no_shuffle_result = sapply(gtex_files, run_perms, custom_func = 'F', plot = T)

#uncomment to extract the ccd values for A and V groups
#simple_ccds = sapply(gtex_files, get_ccd, nPerm = 1000)
