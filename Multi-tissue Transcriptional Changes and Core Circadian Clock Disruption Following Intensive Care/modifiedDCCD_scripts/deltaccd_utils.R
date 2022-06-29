library(data.table)
calcDist = function(r1, r2) sqrt(sum((r1 - r2)^2, na.rm = TRUE))


calcCCDSimple = function(ref, emat, method = 'spearman', scale = FALSE) {
  corVecRef = ref[upper.tri(ref)]
  corMatTest = stats::cor(t(emat), method = method)
  corVecTest = corMatTest[upper.tri(corMatTest)]
  ccd = calcDist(corVecRef, corVecTest)
  
  if (isTRUE(scale)) {
    nPairs = choose(ncol(ref), 2)
    ccd = ccd / nPairs}
  
  return(ccd)}


checkVar = function(emat, groupVec) {
  groupNow = group = variance = NULL
  
  varCheck = foreach (groupNow = sort(unique(groupVec)), .combine = rbind) %do% {
    
    varVec = apply(emat[, groupVec == groupNow], MARGIN = 1, 
                   FUN = stats::var, na.rm = TRUE)
    varDt = data.table::as.data.table(varVec, keep.rownames = 'gene')
    data.table::setnames(varDt, 'varVec', 'variance')
    varDt[, group := groupNow]
    
    zeroVar = varDt[variance == 0]}
  varCheck[, variance := NULL]  
  
  if (nrow(varCheck) > 0) {
    stop('Zero variance in the following gene-group pairs:\n', 
         paste(utils::capture.output(print(varCheck)), collapse = '\n'))}
  
  invisible()}


checkGenes = function(emat, refCor, geneNames = NULL) {
  if (is.null(geneNames)) {
    geneNames = rownames(refCor)[rownames(refCor) %in% rownames(emat)]
  } else { 
    geneNames = rownames(refCor)[rownames(refCor) %in% geneNames]}
  
  if (length(geneNames) < nrow(refCor)) {
    missingGenes = setdiff(rownames(refCor), geneNames)
    stop(paste0('The following gene(s) is/are not in the expression matrix:\n',
                paste0(missingGenes, collapse = '\n')))} 
  
  return(geneNames)}


checkRefCor = function(refCor, refEmat = NULL, geneNames = NULL, method = 'spearman') {
  if (missing(refCor)) {
    if (is.null(refEmat)) {
      stop('Either refCor or refEmat must be supplied.')}
    refCor = stats::cor(t(refEmat[geneNames,]), method = method)
  } else if (any(rownames(refCor) != colnames(refCor)) || !isSymmetric(refCor)) {
    stop('refCor must be a correlation matrix, with identical rownames and colnames.')}
  
  return(refCor)}

calcDeltaCCDSimple = function(ref, emat, idx, method = 'spearman', scale = FALSE) {
  d0 = calcCCDSimple(ref, emat[, !idx], method = method, scale = scale)  #calc CCD normal
  d1 = calcCCDSimple(ref, emat[, idx], method = method, scale = scale)   #calc CCD non-normal
  d = d1 - d0
  return(d)}

get_sub_str = function(x, start = nchar(x), end = nchar(x)){
  return(substr(x, start, end))
}
char_to_binary = function(vect, control = "A", treatment = "V"){
  vect[grepl(control,vect )] = 0
  vect[grepl(treatment, vect)] = 1
  vect = sapply(vect, as.numeric)
  return(vect)
}

binarize_perms = function(df){
  df = apply(df, 2, get_sub_str )
  df = apply(df, 2, char_to_binary)
  return(df)
}


#####debugging function to determine proportion of centers in perms#######
get_perm_centers = function(df){
  df = apply(df, 2, get_sub_str, 1, 2)
  #df = apply(df, 2, char_to_binary)
  return(df)
}

######

binarize_labels_cust = function(mix_vect, start, end){
  ans =unname(sapply(mix_vect, get_sub_str, start, end))
 return(ans)
}
binarize_labels = function(mix_vect){
  ans =unname(sapply(mix_vect, get_sub_str))
  return(ans)
}

makePerms = function(idx, nPerm = 1000, dopar = FALSE, tries = 0, factor = 10) {
  doOp = if (isTRUE(dopar)) `%dorng%` else `%do%`
  if(tries >= 3){
    print(paste("ERROR: idx only has", length(idx), "elements,",nPerm, "permutations not available"))
    return(NA)}
  
  r = doOp(foreach(i = 1:10000, .combine = rbind), {
    sample(idx, length(idx))
    })
  
  unique_perms <- unique(r)
  print(dim(unique_perms)[1])
  if(dim(unique_perms)[1]< nPerm){
    unique_perms = makePerms(idx, nPerm = nPerm, dopar = dopar, tries = tries+1, factor = factor*10)
  }
  if(all(is.na(unique_perms))){
    return(NA)}
  else{
  return(unique_perms[1:nPerm,])
    }
  
  }


