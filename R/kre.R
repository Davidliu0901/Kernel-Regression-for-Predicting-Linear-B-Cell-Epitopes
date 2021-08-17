

library(stringr)
library(ggplot2)
library(ggseqlogo)
library(MASS)
library(Biostrings)



# 'Logistic', 'Sigmoid', 'Gaussian'
kernel_sel<-function(kernel_name){
  
  if(kernel_name=='Logistic'){
    return(function(x) 1/(exp(x)+2+exp(-x)))
  }
  
  if(kernel_name=='Sigmoid'){
    return(function(x) 2/(pi*(exp(x)+exp(-x))))
  }
  
  if(kernel_name=='Gaussian'){
    return(dnorm)
  }
}


self_scaled<-function(result_matrx){
  as.matrix((result_matrx - min(result_matrx))/(max(result_matrx) - min(result_matrx)))}

pred_scaled<-function(result_matrx,tmax,tmin){
  as.matrix((result_matrx - tmin)/(tmax - tmin))
}

get_mode<-function(x){
  as.numeric(names(table(x)))[table(x) == max(table(x))]
}

# data1 row; data2 col
score_matrix<-function(data1,data2,...){
  type = 'local'
  outer(data1,data2, function(x,y) pairwiseAlignment(AAStringSet(x) ,AAStringSet(y),
                                                     substitutionMatrix = "BLOSUM62",
                                                     scoreOnly=T,type = type,...))
}

mNW <- function(dist_matrix, Y, h,cv=F,K = dnorm) {
  Kx <- K(dist_matrix/h)
  # Kx <- apply(dist_matrix, 2, function(x) K(x / h))
  if(cv){
    diag(Kx) = 0
  }
  W <- t(t(Kx) / colSums(Kx)) 
  return(drop(t(Y) %*% W))
}


bw_cv_grid <- function(dist_matrix, Y, original_Y,
                       h.grid = seq(min(dist_matrix),
                                    max(dist_matrix),length.out = 200),
                       method = 'spearman',
                       K = dnorm, plot.cv = FALSE) {
  
  obj <- sapply(h.grid, function(h){
    
    temp = mNW(dist_matrix, Y, h, K = dnorm,cv=T)
    result = cor(temp, original_Y, method = method)
    return(result)})

  h <- h.grid[which.max(obj)]
  if (plot.cv) {
    plot(h.grid, obj, type = "o")
    rug(h.grid)
    abline(v = h, col = 2, lwd = 2)
  }
  h
  
}


# dendrogram<-function(dist_matrix, method = "complete",
#                           members = NULL, ...){
#   dis = as.dist(dist_matrix, diag = FALSE, upper = FALSE)
#   result = hclust(dis, method = method, members = members)
#   plot(result,...)
# }
# 
# draw_dendrogram<-function(data1, ...){
#   self_matrix = score_matrix(data1[,1], data1[,1])
#   dendrogram(- self_matrix,...)
# }


# data1 : dataframe, the first col is pepetide sequence; the second col is score
# use data1 to estimate data2
# use the sequences whose scores are greater than a threshold to estimate
# LOOCV
training<-function(data1, data2 = NULL, threshold = NULL,
                   plot.cv = F,K = dnorm){
  if(is.null(data2)){
    data2 = data1
  }

  # use part of data1 to estimate the whole
  if(!is.null(threshold)){
    # print(1)
    data1 = data1[data1[,2]>threshold,]
  }
  
  self_matrix = score_matrix(data1[,1], data2[,1])
  tmax = max(self_matrix)
  tmin = min(self_matrix)
  hCV = bw_cv_grid(1-self_scaled(self_matrix), data1[,2], data2[,2], plot.cv = plot.cv, K = k)
  kernel_est = mNW(1-self_scaled(self_matrix), Y = data1[,2], h=hCV,cv=T, K = K)
  cv_result = cbind(pep = data2[,1], score = data2[,2], cv_est = kernel_est)
  cv_coe = cor(data2[,2], kernel_est, method = 'spearman')
  config = list(bw = hCV, tmax = tmax, tmin = tmin,K = K)
  return(list(result = cv_result, config = config, threshold = threshold, cv_coe = cv_coe))
}


# fold first, then choose the sequences whose rates are greater than threshold to estmate others
training_fold<-function(data1, k_fold, threshold = NULL, 
                   plot.cv = F, K = dnorm){
  
  n = nrow(data1)
  nk = floor(n/k_fold)
  fold_data = list()
  cur_index = 1:n
  for (k in 1:k_fold) {
    if(k == k_fold){
      temp_index = cur_index
    }else{
      temp_index = sample(cur_index, size = nk, replace = F)
      cur_index =  cur_index[-temp_index]
    }
    fold_data = c(fold_data, list(data1[temp_index,]))
  }
  kk = 1:k_fold
  
  # data1 training; data2 test
  data1 = lapply(1:k_fold,function(x) Reduce(rbind,fold_data[kk[-x]]))
  data2 = fold_data
  # use part of data1 to estimate the whole
  if(!is.null(threshold)){
    # print(1)
    data1 = lapply(data1, function(x) x[x[,2]>threshold,])
  }
  
  self_matrix = lapply(1:k_fold, function(x) score_matrix(data1[[x]][,1], data2[[x]][,1]))
  
  # tmax = mean(sapply(self_matrix,max))
  # tmin = mean(sapply(self_matrix,min))
  tmax = max(sapply(self_matrix,max))
  tmin = min(sapply(self_matrix,min))
  
  h.grid = seq(0.1,1,length.out = 100)
  all_data2 = Reduce(rbind, data2)
  get_coe<-function(h){
    kernel_est = lapply(kk, function(x) 
      mNW(1-self_scaled(self_matrix[[x]]), Y = data1[[x]][,2], h=h,cv=F,K = K))
    est = Reduce(c, kernel_est)
    cv_coe = cor(all_data2[,2], est, method = 'spearman')
    return(cv_coe)
  }
  
  obj <- sapply(h.grid, get_coe)
  
  hCV <- h.grid[which.max(obj)]
  config = list(bw = hCV, tmax = tmax, tmin = tmin,K = K)
  if (plot.cv) {
    plot(h.grid, obj, type = "o")
    rug(h.grid)
    abline(v = hCV, col = 2, lwd = 2)
  }
  kernel_est = lapply(kk, function(x) 
    mNW(1-self_scaled(self_matrix[[x]]), Y = data1[[x]][,2], h=hCV,cv=F,K = K))
  est = Reduce(c, kernel_est)
  cv_result = data.frame(pep = all_data2[,1], score = all_data2[,2], cv_est = est)
  return(list(result = cv_result, config = config, threshold = threshold, cv_coe = max(obj)))
}


# use data1 to pred data2
# data1: dataframe, the first col is pepetide sequence; the second col is score
# data2: a vector contains sequence strings
# use the sequences whose scores are greater than a threshold to predict
pred<-function(data1, data2, k_fold = NULL,threshold = NULL, config=NULL,K = dnorm){
  if(is.null(config)){
    if(is.null(k_fold)){
      temp = training(data1,threshold = threshold,K = K)
    }else{
      temp = training_fold(data1,k_fold =k_fold, threshold = threshold,K = K)
    }
    config = temp$config
    threshold = temp$threshold
    K = temp$K
  }
  if(!is.null(threshold)){
    data1 = data1[data1[,2]>threshold,]
  }
  temp_matrix = score_matrix(data1[,1], data2)
  kernel_est = mNW(1-pred_scaled(temp_matrix,config$tmax,config$tmin), Y = data1[,2], h=config$bw,cv=F,K = K)
  result = data.frame(pep = data2, pred = kernel_est, stringsAsFactors = F)
  return(result)
}


# use data1 to pred data2
# data1: dataframe, the first col is pepetide sequence; the second col is score
# data2: a long AA sequence
# use the sequences whose scores are greater than a threshold to predict
epitope_region<-function(data1, data2, threshold = NULL, config=NULL, k_fold = NULL,
                        main_title=NULL, return_plot = T,
                         region_threshold = NULL,K = dnorm){
  if(is.null(config)){
    if(is.null(k_fold)){
      temp = training(data1,threshold = threshold,K = K)
    }else{
      temp = training_fold(data1,k_fold =k_fold, threshold = threshold,K = K)
    }
    config = temp$config
    threshold = temp$threshold
    K = temp$K
  }
  if(!is.null(threshold)){
    data1 = data1[data1[,2]>threshold,]
  }
  
  # determine the dividing length of data2
  l = get_mode(str_length(data1[,1]))[1]
  # print(paste0('mode=',l))
  # print(l)
  mid = (l+1)%/%2
  
  # divide the long sequence into short ones
  L = str_length(data2)
  seq_array = str_sub(data2,1:L,1:L)
  
  pep = c()
  for (i in 1:(L-l+1)) {
    pep = c(pep,paste0(seq_array[i:(i+l-1)],collapse = ''))
  } 
  
  pred_result = pred(data1, pep, threshold = threshold, config=config)
  prediction = pred_result$pred

  # complete the pep
  if((l+1)%%2==0){
    ttemp = L-mid+1
  }else{
    ttemp = L-mid
  }
  temp1 = c()
  for (i in 1:(mid - 1)) {
    temp1 = c(temp1,paste0(seq_array[1:i],collapse = ''))
  }
  temp2 = c()
  for (i in ttemp:L) {
    temp2 = c(temp2,paste0(seq_array[i:L],collapse = ''))
  }
  pep = c(temp1, pep, temp2)
  prediction = c(rep(mean(head(prediction,mid)),length(temp1)), 
                 prediction, 
                 rep(mean(tail(prediction,mid)),length(temp2)))
  
  final_result = data.frame(pep = pep, pred = prediction, stringsAsFactors = F)
  
  if(return_plot == T){
    axis_theme<-theme(
      legend.title = element_text(size = 30),
      legend.text	= element_text(size = 25),
      legend.spacing= unit(10,"cm"), 
      axis.title=element_text(
        face = "plain", 
        size = 25,
        hjust = .5, 
        vjust = 1, 
        angle = 0 
      ),
      axis.text=element_text(colour="black",size = 17))
    
    if(is.null(region_threshold)){
      region_threshold = median(final_result$pred)
    }
    
    ssize = 1.5
    sssize = 2
    # small_adjust = 0.0001
    result_pred = data.frame(pos=1:length(final_result$pred),pred=final_result$pred)
    
    
    
    f1<-ggplot(result_pred, aes(x=pos, y=pred)) + theme_bw()+axis_theme+
      geom_line(size = ssize,color = '#345d94')+
      geom_line(aes(y=region_threshold),linetype='dashed',color = 'black',size = sssize)+
      geom_ribbon(mapping = aes(ymin = region_threshold, ymax =
                                  ifelse(pred>region_threshold, pred,
                                         region_threshold)
      ), fill = "#8aacdb",outline.type='lower') +
      labs(title = main_title,
           y = 'Prediction Score',
           x = 'Amino Acid Position')+
      theme(plot.title = element_text(hjust = 0.5))
    
    return(list(result=final_result,fig=f1))
  }else{
    return(final_result)
  }
  
}

digit = 3
trunc <- function(data){
  return(signif(data,digits = digit))
}
