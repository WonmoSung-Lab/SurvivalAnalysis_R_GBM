# 0. Import packages
library(randomForestSRC) # : random forest package
library(caret) # : to stratified k-fold sets
library(survival) # : to make survival object
library(mlr) # : to tune hyper parameters etc.
library(ggRandomForests) # : to get variable importance
library(pec) # : to get IBS value
library(pROC)# : to calculate Cindex
library(tidyverse) # : to edit dataframe. 
library(bigmemory) # : for deepcopy. 

# 2. Import customized modules
source("./utils/RF_performance_metrix.R")
source("./utils/rsf_hyperparameter_tuning.R")

# 1. Load train & test set #######
Xy_df = read.csv('dummy_pt_data.csv')
#######

# 2. Settings #######
RandomSeeds = c(10, 42, 50, 100, 250, 567, 750, 1000, 1234, 2022,
                 5678, 8765, 123456, 234567, 345678, 456789, 567890, 678901, 789012, 890123)
outer_K = 5
inner_K = 4
#######

# 3. Tuning. #######
# Empty space to save the results. #####
randomseeds = list()
out_Ks = list()
tunned_ntrees = list()
tunned_mtrys = list()
cindexs = list()
ibss = list()
numFeats = list()
drop_feats_in_order = list()
######

for (numFeat_i in 1:5){ 
  print(colnames(Xy_df)) 
  for (random_i in 1:length(RandomSeeds)){
    # fix random seed.
    set.seed(RandomSeeds[[random_i]])
    
    # split train & test set. 
    outer_strtf_5folds = caret::createFolds(Xy_df$time, k=outer_K, list = TRUE, returnTrain = FALSE)
    for (out_K_i in 1:outer_K){ # outer loop. #outer_K
      train_df = Xy_df[-outer_strtf_5folds[[out_K_i]], ]
      test_df = Xy_df[outer_strtf_5folds[[out_K_i]], ]
      
      # hyperparameter tuning. (inner_K = 4)
      mtry_n = sqrt(ncol(Xy_df)-2) + 2
      params_grid0 = makeParamSet(makeDiscreteParam("mtry", values=seq(1:mtry_n)),makeDiscreteParam("ntree", values = c(50, 100, 150, 200, 250)))#:sqrt(ncol(Xy_df))+1)),makeDiscreteParam("ntree", values = c(50, 100, 150, 200, 250))) #,2)),makeDiscreteParam("ntree", values = c(10))) #
      best_params = get_HPO(train_df0 = train_df, target_vector = c('time', 'event'), cv = inner_K, params_grid = params_grid0, measures_list = list(mlr::cindex, ys_ibrier))
      
      # train a rsf. 
      rsf_i = rfsrc(formula = Surv(time, event)~., data = train_df,
                    mtry = best_params$mtry, 
                    ntree = best_params$ntree,
                    splitrule = "logrank",
                    importance = "permute",
                    split.depth = "all.trees")
      
      rsf_i_copy = rsf_i
      
      # Get predict results. 
      pred_results = predict.rfsrc(rsf_i, test_df)
      
      # Get & Save performance values. 
      cindexs = append(cindexs, get_Cindex(pred_results = pred_results, test_df = test_df))
      ibss = append(ibss, get_IBS(rfsrc.obj = pred_results))

      # Save additional settings.
      randomseeds = append(randomseeds,RandomSeeds[[random_i]])
      out_Ks = append(out_Ks, out_K_i)
      tunned_ntrees = append(tunned_ntrees, best_params$ntree)
      tunned_mtrys = append(tunned_mtrys, best_params$mtry)
      numFeats = append(numFeats, ncol(Xy_df)-2)
      
      # 6. accumulate VIMP(Permutation importance) values #####
      gg_dta = gg_vimp(rsf_i)
      vimp_vars = gg_dta$vars # get feature names
      vimp_value = gg_dta$vimp # get vimp values
      vimp_df_i0 = data.frame("vars"=vimp_vars, "vimp"=vimp_value) # into data.frame
      vimp_df_i = vimp_df_i0[order(-vimp_df_i0$vimp),]
      vimp_df_i$ranking = seq(1:nrow(vimp_df_i))
      vimp_df_i = subset(vimp_df_i, select= -c(vimp))
      if(random_i == 1 && out_K_i==1){
        vimp_results_df <- vimp_df_i # for the first fold
      }else{
        vimp_results_df <- merge(vimp_results_df,vimp_df_i, by='vars')}
      #######
      
    }}
  # Drop the lowest important feature - by: Permutation importance.
  vimp_results_df$mean_ranking = rowMeans(vimp_results_df[, 2:ncol(vimp_results_df)])
  vimp_results_df = vimp_results_df[order(vimp_results_df$mean_ranking),]
  dropFeat = vimp_results_df[nrow(vimp_results_df),1] # the lowest important feature.
  drop_feats_in_order = append(drop_feats_in_order, dropFeat)

  # Save the vimp results.
  fName = paste0('./rev_OS_Vimp_results_numFeat', as.character(ncol(Xy_df)-2), '.csv')
  write.csv(vimp_results_df, file = fName)

  # Prepare for the next round.
  Xy_df = Xy_df[,!(names(Xy_df) %in% dropFeat)]
}


# Save the overall results.
results_df = data.frame(randomseed=unlist(randomseeds), out_K=unlist(out_Ks), numFeat = unlist(numFeats),
                        tunned_ntree= unlist(tunned_ntrees), tunned_mty= unlist(tunned_mtrys), cindex=unlist(cindexs), ibs=unlist(ibss))



