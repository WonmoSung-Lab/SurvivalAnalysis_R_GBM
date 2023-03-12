# Prognosis prediction for patients with glioblastoma

The official code repository of the survival analysis project.

## Paper
For more details, please see our paper [Prognosis prediction for gioblastoma multiforme patients](https://www.naver.com) which has been accepted at Radiotherapy and Oncology in 2023 March. 
If this code is useful for your work, please consider to cite our paper:
```
@inproceedings{
    kim2023--,
    title={Prognosis prediction for glioblastoma multiforme patients using machine learning approaches: development of the clinically applicable model},
    author={Yeseul Kim and --- and Wonmo Sung},
    booktitle={Radiotherapy and Oncology},
    year={2023},
    url={https://--}
}
```


## Requirements
- randomForestSRC
- caret
- survival
- mlr
- ggRandomForests
- pec
- pROC
- tidyverse


## Description

1) ./randomsurvivalforest_backward_feature_selection/
* RSF_BFS_sample_training.R : A backward feature selection code that uses permutation feature importance as a criterion for selecting features for a random survival forest(RSF) model.
* dummy_pt_data.csv : A dummy patient data set with dummy values for the entire input clinical factors.
* utils/RF_performance_metrix.R : Customized functions to compute performance metrics; Harell's concordance index (C-index) and integrated brier score (IBS).
* utils/rsf_hyperparameter_tuning.R : Customized functions to perform multi-criteria (C-index and IBS) hyper parameter tuning.

2) ./webapp/
* gbm_calculator_codes.R : A Rshiny code to develop our web application.
* finalmodels/OS & finalmodels/PFS : The final 100 cross-validated RSF model for predicting overall survival (OS) and progression-free survival (PFS).

## Authors
  - [YeseulKima](https://github.com/YeseulKima) - **Yeseul Kim** - <yeseulkim@catholic.ac.kr>
  - [wonmo](https://github.com/wonmo) - **Wonmo Sung** - <wsung@catholic.ac.kr>
