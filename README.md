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
* RSF_BFS_sample_training.R : 
* dummy_pt_data.csv : 
* utils/RF_performance_metrix.R : 
* utils/rsf_hyperparameter_tuning.R :

2) ./webapp/
* gbm_calculator_codes.R :
* feature_description.csv : 
* finalmodels/OS & finalmodels/PFS :
* www : 

## Authors
  - [YeseulKima](https://github.com/YeseulKima) - **Yeseul Kim** - <yeseulkim@catholic.ac.kr>
  - [wonmo](https://github.com/wonmo) - **Wonmo Sung** - <wsung@catholic.ac.kr>
