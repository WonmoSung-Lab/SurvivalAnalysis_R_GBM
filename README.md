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

- mlr

## Description

* webapp
* randomsurvivalforest

1 .To train `CGS` for porous network problems (Section 6.1)

```console
python pn_train.py
```

2. To train `CGS` for graph value iteration problems (Section 6.2)

```console
python gvi_train.py
```

3. To train `CGS` for graph benchmark problems (Section 6.3)

```console
python benchmark_train.py
```



