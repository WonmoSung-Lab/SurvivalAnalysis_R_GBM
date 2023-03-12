library(pec)
library(survival)
library(mlr)

# Customized function (to fix an error in calculating ibs using mlr package)
ys_ibrier = makeMeasure(
  id = "ibrier", minimize = TRUE, best = 0, worst = 1,
  properties = c("surv", "req.truth", "req.model", "req.task"),
  name = "Integrated brier score using Kaplan-Meier estimator for weighting",
  note = "Only works for methods for which probabilities are provided via pec::predictSurvProb. Currently these are only coxph and randomForestSRC. To set an upper time limit, set argument max.time (defaults to max time in test data). Implemented in pec::pec",
  fun = function(task, model, pred, feats, extra.args) {
    
    targets = getTaskTargets(task)
    tn = getTaskTargetNames(task)
    f = as.formula(sprintf("Surv(%s, %s) ~ 1", tn[1L], tn[2L]))
    newdata = getTaskData(task)[model$subset, ]
    max.time = 3555 #extra.args$max.time %in% max(newdata[[tn[1L]]])
    grid = seq(0, max.time, length.out = extra.args$resolution)
    
    probs = predictSurvProb(getLearnerModel(model, more.unwrap = TRUE), newdata = newdata, times = grid)
    perror = pec(probs, f,
                 data = newdata[, tn], times = grid, exact = FALSE, exactness = 99L,
                 maxtime = max.time, verbose = FALSE, reference = FALSE)
    crps(perror, times = max.time)[[1]]
  },
  extra.args = list(max.time = NULL, resolution = 1000L)
)

# 3. Hyper parameter tuning
get_HPO = function(train_df0, target_vector, cv, params_grid, measures_list){
  traintask = makeSurvTask(data = train_df0, target = target_vector)
  rdesc = makeResampleDesc("CV",iters=cv, stratify = TRUE)
  multi_ctrl <- makeTuneMultiCritControlGrid(resolution = cv)
  rf.lrn <- makeLearner("surv.randomForestSRC")
  tuneMulti <- tuneParamsMultiCrit(learner = rf.lrn
                                   ,task = traintask
                                   , resampling = rdesc
                                   , measures = measures_list
                                   , par.set = params_grid
                                   ,control = multi_ctrl
                                   ,show.info = T)
  
  best_params <- tuneMulti[["x"]][length(tuneMulti[["x"]])][[1]]
  return(best_params)
}
  
