# LogCMs

R codes here include the functions in the alternating iterative algorithm for the Logistic Collaborative Model (LogCM) family -- LogCM, Similarity-regularized LogCM (LogSCM) and Pairwise-fusion LogCM (LogPCM).

`main.R` provide an full example of learning the toy dataset `data_sample.csv` with the LogCM family.

`Algorithms.R` includes LogCM and LogSCM functions, `Algorithms_PCM.R` and `auxs.R` include the modifications for LogPCM cases.

Logistic Collaborative Model (LogCM) is a personalized binary classification model, which can learn a distinct model for each individual in a heterogeneous population. 
It is inspired by the work of collaborative learning where canonical models are used to describe the low-rank structure of the personal models. For more information, please check our papers:

[A Learning Framework for Personalized Random Utility Maximization (RUM) Modeling of User Behavior](http://doi.org/10.1109/TASE.2020.3041411)

[A Collaborative Learning Framework for Estimating Many Individualized Regression Models in a Heterogeneous Population](https://doi.org/10.1109/TR.2017.2767941)
