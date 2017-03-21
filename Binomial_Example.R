rm(list = ls())
source("/home/peter/Dropbox/Data Science/ML Scripts/Log_Functions.R")
set.seed(13)

mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
admission <- as.matrix(mydata$admit)
beta <- as.matrix(rep(0,3))
vars <- as.matrix(mydata[,2:dim(mydata)[2]])

# My two-class algorithm versus the multiple class algorithm:


L_L(y = admission,X = vars,beta = beta)
LL_Multi(Y = admission,X = vars,beta = beta)


# Using Newton's method:
fromGLM <- glm(formula = admission~vars-1,family = binomial())
binom_algo <- newton(function(input){L_L(y = admission,X = vars,beta = input)},initial_values = beta,tol = 1e-13)
multi_algo <- newton(function(input){LL_Multi(Y = admission,X = vars,beta = input)},initial_values = beta,tol = 1e-13)


# 
# 
# 
# GLMcoefficients <- fromGLM$coefficients
# 
# my_results <- L_L(y,binom_algo,x)
# glm_results <- L_L(y,GLMcoefficients,x)
# 
# L_L_Difference <- my_results$f - glm_results$f
# DL_L_Difference <- t(my_results$df - glm_results$df)%*%(my_results$df - glm_results$df)
# print(data.frame(matrix(GLMcoefficients),binom_algo))
# 


# Example of multi-class regression.

set.seed(13)
X <- as.matrix(iris[1:4])
Y <- as.matrix(ifelse(iris[5]=="setosa",0,ifelse(iris[5]=="versicolor",1,2)))
init_betas <- as.matrix(rep(0,2*ncol(X)))



LL_Multi(Y,X,init_betas)
multi_algo2 <- newton(function(input){LL_Multi(Y = Y,X = X,beta = input)},initial_values = init_betas,tol = 1e-10)










