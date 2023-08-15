library(mlbench)
data(BreastCancer)
dim(cleandata)
ncol(BreastCancer)
nrow(BreastCancer)
?BreastCancer
?str

head(BreastCancer)
BC1 <-BreastCancer #not necessary, can work with the existing BreastCancer label
head(BC1)

#all data are in factors. must convert to numeric
for(i in 1:10) { #check to see if you need to convert class
  BC1[,i] <- as.numeric(BC1[,i])
}

#convert categorical variables to numeric using dummy variables: 0,1
levels(Class)
BC1$Class <- ifelse(BC1$Class == "malignant", 1, 0)
head(BC1) #check that class is 0,1

#check if all are now numbers.
str(BC1)
unlist(lapply(BC1, class))
#all numeric using both str and unlist

#check for missing values
sum(is.na(BC1))
colSums(is.na(BC1))
is.na(BC1)
BC2 <- na.omit(BC1)
sum(is.na(BC2)) #check to see if all NA values omitted.
is.na(BC2)

summary(BC2)
ncol(BC2)
nrow(BC2) #rows reduced to 683

#check for duplicate entries
duplicated(BC2) # Boolean representation 
sum(duplicated(BC2)) #how many duplicated

BC2[duplicated(BC2)|duplicated(BC2, fromLast=TRUE),]#extracts all duplicated rows

BC25 <- BC2[!duplicated(BC2), ] #bc25 is dataset without duplicates
sum(duplicated(BC25))
nrow(BC25)

#remove ID column
BC3 <- BC2 [,-1]

#PART 2: EXPLORATORY ANALYSIS

#check how many cases of malignancy and benign
table(BC3$Class)

#pairwise scatterplot
pairs(BC3)
mtext("Pairwise Scatterplot of Cleaned BreastCancer Data", side = 3, line = -2, outer = TRUE)

pairs(BC3[,1:9], col=BC3[,10]+1) 
mtext("Pairwise Scatterplot of BreastCancer Data: Malignancy in Red", side = 3, line = -2, outer = TRUE)

pairs(BC3[,1:9], col=BC3[,10]+9, pch=BC3[,10]+2)
mtext("Pairwise Scatterplot of BreastCancer Data", side = 3, line = -1.5, outer = TRUE)

#pairs(BC3, col=BC2[,7]+2)

#add correlation coefficients through a correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("r = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
pairs(BC3, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

#a nicer, easier to read plot
library(corrplot)
#checking correlation between variables
corrplot(cor(BC3), method = "number", type = "upper", diag = FALSE, col = COL2('PiYG'))
corrplot.mixed(cor(BC3), upper = "circle", lower = "number", number.cex = 1.0, upper.col= COL2("PiYG"), lower.col = COL2("PiYG")) 

#mtext("Correlation Panel: Correlation Coefficients", side = 2, line = -6, outer = TRUE)

#PART 3: LOGISTIC REGRESSION

#store rows and columns
n=nrow(BC3)
p=ncol(BC3)-1
p

install.packages("car")
library(car)
model <- lm(Class~.,data=BC3)
model
vif(model)
#logistic regression on all variables 
lr_fit1 <- glm(Class~., data=BC3, family="binomial")
summary(lr_fit1)


#best subset selection AIC and BIC
library(bestglm)
head(BC3) #verify that response variable is the last column
bss_fit_AIC1 <- bestglm(BC3, family=binomial, IC="AIC")
bss_fit_BIC1 <- bestglm(BC3, family=binomial, IC="BIC")
?bestglm
bss_fit_AIC1$Subsets
bss_fit_BIC1$Subsets

#same for BIC
#to easily find best-fitting models
best_AIC1 <- bss_fit_AIC1$ModelReport$Bestk
best_AIC1 #7
best_BIC1 <- bss_fit_BIC1$ModelReport$Bestk
best_BIC1 #5

#best subset selection k-fold cross validation

## Sample the fold-assignment index
nfolds = 10 #10 folds
fold_index = sample(nfolds, n, replace=TRUE)

logistic_reg_fold_error = function(X, y, test_data) {
  Xy = data.frame(X, y=y)
  if(ncol(Xy)>1) tmp_fit = glm(y ~ ., data=Xy[!test_data,], family="binomial")
  else tmp_fit = glm(y ~ 1, data=Xy[!test_data,,drop=FALSE], family="binomial")
  phat = predict(tmp_fit, Xy[test_data,,drop=FALSE], type="response")
  yhat = ifelse(phat > 0.5, 1, 0) 
  yobs = y[test_data]
  test_error = 1 - mean(yobs == yhat)
  return(test_error)
}

general_cv = function(X, y, fold_ind, fold_error_function) {
  p = ncol(X)
  Xy = cbind(X, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  fold_errors = numeric(nfolds)
  # Compute the test error for each fold
  for(fold in 1:nfolds) {
    fold_errors[fold] = fold_error_function(X, y, fold_ind==fold)
  }
  # Find the fold sizes
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  # Compute the average test error across folds
  test_error = weighted.mean(fold_errors, w=fold_sizes)
  # Return the test error
  return(test_error)
}
logistic_reg_bss_cv = function(X, y, fold_ind) {
  p = ncol(X)
  Xy = data.frame(X, y=y)
  X = as.matrix(X)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  fold_errors = matrix(NA, nfolds, p+1) # p+1 because M_0 included in the comparison
  for(fold in 1:nfolds) {
    # Using all *but* the fold as training data, find the best-fitting models 
    # with 0, 1, ..., p predictors, i.e. identify the predictors in M_0, M_1, ..., M_p
    tmp_fit = bestglm(Xy[fold_ind!=fold,], family=binomial, IC="AIC")
    best_models = as.matrix(tmp_fit$Subsets[,2:(1+p)])
    # Using the fold as test data, find the test error associated with each of 
    # M_0, M_1,..., M_p
    for(k in 1:(p+1)) {
      fold_errors[fold, k] = logistic_reg_fold_error(X[,best_models[k,]], y, fold_ind==fold)
    }
  }
  # Find the fold sizes
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  # For models with 0, 1, ..., p predictors compute the average test error across folds
  test_errors = numeric(p+1)
  for(k in 1:(p+1)) {
    test_errors[k] = weighted.mean(fold_errors[,k], w=fold_sizes)
  }
  # Return the test error for models with 0, 1, ..., p predictors
  return(test_errors)
}

#apply cross-validation for best subset selection
cv_bss1 = logistic_reg_bss_cv(BC3[,1:p], BC3[,p+1], fold_index)
cv_bss1 #10 folds, so 10 CV errors, 1 error for each model.
#Identify the number of predictors in the model which minimises test error
best_cv1 = which.min(cv_bss1 - 1) #7
best_cv1

## Create multi-panel plotting device
par(mfrow=c(2, 2))
## Produce plots, highlighting optimal value of k
plot(0:p, bss_fit_AIC1$Subsets$AIC, xlab="Number of predictors", ylab="AIC", type="b")
points(best_AIC1, bss_fit_AIC1$Subsets$AIC[best_AIC1+1], col="green", pch=16)
plot(0:p, bss_fit_BIC1$Subsets$BIC, xlab="Number of predictors", ylab="BIC", type="b")
points(best_BIC1, bss_fit_BIC1$Subsets$BIC[best_BIC1+1], col="green", pch=16)
plot(0:p, cv_bss1, xlab="Number of predictors", ylab="CV error", type="b")
points(best_cv1, cv_bss1[best_cv1+1], col="green", pch=16)

mtext("Optimal value of predictors: AIC, BIC, cross-validation", side = 3, line = -2, outer = TRUE)

pstar = 7
## Check which predictors are in the 7-predictor model
bss_fit_AIC1$Subsets[pstar+1,]

## Construct a reduced data set containing only the 7 selected predictors
bss_fit_AIC1$Subsets[pstar+1, 2:(p+1)]
indices = which(bss_fit_AIC1$Subsets[pstar+1, 2:(p+1)]==TRUE)
BC4 <- BC3[,c(indices, p+1)]
head(BC4)

## Obtain regression coefficients for this model
logreg1_fit1 = glm(Class ~ ., data=BC4, family="binomial")
summary(logreg1_fit1) #shows that all variables are somewhat useful

#test error of 7-predictor model
test_error1 = general_cv(BC4[,1:pstar], BC4[,pstar+1], fold_index, logistic_reg_fold_error)
test_error1 #0.03513909 but later it turned to 0.03221083

#######TESTING
pstar1 = 8
## Check which predictors are in the 7-predictor model
bss_fit_AIC1$Subsets[pstar+1,]

## Construct a reduced data set containing only the 7 selected predictors
bss_fit_AIC1$Subsets[pstar+1, 2:(p+1)]
indices = which(bss_fit_AIC1$Subsets[pstar+1, 2:(p+1)]==TRUE)
BC45 <- BC3[,c(indices,p+1)]

## Obtain regression coefficients for this model
logreg1_fit15 = glm(Class ~ ., data=BC45, family="binomial")
summary(logreg1_fit15) #shows that all variables are somewhat useful

#TRAINING ERROR FOR 7* MODEL
phat3 <-predict(logreg1_fit1, BC4, type="response")
yhat3 <- as.numeric(ifelse(phat3 > 0.5,1,0))
head(yhat3)

#calculate training error for 7* model
#confusion matrix
confusion2 <- table(Observed=BC4$Class, Predicted=yhat3)
confusion2
lg_TE7 <- 1 - sum(diag(confusion2)/sum(confusion2))
1 - mean(BC4$Class == yhat3) #0.03074671 

#calculate test error for 7* model
test_err_full2 = general_cv(BC4[,1:pstar], BC4[,pstar+1], fold_index, logistic_reg_fold_error)
test_err_full2

#test error of 7-predictor model
test_error15 = general_cv(BC45[,1:pstar], BC45[,pstar+1], fold_index, logistic_reg_fold_error)
test_error15 #0.03513909

#TESTING OVER

#building Bayes classifier for LDA
library(MASS)
#Apply  LDA
lda_fit1 <- lda(Class~., data=BC4)
lda_fit1
plot(lda_fit1)

## Compute predicted values:
lda_predict1 = predict(lda_fit1, BC4)
yhat11 = lda_predict1$class
## Calculate confusion matrix:
(confusion = table(Observed=BC4$Class, Predicted=yhat11))

## Calculate training error of LDA:
1 - mean(BC4$Class == yhat11)

#test error using 10-fold cross-evalidation
lda_fold_error = function(X, y, test_data) {
  Xy = data.frame(X, y=y)
  if(ncol(Xy)>1) tmp_fit = lda(y ~ ., data=Xy[!test_data,])
  tmp_predict = predict(tmp_fit, Xy[test_data,])
    yhat = tmp_predict$class 
  yobs = y[test_data]
  test_error = 1 - mean(yobs == yhat)
  return(test_error)
}

lda_test_error = general_cv(BC4[,1:pstar], BC4[,pstar+1], fold_index, lda_fold_error)
lda_test_error

#Building Bayes classifier for QDA

#Apply  QDA to 7-predictor model
qda_fit1 <- qda(Class~., data=BC4)
qda_fit1

## Compute predicted values:
qda_predict1 = predict(qda_fit1, BC4)
yhat2 = qda_predict1$class
## Calculate confusion matrix:
(confusion = table(Observed=BC4$Class, Predicted=yhat2))

## Calculate training error of QDA:
1 - mean(BC4$Class == yhat2)

#modify the lda_fold_error function (which modified the logistic_reg_fold_error)
#test error of QDA using 10-fold cross-evalidation
qda_fold_error = function(X, y, test_data) {
  Xy = data.frame(X, y=y)
  if(ncol(Xy)>1) tmp_fit = qda(y ~ ., data=Xy[!test_data,])
  tmp_predict = predict(tmp_fit, Xy[test_data,])
  yhat = tmp_predict$class 
  yobs = y[test_data]
  test_error = 1 - mean(yobs == yhat)
  return(test_error)
}

qda_test_error = general_cv(BC4[,1:pstar], BC4[,pstar+1], fold_index, qda_fold_error)
qda_test_error

#may not be fully needed
phat1 <-predict(lr_fit1, BC3, type="response")
yhat1 <- as.numeric(ifelse(phat1 > 0.5,1,0))
head(yhat1)
head(BC3$Class)

#calculate training error for full model
#confusion matrix
confusion1 <- table(Observed=BC3$Class, Predicted=yhat1)
confusion1
lg_trainingE <- 1 - sum(diag(confusion1)/sum(confusion1))
1 - mean(BC3$Class == yhat1) #0.03074671

#calculate test error for full model
test_err_full = general_cv(BC3[,1:p], BC3[,p+1], fold_index, logistic_reg_fold_error)
test_err_full #0.03367496

#test_error1 = general_cv(BC4[,1:pstar], BC4[,pstar+1], fold_index, logistic_reg_fold_error)
#test_error1