#' MARSGWR: a hybrid model that uses the MARS model for important variable selection and the GWR model for prediction at an unknown location based on the selected variables.
#'
#' @param sp_data A dataframe containing the response variables and the predictor variable, as well as the coordinates of the locations. In the dataframe, first column is the response variable (y), last two columns are coordinates i.e., Latitude and Longitudes and in between them is the set of predictor variables(X's).
#' @param bw A numeric value specifying the bandwidth parameter for the GWR model. It can be noted that, optimum bandwidth value can vary depending on the specific dataset and bandwidth parameter depends on the spatial pattern of the data
#' @param deg The degree of interactions to be considered in the MARS model
#' @param sv  Splitting value for dividing the dataset into training and testing set, e.g. 0.8 or 0.7
#' @param gaussian_kernel Spatial weight function of the GWR model, e.g. gaussian_kernel
#' @return A list with the following components:
#'   - `Selected_variables`: The selected variables from the MARS model
#'   - `GWR_y_pred_train`: The GWR predictions at the training locations
#'   - `GWR_y_pred_test`: The GWR predictions at testing locations
#'   - `In_sample_accuracy`: In sample accuracy measures
#'   - `Out_of_sample_accuracy`: Out of sample accuracy measures
#' @examples
#' n<- 100
#' p<- 5
#' m<-sqrt(n)
#' id<-seq(1:n)
#' x<-matrix(runif(n*p), ncol=p)
#' e<-rnorm(n, mean=0, sd=1)
#' xy_grid<-expand.grid(c(1:m),c(1:m))
#' Latitude<-xy_grid[,1]
#' Longitude<-xy_grid[,2]
#' B0<-(Latitude+Longitude)/6
#' B1<-(Latitude/3)
#' B2<-(Longitude/3)
#' B3<-(2*Longitude)
#' B4<-2*(Latitude+Longitude)/6
#' B5<-(4*Longitude/3)
#' y<-B0+(B1*x[,1])+(B2*x[,2])+(B3*x[,3])+(B4*x[,4])+(B5*x[,5])+ e
#' sp_data<-data.frame(y,x,Latitude,Longitude)
#' MARSGWR_gau<-MARSGWR_gaussian(sp_data,5,3,0.7,gaussian_kernel)
#' @references
#' 1. Friedman, J.H. (1991).Multivariate Adaptive Regression Splines. Ann. Statist. 19(1),1-67. <DOI:10.1214/aos/1176347963>.
#' 2. Brunsdon, C., Fotheringham, A.S. and Charlton, M,E. (1996).Geographically weighted regression: a method for exploring spatial non-stationarity. Geogr Anal.28(4),281-298.<DOI:10.1111/j.1538-4632.1996.tb00936.x>.
#' @export
#' @import qpdf
#' @import numbers
#' @import earth
#' @importFrom stats cor dist

MARSGWR_gaussian<- function(sp_data,bw,deg,sv,gaussian_kernel) {

# Step1: Generation of simulated spatial population with spatial coordinates
  sp_data<-as.data.frame(sp_data)
# Step2: Split data into training and testing sets
train_index <- sample(nrow(sp_data), sv * nrow(sp_data)) # Indices of training observations
train_data <- sp_data[train_index,]
nrow(train_data)
test_data <- sp_data[-train_index, ]
nrow(test_data)
x_train <-subset(train_data[,-c(1,ncol(train_data)-1,ncol(train_data))])
y_train <-train_data[,1]
x_test <- subset(test_data,select = -c(1,ncol(test_data)-1,ncol(test_data)))
y_test <- test_data[,1]
# Step3: Important variable selection using MARS model
train_data_MARS<-subset(train_data,select = -c(ncol(train_data)-1,ncol(train_data)))
MARS_fit_train<-earth(train_data_MARS[,1]~.,data=train_data_MARS[,-1],degree= deg)

# Step4 : Select the important variables
imp_vars_MARS<-evimp(MARS_fit_train)
imp_vars_MARS

## Important variable selection for GWR model fitting

nos_imp_var<-length(imp_vars_MARS[,1])
nos_imp_var

for(i in 1:nos_imp_var)
{
  mars_imp_X_vars<-x_train[,imp_vars_MARS[,1]]
}
mars_imp_X_vars

# Step5 : Fitting of the GWR model on training data


coords_train<-cbind(train_data[,ncol(train_data)-1],train_data[,ncol(train_data)])
dists<- as.matrix(dist(coords_train))

# Define the Gaussian kernel function
gaussian_kernel <- function(dists, bw) {
  exp(-0.5 * (dists/bw)^2)
}
weights <- gaussian_kernel(dists, bw)
dim(weights)
Xt <- as.matrix(cbind(rep(1, nrow(x_train)),mars_imp_X_vars))
dim(Xt)
y_tr<- matrix(y_train)
dim(y_tr)
W.train<-list()
for (i in 1:ncol(weights)){
  t <- diag(weights[,i],nrow=nrow(x_train), ncol=nrow(x_train))
  W.train[[i]]<-t
}
W.train

Beta.train<-list()
for (i in 1:length(W.train)){
  lm<- solve((t(Xt)%*%W.train[[i]]%*%Xt))%*%t(Xt)%*%W.train[[i]]%*%y_tr
  Beta.train[[i]]<-lm
}
Beta.train

X.tr_row<-list()
for(i in 1:nrow(Xt)){

  X.tr_row[[i]]<-(Xt[i,])
}


y_hat.train<-mapply("%*%", X.tr_row,Beta.train)
y_hat.train


# Step6 : Make predictions at the test locations

for(i in 1:nos_imp_var)
{
  mars_imp_X_test<-x_test[,imp_vars_MARS[,1]]
}
mars_imp_X_test

coords_test<-cbind(test_data[,ncol(test_data)-1],test_data[,ncol(train_data)])
dists_test <- as.matrix(dist(coords_test))
weights_test <- gaussian_kernel(dists_test, bw)
dim(weights_test)

Xtest <- as.matrix(cbind(rep(1, nrow(x_test)),mars_imp_X_test))
dim(Xtest)
ytest<- matrix(y_test)
dim(ytest)

W.test<-list()
for (i in 1:ncol(weights_test)){
  test <- diag(weights_test[,i],nrow=nrow(x_test), ncol=nrow(x_test))
  W.test[[i]]<-test
}
W.test

Beta.test<-list()
for (i in 1:length(W.test)){
  lm_test<- solve((t(Xtest)%*%W.test[[i]]%*%Xtest))%*%t(Xtest)%*%W.test[[i]]%*%ytest
  Beta.test[[i]]<-lm_test
}
Beta.test

X.test_row<-list()
for(i in 1:nrow(Xtest)){

  X.test_row[[i]]<-(Xtest[i,])
}


y_hat.test<-mapply("%*%", X.test_row,Beta.test)
y_hat.test

# Step7 : Create summary output

# Compute in-sample accuracy measures

rrmse_train <- sqrt(mean((y_hat.train - y_train)^2)) / mean(y_train)
mae_train <- mean(abs(y_hat.train - y_train))
mse_train <- mean((y_hat.train - y_train)^2)
r2_train <- cor(y_hat.train, y_train)^2

# Compute out-of-sample accuracy measures
rrmse_test <- sqrt(mean((y_hat.test - y_test)^2)) / mean(y_test)
mae_test <- mean(abs(y_hat.test - y_test))
mse_test <- mean((y_hat.test - y_test)^2)
r2_test <- cor(y_hat.test, y_test)^2

# Summary output includes selected variables, GWR based predicted values at both training and testing locations and accuracy measures
Summary_MARSGWR_output <-list(Selected_variables = imp_vars_MARS[,1], GWR_y_pred_train =y_hat.train, GWR_y_pred_test = y_hat.test,
                              In_sample_accuracy = c(RRMSE = round(rrmse_test, 4),
                                                     MAE = round(mae_train, 4),
                                                     MSE = round(mse_train, 4),
                                                     R2 = round(r2_train, 4)),
                              Out_of_sample_accuracy = c(RRMSE = round(rrmse_test, 4),
                                                         MAE =  round(mae_test, 4),
                                                         MSE = round(mse_test, 4),
                                                         R2 = round(r2_test, 4)))
return(Summary_MARSGWR_output)

}
