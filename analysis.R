
# Load function definitions
source("functions.R")

# Set the random seed to a fixed value to allow reproducible results.
set.seed(12345L)

##############################initializations##################################
x <- seq(-10,10,len=1000)
dx <- (x[length(x)]-x[1])/length(x)
y <- f(x)
y_hat <- rep(NA, length(x))

# Initial mean
mu <- rep(0,length(x))
mu_new <-mu

# Initial covariance:
Sigma <- Matern_k(x,x)
Sigma_new <- Sigma

save(x,dx,y,y_hat,mu,mu_new,Sigma,Sigma_new, file="initializations.RData")

##############################table 1#################################################
initial_design <- c(1,length(x))
iterations <- seq(10,40,len=11)


data <- matrix(NA,length(iterations),13)
data[,1] <- iterations
for (i in  1:length(iterations)){
  
  obs_sf <- rep(FALSE, length(x))
  obs_sf[seq(1,length(x),len=iterations[i])] <- TRUE
  t0 <- Sys.time()
  c0 <- update_estimate(mu= rep(0,length(y)), Sigma=Matern_k(x,x), sigma2_noise = 0.01,
                        y = y[obs_sf], obs = obs_sf)
  data[i,11] <- Sys.time()-t0
  data[i,12] <- mean(SE_score(y_pred=c0$mu[!obs_sf],y_true=y[!obs_sf]))
  data[i,13] <- mean(DS_score(y_pred=c0$mu[!obs_sf],y_true=y[!obs_sf],
                              sigma2=pmax(0,(diag(c0$Sigma))[!obs_sf])))
  
  
  t0 <- Sys.time()
  c1 <- seq_DoE(x=x,initial_design=initial_design, y=y,  sigma2_noise = 0.01,
                criterion=criterion_MMSE, max_iterations=iterations[i]-2)
  data[i,8] <- Sys.time()-t0
  data[i,2] <- mean(SE_score(y_pred=c1$y_hat[!c1$obs],y_true=c1$y[!c1$obs]))
  data[i,5] <- mean(DS_score(y_pred=c1$y_hat[!c1$obs],y_true=c1$y[!c1$obs],
                             sigma2=c1$var[!c1$obs]))
  
  t0 <- Sys.time()
  c2 <- seq_DoE(x=x,initial_design=initial_design, y=y, sigma2_noise = 0.01,
                criterion=criterion_IMSE, max_iterations=iterations[i]-2)
  data[i,9] <- Sys.time()-t0
  data[i,3] <- mean(SE_score(y_pred=c2$y_hat[!c2$obs],y_true=c2$y[!c2$obs]))
  data[i,6] <- mean(DS_score(y_pred=c2$y_hat[!c2$obs],y_true=c2$y[!c2$obs],
                             sigma2=c2$var[!c2$obs]))
  
  t0 <- Sys.time()
  c3 <- seq_DoE(x=x,initial_design=initial_design,y=y, sigma2_noise = 0.01,
                criterion=criterion_IMDS, max_iterations=iterations[i]-2)
  data[i,10] <- Sys.time()-t0
  data[i,4] <- mean(SE_score(y_pred=c3$y_hat[!c3$obs],y_true=c3$y[!c3$obs]))
  data[i,7] <- mean(DS_score(y_pred=c3$y_hat[!c3$obs],y_true=c3$y[!c3$obs],
                             sigma2=c3$var[!c3$obs]))
}

data <- data.frame(max_iterations=data[,1],
                   SE_MMSE=data[,2],
                   SE_IMSE=data[,3],
                   SE_IMDS=data[,4],
                   SE_USF=data[,12],
                   DS_MMSE=data[,5],
                   DS_IMSE=data[,6],
                   DS_IMDS=data[,7],
                   DS_USF=data[,13],
                   time_cost_MMSE=data[,8],
                   time_cost_IMSE=data[,9],
                   time_cost_IMDS=data[,10],
                   time_cost_USF=data[,11])

saveRDS(data, file = "data/table1.rds")

############################table 2###################################################

initial_design <- c(1,length(x))
iterations <- c(5,10,20,30)

data <- matrix(NA,4,4)
data[,1] <- iterations
for (i in  1:length(iterations)){
  c1 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_MMSE,sigma2_noise=0.1^2, max_iterations=iterations[i])
  data[i,2] <- sum(c1$var)
  
  c2 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_IMSE, sigma2_noise=0.1^2,max_iterations=iterations[i])
  data[i,3] <- sum(c2$var)
  
  c3 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_IMDS, sigma2_noise=0.1^2,max_iterations=iterations[i])
  data[i,4] <- sum(c3$var)
}

data <- data.frame(max_iterations=data[,1],
                   MMSE=dx*data[,2],
                   IMSE=dx*data[,3],
                   IMDS=dx*data[,4])

saveRDS(data, file = "data/table2.rds")

##########################table 3####################################################

initial_design <- c(length(x)/2)
iterations <- c(3,5,30,50)

data <- matrix(NA,4,3)
data[,1] <- iterations
for (i in  1:length(iterations)){
  c1 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_1, max_iterations=iterations[i])
  data[i,2] <- sum(c1$var)
  
  c2 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_2, max_iterations=iterations[i])
  data[i,3] <- sum(c2$var)
}

data <- data.frame(max_iterations=data[,1],
                   criterion1=dx*data[,2],
                   criterion2=dx*data[,3])

saveRDS(data, file = "data/table3.rds")

###########################table 4##################################################

data <- matrix(NA,51,5)

for (i in seq_len(nrow(data))){
  n <- floor((i-1)/(nrow(data)-1)*(length(x)-1)+1)
  initial_design <- c(n)
  data[i,1] <- n
  c1 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_1, max_iterations=20)
  data[i,2] <- dx*sum(c1$var)
  data[i,3] <- max(c1$va)
  c2 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_2, max_iterations=20)
  data[i,4] <- dx*sum(c2$var)
  data[i,5] <- max(c2$var)
  
}


data <- as.data.frame(data)
colnames(data) <- c("initial_obs","max_int", "max_max", "int_int", "int_max")

saveRDS(data, file = "data/table4.rds")



############################################13/7#######################################

initial_design <- c(100)
var_dat1 <- data.frame(x=x)
var_dat2 <- data.frame(x=x)
var_datDS <- data.frame(x=x)
for (i in 1:10){
  df1 <- seq_DoE(x,y,initial_design,criterion_MMSE,sigma2_noise=0.1^2, max_iterations = i)
  var_dat1 <- cbind(var_dat1, df1$var)

  df2 <- seq_DoE(x,y,initial_design,criterion_IMSE,sigma2_noise=0.1^2,max_iterations = i)
  var_dat2 <- cbind(var_dat2, df2$var)
  
  df3 <- seq_DoE(x,y,initial_design,criterion_IMDS,sigma2_noise=0.1^2,max_iterations = i)
  var_datDS <- cbind(var_datDS, df3$var)
  
}
colnames(var_dat1) <- c("x",1:10)
var_dat1 <- pivot_longer(var_dat1, cols = 2:11, names_to="iterations",values_to='var')
colnames(var_dat2) <- c("x",1:10)
var_dat2 <- pivot_longer(var_dat2, cols = 2:11 ,names_to="iterations",values_to='var')
colnames(var_datDS) <- c("x",1:10)
var_datDS <- pivot_longer(var_datDS, cols = 2:11, names_to="iterations",values_to='var')


saveRDS(var_dat1, file = "data/var_dat1.rds")
saveRDS(var_dat2, file = "data/var_dat2.rds")
saveRDS(var_datDS, file = "data/var_datDS.rds")






#########################################################13/7#####################################

initial_design <- c(1,length(x))
var_dat3 <- data.frame(x=x)
var_dat4 <- data.frame(x=x)
var_datDS2 <- data.frame(x=x)
for (i in 1:10){
  df1 <- seq_DoE(x,y,initial_design,criterion_MMSE,sigma2_noise=0.1^2, max_iterations = i)
  var_dat3 <- cbind(var_dat3, df1$var)
  
  df2 <-seq_DoE(x,y,initial_design,criterion_IMSE,sigma2_noise=0.1^2,max_iterations = i)
  var_dat4 <- cbind(var_dat4, df2$var)
  
  df3 <- seq_DoE(x,y,initial_design,criterion_IMDS,sigma2_noise=0.1^2,max_iterations = i)
  var_datDS2 <- cbind(var_datDS2, df3$var)
  
}
colnames(var_dat3) <- c("x",1:10)
var_dat3 <- pivot_longer(var_dat3, cols = 2:11, names_to="iterations",values_to='var')
colnames(var_dat4) <- c("x",1:10)
var_dat4 <- pivot_longer(var_dat4, cols = 2:11 ,names_to="iterations",values_to='var')
colnames(var_datDS2) <- c("x",1:10)
var_datDS2 <- pivot_longer(var_datDS2, cols = 2:11, names_to="iterations",values_to='var')

var_dat3 %>%
  ggplot() +
  geom_line(aes(x,var)) +
  xlab("x") +
  ylab("variance") +
  ggtitle("Prediction variance plot (MMSE)") +
  facet_grid(rows = vars(iterations),scales = "free_y")

var_dat4 %>%
  ggplot() +
  geom_line(aes(x,var)) +
  xlab("x") +
  ylab("variance") +
  ggtitle("Prediction variance plot (IMSE") +
  facet_grid(rows = vars(iterations),scales = "free_y")


saveRDS(var_dat3, file = "data/var_dat3.rds")
saveRDS(var_dat4, file = "data/var_dat4.rds")
