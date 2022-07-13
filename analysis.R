
# Load function definitions
source("functions.R")

# Set the random seed to a fixed value to allow reproducible results.
set.seed(12345L)

##############################initializations##################################
x <- seq(-5,5,len=1000)
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
initial_design <- sample.int(length(x),1)
iterations <- c(3, 5, 10, 20, 30, 40, 50)


data <- matrix(NA,length(iterations),5)
data[,1] <- iterations
for (i in  1:length(iterations)){
  t0 <- Sys.time()
  c1 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_1, max_iterations=iterations[i])
  data[i,4] <- Sys.time()-t0
  data[i,2] <- sum(c1$var)
  
  t0 <- Sys.time()
  c2 <- seq_DoE(x=x,initial_design=initial_design, y=y,
                criterion=criterion_2, max_iterations=iterations[i])
  data[i,5] <- Sys.time()-t0
  data[i,3] <- sum(c2$var)
}

data <- data.frame(max_iterations=data[,1],
                   variances_1=dx*data[,2],
                   variances_2=dx*data[,3],
                   time_cost_1=data[,4],
                   time_cost_2=data[,5])

ggplot(data) +
  geom_line(aes(x=max_iterations, y=variances_1, col="Variances (criterion 1)")) +
  geom_line(aes(x=max_iterations, y=variances_2, col="Variances (criterion 2)")) +
  geom_line(aes(x=max_iterations, y=time_cost_1, col="Time cost (criterion 1)")) +
  geom_line(aes(x=max_iterations, y=time_cost_2, col="Time cost (criterion 2)")) +
  xlab("iterations") +
  ylab("Variances/Time cost")

saveRDS(data, file = "data/table1.rds")

############################table 2###################################################

initial_design <- c(1,length(x))
iterations <- c(5,10,20,30)

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
                   MMSE=dx*data[,2],
                   IMSE=dx*data[,3])

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
for (i in 1:10){
  df1 <- seq_DoE(x,y,initial_design,criterion_1,sigma2_noise=0.1^2, max_iterations = i)
  var_dat1 <- cbind(var_dat1, df1$var)

  df2 <-seq_DoE(x,y,initial_design,criterion_2,sigma2_noise=0.1^2,max_iterations = i)
  var_dat2 <- cbind(var_dat2, df2$var)
  
}
colnames(var_dat1) <- c("x",1:10)
var_dat1 <- pivot_longer(var_dat1, cols = 2:11, names_to="iterations",values_to='var')
colnames(var_dat2) <- c("x",1:10)
var_dat2 <- pivot_longer(var_dat2, cols = 2:11 ,names_to="iterations",values_to='var')


var_dat1 %>%
  ggplot() +
  geom_line(aes(x,var)) +
  xlab("x") +
  ylab("variance") +
  ggtitle("Prediction variance plot (MMSE)") +
  facet_grid(rows = vars(iterations),scales = "free_y")

var_dat2 %>%
  ggplot() +
  geom_line(aes(x,var)) +
  xlab("x") +
  ylab("variance") +
  ggtitle("Prediction variance plot (IMSE") +
  facet_grid(rows = vars(iterations),scales = "free_y")


saveRDS(var_dat1, file = "data/var_dat1.rds")
saveRDS(var_dat2, file = "data/var_dat2.rds")







#########################################################13/7#####################################

initial_design <- c(1,length(x))
var_dat3 <- data.frame(x=x)
var_dat4 <- data.frame(x=x)
for (i in 1:10){
  df1 <- seq_DoE(x,y,initial_design,criterion_1,sigma2_noise=0.1^2, max_iterations = i)
  var_dat3 <- cbind(var_dat3, df1$var)
  
  df2 <-seq_DoE(x,y,initial_design,criterion_2,sigma2_noise=0.1^2,max_iterations = i)
  var_dat4 <- cbind(var_dat4, df2$var)
  
}
colnames(var_dat3) <- c("x",1:10)
var_dat3 <- pivot_longer(var_dat3, cols = 2:11, names_to="iterations",values_to='var')
colnames(var_dat4) <- c("x",1:10)
var_dat4 <- pivot_longer(var_dat4, cols = 2:11 ,names_to="iterations",values_to='var')


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
