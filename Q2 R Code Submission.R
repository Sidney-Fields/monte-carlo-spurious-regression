#EC2062 Programming Project

set.seed(12345)

library(ggplot2)

T<- 500
r <- 15000

phi2 <- 0.7
phi1_vals <- c(0.5,0.6,0.7,0.8,0.9,0.95,0.97,0.985)
phi_valslength<- length(phi1_vals)
inpt<- rep(1,T)
e<- rnorm(T)
x<- matrix(0, nrow = T, ncol = 1)
for (i in 2:T) {
  x[i]<- phi2 * x[i-1]+e[i]  
}

phi1<- 0.5
e1<- rnorm(T)
y<- matrix(0, nrow = T, ncol = 1)
for (i in 2:T) {
  y[i]<- phi1 * y[i-1] + e1[i]  
}

X<- cbind(inpt, x)
k <- 2

M<- solve(t(X)%*%X)
b <- as.numeric(M %*% t(X) %*% y)

uhat<- y - X %*% b
s2<- t(uhat) %*% uhat / (T-k)

v_b <- s2 * M[k,k] 
t_ratio <- b[k] / sqrt(v_b)

results<- matrix(0, nrow = phi_valslength, ncol = 3)
colnames(results) <- c("Mean_beta", "Var_Beta", "Rejection_Rate")

t_05 <- matrix(0, nrow = r, ncol = 1)
t_0985 <- matrix(0, nrow = r, ncol = 1)

for (p in 1:phi_valslength)
{
  phi1<- phi1_vals[p]
  
  beta_store <- matrix(0, nrow = r, ncol = 1)
  t_store <- matrix(0, nrow = r, ncol = 1)
  reject <- matrix(0, nrow = r, ncol = 1)
  
  for (j in 1:r) {
    e<- rnorm(T)
    x<- matrix(0, nrow = T, ncol = 1)
    for (i in 2:T) 
    {
      x[i] <- phi2 * x[i-1] + e[i]
    }
    e1<- rnorm(T)
    y<- matrix(0, nrow = T, ncol = 1)
    for (i in 2:T) 
    {
      y[i]<- phi1 * y[i-1] + e1[i]
    }
    X <- cbind(inpt, x)
    M <- solve(t(X)%*%X)
    b <- as.numeric(M %*% t(X) %*% y)
    
    uhat<- y - X%*% b
    s2 <- t(uhat) %*% uhat / (T-k)
    
    v_b<- s2 * M[k,k]
    t_ratio <- as.numeric(b[k] / sqrt(v_b))
    
    beta_store[j]<- b[k]
    t_store[j]<- t_ratio
    
    if (abs(t_ratio) > 1.96)
    {
      reject[j]<-1 
    }
  }
  results[p,1] <- mean(as.vector(beta_store))
  results[p,2] <- var(as.vector(beta_store))
  results[p,3] <- mean(as.vector(reject))
  
  if (phi1 == 0.5)
  {
    t_05<-t_store
  }
  
  if (phi1 == 0.985)
  {
    t_0985<- t_store
  }
}

 t_05_vec<- as.vector(t_05)
 t_0985_vec<- as.vector(t_0985) 
 t_all <- c(t_05_vec, t_0985_vec)
 group <- c(rep("phi1 = 0.5", r),rep("phi1 = 0.985", r))
 t_data <-data.frame(t_value = t_all, group = group)
 
 print(ggplot(t_data, aes(x = t_value, fill = group)) +
   geom_density(alpha = 0.5) +
   labs(title = "Distribution of t-ratios under different persistence levels",
     x = "t-ratio",
     y = "Density"))
 