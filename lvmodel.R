#LV model
library(tidyverse)         
library(doSNOW)
library(doParallel) 
library(data.table)
library(deSolve)
library(foreach)
library(tibble)

#the equation of the model:
# dR/dt = rR – αN*R
# dN/dt = fαNR – qN – eNP
# dP/dt = zeNP – bP

lv3 <- function(t, start, parms) {
  #set the start values for each of the species
  R <- start["R"]                          # basal prey species
  
  N <- start["N"]                          # primary consumer/intermediate predator
  
  P <- start["P"]                        # top predators
  
  # allow R to look within the parameters in the Lv equation.
  with(as.list(parms), {
    
    dR <- r*R - a*N*R                 # dynamics of the resources(dR)
    
    dN <- f*a*N*R - q*N - e*N*P      # dynamics of the intermediate predators(dN)
    
    dP <- z*e*N*P - b*P         # dynamics of the top predators(dP)
    
    list(c(dR, dN, dP))  # returns a list of the abundances each species
  })
}

# set abundance of each species. if abundance falls to 0, species 
#would be considered as extinct


eventFun <- function(t,y,p){
  y <- c(0, 0, 0)         
  return(y)
}

rootfun <- function(t,y,p){
  if(min(y)<1){
    y <- c(0, 0, 0)
  }
  return(y)
}

# place function in a wrapper function to get the dynamics of each species
lv3_lsoda <- function(start, parms, time){
  ##run the simulation
  sim <- as_tibble(lsoda(y = start,
                         times = time,
                         func = lv3,
                         parms = parms,
                         events = list(func = eventFun,
                                       root = TRUE,
                                       terminalroot = 1),
                         rootfun = rootfun))
  
  # reshape data into a long format
  longer_sim <- sim %>% pivot_longer(-time, 
                                     values_to = "abundance",
                                     names_to="species")
  # set the number of years the time series would run for:
  longer_sim$sim_length<-max(longer_sim$time)
  
  # split the results using columns
  longer_sim$parms <- paste(names(parms), parms, sep = "=", collapse = ",")
  
  # a list of the simulation results and the parameters are made and saved as res
  res <- list("sim" = longer_sim,
              "parameters" = parms)
  
  return(res)

}


# assign values to the list of parameters and save as parms
parms <- c(r = 0.5,              # growth rate of prey
           
           a = 0.01,             # feeding rate of N on R
           
           f = 0.001,           # conversion efficiency of R's to new N's
           
           q = 0.00001,         # death rate of N
           
           z = 0.00001,        # conversion efficiency of N's to new P's
           
           b = 0.01,           # decline rate of P
           
           e = 0.001           # feeding rate of P on N
)   

# start values of each species
start <- c(R = 1000,
           
           N = 100,
           
           P = 100)

#set the times for simulation to 100    where t is in years
time <- seq(0, 100, 1)

#run the simulation by passing each parameters to the function
dd <- lv3_lsoda(start, parms, time)

head(dd)  
#plot the output of the Model
lv_mod <- ggplot(data = dd$sim, aes(x=time, y=abundance, col=species))+
  geom_line(size = 3)+
  scale_y_log10() +
  theme_minimal()+
  ggtitle("Lotkha volterra three species Model")
lv_mod

# make a list of sequence from 0.01 to 0.5 which would be replicated for each parameter
parms_list <- rep(list(seq(0.01, 0.5, length.out = 3)), 
                  length(parms))
head(parms_list)
# use expand.grid() to make a data.frame of every possible 
# combination of the 7 parameters
all_var <- expand.grid(parms_list)


head(all_var)     

# view the dimension of all_var
dim(all_var)                  
# there is a total of 2187 simulations from our sequence

# convert list into a list of 2187 simulations
all_var_list <- as.list(as.data.frame(t(all_var)))

# find number of cores
n_cores <- detectCores()      
n_cores 

# define the number of cores to use
clus  <- makeCluster(detectCores() - 1)

# register the cluster in other to have a backend to run parallel
registerDoSNOW(clus)      

#set the maximum number for txtProgressBar
pb <- txtProgressBar(max = length(all_var_list), style = 3)    

#checks the progress of simulations
progress <- function(n) setTxtProgressBar(pb, n) 

#set a progressbar
opts <- list(progress = progress)

# iterate over all_var_List
result <- foreach(i = all_var_list,  
                  .combine = "rbind",
                  .packages = c("dplyr","foreach", "deSolve", "parallel" ,"tidyverse", "doSNOW"),
                  .options.snow = opts,   # to display progress bar
                  .errorhandling = "pass") %dopar% {
                    name_vars = names(parms)  # rename parms new list
                    names(i) <- name_vars
                    lv3_res <- lv3_lsoda(start = start, parms = i, time = time)
                    return(lv3_res)
                  }

head(result) 

# convert simulation into a dataframe in other to plot the model
result <- as.data.frame(t(result))

# create a dataframe using rbindlist to visualise the model further
extract_res_sim <-rbindlist(lapply(result, function(h){
  return(h$sim)
}))

#plot
ggplot(extract_res_sim, aes(x = time,
                            y = abundance,
                            group = parms))+
  geom_line(col = "white", alpha =0.3) +
  facet_wrap(~species)+
  #log of the y axis
  scale_y_log10() +
  
  theme_dark()

#stop cluster
stopCluster(clus)
