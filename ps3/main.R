# Federico Rodari 2022
# ECON8825 Problem Set 3

### -- Initialize Session ---------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(hrbrthemes))
suppressMessages(library(viridis))
suppressMessages(library(kableExtra))
suppressMessages(library(stargazer))
suppressMessages(library(glmnet))




### -- Setup Directories  ------------------------------------------------------------------

root <- getwd()
script_name <- "het"

# Results
results.dir <- paste0(root,"/_results")


# figures (child)
figures.dir <- paste0(results.dir,
                      "/figures")

dir.create(figures.dir,
           showWarnings = FALSE, recursive = TRUE)

# data (child)
data.dir <- paste0(results.dir,
                   "/data")

dir.create(data.dir,
           showWarnings = FALSE, recursive = TRUE)


digits <- 2


# Load functions
message('\nLoading functions...')

funcs <- list.files(path = paste0(getwd(),"/functions"))

for (f in 1:length(funcs)) {
  
  source(paste0(getwd(), "/functions/", funcs[f]))
  
}

rm(funcs)


# Load Data
names_df <- list.files(path = paste0(getwd(),"/_data"))
names_df

# Data
df_bribes <- read_csv(paste0("_data/",names_df[1]))
cols_bribes <- colnames(df_bribes)

ols <- lm(lba~ tariff_change_post2008 + tariff_change_2008 +  + tariff2007 + differentiated + 
                    agri + perishable + dfs + clear_agent + lvalue_tonnage + day_w_arrival + 
                    psi + monitor + post_2008 + hc_group + hc_4digits + rsa + term +  + diff_post_2008 +
                    agri_post_2008 + lvalue_ton_post_2008 +  + perishable_post_2008 + day_w_arrival_post2008 +
                    dfs_post_2008 + psi_post_2008 + tariff2007_post_2008, 
                    data = df_bribes)
summary(ols)

stargazer(ols, title="Regression Results",
          align=TRUE, dep.var.labels=c("Bribes"), no.space=TRUE)




df_bribes %>%
  mutate(psi = ifelse(psi == "inspected ", 1, ifelse(is.na(psi) == TRUE,NA,0)),
         monitor = ifelse(monitor == "yes ", 1, ifelse(is.na(monitor) == TRUE,NA,0))) %>%
  mutate(inter_diff=tariff_change_post2008*differentiated,
         inter_agri=tariff_change_post2008*agri,
         inter_psi=tariff_change_post2008*psi,
         inter_perish=tariff_change_post2008*perishable,
         inter_dfs=tariff_change_post2008*dfs,
         inter_lvalue=tariff_change_post2008*lvalue_tonnage,
         #inter_dwa=tariff_change_post2008*day_w_arrival,
         inter_rsa=tariff_change_post2008*rsa,
         inter_tariff=tariff_change_post2008*tariff2007,
         inter_monitor=tariff_change_post2008*monitor
        ) -> df_bribes


ols.het <- lm(lba ~ tariff_change_post2008 + inter_diff + inter_agri + inter_lvalue + inter_perish +
                    inter_dfs + inter_psi + inter_rsa + inter_tariff + inter_monitor +
                    tariff_change_2008 + tariff2007 + differentiated + agri + perishable + dfs + 
                    clear_agent +  lvalue_tonnage + day_w_arrival + psi + monitor + post_2008 + 
                    hc_group + hc_4digits + rsa + term + diff_post_2008 + agri_post_2008 + lvalue_ton_post_2008 +
                    perishable_post_2008 + dfs_post_2008 + psi_post_2008 + tariff2007_post_2008 + day_w_arrival_post2008, 
          data = df_bribes)
summary(ols.het)


####
Y = df_bribes$lba








X=cbind(poly(data$age, 6, raw=TRUE), poly(data$inc, 8, raw=TRUE), poly(data$educ, 4, raw=TRUE), poly(data$fsize, 2, raw=TRUE), data$marr, data$twoearn, data$db, data$pira, data$hown)

X=cbind(data$age, data$inc, data$educ, data$fsize, data$marr, data$twoearn, data$db, data$pira, data$hown)


Y= data$net_tfa
D=data$e401

index=t(combn(length(X[1,]),2))
p=X

DMLML(Y,D,p,6)


DMLML=function(Y,D,p,k){
  N=length(Y)
  set.seed(123)
  
  # Create the folds
  dats <- data.frame(p)
  folds <- list() # flexible object for storing folds
  fold.size <- nrow(dats)/k
  remain <- 1:nrow(dats) # all obs are in
  
  for (i in 1:k){
    select <- sample(remain, fold.size, replace = FALSE)
    #randomly sample “fold_size” from the ‘remaining observations’
    
    folds[[i]] <- select # store indices
    #write a special statement for the last fold — if there are ‘leftover points’
    
    if (i == k){
      folds[[i]] <- remain
    }
    
    #update remaining indices to reflect what was taken out
    remain <- setdiff(remain, select)
    remain
  }
  
  
  
  thetaDML=c(0)
  
  for (q in 1:k){
    
    idx <- folds[[q]] #unpack into a vector
    
    ##Fit P(D = 1 | X) using LASSO
    set.seed(333)
    CV=cv.glmnet(p[-idx,],D[-idx],family="binomial",alpha=1)
    fit=glmnet(p[-idx,],D[-idx],family="binomial",alpha=1,lambda=CV$lambda.1se)
    beta1hat=fit$beta
    beta1hat <- as.numeric(as.character(beta1hat))
    
    # predict using new data
    mhat=1/(1+exp(-p[idx,]%*%beta1hat))
    
    #index1=q[which(mhat<0.97 & mhat>0.03)]
    
    ##Estimation of P(D = 1 | X)
    #mhat=1/(1+exp(-p[index1,]%*%beta1hat))
    
    
    # Estimation of g(D,X)
    YY=Y[-idx]
    DD=D[-idx]
    XX=cbind(DD,p[-idx,])
    
    set.seed(333)
    CV=cv.glmnet(XX,YY,family="gaussian",alpha=1)
    fit=glmnet(XX,YY,family="gaussian",alpha=1,lambda=CV$lambda.1se)
    beta2hat=fit$beta
    beta2hat <- as.numeric(as.character(beta2hat))
    
    # predict using new data
    ghat=cbind(D[idx],p[idx,])%*%beta2hat
    
    s = (D[idx]*(ghat) - (1-D[idx])*ghat) + D[idx]*(Y[idx] - ghat)/mhat - (1-D[idx])*(Y[idx] - ghat)/(1-mhat)
    #s=s[which(s<abs(min(s)))]
    
    thetaDML[q]=mean(s)
    
  }
  finaltheta=mean(thetaDML)
  finaltheta
}






