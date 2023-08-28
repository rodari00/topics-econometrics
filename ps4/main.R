# Federico Rodari 2022
# ECON8825 Problem Set 4

### -- Initialize Session ---------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(hrbrthemes))
suppressMessages(library(viridis))
suppressMessages(library(kableExtra))
suppressMessages(library(rdrobust))
suppressMessages(library(ggthemes))


# Settings:
# use mimicking variance criterion for the number of bins


### -- Setup Directories  ------------------------------------------------------------------

root <- setwd("C:/Users/feder/Dropbox/Github/topics-econometrics/ps4")
script_name <- "ps4"

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

# Load functions
message('\nLoading functions...')

funcs <- list.files(path = paste0(getwd(),"/functions"))

for (f in 1:length(funcs)) {
  
  source(paste0(getwd(), "/functions/", funcs[f]))
  
}

rm(funcs)



digits <- 2


# Load Data
names_df <- list.files(path = paste0(root,"/_data"))

# Probation
df <- read_csv(paste0("_data/",names_df[1]))
cols <- colnames(df)


### -- Descriptive analysis ----------------------------

# Plot left_school

df %>%
  mutate(dist_from_cut = round(dist_from_cut, 2)) %>%
  group_by(dist_from_cut) %>%
  summarise(mean_left_school = mean(left_school)) %>%
  #filter(dist_from_cut <= 1.5) %>%
  ggplot(aes(x = dist_from_cut, y = mean_left_school)) +
  geom_point(shape = 21, colour = "#494949", fill = "white", size = 2, stroke = 1) +
  geom_vline(xintercept = 0, color= "darkgrey", linetype = "dotted", linewidth =1.1)  +
  theme_minimal() +
  scale_x_continuous(breaks=seq(-1.5,1.5,0.5)) +
  labs(x = " First year GPA minus probation cutoff",
       y = "Mean Left School") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
saveFigure("desc_left_school.pdf","pdf",9,6,figures.dir)  

# Plot nextGPA

df %>%
  mutate(dist_from_cut = round(dist_from_cut, 2)) %>%
  group_by(dist_from_cut) %>%
  summarise(mean_nextGPA = mean(nextGPA, na.rm = TRUE)) %>%
  #filter(dist_from_cut <= 1.5) %>%
  ggplot(aes(x = dist_from_cut, y = mean_nextGPA)) +
  geom_point(shape = 21, colour = "#494949", fill = "white", size = 2, stroke = 1) +
  geom_vline(xintercept = 0, color= "darkgrey", linetype = "dotted", linewidth = 1.1) +
  theme_minimal() +
  scale_x_continuous(breaks=seq(-1.5,1.5,0.5)) +
  labs(x = " First year GPA minus probation cutoff",
       y = "Mean Next GPA") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

saveFigure("desc_nextGPA.pdf","pdf",9,6,figures.dir)  

### -- Replicate Lindo (2010) using their specification ----------------------------


# 1) left_school -----------------------------------------
Y = df$left_school
X = df$dist_from_cut

# 1.1) Uniform Kernel
T = as.integer( X<= 0) # treatment
T_X = X*T
lm_lindo10_left_school <- lm(Y[X >= -0.6 & X <= 0.6] ~ X[X >= -0.6 & X <= 0.6]+ T[X >= -0.6 & X <= 0.6]+
                                    T_X[X >= -0.6 & X <= 0.6])

lindo10_left_school <- data.frame("First year GPA < cutoff" = c(as.numeric(lm_lindo10_left_school$coefficients[3]),
                                                                as.numeric(sqrt(diag(vcov(lm_lindo10_left_school)))[3])),
                                  "Constant (control mean)" = c(as.numeric(lm_lindo10_left_school$coefficients[2]),
                                                                as.numeric(sqrt(diag(vcov(lm_lindo10_left_school)))[2]))
)

row.names(lindo10_left_school) <- c("Estimate", "SE")

lindo10_left_school %>%
  kbl(caption="Estimated Discontinuity in Left School",
      format="latex",
      col.names = c("First year GPA < cutoff","Constant (control mean)"),
      align="l", digits = 3) %>%
  kable_minimal(full_width = F) %>%
  save_kable(file = paste0(data.dir,"/lindo2010_left_school.tex"))

# # Triangular Kernel
# w = NA
# w[X <= 0 & X >= -0.6] = 1 - abs(X[X <= 0 & X >= -0.6]/0.6) # Weights for treatment
# w[X > 0 & X <= 0.6] = 1 - abs(X[X > 0 & X <= 0.6]/0.6)
# 
# 
# # Control
# out <- lm(Y[X > 0]~ X[X > 0], weights = w[X > 0])
# right_intercept = out$coefficients[1]
# # Treatment
# out <- lm(Y[X <= 0]~ X[X <= 0], weights = w[X <= 0])
# left_intercept = out$coefficients[1]
# 
# left_intercept - right_intercept


# 2) nextGPA --------------------------------
Y = df$nextGPA
X = df$dist_from_cut

# 2.1) Uniform Kernel

T = as.integer(X< 0) # treatment
T_X = X*T
lm_lindo10_nextGPA <- lm(Y[X >= -0.6 & X <= 0.6] ~ X[X >= -0.6 & X <= 0.6]+ T[X >= -0.6 & X <= 0.6]+
                               T_X[X >= -0.6 & X <= 0.6])

lindo10_nextGPA <- data.frame("First year GPA < cutoff" = c(as.numeric(lm_lindo10_nextGPA$coefficients[3]),
                                                                as.numeric(sqrt(diag(vcov(lm_lindo10_nextGPA)))[3])),
                                  "Constant (control mean)" = c(as.numeric(lm_lindo10_nextGPA$coefficients[2]),
                                                                as.numeric(sqrt(diag(vcov(lm_lindo10_nextGPA)))[2]))
)

row.names(lindo10_nextGPA) <- c("Estimate", "SE")

lindo10_nextGPA %>%
  kbl(caption="Estimated Discontinuity in Next GPA",
      format="latex",
      col.names = c("First year GPA < cutoff","Constant (control mean)"),
      align="l", digits = 3) %>%
  kable_minimal(full_width = F) %>%
  save_kable(file = paste0(data.dir,"/lindo2010_nextGPA.tex"))


 

### -- Replicate Lindo (2010) with rdrobust + alternatives  ---------------------

# Using RDROBUST (need to flip X to let the code interpret the standard X >= c)
# kernel <- "uniform"
# binwidth <- 0.6
# poly <- 1
# message(paste0("Estimating LATE", 
#                "\nkernel = " , kernel,
#                "\nbinwidth = ", binwidth,
#                "\npoly = ", poly))
# 
# out <- rdrobust( Y,-X,
#                  kernel = kernel,
#                  p = poly,
#                  h = binwidth)
# 
# toprint <- data.frame("beta_left" = out$beta_Y_p_l,
#            "beta_right" = out$beta_Y_p_r,
#            "tau" = c(out$Estimate[1], out$Estimate[3]))
# row.names(toprint) <- c("Estimate", "SE")
# 
# print(toprint)


# We can approximate the two CE functions for potential outcomes using a linear
# local polynomial approximation. 

# Step 1: Plot RD effect as in Lindo (2010)

# Variable: left_school

Y = df$left_school

left_school_rdplot <- rdplot(Y, X, p = 1, kernel = "uniform", h = 0.6, binselect = "esmv", col.dots = "darkgrey", masspoints	
                              = "off")

left_school_rdplot$rdplot +
  theme_minimal() +
  xlim(-0.6,0.6) +
  xlab("First Year GPA Minus Probation Cutoff") +
  ylab("Left University Voluntarily") +
  labs(title = "") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

saveFigure("rd_lindo10_left_school.pdf","pdf",9,6,figures.dir)  

# Variable: nextGPA
Y = df$nextGPA

nextGPA_rdplot <- rdplot(Y, X, p = 1, h = 0.6, kernel = "uniform", binselect = "esmv", col.dots = "darkgrey", masspoints	
                         = "off")

nextGPA_rdplot$rdplot +
  theme_minimal() +
  xlim(-0.6,0.6) +
  xlab("First Year GPA Minus Probation Cutoff") +
  ylab("Next Year GPA") +
  labs(title = "") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

saveFigure("rd_lindo10_nextGPA.pdf","pdf",9,6,figures.dir) 


# Step 2: Estimation of both specifications
varlist <- c("left_school",
             "nextGPA")

varnames <- c("Left Voluntarily School",
              "Next Year GPA")

if (file.exists(paste0(data.dir,"/rd-estim.txt"))){
  file.remove(paste0(data.dir,"/rd-estim.txt"))
}


for(i in 1:length(varlist)){
  
  for(k in c("lindo", "optimal")){
    
    message(paste0("Running robustness check for ", varlist[i],"...\n"))
    
    Y = df[[varlist[i]]] 
    if(k == "lindo"){
      rd_output <- rdrobust(Y, -X, p = 1, h = 0.6, kernel = "uniform", masspoints	= "off")
    } else{
      rd_output <- rdrobust(Y, -X, masspoints	= "off")
    }
    
    
    # Create LaTeX table to print
    if(i == 1){
      cat(
        "\\begin{table}[]",
        "\\begin{tabular}{@{}lccccc@{}}",
        "\\midrule\\midrule\\\\",
        " Variable                          & MSE-Optimal &\\multicolumn{1}{c}{RD} & \\multicolumn{2}{c}{\\underline{Robust Inference}} & Eff. Number  \\\\",
        "                                   &  Bandwitdh  &  Estimator    &       p-value          & Conf.Int     & Observations \\\\\\midrule",
        paste(varnames[i], " & ", round(as.numeric(rd_output$bws[1,1]),3) ,  " & ", round(as.numeric(rd_output$coef[1]),3), " & ", round(as.numeric(rd_output$pv[3]),3), " & ", "[",round(as.numeric(rd_output$ci[3,1]),3),", ",round(as.numeric(rd_output$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_output$N_h[1]),", ",as.numeric(rd_output$N_h[2]),") ", " \\\\" ,sep = ""),
        sep = "\n", file=paste0(data.dir,"/rd-estim.txt"), append=TRUE
      )
    } else if(i >1 & i < length(varlist)){
      
      cat(
        paste(varnames[i] , " & ", round(as.numeric(rd_output$bws[1,1]),3) ,  " & ", round(as.numeric(rd_output$coef[1]),3), " & ", round(as.numeric(rd_output$pv[3]),3), " & ", "[",round(as.numeric(rd_output$ci[3,1]),3),", ",round(as.numeric(rd_output$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_output$N_h[1]),", ",as.numeric(rd_output$N_h[2]),") ", " \\\\" ,sep = ""),
        sep = "\n", file=paste0(data.dir,"/rd-estim.txt"), append=TRUE
      )  
      
    } else{
      cat(
        paste(varnames[i], " & ", round(as.numeric(rd_output$bws[1,1]),3) ,  " & ", round(as.numeric(rd_output$coef[1]),3), " & ", round(as.numeric(rd_output$pv[3]),3), " & ", "[",round(as.numeric(rd_output$ci[3,1]),3),", ",round(as.numeric(rd_output$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_output$N_h[1]),", ",as.numeric(rd_output$N_h[2]),") ", " \\\\\\bottomrule" ,sep = ""),
        "\\end{tabular}",
        "\\end{table}",
        sep = "\n", file=paste0(data.dir,"/rd-estim.txt"), append=TRUE
      )  
      
    }
    
  } # specifications 
} # variables



### -- Robustness checks -----------------------------
# Need continuity argument.
# Distribution of the running variable, should not observe some treatment manipulation

# TEMP: Predetermined covariates

# vars <- c("hsgrade_pct", "totcredits_year1", "age_at_entry",
#           "male", "bpl_north_america", "english" )
# 
# for(i in 1:length(vars)){
# 
#   message(paste0("Running robustness check for ", vars[i],"...\n"))
#   
#   Y = df[[vars[i]]] 
#   
#   lm_robcheck <- lm(Y[X >= -0.6 & X <= 0.6] ~ X[X >= -0.6 & X <= 0.6]+ T[X >= -0.6 & X <= 0.6]+
#                       T_X[X >= -0.6 & X <= 0.6])
#   
#   robcheck <- data.frame("First year GPA < cutoff" = c(as.numeric(lm_robcheck$coefficients[3]),
#                                                        as.numeric(sqrt(diag(vcov(lm_robcheck)))[3]),
#                                                        as.numeric(summary(lm_robcheck)$coefficients[,4][3])),
#                          "Constant (control mean)" = c(as.numeric(lm_robcheck$coefficients[2]),
#                                                        as.numeric(sqrt(diag(vcov(lm_robcheck))[2])),
#                                                        as.numeric(summary(lm_robcheck)$coefficients[,4][2]))
#   )
#   
#   row.names(robcheck) <- c("Estimate", "SE", "p-value")
#   
#   print(robcheck)
#   
# 
# }

# 1) Running variable check ------------------------
df %>%
  count(dist_from_cut = round(dist_from_cut, 1)) %>%
  filter(dist_from_cut <= 1.5) %>%
  ggplot(aes(x = dist_from_cut, y = n )) +
  geom_point(shape = 21, colour = "#494949", fill = "white", size = 2, stroke = 1) +
  geom_vline(xintercept = 0, color= "darkgrey", linetype = "dotted") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(-1.5,1.5,0.5)) +
  labs(x = " First year GPA minus probation cutoff",
       y = "Frequency count")

df %>%
  mutate(dist_from_cut = round(dist_from_cut, 1)) %>%
  group_by(dist_from_cut) %>%
  summarise(mean_left_school = mean(left_school)) %>%
  #filter(dist_from_cut <= 1.5) %>%
  ggplot(aes(x = dist_from_cut, y = mean_left_school)) +
  geom_point(shape = 21, colour = "#494949", fill = "white", size = 2, stroke = 1) +
  geom_vline(xintercept = 0, color= "darkgrey", linetype = "dotted") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(-1.5,1.5,0.5)) +
  labs(x = " First year GPA minus probation cutoff",
       y = "Mean Let School")





# 2) Predetermined/placebo tests --------------------------

vars <- c("hsgrade_pct", "totcredits_year1", "age_at_entry",
          "male", "bpl_north_america", "english" )

varnames <- c("HS grade percentile ranking",
              "Credits attempted in 1st Year",
              "Age at Entry",
              "Male",
              "Born in North America",
              "English is 1st language")

# 2.1) Replication as in Lindo (2010): predetermined variables

if (file.exists(paste0(data.dir,"/rob-check-lindo10.txt"))){
  file.remove(paste0(data.dir,"/rob-check-lindo10.txt"))
}

for(i in 1:length(vars)){
  
  message(paste0("Running robustness check for ", vars[i],"...\n"))
  
  Y = df[[vars[i]]] 
  
  rd_robcheck <- rdrobust(Y, -X, p = 1, h = 0.6, kernel = "uniform", masspoints	= "off")
  #rd_robcheck$bws[1,1]
  
  # Create LaTeX table to print
  if(i == 1){
  cat(
    "\\begin{table}[]",
    "\\begin{tabular}{@{}lccccc@{}}",
    "\\midrule\\midrule\\\\",
    " Variable                          & MSE-Optimal &\\multicolumn{1}{c}{RD} & \\multicolumn{2}{c}{\\underline{Robust Inference}} & Eff. Number  \\\\",
    "                                   &  Bandwitdh  &  Estimator    &       p-value          & Conf.Int     & Observations \\\\\\midrule",
    paste(varnames[i], " & ", round(as.numeric(rd_robcheck$bws[1,1]),3) ,  " & ", round(as.numeric(rd_robcheck$coef[1]),3), " & ", round(as.numeric(rd_robcheck$pv[3]),3), " & ", "[",round(as.numeric(rd_robcheck$ci[3,1]),3),", ",round(as.numeric(rd_robcheck$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_robcheck$N_h[1]),", ",as.numeric(rd_robcheck$N_h[2]),") ", " \\\\" ,sep = ""),
    sep = "\n", file=paste0(data.dir,"/rob-check-lindo10.txt"), append=TRUE
    )
  } else if(i >1 & i < length(vars)){
    
  cat(
      paste(varnames[i] , " & ", round(as.numeric(rd_robcheck$bws[1,1]),3) ,  " & ", round(as.numeric(rd_robcheck$coef[1]),3), " & ", round(as.numeric(rd_robcheck$pv[3]),3), " & ", "[",round(as.numeric(rd_robcheck$ci[3,1]),3),", ",round(as.numeric(rd_robcheck$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_robcheck$N_h[1]),", ",as.numeric(rd_robcheck$N_h[2]),") ", " \\\\" ,sep = ""),
      sep = "\n", file=paste0(data.dir,"/rob-check-lindo10.txt"), append=TRUE
      )  
  
  } else{
  cat(
      paste(varnames[i], " & ", round(as.numeric(rd_robcheck$bws[1,1]),3) ,  " & ", round(as.numeric(rd_robcheck$coef[1]),3), " & ", round(as.numeric(rd_robcheck$pv[3]),3), " & ", "[",round(as.numeric(rd_robcheck$ci[3,1]),3),", ",round(as.numeric(rd_robcheck$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_robcheck$N_h[1]),", ",as.numeric(rd_robcheck$N_h[2]),") ", " \\\\\\bottomrule" ,sep = ""),
      "\\end{tabular}",
      "\\end{table}",
      sep = "\n", file=paste0(data.dir,"/rob-check-lindo10.txt"), append=TRUE
      )  
    
  }
  
}

# 2.2) Robust RDD: predetermined variables

if (file.exists(paste0(data.dir,"/rob-check.txt"))){
  file.remove(paste0(data.dir,"/rob-check.txt"))
}

for(i in 1:length(vars)){
  
  message(paste0("Running robustness check for ", vars[i],"...\n"))
  
  Y = df[[vars[i]]] 
  
  rd_robcheck <- rdrobust(Y, -X, masspoints	= "off")
  rd_robcheck$bws[1,1]
  
  # Create LaTeX table to print
  if(i == 1){
    cat(
      "\\begin{table}[]",
      "\\begin{tabular}{@{}lccccc@{}}",
      "\\midrule\\midrule\\\\",
      " Variable                          & MSE-Optimal &\\multicolumn{1}{c}{RD} & \\multicolumn{2}{c}{\\underline{Robust Inference}} & Eff. Number  \\\\",
      "                                   &  Bandwitdh  &  Estimator    &       p-value          & Conf.Int     & Observations \\\\\\midrule",
      paste(varnames[i], " & ", round(as.numeric(rd_robcheck$bws[1,1]),3) ,  " & ", round(as.numeric(rd_robcheck$coef[1]),3), " & ", round(as.numeric(rd_robcheck$pv[3]),3), " & ", "[",round(as.numeric(rd_robcheck$ci[3,1]),3),", ",round(as.numeric(rd_robcheck$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_robcheck$N_h[1]),", ",as.numeric(rd_robcheck$N_h[2]),") ", " \\\\" ,sep = ""),
      sep = "\n", file=paste0(data.dir,"/rob-check.txt"), append=TRUE
    )
  } else if(i >1 & i < length(vars)){
    
    cat(
      paste(varnames[i] , " & ", round(as.numeric(rd_robcheck$bws[1,1]),3) ,  " & ", round(as.numeric(rd_robcheck$coef[1]),3), " & ", round(as.numeric(rd_robcheck$pv[3]),3), " & ", "[",round(as.numeric(rd_robcheck$ci[3,1]),3),", ",round(as.numeric(rd_robcheck$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_robcheck$N_h[1]),", ",as.numeric(rd_robcheck$N_h[2]),") ", " \\\\" ,sep = ""),
      sep = "\n", file=paste0(data.dir,"/rob-check.txt"), append=TRUE
    )  
    
  } else{
    cat(
      paste(varnames[i], " & ", round(as.numeric(rd_robcheck$bws[1,1]),3) ,  " & ", round(as.numeric(rd_robcheck$coef[1]),3), " & ", round(as.numeric(rd_robcheck$pv[3]),3), " & ", "[",round(as.numeric(rd_robcheck$ci[3,1]),3),", ",round(as.numeric(rd_robcheck$ci[3,2]),3),"] ", " & ", "(",as.numeric(rd_robcheck$N_h[1]),", ",as.numeric(rd_robcheck$N_h[2]),") ", " \\\\\\bottomrule" ,sep = ""),
      "\\end{tabular}",
      "\\end{table}",
      sep = "\n", file=paste0(data.dir,"/rob-check.txt"), append=TRUE
    )  
    
  }
  
}














