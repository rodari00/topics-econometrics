# ECON8852 Rodari - Problem Set 3 ---------------------------------------------

# Load libraries
library(hdm) 
library(DoubleML)
library(mlr3learners)
library(mlr3)
library(data.table)
library(randomForest)

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


# -- Load data --------------------------------------------------------------------
data(pension)
help(pension)
data <- pension 
# Constructing the data (as DoubleMLData)
data['age_demean'] <- data$age - mean(data$age)
data['inc_demean'] <- data$inc - mean(data$inc)
data['educ_demean'] <- data$educ - mean(data$educ)
data['fsize_demean'] <- data$fsize - mean(data$fsize)

# Main OLS specifications
formula_flex = "net_tfa ~ e401 + poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + marr + twoearn + db + pira + hown"
formula_flex_het = "net_tfa ~ e401 + poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) +
                                 e401*age_demean +e401*inc_demean + e401*educ_demean + e401*fsize_demean +
                                 e401*marr + e401*twoearn + e401*db + e401*pira + e401*hown"

# Estimate OLS
ols <- lm(formula_flex, data = data)
summary(ols)$coefficients["e401","Estimate"]
summary(ols)$coefficients["e401", "Pr(>|t|)"]

# Estimate HetOLS
ols.het <- lm(formula_flex_het, data = data)
summary(ols.het)$coefficients["e401","Estimate"]
summary(ols.het)$coefficients["e401", "Pr(>|t|)"]

# Save results
fileConn<-file(paste0(results.dir ,"/ols.txt"))
writeLines(stargazer(ols, ols.het, title="Regression Results",
                     align=TRUE, dep.var.labels=c("Financial Wealth"), no.space=TRUE),
           fileConn)
close(fileConn)



# -- DML -----------------------------------------------------------------------

model_flex = as.data.table(model.frame(formula_flex, data))
x_cols = colnames(model_flex)[-c(1,2)]
data_ml = DoubleMLData$new(model_flex, y_col = "net_tfa", d_cols = "e401", x_cols=x_cols)


p <- dim(model_flex)[2]-2
p

# Double-Machine Learning
set.seed(123)
lgr::get_logger("mlr3")$set_threshold("warn") 
lasso <- lrn("regr.cv_glmnet",nfolds = 5, s = "lambda.min")
lasso_class <- lrn("classif.cv_glmnet", nfolds = 5, s = "lambda.min")
dml_irm = DoubleMLIRM$new(data_ml, ml_g = lasso, 
                          ml_m = lasso_class, 
                          trimming_threshold = 0.01, n_folds=5)
dml_irm$fit(store_predictions=TRUE)
dml_irm$summary()
lasso_irm <- dml_irm$coef
lasso_std_irm <- dml_irm$se


# -- Heterogeneous Treatment Effects -------------------------------------------------------

# 1) Twoearn -------------------------------------------------------------------------------

formula_flex = "net_tfa ~ e401 + poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + marr + db + pira + hown"


# Twoearn == 1
data.twoearn1 <- data[data$twoearn == 1,]


model_flex = as.data.table(model.frame(formula_flex, data.twoearn1))
x_cols = colnames(model_flex)[-c(1,2)]
data_ml = DoubleMLData$new(model_flex, y_col = "net_tfa", d_cols = "e401", x_cols=x_cols)

# Double-Machine Learning
set.seed(123)
lgr::get_logger("mlr3")$set_threshold("warn") 
lasso <- lrn("regr.cv_glmnet",nfolds = 5, s = "lambda.min")
lasso_class <- lrn("classif.cv_glmnet", nfolds = 5, s = "lambda.min")
dml_irm = DoubleMLIRM$new(data_ml, ml_g = lasso, 
                          ml_m = lasso_class, 
                          trimming_threshold = 0.01, n_folds=5)
dml_irm$fit(store_predictions=TRUE)
dml_irm$summary()
lasso_irm_twoearn1 <- dml_irm$coef
lasso_std_irm_twoearn1 <- dml_irm$se


# Twoearn == 0
data.twoearn0 <- data[data$twoearn == 0,]


model_flex = as.data.table(model.frame(formula_flex, data.twoearn0))
x_cols = colnames(model_flex)[-c(1,2)]
data_ml = DoubleMLData$new(model_flex, y_col = "net_tfa", d_cols = "e401", x_cols=x_cols)

# Double-Machine Learning
set.seed(123)
lgr::get_logger("mlr3")$set_threshold("warn") 
lasso <- lrn("regr.cv_glmnet",nfolds = 5, s = "lambda.min")
lasso_class <- lrn("classif.cv_glmnet", nfolds = 5, s = "lambda.min")
dml_irm = DoubleMLIRM$new(data_ml, ml_g = lasso, 
                          ml_m = lasso_class, 
                          trimming_threshold = 0.01, n_folds=5)
dml_irm$fit(store_predictions=TRUE)
dml_irm$summary()
lasso_irm_twoearn0 <- dml_irm$coef
lasso_std_irm_twoearn0 <- dml_irm$se



# 2) Income --------------------------------------------------------------------------------

formula_flex = "net_tfa ~ e401 + poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + twoearn + marr + db + pira + hown"

quantile_income <- quantile(data$inc, probs = seq(0, 1, 0.25), na.rm = FALSE)
quantile_income <- quantile_income[2:5]
tau_inc <- c()
std_inc <- c()
for (q in c(1,2,3,4)){
  
  
  
  if (q == 1){
    data_tmp <- data[data$inc <= quantile_income[q],]
  } else{
    data_tmp <- data[(quantile_income[q-1]<= data$inc) & (data$inc<= quantile_income[q]),]
  }
  
  model_flex = as.data.table(model.frame(formula_flex, data_tmp))
  x_cols = colnames(model_flex)[-c(1,2)]
  data_ml = DoubleMLData$new(model_flex, y_col = "net_tfa", d_cols = "e401", x_cols=x_cols)
  
  # Double-Machine Learning
  set.seed(123)
  dml_irm = DoubleMLIRM$new(data_ml, ml_g = lasso, 
                            ml_m = lasso_class, 
                            trimming_threshold = 0.01, n_folds=5)
  dml_irm$fit(store_predictions=TRUE)
  #dml_irm$summary()
  lasso_irm_inc <- dml_irm$coef
  tau_inc[q] <- lasso_irm_inc
  lasso_std_irm_inc <- dml_irm$se
  std_inc[q] <- lasso_std_irm_inc
  
  
  
  message(paste0("Quantile ", as.character(q), ' gives ATE = ', as.character(lasso_irm_inc), ' with SE = ',lasso_std_irm_inc,'\n'))
  
}

data.frame(tau_inc,std_inc)



