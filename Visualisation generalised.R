library(lme4)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(haven)
fracpol <- read_sas("C:/Users/u0118563/OneDrive - KU Leuven/Projecten/Rare diseases/Data/fracpol2.sas7bdat", 
                    NULL)
#recode the responses
recode=function(x){
  ifelse(x=='Yes',1,0)
}
fracpol$gen_01=sapply(as.vector(fracpol$gen_bin),FUN=recode)
fracpol$uni_01=sapply(as.vector(fracpol$uni_bin),FUN=recode)
fracpol$foc_01=sapply(as.vector(fracpol$foc_bin),FUN=recode)

#Sex variable
recode=function(x){
  ifelse(x=='F','Female','Male')
}
fracpol$Sex=sapply(as.vector(fracpol$sex),FUN=recode)

# Define the model formula
fracpol$log_plushalf=fracpol$log*fracpol$plushalf
response <- "gen_01"
poww1 <- "plushalf"
poww2 <- "log_plushalf"
formula <- as.formula(paste('as.factor(', response, ")~", "+ baseline_age + Male + mut_Other + mutPCDH19 +", poww1, "+",
                            paste(poww1, "*", c("baseline_age", "Male"), collapse = "+"), "+",
                            poww2, "+", paste(poww2, "*", c("baseline_age", "Male"), collapse = "+"), '+(1|id)'))

# Fit the initial model with less strict convergence criterion
initial_control_settings <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000), tolPwrss = 1e-3)
initial_model <- glmer(formula, data = fracpol, family = binomial(link = "probit"), nAGQ = 1, control = initial_control_settings)

# Extract the fixed effect coefficients and random effect variances
#initial_model=final_model
start_vals <- list(fixef = fixef(initial_model), theta = getME(initial_model, "theta"))

# Fit the final model with stricter convergence criterion and using initial parameter estimates
final_control_settings <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
final_model <- glmer(formula, data = fracpol, family = binomial(link = "probit"), nAGQ = 5, start = start_vals, control = final_control_settings)

fracpol$age=fracpol$baseline_age+fracpol$time
fracpol<-fracpol[!(is.na(fracpol$mutPCDH19)),]


# Function to create visualizations
visualize_results <- function(model, data=fracpol, response, age_var = "age", gender_var = "Sex", id_var = "id") {
  
  # Create a new data frame with the original data and predicted values
  data_with_preds <- data %>%
    mutate(predicted = predict(model, type = "response"),
           fitted = predict(model, re.form = NA, type = "response"))
  
  # Spaghetti plot of predicted values
  spaghetti_plot <- ggplot(data_with_preds, aes_string(x = age_var, y = "predicted", group = id_var)) +
    geom_line(alpha = 0.3) +ylim(0,1)+
    labs(title = "",
         x = "Age",
         y = "Predicted Probability") +
    theme_minimal()
  
  # Loess plot of fitted values by gender
  loess_plot <- ggplot(data_with_preds, aes_string(x = age_var, y = "fitted", color = gender_var)) +
    geom_smooth(method = "loess", se = TRUE) +ylim(0,1)+
    labs(title = "",
         x = "Age",
         y = "Predicted Probability") +
    theme_minimal()
  
  return(list(spaghetti_plot = spaghetti_plot, loess_plot = loess_plot))
}

# Visualize the results
plots <- visualize_results(final_model, fracpol, response='gen_01', age_var='age', gender_var='Sex', id_var="id")

# Display the plots
ggarrange(plots$spaghetti_plot,plots$loess_plot,ncol=2)
