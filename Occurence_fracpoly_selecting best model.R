# Load necessary libraries
library(lme4)
library(haven)
364+111
# Create the data frame
364/475
View(fracpol)
MyNames <- data.frame(FirstName = c("minushalf", "minus2", "minus1", "log", "plushalf", "plus1", "plus2", "plus3"))
fracpol <- read_sas("C:/Users/u0118563/OneDrive - KU Leuven/Projecten/Rare diseases/Data/fracpol2.sas7bdat", 
                   NULL)
names(fracpol)
length(unique(fracpol$id))
#recode the responses
recode=function(x){
  ifelse(x=='Yes',1,0)
}
fracpol$gen_01=sapply(as.vector(fracpol$gen_bin),FUN=recode)
fracpol$uni_01=sapply(as.vector(fracpol$uni_bin),FUN=recode)
fracpol$foc_01=sapply(as.vector(fracpol$foc_bin),FUN=recode)

# Fractional polynomial with one polynomial or two different polynomials 
doublepow <- function(response, poww1, poww2, optimizer = "nloptwrap",tolPwrss=1e-3) {
  # Assuming fracpol is a pre-existing data frame
  pol_formula <- as.formula(paste('as.factor(',response, ")~", "+ baseline_age + Male + mut_Other + mutPCDH19 +"
                              ,poww1,"+",
                              paste(poww1, "*", c("baseline_age", "Male"#,"mut_Other","mutPCDH19"
                                                  ), collapse = "+"), "+",
                              poww2, "+", paste(poww2, "*", c("baseline_age", "Male"#,"mut_Other","mutPCDH19"
                                                              ), collapse = "+"),'+(1|id)'))
   # Control settings
  initial_control_settings <- glmerControl(optimizer = optimizer, optCtrl = list(maxfun = 100000), tolPwrss = tolPwrss)
  final_control_settings <- glmerControl(optimizer = optimizer, optCtrl = list(maxfun = 1000000))
  
  # Try fitting the initial model with less quadrature points
  initial_model <- tryCatch({
    glmer(pol_formula, data = fracpol, family = binomial(link = "probit"), nAGQ = 0, control = initial_control_settings)
  }, warning = function(w) {
    warning("Initial model convergence warning: ", conditionMessage(w))
    NULL
  }, error = function(e) {
    message("Initial model did not converge: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(initial_model)) {
     return(data.frame(Power1 = poww1, Power2 = poww2, resp = response, logLik = 9999))
  }
  
  # Extract the fixed effect coefficients and random effect variances
   start_vals <- list(fixef = fixef(initial_model), theta = getME(initial_model, "theta"))
  
  # Try fitting the final model with stricter convergence criterion and using initial parameter estimates
  final_model <- tryCatch({
     glmer(pol_formula, data = fracpol, family = binomial(link = "probit"), nAGQ = 5, start = start_vals, control = final_control_settings)
  }, warning = function(w) {
    warning("Final model convergence warning: ", conditionMessage(w))
    NULL
  }, error = function(e) {
    message("Final model did not converge: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(final_model)) {
    return(data.frame(Power1 = poww1, Power2 = poww2, resp = response, logLik = 9999))
  }
  
  #extract the likelihood
  logLik_val <- logLik(final_model)
  lrt <- data.frame(Power1 = poww1, Power2 = poww2, resp = response, logLik = round(-2*as.numeric(logLik_val),2))
  
  return(lrt)}

# Fractional polynomial with the same polynomial (log is added)
samepow <- function(response, poww1, optimizer = "nloptwrap",tolPwrss=1e-3) {
  # Assuming fracpol is a pre-existing data frame
  pol_formula <- as.formula(paste('as.factor(',response, ")~", poww1, "+ baseline_age + Male + mut_Other + mutPCDH19 +",
                              paste(poww1, "*", c("baseline_age", "Male"#,"mut_Other","mutPCDH19"
                                                  ), collapse = "+"), "+ log +",
                              paste(poww1, "*log*", c("baseline_age", "Male"#,"mut_Other","mutPCDH19"
                                                      , collapse = "+"),'+(1|id)')))
  
  # Control settings
  initial_control_settings <- glmerControl(optimizer = optimizer, optCtrl = list(maxfun = 100000), tolPwrss = tolPwrss)
  final_control_settings <- glmerControl(optimizer = optimizer, optCtrl = list(maxfun = 1000000))
  
  # Try fitting the initial model with less quadrature points
  initial_model <- tryCatch({
    glmer(pol_formula, data = fracpol, family = binomial(link = "probit"), nAGQ = 0, control = initial_control_settings)
  }, warning = function(w) {
    warning("Initial model convergence warning: ", conditionMessage(w))
    NULL
  }, error = function(e) {
    message("Initial model did not converge: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(initial_model)) {
    return(data.frame(Power1 = poww1, Power2 = poww1, resp = response, logLik = 9999))
  }
  
  # Extract the fixed effect coefficients and random effect variances
  start_vals <- list(fixef = fixef(initial_model), theta = getME(initial_model, "theta"))
  
  
  # Try fitting the final model with stricter convergence criterion and using initial parameter estimates
  final_model <- tryCatch({
    glmer(pol_formula, data = fracpol, family = binomial(link = "probit"), nAGQ = 5, start = start_vals, control = final_control_settings)
  }, warning = function(w) {
    warning("Final model convergence warning: ", conditionMessage(w))
    NULL
  }, error = function(e) {
    message("Final model did not converge: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(final_model)) {
    return(data.frame(Power1 = poww1, Power2 = poww1, resp = response, logLik = 9999))
  }
  
  #extract the likelihood
  logLik_val <- logLik(final_model)
  lrt <- data.frame(Power1 = poww1, Power2 = poww1, resp = response, logLik = round(-2*as.numeric(logLik_val),2))
  
  return(lrt)
}


######################################################
#######first response: unilateral seizures############
######################################################
# Initialize data frames to store results: uni_bin
lrt2 <- data.frame(Power1 = character(), Power2 = character(), resp = character(), logLik = numeric())
lrt3 <- data.frame(Power1 = character(), Power2 = character(), resp = character(), logLik = numeric())

# Loop through the names and call the doublepow function: uni_bin
for (i in 1:7) {
  for(j in i:8){
 lrt2 <- rbind(lrt2, doublepow(response = "uni_01", poww1 = MyNames$FirstName[i], poww2 = MyNames$FirstName[j],optimizer='bobyqa'))
 print(c(i,j))
  }
}

# Update lrt3 with lrt2
lrt3 <- lrt2
lrt3$Power2 <- ifelse(lrt3$Power1 == lrt3$Power2, "/", lrt3$Power2)

# Loop through the names and call the samepow function: uni_bin
for (name in MyNames$FirstName) {
  lrt3 <- rbind(lrt3, samepow(response = "uni_01", poww1 = name,optimizer='bobyqa'))
  print(name)
}
results_uni<-lrt3

results_uni[order(results_uni$logLik),]

######################################################
#######second response: generalized seizures##########
######################################################
# Initialize data frames to store results: gen_01
lrt2 <- data.frame(Power1 = character(), Power2 = character(), resp = character(), logLik = numeric())
lrt3 <- data.frame(Power1 = character(), Power2 = character(), resp = character(), logLik = numeric())


# Loop through the names and call the doublepow function: gen_01
for (i in 1:7) {
  for(j in i:8){
    lrt2 <- rbind(lrt2, doublepow(response = "gen_01", poww1 = MyNames$FirstName[i], poww2 = MyNames$FirstName[j],optimizer='bobyqa'))
    print(c(i,j))
  }
}

# Update lrt3 with lrt2
lrt3 <- lrt2
lrt3$Power2 <- ifelse(lrt3$Power1 == lrt3$Power2, "/", lrt3$Power2)


# Loop through the names and call the samepow function
for (name in MyNames$FirstName) {
  lrt3 <- rbind(lrt3, samepow(response = "gen_01", poww1 = name,optimizer='bobyqa'))
  print(name)
}

results_gen=lrt3

###############################################
#######third response: focal seizures##########
###############################################
# Initialize data frames to store results: foc_01
lrt2 <- data.frame(Power1 = character(), Power2 = character(), resp = character(), logLik = numeric())
lrt3 <- data.frame(Power1 = character(), Power2 = character(), resp = character(), logLik = numeric())


# Loop through the names and call the doublepow function: foc_01
for (i in 1:7) {
  for(j in i:8){
    lrt2 <- rbind(lrt2, doublepow(response = "foc_01", poww1 = MyNames$FirstName[i], poww2 = MyNames$FirstName[j],optimizer='bobyqa'))
    print(c(i,j))
  }
}

# Update lrt3 with lrt2
lrt3 <- lrt2
lrt3$Power2 <- ifelse(lrt3$Power1 == lrt3$Power2, "/", lrt3$Power2)

# Loop through the names and call the samepow function for foc_01
for (name in MyNames$FirstName) {
  lrt3 <- rbind(lrt3, samepow(response = "foc_01", poww1 = name,optimizer='bobyqa'))
  print(name)
}
results_foc=lrt2
results_foc[order(results_foc$logLik),]
