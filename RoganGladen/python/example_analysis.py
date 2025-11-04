######################################################################################################################
# The Rogan-Gladen estimator for outcome misclassification
#
# Paul Zivich (2025/11/03)
######################################################################################################################

##########################################
# Loading dependencies
import numpy as np 
import pandas as pd
from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_regression, ee_rogan_gladen
from delicatessen.utilities import inverse_logit
from formulaic import model_matrix


##########################################
# Data Formatting

# Loading EACBIHS data
d1 = pd.read_csv("data/rg_dat.csv")                    # Loading csv data file
d1['V'] = 0                                       # Indicating not the validation data

# Loading Chetty validation data from the Table
d0 = pd.DataFrame()                               # Empty data object
d0['ystar'] = [1]*(208+1) + [0]*(13+31)           # Copy columns from table
d0['y_vld'] = [1]*208 + [0]*1 + [1]*13 + [0]*31   # Dividing into row categories
d0['V'] = 1                                       # Indicating the validation data

# Stacking interest and validation data sets together
d = pd.concat([d1, d0], ignore_index=True)

##########################################
# Uncorrected Prevalence

# Setup
dcc = d.dropna(subset='ystar')                             # Subsetting to non-missing outcome observations
y_star = np.asarray(dcc.loc[dcc['V'] == 0, 'ystar'])       # Mismeasured outcomes in main study data

# By-hand
print("Uncorrected -- By-Hand")
print("Prevalence:", np.round(np.mean(y_star), 3))


# M-estimator
def ef_uncorrected(theta):
    return y_star - theta[0]   # Simple estimating function for the mean


estr = MEstimator(ef_uncorrected, init=[0.5, ])   # Specifying the M-estimator
estr.estimate()                                   # Fitting the M-estimator

print("Uncorrected -- M-estimator")
print("Prevalence:", np.round(estr.theta[0], 3))
print("Stand Err: ", np.round(estr.variance[0, 0]**0.5, 3))
print("95% CI:    ", np.round(estr.confidence_intervals()[0], 3))
print("")

##########################################
# Rogan-Gladen Corrected Prevalence

# Setup
dcc = d.dropna(subset='ystar')                     # Subset to non-missing main study data
y_star = np.asarray(dcc['ystar'])                  # Mismeasured outcome
y_no_nan = np.asarray(dcc['y_vld'].fillna(-99))    # True outcome (only in validation data)
v = np.asarray(1 - dcc['V'])                       # Indicator for main study data

# By-hand
rho_star = np.mean(y_star[v == 1])                        # Prevalence of mismeasured outcome
sens = np.mean(y_star[(v == 0) & (y_no_nan == 1)])        # Sensitivity from validation data
spec = 1 - np.mean(y_star[(v == 0) & (y_no_nan == 0)])    # Specificity from validation data
rho = (rho_star + spec - 1) / (sens + spec - 1)           # Rogan-Gladen estimator
n = dcc.loc[v == 1].shape[0]                              # Observations in main study
m = dcc.loc[(v == 0) & (y_no_nan == 1)].shape[0]          # Observations for sensitivity calc
r = dcc.loc[(v == 0) & (y_no_nan == 0)].shape[0]          # Observations for specificity calc
var = ((rho_star * (1 - rho_star)) / (n * (sens + spec - 1)**2)
       + (sens*(1 - sens)*rho**2)/(m * ((sens + spec - 1)**2))
       + (spec * (1 - spec) * (1 - rho)**2)/(r * ((sens + spec - 1)**2)))
print("RG-Corrected -- By-Hand")
print("Sensitivity:", np.round(sens, 3))
print("Specificity:", np.round(spec, 3))
print("Prevalence:", np.round(rho, 3))
print("Stand Err: ", np.round(var**0.5, 3))


# M-estimator
def ef_rg(theta):
    return ee_rogan_gladen(theta=theta,            # Built-in Rogan-Gladen estimating equation
                           y=y_no_nan,             # ... with validated outcomes
                           y_star=y_star,          # ... mismeasured data
                           r=v)                    # ... and indicator for main study


estr = MEstimator(ef_rg, init=[0.5, 0.5, 0.75, 0.75])  # Specify M-estimator
estr.estimate()                                        # Fit M-estimator

print("RG-Corrected -- M-estimator")
print("Sensitivity:", np.round(estr.theta[2], 3))
print("95% CI:    ", np.round(estr.confidence_intervals()[2, :], 3))
print("Specificity:", np.round(estr.theta[3], 3))
print("95% CI:    ", np.round(estr.confidence_intervals()[3, :], 3))
print("Prevalence:", np.round(estr.theta[0], 3))
print("Stand Err: ", np.round(estr.variance[0, 0]**0.5, 3))
print("95% CI:    ", np.round(estr.confidence_intervals()[0, :], 3))
print("")

##########################################
# Corrected Prevalence with IPMW

# Setup
d_no_nan = d.fillna(-99)                                                   # Filling missing data with placeholders
y_star = np.asarray(d_no_nan['ystar'])                                     # Mismeasured outcome
y_no_nan = np.asarray(d_no_nan['y_vld'])                                   # Validated outcomes from validation data
v = np.asarray(1 - d_no_nan['V'])                                          # Indicator if in main study data
r = np.asarray(np.where(y_star != -99, 1, 0))                              # Missing indicator for mismeasured outcome
W = model_matrix("age_g30 + C(edu) + drinking + selfreportsw", d_no_nan)   # Missingness score model covariates


def ef_rg_ipmw(theta):
    # Subset parameters for estimating functions
    alpha = theta[:4]
    eta = theta[4:]
    # Inverse Probability of Missingness weight calculations
    ee_log = ee_regression(theta=eta, X=W, y=r, model='logistic') * v     # Fitting with only main study data
    ipmw = (v*r) / inverse_logit(np.dot(W, eta)) + (1-v)                  # Assign IPMW, where validation set to 1
    # Weighted Rogan-Gladen estimator
    ee_rg = ee_rogan_gladen(theta=alpha,                                  # Built-in Rogan-Gladen estimating equation
                            y=y_no_nan,                                   # ... with validated outcomes
                            y_star=y_star,                                # ... mismeasured data
                            r=v,                                          # ... indicator for main study
                            weights=ipmw)                                 # ... and estimated weights
    # Return the stacked estimating functions
    return np.vstack([ee_rg, ee_log])


init_vals = [0.5, 0.5, 0.75, 0.75] + [0., ]*W.shape[1]    # Starting values for root-finding procedure
estr = MEstimator(ef_rg_ipmw, init=init_vals)             # Specify M-estimator
estr.estimate()                                           # Fit M-estimator

print("RG + IPMW")
print("Sensitivity:", np.round(estr.theta[2], 3))
print("95% CI:    ", np.round(estr.confidence_intervals()[2, :], 3))
print("Specificity:", np.round(estr.theta[3], 3))
print("95% CI:    ", np.round(estr.confidence_intervals()[3, :], 3))
print("Prevalence:", np.round(estr.theta[0], 3))
print("Stand Err: ", np.round(estr.variance[0, 0]**0.5, 3))
print("95% CI:    ", np.round(estr.confidence_intervals()[0, :], 3))

# ------------------------------------------------------------------------------------------------------------
# OUTPUT (2025/11/03)
# ------------------------------------------------------------------------------------------------------------
#
# Uncorrected -- By-Hand
# Prevalence: 0.053
# Uncorrected -- M-estimator
# Prevalence: 0.053
# Stand Err:  0.015
# 95% CI:     [0.024 0.083]
#
# RG-Corrected -- By-Hand
# Sensitivity: 0.941
# Specificity: 0.969
# Prevalence: 0.024
# Stand Err:  0.037
# RG-Corrected -- M-estimator
# Sensitivity: 0.941
# 95% CI:     [0.91  0.972]
# Specificity: 0.969
# 95% CI:     [0.908 1.029]
# Prevalence: 0.024
# Stand Err:  0.037
# 95% CI:     [-0.048  0.097]
#
# RG + IPMW
# Sensitivity: 0.941
# 95% CI:     [0.91  0.972]
# Specificity: 0.969
# 95% CI:     [0.908 1.029]
# Prevalence: 0.025
# Stand Err:  0.037
# 95% CI:     [-0.048  0.097]
