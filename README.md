# MoE-SMN-CR
Robust Bayesian inference for the censored mixture of experts model using heavy-tailed distributions

By: Elham Mirfarah · Mehrdad Naderi · Tsung-I Lin · Wan-Lun Wang

contains 
Gibbs.MoE.Cr: Gibbs sampling for the family of scale mixture of normal.

dLaplace: pdf of laplace distribution;
dT.c: Likelihood function of censored t student;
MHnu2 and MHnu: Metropolis- Hasting for sampling for scale parametr of t student;
dSlash: pdf of Slash distribution;
dVG: pdf of Variance-Gamma distribution; 
dVG.c: Likelihood function of censored Variance-Gamma;
PVG.c and PVG.c.prim: cdf and survival function of Variance-Gamma distribution;
MHnu2.VG and MHnu.Vg: Metropolis- Hasting for sampling for scale parametr of t Variance-Gamma distribution;
dCN: pdf of Contaminated Normal distribution;
PTIN and PTIN.prim: cdf and survival function of Tail inflated normal distribution(TIN);
dTIN: pdf of TIN distrbition;
d.TIN.c: Likelihood function of censored TIN;
MHnuTIN: Metropolis- Hasting for sampling for scale parametr of TIN;
PI.Exp: softmax (multinomial logistic) model.

Notes:
The "UPG", "truncnorm", "zipfR", "GIGrvg" R packages are required. 

