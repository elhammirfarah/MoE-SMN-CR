PI.Exp = function(r, alpha){
  alpha = as.matrix(alpha)
  G = ncol(alpha) + 1
  alpha0 = exp(cbind(r %*% alpha, 0))
  Out = alpha0 / rowSums(alpha0)
  return(Out)
}

dLaplace = function(x, mu = 0, sigma = 1, log = T){
  PDF = -abs(x - mu)/sigma - log(2 * sigma)
  ifelse(log == T, return(PDF), return(exp(PDF)))
}

dT.c = function(cc,y,mu,sigma2 = 1,nu = 4,log = F,Cens.type = c("Left", "Right")){
  PDF = vector(mode = "numeric", length = length(y))
  y = (y - mu) / sqrt(sigma2)
  PDF = dt(y, nu, log = T) - 0.5 * log(sigma2)
  if(Cens.type == "Left") PDF[cc == 1] = pt(y[cc == 1], nu, log.p = T) 
  if(Cens.type == "Right") PDF[cc == 1] = pt(y[cc == 1], nu, lower.tail = F, log.p = T)
  ifelse(log == T, return(PDF), return(exp(PDF)))
}

MHnu2 = function(last, U, nuprior, hypernu0){
  n <- length(U)
  
  if (nuprior == "JandS"){
    d = 9 / (1 + sqrt(2))
    gSylvia = function(nu, U, d) {
      n = length(U)
      ff = log(nu - 1) - 3 * log(nu - 1 + d) +  
        0.5 * n * nu * log(nu/2) - (0.5 * nu) * sum(U -log(U)) - n * lgamma(nu/2)
      ifelse((nu > 1), return(ff), return(0))
    }
    Fonseca1 = deriv(~log(nu - 1) - 3 * log(nu - 1 + d) + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2), c("nu"), function(nu) { }, hessian = TRUE)
    aux1 = Fonseca1(last)
    aux2 = -0.5 * sum(U-log(U))
    q1 = attr(aux1, "gradient")[1]  + aux2
    q2 = attr(aux1, "hessian")[1] 
    bw = max(0.0001, -1/q2) # # -q2 #  
    aw = last + q1 * bw
    
    cand = truncnorm :: rtruncnorm(1, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    aux11 = Fonseca1(cand)
    aux21 = -0.5 * sum(U-log(U))
    q11 = attr(aux11, "gradient")[1]  + aux21
    q21 = attr(aux11, "hessian")[1] 
    bw1 = max(0.0001, -1/q21) #  -q21 #    
    aw1 = cand + q11 * bw1
    
    alfa1 = truncnorm :: dtruncnorm(last, a = 2, b = 40, mean = aw1, sd = sqrt(bw1))
    alfa2 = truncnorm :: dtruncnorm(cand, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    alfa = exp(gSylvia(cand, U, d)-gSylvia(last, U, d)) * alfa1/ alfa2
  }
  
  if(nuprior == "Jeffreys") {
    gJeffreys = function(nu, U) {
      n = length(U)
      ff = log(sqrt(nu/(nu + 3)) * sqrt(trigamma(nu/2) - trigamma((nu + 1)/2) - 
                                          2 * (nu + 3)/((nu) * (nu + 1)^2))) + 
        0.5 * n * nu * log(nu/2) - (0.5 * nu) * sum(U -log(U)) - n * lgamma(nu/2)
      return(ff)
    }
    Fonseca1 = deriv(~log(sqrt(nu/(nu + 3)) * sqrt(trigamma(nu/2) - trigamma((nu + 1)/2) - 2 * (nu + 3)/((nu) * (nu + 1)^2))) + 
                       0.5 * n * nu * log(nu/2) - n * lgamma(nu/2), c("nu"), function(nu) { }, hessian = TRUE)
    aux1 = Fonseca1(last)
    aux2 = -0.5 * sum(U-log(U))
    q1 = attr(aux1, "gradient")[1]  + aux2
    q2 = attr(aux1, "hessian")[1] 
    bw = max(0.0001, -1/q2) # # -q2 #  
    aw = last + q1 * bw
    
    cand = truncnorm :: rtruncnorm(1, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    aux11 = Fonseca1(cand)
    aux21 = -0.5 * sum(U-log(U))
    q11 = attr(aux11, "gradient")[1]  + aux21
    q21 = attr(aux11, "hessian")[1] 
    bw1 = max(0.0001, -1/q21) #  -q21 #    
    aw1 = cand + q11 * bw1
    
    alfa1 = truncnorm :: dtruncnorm(last, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    alfa2 = truncnorm :: dtruncnorm(cand, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    alfa = exp(gJeffreys(cand, U)-gJeffreys(last, U)) * alfa1/ alfa2
  }
  
  if (nuprior == "Pareto"){
    gPareto = function(nu, U, hypernu0) {
      n = length(U)
      ff = (- (hypernu0+1) * log(nu)) +  0.5 * n * nu * log(nu/2) - (0.5 * nu) * sum(U -log(U)) - n * lgamma(nu/2)
      ifelse((nu > 2), return(ff), return(0))
    }
    Fonseca1 = deriv(~ (- (hypernu0+1) * log(nu)) + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2), c("nu"), function(nu) { }, hessian = TRUE)
    aux1 = Fonseca1(last)
    aux2 = -0.5 * sum(U-log(U))
    q1 = attr(aux1, "gradient")[1]  + aux2
    q2 = attr(aux1, "hessian")[1] 
    bw = max(0.0001, -1/q2) # # -q2 #  
    aw = last + q1 * bw
    
    cand = truncnorm :: rtruncnorm(1, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    aux11 = Fonseca1(cand)
    aux21 = -0.5 * sum(U-log(U))
    q11 = attr(aux11, "gradient")[1]  + aux21
    q21 = attr(aux11, "hessian")[1] 
    bw1 = max(0.0001, -1/q21) #  -q21 #    
    aw1 = cand + q11 * bw1
    
    alfa1 = truncnorm :: dtruncnorm(last, a = 2, b = 40, mean = aw1, sd = sqrt(bw1))
    alfa2 = truncnorm :: dtruncnorm(cand, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    alfa = exp(gPareto(cand, U, hypernu0)-gPareto(last, U, hypernu0)) * alfa1/ alfa2
    
  }
  
  if (nuprior == "Gama"){
    hyper1 = hypernu0[1]
    hyper2 = hypernu0[2]
    gGama = function(nu, U, hyper1, hyper2) {
      n = length(U)
      ff = (hyper1-1) * log(nu) - hyper2 * nu + 0.5 * n * nu * log(nu/2) - (0.5 * nu) * sum(U - log(U)) - n * lgamma(nu/2)
      return(ff)
    }
    Fonseca1 = deriv(~(hyper1-1) * log(nu) - hyper2 * nu + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2), c("nu"), function(nu) { }, hessian = TRUE)
    aux1 = Fonseca1(last)
    aux2 = -0.5 * sum(U-log(U))
    q1 = attr(aux1, "gradient")[1]  + aux2
    q2 = attr(aux1, "hessian")[1] 
    bw = max(0.0001, -1/q2) # # -q2 #  
    aw = last + q1 * bw
    
    cand = truncnorm :: rtruncnorm(1, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    aux11 = Fonseca1(cand)
    aux21 = -0.5 * sum(U-log(U))
    q11 = attr(aux11, "gradient")[1]  + aux21
    q21 = attr(aux11, "hessian")[1] 
    bw1 = max(0.0001, -1/q21) #  -q21 #    
    aw1 = cand + q11 * bw1
    
    alfa1 = truncnorm :: dtruncnorm(last, a = 2, b = 40, mean = aw1, sd = sqrt(bw1))
    alfa2 = truncnorm :: dtruncnorm(cand, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    
    alfa = exp(gGama(cand, U, hyper1, hyper2)-gGama(last, U, hyper1, hyper2)) * alfa1/ alfa2
  }
  
  last = ifelse(runif(1) < min(alfa, 1),  cand, last)
  return(last)
  
}



MHnu = function(cc, y, r, mu, sigma2, alpha, nu,
                Cens.type = c("Left", "Right"),
                cnu = cnu, g, nuprior, hypernu0) {
  
  d.T.Exp = function(cc, y, r, mu, sigma2, alpha, nu,
                     Cens.type = c("Left", "Right")){
    g = length(sigma2)
    pi1 = PI.Exp(r = r, alpha = alpha)
    PDF = 0
    for (j in 1:g) PDF = PDF + exp(log(pi1[, j]) +
                                     dT.c(cc, y, mu[, j], sigma2[j], nu[j], log = T,
                                          Cens.type = Cens.type))
    return(PDF)
  }
  
  if (nuprior == "JandS") {
    d = hypernu0   
    f.pr = function(nu, d){
      f1 = (2 * d * (nu - 1)) / (nu - 1 + d)^3
      ifelse((nu > 2)&& (nu < 40) , return(f1), return(0))  
    }
    last = nu[g]
    F.a = plnorm(2, log(last), sqrt(cnu))
    F.b = plnorm(40, log(last), sqrt(cnu))
    vi = runif(1, F.a, F.b)
    cand = qlnorm(vi, log(last), sqrt(cnu))
    nu2 = nu; nu2[g] = cand
    alfa1 = (last) * f.pr(last, d) 
    alfa2 = (cand) * f.pr(cand, d)
    }
  
  if (nuprior == "Gama"){ 
    al = hypernu0[1] 
    be = hypernu0[2] 
    f.pr = function(nu, al, be){
      f1 = ((be^al)/gamma(be)) * (nu^(al-1)) * exp(-nu*be)
      ifelse((nu > 2)  && (nu < 40) , return(f1), return(0)) 
    }
    last = nu[g]
    F.a = plnorm(2, log(last), sqrt(cnu))
    F.b = plnorm(40, log(last), sqrt(cnu))
    vi = runif(1, F.a, F.b)
    cand = qlnorm(vi, log(last), sqrt(cnu))
    nu2 = nu; nu2[g] = cand
    alfa1 = (last) * f.pr(last, al, be) 
    alfa2 = (cand) * f.pr(cand, al, be)
    }
  if (nuprior == "Jeffreys"){
    f.pr = function(nu){
      f1 = sqrt(nu/(nu + 3)) * sqrt(trigamma(nu/2) - trigamma((nu + 1)/2) - 
                                      2 * (nu + 3)/((nu) * (nu + 1)^2))
      ifelse((nu > 2)&& (nu < 40) , return(f1), return(0)) # 
    }
    last = nu[g]
    F.a = plnorm(2, log(last), sqrt(cnu))
    F.b = plnorm(40, log(last), sqrt(cnu))
    vi = runif(1, F.a, F.b)
    cand = qlnorm(vi, log(last), sqrt(cnu))
    nu2 = nu; nu2[g] = cand
    alfa1 = (last) * f.pr(last) 
    alfa2 = (cand) * f.pr(cand)
    }
  if (nuprior == "Pareto"){
    d = hypernu0  
    f.pr = function(nu, d){
      f1 = (d * 2^d)/ (nu^(d+1))
      ifelse((nu > 2)  && (nu < 40), return(f1), return(0)) #  
    }
   last = nu[g]
   F.a = plnorm(2, log(last), sqrt(cnu))
   F.b = plnorm(40, log(last), sqrt(cnu))
   vi = runif(1, F.a, F.b)
   cand = qlnorm(vi, log(last), sqrt(cnu))
   nu2 = nu; nu2[g] = cand
   alfa1 = (last) * f.pr(last, d) 
   alfa2 = (cand) * f.pr(cand, d)
   } 
  alfa0 = exp(sum(log(d.T.Exp(cc, y, r, mu, sigma2, alpha, nu2, Cens.type = Cens.type)))-
                sum(log(d.T.Exp(cc, y, r, mu, sigma2, alpha, nu, Cens.type = Cens.type))))
  alfa = alfa0 * (alfa2 / alfa1)
  out = ifelse(runif(1) < min(alfa, 1),  cand, last)
  return(out)
}

dSlash = function(y, mu, sigma2, nu, log = F)
{
  z = (y - mu) / sqrt(sigma2)
  PDF = log(nu/sqrt(2 * pi * sigma2)) + (nu + 0.5) * log(2/z^2) +
    pgamma(1, nu + 0.5, rate = 0.5 * (z^2), log.p = T) +
    lgamma(nu + 0.5)
  ifelse(log == F, return(exp(PDF)), return(PDF))
}

BesselK = function(x, kappa, log.val = T){
  F1 = log( besselK(x, kappa, expon.scaled = TRUE) ) - x
  F2 = besselK(x, kappa)
  ifelse(log.val == T, return(F1), return(F2))
}

PVG.c = function(y, mu, sigma2, nu)
{
  z = (y - mu)/sqrt(sigma2)
  f1 = function(u) exp(dgamma(u, nu/2, rate = nu/2, log = T) + pnorm(z/sqrt(u), 0, 1, log.p = T)) 
  Out = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                  stop.on.error = FALSE)$value
  Out[Out > 1] = 1
  return(Out)
}

PVG.c.prim = function(y, mu, sigma2, nu)
{
  z = (y - mu)/sqrt(sigma2)
  f1 = function(u) exp(dgamma(u, nu/2, rate = nu/2, log = T) + pnorm(z/sqrt(u), 0, 1, lower.tail = F,  log.p = T)) 
  Out = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                  stop.on.error = FALSE)$value
  Out[Out > 1] = 1
  return(Out)
}

dVG = function(y, mu, sigma2 = 1, nu = 4, log = F){
  PDF = vector(mode = "numeric", length = length(y))
  y = (y - mu) / sqrt(sigma2)
  PDF = ((nu-1)/2) * log(abs(y)) + ((nu+1)/4) * log(nu) + ((1-nu)/2) * log(2) -
    lgamma(nu/2) - 0.5 * log(pi) - 0.5 * log(sigma2) + BesselK(sqrt(nu) * abs(y), (nu-1)/2, log.val = T)
  
  ifelse(log == T, return(PDF), return(exp(PDF)))
}

dVG.c = function(cc, y, mu, sigma2 = 1, nu = 4, log = F,
                 Cens.type = c("Left", "Right")){
  PDF = vector(mode = "numeric", length = length(y))
  y = (y - mu) / sqrt(sigma2)
  PDF = ((nu-1)/2) * log(abs(y)) + ((nu+1)/4) * log(nu) + ((1-nu)/2) * log(2) -
    lgamma(nu/2) - 0.5 * log(pi) - 0.5 * log(sigma2) + BesselK(sqrt(nu) * abs(y), (nu-1)/2, log.val = T)
  IND = which(cc == 1)
  for(j in IND){
    if(Cens.type == "Left") PDF[j] = log(PVG.c(y[j], 0, 1, nu)) 
    if(Cens.type == "Right") PDF[j] = log(PVG.c.prim(y[j], 0, 1, nu))
  }
  ifelse(log == T, return(PDF), return(exp(PDF)))
}


MHnu2.VG = function(last, U, nuprior, hypernu0) {
  n <- length(U)
  
  if (nuprior == "Pareto"){
    gPareto = function(nu, U, hypernu0) {
      n = length(U)
      ff = (- (hypernu0+1) * log(nu)) +  0.5 * n * nu * log(nu/2) - (0.5 * nu) * sum(1/U +log(U)) - n * lgamma(nu/2)
      ifelse((nu > 2), return(ff), return(0))
    }
    Fonseca1 = deriv(~ (- (hypernu0+1) * log(nu)) + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2), c("nu"), function(nu) { }, hessian = TRUE)
    aux1 = Fonseca1(last)
    aux2 = -0.5 * sum(1/U +log(U))
    q1 = attr(aux1, "gradient")[1]  + aux2
    q2 = attr(aux1, "hessian")[1] 
    bw = max(0.0001, -1/q2)  
    aw = last + q1 * bw
    cand = truncnorm :: rtruncnorm(1, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    aux11 = Fonseca1(cand)
    aux21 = -0.5 * sum(1/U +log(U))
    q11 = attr(aux11, "gradient")[1]  + aux21
    q21 = attr(aux11, "hessian")[1] 
    bw1 = max(0.0001, -1/q21) #  -q21 #    
    aw1 = cand + q11 * bw1
    alfa1 = truncnorm :: dtruncnorm(last, a = 2, b = 40, mean = aw1, sd = sqrt(bw1))
    alfa2 = truncnorm :: dtruncnorm(cand, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    alfa = exp(gPareto(cand, U, hypernu0)-gPareto(last, U, hypernu0)) * alfa1/ alfa2
    }
  
  if (nuprior == "Gama"){
    hyper1 = hypernu0[1]
    hyper2 = hypernu0[2]
    gGama = function(nu, U, hyper1, hyper2) {
      n = length(U)
      ff = (hyper1-1) * log(nu) - hyper2 * nu + 0.5 * n * nu * log(nu/2) - (0.5 * nu) * sum(1/U +log(U)) - n * lgamma(nu/2)
      return(ff)
    }
    Fonseca1 = deriv(~(hyper1-1) * log(nu) - hyper2 * nu + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2), c("nu"), function(nu) { }, hessian = TRUE)
    aux1 = Fonseca1(last)
    aux2 = -0.5 * sum(1/U +log(U))
    q1 = attr(aux1, "gradient")[1]  + aux2
    q2 = attr(aux1, "hessian")[1] 
    bw = max(0.0001, -1/q2) 
    aw = last + q1 * bw
    cand = truncnorm :: rtruncnorm(1, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    aux11 = Fonseca1(cand)
    aux21 = -0.5 * sum(1/U +log(U))
    q11 = attr(aux11, "gradient")[1]  + aux21
    q21 = attr(aux11, "hessian")[1] 
    bw1 = max(0.0001, -1/q21)   
    aw1 = cand + q11 * bw1
    alfa1 = truncnorm :: dtruncnorm(last, a = 2, b = 40, mean = aw1, sd = sqrt(bw1))
    alfa2 = truncnorm :: dtruncnorm(cand, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    alfa = exp(gGama(cand, U, hyper1, hyper2)-gGama(last, U, hyper1, hyper2)) * alfa1/ alfa2
  }
  
  if (nuprior == "JandS"){
    d = hypernu0
    gSylvia = function(nu, U, d) {
      n = length(U)
      ff = log(nu - 1) - 3 * log(nu - 1 + d) +  
        0.5 * n * nu * log(nu/2) - (0.5 * nu) * sum(1/U +log(U)) - n * lgamma(nu/2)
      return(ff)
    }
    Fonseca1 = deriv(~log(nu - 1) - 3 * log(nu - 1 + d) + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2), c("nu"), function(nu) { }, hessian = TRUE)
    aux1 = Fonseca1(last)
    aux2 = -0.5 * sum(1/U +log(U))
    q1 = attr(aux1, "gradient")[1]  + aux2
    q2 = attr(aux1, "hessian")[1] 
    bw = max(0.0001, -1/q2)   
    aw = last + q1 * bw
    cand = truncnorm :: rtruncnorm(1, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    aux11 = Fonseca1(cand)
    aux21 = -0.5 * sum(1/U +log(U))
    q11 = attr(aux11, "gradient")[1]  + aux21
    q21 = attr(aux11, "hessian")[1] 
    bw1 = max(0.0001, -1/q21)     
    aw1 = cand + q11 * bw1
    alfa1 = truncnorm :: dtruncnorm(last, a = 2, b = 40, mean = aw1, sd = sqrt(bw1))
    alfa2 = truncnorm :: dtruncnorm(cand, a = 2, b = 40, mean = aw, sd = sqrt(bw))
    alfa = exp(gSylvia(cand, U, d)-gSylvia(last, U, d)) * alfa1/ alfa2
  }
  last = ifelse(runif(1) < min(alfa, 1),  cand, last)
  return(last)
  
}

MHnu.Vg = function(cc, y, r, mu, sigma2, alpha, nu,
                Cens.type = c("Left", "Right"),
                cnu = cnu, g, nuprior, hypernu0) {
  d.VG.Exp = function(cc, y, r, mu, sigma2, alpha, nu,
                      Cens.type = c("Left", "Right")) {
    g = length(sigma2)
    pi1 = PI.Exp(r = r, alpha = alpha)
    PDF = 0
    for (j in 1:g) PDF = PDF + exp(log(pi1[, j]) +
                                     dVG.c(cc, y, mu[, j], sigma2[j], nu[j], log = T,
                                           Cens.type = Cens.type))
    return(PDF)
  }
  if (nuprior == "JandS"){
    d = hypernu0  
    f.pr = function(nu, d){
      f1 = (nu - 1) / (nu - 1 + d)^3
      ifelse((nu > 2) && (nu < 40) , return(f1), return(0)) 
    }
    last = nu[g]
    F.a = plnorm(2, log(last), sqrt(cnu))
    F.b = plnorm(40, log(last), sqrt(cnu))
    vi = runif(1, F.a, F.b)
    cand = qlnorm(vi, log(last), sqrt(cnu))
    nu2 = nu; nu2[g] = cand
    alfa1 = (last) * f.pr(last, d) 
    }
  if (nuprior == "Pareto"){
    d = hypernu0   
    f.pr = function(nu, d){
      f1 = (d * 2^d)/ (nu^(d+1))
      ifelse((nu > 2), return(f1), return(0)) #  
    }
   last = nu[g]
   F.a = plnorm(2, log(last), sqrt(cnu))
   F.b = plnorm(40, log(last), sqrt(cnu))
   vi = runif(1, F.a, F.b)
   cand = qlnorm(vi, log(last), sqrt(cnu))
   nu2 = nu; nu2[g] = cand
   alfa1 = (last) * f.pr(last, d) 
   alfa2 = (cand) * f.pr(cand, d)
   } 
  if (nuprior == "Gama"){ 
    al = hypernu0[1]   
    be = hypernu0[2]   
    f.pr = function(nu, al, be){
      f1 = ((be^al)/gamma(be)) * (nu^(al-1)) * exp(-nu*be)
      ifelse((nu > 2)  && (nu < 40) , return(f1), return(0)) 
    }
   last = nu[g]
   F.a = plnorm(2, log(last), sqrt(cnu))
   F.b = plnorm(40, log(last), sqrt(cnu))
   vi = runif(1, F.a, F.b)
   cand = qlnorm(vi, log(last), sqrt(cnu))
   nu2 = nu; nu2[g] = cand
   alfa1 = (last) * f.pr(last, al, be) 
   alfa2 = (cand) * f.pr(cand, al, be)
   }
  alfa0 = exp(sum(log(d.VG.Exp(cc, y, r, mu, sigma2, alpha, nu2, Cens.type = Cens.type)))-
                sum(log(d.VG.Exp(cc, y, r, mu, sigma2, alpha, nu, Cens.type = Cens.type))))
  alfa = alfa0 * (alfa2 / alfa1)
  out = ifelse(runif(1) < min(alfa, 1),  cand, last)
  return(out)
}

dCN  = function(y, mu, sigma2, nu, gama, log = F){
  PDF = vector(mode = "numeric", length = length(y))
  D1 = dnorm(y, mu, sqrt(sigma2 / gama), log = T)
  D2 = exp(dnorm(y, mu, sqrt(sigma2), log = T) -  D1)
  PDF = D1 + log(nu + (1 - nu) * D2)
  ifelse(log == T, return(PDF), return(exp(PDF)))
}

PTIN = function(y, mu, sigma2, nu)
{
  z = (y - mu)/sqrt(sigma2)
  f1 = function(u) exp(-log(nu)+ pnorm(z * sqrt(u), log.p = T))
  Out = integrate(f1, 1-nu, 1, rel.tol = .Machine$double.eps^0.25, stop.on.error = FALSE)$value
  Out[Out > 1] = 1
  return(Out)
}

PTIN.prim = function(y, mu, sigma2, nu)
{
  z = (y - mu)/sqrt(sigma2)
  f1 = function(u) exp(-log(nu) + pnorm(z * sqrt(u), log.p = T, lower.tail = F))
  Out = integrate(f1, 1-nu, 1, rel.tol = .Machine$double.eps^0.25, stop.on.error = FALSE)$value
  Out[Out > 1] = 1
  return(Out)
}

dTIN = function(y, mu, sigma2, nu, log = F)
{
  z = (y - mu)^2 / (2 * sigma2)
  PDF = - 0.5 * log(2 * pi * sigma2) - log(nu) - 1.5 * log(z) +
    zipfR::Igamma(1.5, (1-nu) * z, lower = FALSE, log = TRUE) + 
    log(1 - exp(zipfR::Igamma(1.5, z, lower = FALSE, log = TRUE) - 
                  zipfR::Igamma(1.5, (1-nu) * z, lower = FALSE, log = TRUE)))
  for (j in which(PDF == -Inf)) PDF[j] = log(.Machine$double.eps)
  ifelse(log == F, return(exp(PDF)), return(PDF))
}

d.TIN.c = function(cc, y, mu, sigma2 = 1, nu = 2, log = F,
                   Cens.type = c("Left", "Right"))
{
  PDF = vector(mode = "numeric", length = length(y))
  delta = (y - mu) / sqrt(sigma2)
  PDF = dTIN(delta, 0, 1, nu, log = T) - 0.5 * log(sigma2)
  IND = which(cc == 1)
  for(j in IND){
    if(Cens.type == "Left") PDF[j] = log(PTIN(delta[j], 0, 1, nu)) 
    if(Cens.type == "Right") PDF[j] = log(PTIN.prim(delta[j], 0, 1, nu)) 
  }
  ifelse(log == T, return(PDF), return(exp(PDF)))
}

MHnuTIN <- function(cc, y, r, mu, sigma2, alpha, nu,
                    Cens.type = c("Left", "Right"),
                    cnu = 0.001, hypernu0, g){
  
  d.TIN.Exp = function(cc, y, r, mu, sigma2, alpha, nu, Cens.type = c("Left", "Right")){
    G = length(sigma2)
    pi1 = PI.Exp(r = r, alpha = alpha)
    PDF = 0
    for (j in 1:G) PDF = PDF + exp(log(pi1[, j]) + d.TIN.c(cc, y, mu[, j], sigma2[j], nu[j], log = T,
                                                           Cens.type = Cens.type))
    return(PDF)
  }
  
  s0 = hypernu0[1]
  s1 = hypernu0[2]
  f.pr = function(nu){
    nu1 = nu/(1 + nu)
    f1 = (nu1)^(s0 - 1) * (1 - nu1)^(s1 - 1) * (1/(1 + nu)^2)
    return(f1)
  }
  
  last = nu[g]
  
  last1 = last/(1 - last)
  cand = rlnorm(1, meanlog = log(last1), sdlog = sqrt(cnu)) 
  cand1 = cand / (cand + 1)
  
  nu2 = nu; nu2[g] = cand1
  
  alfa1 = last1 * f.pr(last1);  alfa2 = cand * f.pr(cand)
  
  alfa0 = exp(sum(log(d.TIN.Exp(cc, y, r, mu, sigma2, alpha, nu2, Cens.type = Cens.type)))-
                sum(log(d.TIN.Exp(cc, y, r, mu, sigma2, alpha, nu, Cens.type = Cens.type))))
  
  alfa = alfa0 * (alfa2 / alfa1)
  out = ifelse(runif(1) < min(alfa, 1),  cand1, last)
  return(out)
}

dSEN = function(y, mu, sigma2, nu, log = F)
{
  PDF = vector(mode = "numeric", length = length(y))
  z = (y - mu) / sqrt(sigma2)
  PDF = nu + log(nu/sqrt(2 * pi * sigma2)) - 1.5 * log(nu + (z^2)/2) +
    zipfR::Igamma(1.5, nu + (z^2)/2, lower = FALSE, log = TRUE)
  for (j in which(PDF == -Inf)) PDF[j] = log(.Machine$double.eps)
  ifelse(log == F, return(exp(PDF)), return(PDF))
}
