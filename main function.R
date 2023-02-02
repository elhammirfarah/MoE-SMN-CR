Gibbs.MoE.Cr <- function(cc, y, r, x, Cens.type = c("Non", "Left", "Right"), G = NULL, type,  
                           mu0, Sigma0, alfa1, alfa2, nu01 = NULL, nu02 = NULL,
                           numethod = NULL, nuprior = NULL, hypernu0 = NULL,
                           nsim, burnin, pivotal.z = NULL)
{
  n = length(y) 
  Yobs = y
  type.upg = "mnl"
  if(G == 2) {type.upg = "logit"}
  if(Cens.type == "Left") {Lim1 =  rep(-Inf, n); Lim2 = y}
  if(Cens.type == "Right") {Lim2 =  rep(Inf, n); Lim1 = y}
  x = as.matrix(x)
  r = as.matrix(r)
  p = ncol(x) 
  q1 = ncol(r) 
  S0 = solve(Sigma0)
  z = pivotal.z
  beta.list = vector("list",G) 
  tau = matrix(0,nrow = q1, ncol = G)  
  sig2.list = c()
  nu.list = rep(0.05, G) 
  if((type == "T")||(type == "Slash")||(type == "VG")||(type == "SEN")){
  nu.list = rep(10, G)
  }
  gamma.list = rep(0.5, G)
  
  for(k in 1:G){
    ll = lm(y[z == k] ~ x[z == k, 2:p])
    beta.list[[k]] = as.vector(coefficients(ll), mode = "numeric")
    sig2.list[k] = mean(ll$residuals^2) 
  }
  n.iter <- nsim - burnin 
  Z = matrix(0,nrow = n.iter,ncol = n) 
  BETA.list = vector("list",G) 
  SIGMA.list = c()
  NU.list = c()
  TAU = vector("list", G) 
  
  for(i in 1:nsim){
    
    mu.j = c()
    mu.j = apply(cbind(x %*% matrix(unlist(beta.list), p, G), z), 1, function(x) x[x[G + 1]])
    if(type == "La"){
      a2 = (((y - mu.j) ^ 2) / sig2.list[z])
      ui = apply(cbind(1, a2), 1, function(x) GIGrvg :: rgig(1, lambda = -0.5, chi = 1, psi = x[2] ))
    }
    if(type == "T"){
      a1 = (nu.list[z] + 1) / 2
      a2 = ((((y - mu.j) ^ 2) / sig2.list[z]) + nu.list[z]) / 2
      ui = apply(cbind(a1, a2), 1, function(x) rgamma(1, x[1], rate = x[2]))
    }
    if(type == "Slash"){
      a1 = nu.list[z]+ 0.5
      a2 = ((y-mu.j)^2)/(2*sig2.list[z])
      ui = apply(cbind(a1, a2), 1, 
                 function(x) qgamma(runif(1, 0, pgamma(1, shape = x[1], rate = x[2])), 
                                    shape = x[1], rate = x[2]))
    }
    if(type == "VG"){
      a1 = (1 - nu.list[z])/2
      a2 = (((y - mu.j) ^ 2) / sig2.list[z])
      a3 = nu.list[z]
      ui = apply(cbind(a1, a2, a3), 1, function(x) GIGrvg :: rgig(1, lambda = x[1], chi = x[3], psi = x[2] ))
    }
    if(type == "CN"){
      meany2 = ((y - mu.j) ^ 2) / (2 * sig2.list[z])
      nuz    = nu.list[z]
      gammaz = gamma.list[z]
      pp1 = exp((gammaz - 1) * meany2)
      pn  = nuz * sqrt(gammaz) / ((1 - nuz) * pp1 + nuz * sqrt(gammaz))
      ui  = apply(cbind(gammaz, pn), 1, function(x) sample(c(1, x[1]), size = 1, prob = c(1-x[2], x[2])))
    }
    if(type == "TIN"){
      a1 = 1 - nu.list[z]
      a2 = ((y-mu.j)^2)/(2*sig2.list[z])
      ui = apply(cbind(a1, a2), 1, 
                 function(x) qgamma(runif(1, pgamma(x[1], shape = 1.5, rate = x[2]), pgamma(1, shape = 1.5, rate = x[2])), 
                        shape = 1.5, rate = x[2]))
    }
    if(type == "SEN"){
      a2 = nu.list[z] + ((y-mu.j)^2)/(2 * sig2.list[z])
      ui = apply(cbind(1.5, a2), 1, 
                function(x) qgamma(runif(1, pgamma(1, shape = 1.5, rate = x[2]), 1), 
                                   shape = 1.5, rate = x[2]))
    }
    if ( Cens.type != "Non"){
      y[cc == 1] = truncnorm::rtruncnorm(
        1,
        a = Lim1[cc == 1],
        b = Lim2[cc == 1],
        mean = mu.j[cc == 1],
        sd = sqrt(sig2.list[z[cc == 1]] )
      )
    }
    
    for(l in 1:G) 
    {
      X.l = as.matrix(x[z == l,]) 
      Y.l = y[z == l] 
      n.l = sum(z == l) 
      mu.j.l = mu.j[z == l]
      if(type != "N"){
        ui.l = ui[z == l]
        if(type != "La"){
        nu.l = nu.list[l]
        }}
      Vn.l = (alfa2 + t(Y.l - mu.j.l) %*% (Y.l - mu.j.l))/2
      nun.l = (alfa1 + n.l)/2
      sig2.l = 1/rgamma(1, nun.l, rate = Vn.l)
      sig2.list[l] = sig2.l
      SigmaA = solve(S0 + t(X.l) %*% (X.l)/sig2.l)
      Bn.l = SigmaA %*% (S0 %*% mu0 + t(X.l) %*% Y.l/sig2.l)
      beta.list[[l]] = matrix(mvtnorm :: rmvnorm(1, Bn.l, SigmaA), nrow = 1, ncol = p)
      
      if((type == "T") && (numethod == "Hierar")){
        nu.list[l] = MHnu2(nu.l, ui.l, nuprior, hypernu0) 
      }
      if(type == "Slash"){
        nu01 = hypernu0[1]
        nu02 = hypernu0[2]
        nu.list[l] = qgamma(runif(1,pgamma(1, shape = nu01+n.l, rate = nu02-sum(log(ui.l))),
                                  pgamma(20, shape = nu01+n.l, rate = nu02-sum(log(ui.l)))),
                            shape = nu01+n.l, rate = nu02-sum(log(ui.l)))
      }
      if((type == "VG") && (numethod == "Hierar")){
        nu.list[l] = MHnu2.VG(nu.l, ui.l, nuprior, hypernu0)
        }
    if(type == "CN"){
      m1 = n.gamma 
      m2 = n.l - n.gamma 
      nu.list[l] = qbeta(runif(1, pbeta(0.01, nu01 + m1, nu02 + m2), pbeta(0.5, nu01 + m1, nu02 + m2)),
                         nu01 + m1, nu02 + m2)
      Y.lb = Y.l[ui.l == gamma.l]
      X.lb = X.l[ui.l == gamma.l,]
      meany.l = t(Y.lb - X.lb %*% Bn.l) %*% (Y.lb - X.lb %*% Bn.l) / (2 * sig2.l)
      gamma.list[l] = qgamma(runif(1, pgamma(0.001, shape = 1 + 0.5 * n.gamma, rate = meany.l), 
                             pgamma(1, shape = 1 + 0.5 * n.gamma, rate = meany.l)),
                       shape = 1 + 0.5 * n.gamma, rate = meany.l)
    }
      if(type == "SEN"){
      nu.list[l] = rgamma(1, shape = n.l + 7.5, rate = 1 + sum(ui.l - 1))
      }
    }
    
    mu = x %*% matrix(unlist(beta.list), p, G)
    
    if ((type == "T") && (numethod == "NonHierar")){
      old.nu.list = nu.list
      for (l in 1:G) {
        nu.list[l] = MHnu(cc, Yobs, r, mu, sig2.list, tau[, -G], nu = old.nu.list, 
                          Cens.type = Cens.type, cnu = 2.4, 
                          g = l, nuprior, hypernu0)
      }}
    if ((type == "VG") && (numethod == "NonHierar")){
      old.nu.list = nu.list
      for (l in 1:G){
        nu.list[l] = MHnu.VG(cc, Yobs, r, mu, sig2.list, tau[, -G], nu = old.nu.list,
                          Cens.type = Cens.type, 
                          cnu = 2.4 , g = l, nuprior, hypernu0) 
      }}
    if(type == "TIN"){
      for (l in 1:G) {
        nu.list[l] = MHnuTIN(cc, Yobs, r, mu, sig2.list, tau, nu = nu.list, 
                             Cens.type = Cens.type, cnu = 0.005, hypernu0, g = l) 
      }
    }  
    eta0 = exp(r %*% tau) 
    pii = eta0 / rowSums(eta0) 
    p.z = PLOG = matrix(NA, nrow = n, ncol = G)
    for (l in 1:G) {
      mu.l = c(x %*% c(beta.list[[l]]))
      delta = (y - mu.l)/ sqrt(sig2.list[l])
      if(type == "N"){
      PLOG[, l] = log(pii[, l]) + dnorm(y, mu.l, sqrt(sig2.list[l]), log = T)}
      if(type == "La"){
        PLOG[, l] = log(pii[, l]) + dLaplace(delta, mu = 0, sigma = 1, log = T) - 0.5 * log(sig2.list[l]) 
      }
     if(type == "T"){
     PLOG[, l] = log(pii[, l]) + dt(delta, nu.list[l], log = TRUE) - 0.5 * log(sig2.list[l]) 
     }
      if(type == "Slash"){
        PLOG[, l] = log(pii[, l]) + dSlash(y, mu.l, sig2.list[l], nu.list[l], log = T)
      }
      if(type == "VG"){
        PLOG[, l] = log(pii[, l]) + dVG(y, mu.l, sig2.list[l], nu.list[l], log = T) 
      }
      if(type == "CN"){
        PLOG[, l] = log(pii[, l]) + dCN(y, mu.l, sig2.list[l], nu.list[l], gamma.list[l], log = T)
      }
      if(type == "TIN"){
        PLOG[, l] = log(pii[, l]) + dTIN(y, mu.l, sig2.list[l], nu.list[l], log = T)
      }
      if(type == "SEN"){
        PLOG[, l] = log(pii[, l]) + dSEN(y.j, mu.j, sig2.list[l], nu.list[l], log = T)
      }
    }
    
    if(G != 2){
      for (l in 1:G){
        p.z[, l] = 1 / rowSums(exp(PLOG - PLOG[, l]))
      }}else{
        p.z[,1] = 1 /(1 + exp(PLOG[,2] - PLOG[, 1]))
        p.z[,2] = 1 - p.z[,1] 
      } 
    
    z = apply(p.z, 1, function(x) which(rmultinom(1, 1, x) == 1))
    z = factor(z,levels = 1:G)
    n.z = as.vector(unname(table(z))) 
    if (any(n.z == 0)){
      print("MCMC terminated due to an empty class.")
      break}
    if (any(n.z == 1)){
      print("MCMC terminated due to a singleton class.")
      break}
    if(G == 2){
      z = as.numeric(z)
      z[z == 2] = 0
      tau.old = tau[,1]
      result.upg = UPG :: UPG(z, r, type.upg, Ni = NULL, draws = 1, burnin = 0,
                              A0 = 2000, d0 = 0.7, D0 = 0.7, G0 = 99, verbose  = FALSE,  
                              beta.start = tau.old)
      z[z == 0] = 2
      z = factor(z,levels = 1:G)
      tau = cbind(result.upg$posterior$beta, rep(0,q1))
    }else{
      tau.old = tau 
      result.upg = UPG :: UPG(z, r, type.upg, Ni = NULL, baseline = G, draws = 1, burnin = 0,
                              A0 = 2000, d0 = 0.7, D0 = 0.7, G0 = 99, verbose  = FALSE, 
                              beta.start = tau.old)
      tau = matrix(result.upg$posterior$beta, nrow= q1, ncol = G)
    }
    
    if (i > burnin)
    {
      j = i-burnin
      Z[j,] = z
      for(l in 1:G)
      {
        TAU[[l]]= rbind(TAU[[l]],tau[,l])
        BETA.list[[l]] = rbind(BETA.list[[l]],c(beta.list[[l]]))
      }
      SIGMA.list = rbind(SIGMA.list,c(sig2.list))
      if((type != "N") && (type != "La")){
        NU.list <- rbind(NU.list, c(nu.list))
      }
    }
    
   
  }
  ret_list = list(BETA = BETA.list,
                  TAU = TAU,
                  SIGMA = SIGMA.list,
                  NU = NU.list,
                  Z = Z,
                  )
  return(ret_list)
} 

