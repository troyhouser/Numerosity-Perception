KL = function(p,q){
  return(rowSums(q*(log2(q)-log2(p))))
}

compute_q_nk = function(ns,ks,p_n,p_k,lambda,sm=1e-10){
  #compute posterior Q for given lambda (lagrange multiplier)
  lambda = matrix(lambda,length(lambda),1)
  p_n = matrix(p_n,length(p_n),1)
  q_nk = matrix(0,15,24)
  m = sweep(ns,2,ks)**2
  for(k in 1:15){
    q_nk[k,] = p_n[k] * m[k,]
    q_nk[k,] = q_nk[k,] / lambda[k]
  }
  q_nk = exp(-q_nk)
  q_nk = sweep(q_nk,2,p_k,"*")
  q_nk = q_nk / rowSums(q_nk)
  q_nk = q_nk + sm
  q_nk = q_nk/rowSums(q_nk)
  return(q_nk)
}

find_q_nk = function(ns,ks,p_n,p_k,info_bound,n_steps=1500){
  lams = rep(1,nrow(ns))*0.5
  q_nk = compute_q_nk(ns,ks,p_n,p_k,lams)
  ents = KL(p_k,q_nk)
  
  for(i in 1:n_steps){
    diffs = ents - info_bound
    deltas = diffs * 0.025
    lams = exp(log(lams) + deltas)
    q_nk = compute_q_nk(ns,ks,p_n,p_k,lams)
    ents = KL(p_k,q_nk)
  }
  return(q_nk)
}

PP = function(x,a){
  p = 1/x**a
  return(p/sum(p))
}

step_n=1;step_k=1;alpha=1
min_n=1;max_n=15;min_k=step_k;max_k=24
ks = seq(min_k,max_k,by=step_k)
ns = seq(min_n,max_n,by=step_n)
ns = matrix(ns,length(ns),1)
p_ks = PP(ks,alpha)
p_ns = PP(ns,alpha)
info_bound=3
ns = matrix(ns,length(ns),length(ks))
Q = find_q_nk(ns,ks,p_ns,p_ks,3)

matplot(1:24,t(Q),"l")
plot(Q)
