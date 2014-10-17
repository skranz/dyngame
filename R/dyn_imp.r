# Solve stochastic games with imperfect monitoring. 
# For the moment, I assume that only states can be observed

# In the static problem, we optimize over a payment vector of size nx 

# Storage groups
#
# ax:   one entry per action profile / state combination
#       indexed by VectorListInd
#
# ax_i: one entry per state and action profile of other players than player i
#       store replies on this level and use to calculate cheating payoffs
#       ind.ax_i: assigns to each ax the corresponding ax_i  

# Problem we need a notation that makes clear to which list indices refer
# nax : number of ax profiles
# nad : number of admissible action profiles



# Description of some used variables:
#
# nad          : NROW(m$admiss)
# m$admiss     : nad*1; ax index of admissible action profiles
# ci.mat       : nad*n; for each player i his cheating payoffs for any admissible action profile
# EUe          : nad*1; the joint payoffs for any given action profile using Ue.x as continuation payoff 
# infeasible   : nad*1; boolean whether a previously admissible profile becomes infeasible
# infeasible.ax: nad*1; ax index of those infeasible profiles

# Ue.x         : nx * 1; Joint continuation payoff in state x
# V.x          : nx * 1; sum of punishment payoffs in state x



solve.LP.e = function(m,ax,U,v,mod.file = "D:/libraries/dyngame/LPe.mod",dat.file = "D:/libraries/dyngame/LPk.dat") {
  make.lp.k.constraints.GMPL(m,ax,U,v)
  res = gmpl.solve(mod.file,dat.file)
  print(res)
  return(res)
}

solve.LP.i = function(m,i=i,ax,U,v,mod.file = "D:/libraries/dyngame/LPe.mod",dat.file = "D:/libraries/dyngame/LPk.dat") {
  make.lp.k.constraints.GMPL(m,ax,U,v)
  res = gmpl.solve(mod.file,dat.file)
  print(res)
  return(res)
}



make.lp.k.constraints.GMPL = function(m,ax,U,v,delta=m$delta, mod.file = "D:/libraries/dyngame/LPe.mod",dat.file = "D:/libraries/dyngame/LPk.dat") {
  library(rgmpl)
  
  n = m$n
  nx = m$nx
  # Give all excess liquidity to player 1
  u = v
  u[,1] = U - colSums(v)
  
  REPLIES = list()
  for (i in 1:n) {
    REPLIES[[i]] = m$ax.lab[get.replies(m=m,i=i,ax=ax,keep.ax=FALSE)]
  }
  sets = list(N=1:n,X=m$x.lab,Ax=m$ax.lab,REPLIES =REPLIES) 
  param = list(pi = m$g,tau=m$tau,u=u,v=v,delta=delta,ax=m$ax.lab[ax])
  
  gmpl.make.dat.file(mod.file=mod.file,dat.file=dat.file,sets = sets, param=param)
}

# We optimize over a (nx * n) x 1 payment vector
# The first nx elements correspond to the payments of player 1 etc
make.lp.k.constraints = function(m,ax,U,v,delta) {
  nx = m$nx
  n = m$n
  
  # Give all excess liquidity to player 1
  u = v
  u[,1] = U - colSums(v)
  
  ##############################################################
  # Payment constraints
  #
  # (1-delta)p_i(x',y)<= u_i(x',y)-v_i(x')
  ##############################################################
  nPC =  nx*n # Number of non-zero entries

  rowPC = 1:nPC # Rows
  colPC = 1:nPC # Columns
  valPC = rep(1-delta,nPC) # Values 

  dirPC = rep("<=",nPC) # Dir
  rhsPC = as.numeric(u-v)
  
  ##############################################################
  # Budget constraints
  #
  # sum(p_{i}(x',y)) >= 0
  ##############################################################
  nBC = nx * n  # We have nx constraints but nx*n non-zeros
  
  rowBC = (1:nBC)+nPC # Rows
  colBC = rep(1:n,times=nx) + rep(1:nx,each=n) - 1  # Columns
  valBC = rep(1,nBC) # Values  
  dirBC = rep(">=",nBC) # Dir
  rhsBC = rep(0,nBC)
  
  ##############################################################
  # Action constraints
  #
  ##############################################################

  # We first need a list of possible replies of player i
  make.ac.i = function(m,i,ax,u,v,delta) {
    ax.hat = get.replies(m,i,ax,keep.ax=FALSE) # Possible replies    
    nrep = length(ax.hat)
    
    nAC = nx * nrep # A payment of player i for every reply and every state

    tau.ax = rep(m$tau[ax,],times=nrep)
    tau.ax.hat = as.numeric(t(m$tau[ax.hat,]))

    # For a single ax.hat
    # ri  = 1:nx
    # ci  = 1:nx + (i-1)*nx
    # val = - delta * (1-delta) * (tau.ax - tau.ax.hat) 
    # rhs = - ((1-delta) * (m$g[ax,i]-m$g[ax.hat,i]) + delta * (tau.ax - tau.ax.hat) * u[,i])

    row  = 1:nAC
    col  = rep(1:nx + (i-1)*nx,nrep) # p_i is repeated nrep times
    val = - delta * (1-delta) * (tau.ax - tau.ax.hat)
    rhs = - ((1-delta) * (m$g[ax,i]-m$g[ax.hat,i]) + delta * (tau.ax - tau.ax.hat) * u[,i])
    dir = rep(">=",nAC)
    return(list)
  } 
  
  
  
  
  
}
