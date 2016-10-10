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


print.dyngame = function(x,...) {
  cat("\nA dynamic game:\n")
  fields = ls(x)
  for (field in fields) {
    try(cat(paste0("$",field, ": ",paste0(str(x[[field]]), collapse="\n"))), silent=TRUE)
  }
}


mywarning = function(txt) {
	warning(txt)
	print(txt)
}	

adapt.for.multistage = function(m,delta) {
 
  restore.point("adapt.for.multistage")
  
  m$g.org = m$g
  m$G.org = m$G
  m$delta.org = delta
  
  if (m$n.stage==1)
    return()
  
  delta.stage = delta^(m$n.stage)
    
  discount.vec = delta.stage^(m$n.stage-m$x.stage)
  
  discount.vec = discount.vec[m$ax[,"x"]]
  discount.mat = matrix(discount.vec,NROW(discount.vec),m$n)  

   
  #m$g.disc = m$g = m$g * discount.mat
  #m$G.disc = m$G = m$G * discount.vec

  m$g.disc = m$g = m$g * discount.mat #* m$n.stage
  m$G.disc = m$G = m$G * discount.vec #* m$n.stage
}

reset.from.multistage = function(m,ms=NULL) {
  
  restore.point("reset.from.multistage")
  
  if (m$n.stage==1)
    return()
  
  delta = m$delta.org
  delta.stage = delta^(m$n.stage)
  if (is.null(ms)) {
    m$delta = m$delta.org
    m$g = m$g.org
    m$G = m$G.org
    return()
  }
  
    
  #######################################################################
  # Adapt Equilibrium State Payoffs
  # There will be no discounting within a period
  #######################################################################

  # Some shortcuts
  g = m$g.org
  G = m$G.org
  sol = ms$sol.mat
  
  # Transition probabilities given optimal action rule
  p = ms$sol.mat[,"ae"]
	
	
  T  = m$tau.sub(m$tau[p,,drop=FALSE])
  # TW will contain transition probabilities in period. Non direct transitions
  # are included, i.e. we also describe for states in stage 1 what will be the
  # the final probability to end up in a particular state say in state 3.
  TW = matrix(0,m$nx,m$nx)
  # W will contain the excpected continuation payoffs only within the actual period
  # there will be no discounting between stages
  # UW will denote the corresponding continuation payoff that accounts for the
  # fact that we do not have discounting between stages  
  UW = W = rep(0,m$nx)      
  
  rownames(T) = rownames(TW) = names(UW) = names(W) = m$ax.lab[p]
  colnames(T) = colnames(TW) = m$x.lab
  
  stage = m$n.stage
  xs = which(m$x.stage == stage)
  W[xs] = G[p[xs]]
  TW[xs,] = T[xs,]
  for (stage in (m$n.stage-1):1) {
    xs = which(m$x.stage == stage)
    xn = which(m$x.stage == stage+1)
    W[xs]   = G[p[xs]]  + T[xs,] %*% W
    TW[xs,] = T[xs,] + T[xs,] %*% TW 
  }
  
  # Continuation payoffs in stage 1
  UW[xs] = solve(diag(NROW(xs))-delta*TW[xs,xs], (1-delta)*W[xs]) 
  
  stage = m$n.stage
  xs = which(m$x.stage == stage)
  UW[xs] = (1-delta)*W[xs] + delta * T[xs,] %*% UW  
  
  if (m$n.stage > 2) {
    for (stage in (m$n.stage-1):2) {
      xs = which(m$x.stage == stage)
      UW[xs] = (1-delta)* G[p[xs]] + T[xs,] %*% UW 
    }
  }
  ms$sol.mat[,"U"] = UW
  
  #######################################################################
  # Adapt Punishment Payoffs
  #######################################################################
  
  for (i in 1:m$n) {
    ai.col = paste("a",i,sep="")
    ax_i = m$ind.ax.to.ax_i[[i]][ms$sol.mat[,ai.col]]
    ai = get.full.dyn.vi(m,i,ax_i)$ax
    p = ai
    
    T  = m$tau.sub(m$tau[p,,drop=FALSE])
    Tw = matrix(0,m$nx,m$nx)
    vw = w = rep(0,m$nx) 
    
    rownames(T) = rownames(Tw) = names(vw) = names(w) = m$ax.lab[p]
    colnames(T) = colnames(Tw) = m$x.lab
    
    stage = m$n.stage
    xs = which(m$x.stage == stage)
    w[xs] = g[p[xs],i]
    Tw[xs,] = T[xs,]
    for (stage in (m$n.stage-1):1) {
      xs = which(m$x.stage == stage)
      xn = which(m$x.stage == stage+1)
      w[xs]   = g[p[xs],i]  + T[xs,] %*% w
      Tw[xs,] = T[xs,] + T[xs,] %*% Tw 
    }
    
    # Continuation punishment payoffs in stage 1
    vw[xs] = solve(diag(NROW(xs))-delta*Tw[xs,xs], (1-delta)*w[xs]) 
    
    stage = m$n.stage
    xs = which(m$x.stage == stage)
    vw[xs] = (1-delta)*w[xs] + delta * T[xs,] %*% vw  
    
    if (m$n.stage > 2) {
      for (stage in (m$n.stage-1):2) {
        xs = which(m$x.stage == stage)
        vw[xs] = (1-delta)* g[p[xs],i] + T[xs,] %*% vw 
      }
    }
  
    ms$sol.mat[,paste("v",i,sep="")] = vw
  }      
  
  # Don't reset delta, g and G earlier since we have a call 
  # get.full.dyn.vi(m,i,ax_i) above, which would not work properly
  m$delta = m$delta.org
  m$g = m$g.org
  m$G = m$G.org
}

solve.integrated = function(m,delta,tol.feasible = 1e-8) {
  
  restore.point("solve.integrated")

  stopifnot(m$n.stage==1)

  
  ms = new.dynsol(m)
  ms$name = paste("FB delta=",delta, m$name)
  
  m$delta = delta
  ms$delta = delta

  # admissible ax action profiles
  m$admiss = 1:m$nax
  m$admiss.ind = VectorListInd(m$na)
  names(m$admiss) = m$ax.lab
  
  ret.i = list()
  Ue.x = rep(NA,m$nx)
  
  print("")
  print("*************************************************************")
  print(paste("Integrated (First Best) Solution"))
  ret.e = get.highest.joint.payoff(m)
  names(ret.e$ax) = m$ax.lab[ret.e$ax]
  Ue.x = ret.e$Ue

  
  plot(1:m$nx,Ue.x,main=paste(ms$name),col="blue",pch=3)
  
  mat = matrix(NA,m$nx,1+(m$n+1)*2)
  colnames(mat)=c("x","ae",paste("a",1:m$n,sep=""),"U",paste("v",1:m$n,sep=""))
  mat[,"x"]  =  1:m$nx
  mat[,"ae"] =  ret.e$ax
  mat[,"U"]  =  ret.e$Ue
  rownames(mat) = m$ax.lab[ret.e$ax]
  ms$sol.mat = mat
  #reset.from.multistage(m,ms)
  print(paste("Integrated (First best) successfully solved for delta = ",delta))
  ms$adprob = get.average.discounted.prob(ms)
  make.extra.sol(ms)
  #ms$sol.mat = cbind(ms$sol.mat,ms$extra.sol)
  return(ms)
}    

  
  
solve.game = function(m,delta=m$delta,tol.feasible = 1e-8, verbose = interactive(), plots=verbose) {
  
  restore.point("solve.game")

  if (m$integrated) {
    return(solve.integrated(m,delta,tol.feasible = 1e-8))
  }   
  ms = new.dynsol(m)
  ms$name = m$name
  adapt.for.multistage(m,delta)
  delta.stage = delta^(m$n.stage)
  delta = delta.stage  
  
  m$delta = delta
  ms$delta = delta

  # admissible ax action profiles
  m$admiss = 1:m$nax
  m$admiss.ind = VectorListInd(m$na)
  names(m$admiss) = m$ax.lab
  
  ret.i = list()
  Ue.x = rep(NA,m$nx)
  
  ci.mat = matrix(NA,m$nax,m$n)

  n.sym = ifelse(m$symmetric,1,m$n)
  infeas.k = rep(TRUE,n.sym+1)
	
	
  iter = 0
  while(TRUE) {
    iter = iter+1
    
    if (verbose) {
      cat("\n")
      cat("\n*************************************************************")
      cat(paste("\nRound",iter))
    }
    
    #if (iter == 4)
    #  stop()
    
    # Calculate optimal equilibrium state actions
    # if some of the previous action profiles became infeasible
    if (infeas.k[1]) {
      if (verbose) {
        print("get.highest.joint.payoff")
      }
      ret.e = get.highest.joint.payoff(m)
      names(ret.e$ax) = m$ax.lab[ret.e$ax]
      Ue.x = ret.e$Ue
    }
    
    V.x = rep(0,m$nx)
		
		
		for (i in 1:m$n) {
			# Calculate optimal punishment profiles for player i
			# if some of the previous optimal punishment profiles were infeasible			
			if (infeas.k[min(i+1,n.sym+1)]) {	
				if (i==1 | (!m$symmetric)) {
          if (verbose) {
					  print(paste("get.harshest.punishment",i))
          }
					ret.i[[i]] = get.harshest.punishment(m,i)
					names(ret.i[[i]]$ax) = m$ax.lab[ret.i[[i]]$ax]				
					v1 = ret.i[[i]]$vi
				} else {
					ret.i[[i]] = list()
					ret.i[[i]]$vi=v1[m$x.perm[,i]]
				}
				vi = ret.i[[i]]$vi
				#Cheating payoffs for all admisisble ax given the just
				#calculated punishment payoff vefor all states x
				ci.mat[,i] = get.cheating.payoffs.ax(m,i,delta=delta,v = vi,
																						 admiss=m$admiss)
			}
			V.x = V.x + vi
		}
		
    EUe = as.numeric((1-delta)*m$G[m$admiss] + delta * (m$tau[m$admiss,,drop=FALSE] %*% Ue.x))
    infeasible = which(EUe-rowSums(ci.mat) < -tol.feasible)
    infeasible.ax = m$admiss[infeasible]             
    
    #Give names
    rownames(ci.mat) = names(EUe)        = m$ax.lab[m$admiss]
    names(infeasible) = names(infeasible.ax) =  m$ax.lab[infeasible.ax]

#     print("Ue.x:")
#     print(Ue.x)
#     print("V.x:")
#     print(V.x)
#     print("ret.i[[1]]$ax:")
#     print(ret.i[[1]]$ax)
#     print("ret.i[[2]]$ax:")
#     print(ret.i[[2]]$ax)
#     print("ret.e$ax:")
#     print(ret.e$ax)
#     print("infeasible:")
#     print(infeasible)
#     
    # Check whether all optimal action profiles are feasible
    # Note that policies are indixed on ax (not on admiss)
    infeas.k[1] = any(ret.e$ax %in% infeasible.ax)
		for (i in 1:n.sym) {
			infeas.k[i+1] = any(ret.i[[i]]$ax %in% infeasible.ax)
		}
		# If none of the optimal action profiles is infeasible, we can stop
    if (!any(infeas.k)) break;
    
    # Remove infeasible action profiles from the set of admissible action profiles 
    m$admiss = m$admiss[-infeasible]
    m$admiss.ind = remove.elements(m$admiss.ind,infeasible)
    
#     print("admiss (end of round):")
#     print(m$admiss)

    if (length(infeasible)>0) {
      ci.mat = ci.mat[-infeasible,]
    }
    #browser()
    if (min(m$admiss.ind$ncols) <= 0) {
      warning("There does not exist a subgame perfect equilibrium")
      reset.from.multistage(m)
      ms$sol.exists = FALSE
      return(ms)
    }
    if (plots) {
      plot(1:m$nx,Ue.x,ylim=range(c(Ue.x,V.x)),main=paste(iter,ms$name),col="blue",pch=3)
      points(1:m$nx,V.x,col="red")
    }
  }
  
  
  if (plots) {
    plot(1:m$nx,Ue.x,ylim=range(c(Ue.x,V.x)),main=paste("*",iter,ms$name),col="blue",pch=3)
    points(1:m$nx,V.x,col="red")
  }

  ax.i = vi = list()
  for (i in 1:m$n) {
		if (i<=n.sym) {
			ax.i[[i]] = ret.i[[i]]$ax
		} else {
			ax.i[[i]] = rep(NA,m$nx)
		}
    vi[[i]] = ret.i[[i]]$vi
  }
  
  sol = list(x=1:m$nx,ax.e = ret.e$ax, ax.i = ax.i, Ue = ret.e$Ue, vi = vi)
  mat = matrix(unlist(sol),m$nx,1+(m$n+1)*2)
  colnames(mat)=c("x","ae",paste("a",1:m$n,sep=""),"U",paste("v",1:m$n,sep=""))
  rownames(mat) = m$x.lab
  ms$sol.mat = mat
  #m$dist.e = get.eq.dist(m)
  reset.from.multistage(m,ms)

  cat(paste("\nGame successfully solved for delta = ",delta))

  ms$adprob = get.average.discounted.prob(delta=ms$delta, ax=ms$sol.mat[,"ae"], tau=m$tau)
  make.extra.sol(ms)
  #ms$sol.mat = cbind(ms$sol.mat,ms$extra.sol)
  return(ms)
}    


# Calls the users' extra sol fun
make.extra.sol = function(m,ax=m$sol.mat[,"ae"]) {
  
 	restore.point("set.g.and.tau")
 	
 	if (is.null(m$extra.sol.fun)) {
   	m$extra.sol.cur  = m$extra.sol.ad= NULL
   	return()
 	}
 	avm = ax.to.av(m,ax)
 	# Yields values for every value of x
 	m$extra.sol.cur = m$extra.sol.fun(avm,m$xv.mat,m=m)
 	
 	# Average discounted values
 	m$extra.sol.ad = as.matrix(m$adprob %*% m$extra.sol.cur)
 	colnames(m$extra.sol.ad) = paste(colnames(m$extra.sol.cur),"ad",sep="")
}


get.eq.dist = function(m, iter = 10) {
  aex = m$sol$ax.e
  tau.e = m$tau.sub(m$tau[aex,,drop=FALSE])
  ret = NULL
  if (!is.finite(iter)) {
    ret = try(solve(diag(m$nx)-tau.e))
  } else {
    ret = tau.e
    for (i in 1:iter) {
      ret = tau.e %*% ret
    }
  }
  return(ret)
}



#' Calculates the highest joint payoff
#' The returned policy are indexed on ax (not on admiss)
get.highest.joint.payoff = function(m) {
  
  restore.point("get.highest.joint.payoff")
  
  T = m$tau[m$admiss,,drop=FALSE]
  R = m$G[m$admiss]
  
  ret = policy_iteration(T=T,R=R,delta=m$delta,na=m$admiss.ind)
  names(ret$p) = names(ret$V) =  m$x.lab
  #print(rbind(label.a(m,opt.a), ret$V))

  return(list(Ue=ret$V, ax = m$admiss[ret$p]))
}



# Counters
# nax_i        : Number of ax_i, i.e. action profiles of other players summed over all states
# nad_i        : Number of admissible ax_i
#
# Description of some used variables:
#
# admiss.ax_i    : nad_i      ; ax_i index of admissible ax_i
# not.admiss.ax_i: nax_i-nad_i; ax_i index of non-admissible ax_i
#
# ind.ad.ax_i.by.x: nad_i; VectorListInd that maps admissible ax_i to x
 
# ci.mat       : nad*n; for each player i his cheating payoffs for any admissible action profile
# EUe          : nad*1; the joint payoffs for any given action profile using Ue.x as continuation payoff 
# infeasible   : nad*1; boolean whether a previously admissible profile becomes infeasible
# infeasible.ax: nad*1; ax index of those infeasible profiles

# Ue.x         : nx * 1; Joint continuation payoff in state x
# V.x          : nx * 1; sum of punishment payoffs in state x

get.harshest.punishment = function(m,i, tol=1e-8) {
  
  restore.point("get.harshest.punishment")
  
  
  tab = tabulate(m$ind.ax.to.ax_i[[i]][m$admiss],nbin=length(m$ind.ax_i.by.x[[i]]))
  admiss.ax_i = which(tab > 0) 
  not.admiss.ax_i = which(tab==0)
  
  # Only for labeling purposes
  ax.of.ax_i = match(1:max(m$ind.ax.to.ax_i[[i]]),m$ind.ax.to.ax_i[[i]])
  ax_i.lab = m$ax.lab[ax.of.ax_i]
  names(admiss.ax_i)=ax_i.lab[admiss.ax_i]
  names(not.admiss.ax_i)=ax_i.lab[not.admiss.ax_i]
  
  admiss.ax_i
  not.admiss.ax_i
  
  # Create a VectorListIndex on the admissible ax_i
  ind.ad.ax_i.by.x = remove.elements(m$ind.ax_i.by.x[[i]], not.admiss.ax_i)
   
  # Start with action profiles that minimize player i's static cheating payoff
  # returns a cheating payoff for every admissible ax_i profile
  cheat.pay = get.cheating.payoffs(m,i,delta = 0, admiss = admiss.ax_i)
  names(cheat.pay) = names(admiss.ax_i)
  cheat.pay
  
  # Get for every state that admissible action profile that minimizes
  # player i's cheating payoff
  # act.ax_i indexes on ax_i (not on admissible ax_i)
  act.ax_i = admiss.ax_i[which.RowMins(ind.ad.ax_i.by.x,cheat.pay)]
  
  # Get corresponding cheating payoffs
  # We solve a MDP for player i
  v = get.full.dyn.vi(m,i,act.ax_i)$vi
  v
  old.cheat.pay = cheat.pay
  counter = 0
  while (TRUE) {
    counter = counter+1
    cheat.pay = get.cheating.payoffs(m,i,delta = m$delta, v=v, admiss = admiss.ax_i)
    
    #cheat.pay
    #old.cheat.pay
    # Stop if player i cannot improve his cheating payoff in any state
    if (approxeq(cheat.pay,old.cheat.pay,tol)) {
      break;
    }
    
    act.ax_i = admiss.ax_i[which.RowMins(ind.ad.ax_i.by.x,cheat.pay)]
  
    # Get corresponding cheating payoffs
    # We solve a MDP for player i
    v = get.full.dyn.vi(m,i,act.ax_i)$vi
    
    old.cheat.pay = cheat.pay
    
    #cheat.pay
    #v
    #act.ax_i
    
  }
  #print(rbind(label.ax(m,a.to.ax(m,act.a)), v))
  cat("\n\nget.harshest.punishment iterations: ", counter,"\n")
  
  # Translate ax_i to some admissible ax profile
  act.ax = m$admiss[match(act.ax_i,m$ind.ax.to.ax_i[[i]][m$admiss])] 
  return(list(vi=v,ax = act.ax))
}    

 

# Calculates the cheating payoffs from all admissible action profiles
# admiss.ax_i    : nad_i      ; ax_i index of admissible ax_i
# m$replies[[i]] : a vector based on a vector list of dimensions nax_i * na_i[x]
# g.replies      : a vector based on a vector list of dimensions nax_i * na_i[x]
# cheat.payoff   : nax_i ; later reduced to nad_i

# Returns player i's cheating payoff for every admissible ax_i
get.cheating.payoffs = function(m, i, delta = 0,v=rep(0,m$nx), admiss.ax_i = NULL) {
  
  restore.point("get.cheating.payoffs")
    
  g.replies = m$g[m$replies[[i]],i] #Need to transpose, since by default columns instead of rows are stacked
    
  if (delta == 0) {
    cheat.payoff = RowMaxs(m$ind.replies.by.ax_i[[i]],g.replies)
  } else {
    V.ax = m$tau %*% v
    V.replies = V.ax[m$replies[[i]]]    
    cheat.payoff = RowMaxs(m$ind.replies.by.ax_i[[i]],
                    (1-delta)*g.replies + delta*V.replies)
  }
  # Map cheating payoff to all admissible action profiles
  if (!is.null(admiss.ax_i)) {
    return(cheat.payoff[admiss.ax_i])
  } else {
    return(cheat.payoff)
  }   
}    

# Returns player i's cheating payoff for every admissible ax profile
# This means we have duplication as several ax profiles correspond to one ax_i profile
get.cheating.payoffs.ax = function(m, i, delta = 0,v=rep(0,m$nx),admiss=NULL) {
  
  restore.point("get.cheating.payoffs.ax")

  # ax_i index
  cheat.payoff = get.cheating.payoffs(m,i,delta,v)

  if (is.null(admiss)) {
    return(cheat.payoff[m$ind.ax.to.ax_i[[i]] ])
  } else {
    return(cheat.payoff[m$ind.ax.to.ax_i[[i]][admiss]])
  }    
}    


# ax_i : nx * 1; describes the ax_i indeces of ax_i profiles played in every state x
get.full.dyn.vi = function(m,i,ax_i) {
  
  restore.point("get.full.dyn.vi")

  replies = m$replies[[i]][rows.to.v.ind(m$ind.replies.by.ax_i[[i]],ax_i)]
  ind.replies = make.child.VectorListInd(m$ind.replies.by.ax_i[[i]],rows = ax_i)  

  #names(replies) = m$ax.lab[replies]
  
  # Reward function
  R = m$g[replies,i]  
  # Transition function between states
  T = m$tau[replies,,drop=FALSE]

  ret = policy_iteration(T,R,delta=m$delta,na=ind.replies)
  names(ret$p) = names(ret$V) =  m$x.lab
  opt.ax = replies[ret$p]
  return(list(vi=ret$V, ax=opt.ax))
}  







