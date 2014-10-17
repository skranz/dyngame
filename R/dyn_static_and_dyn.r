# # Quickly solving the game with static, separate additive activities

# Idea: Consider a game where action profiles can be split in dynamic and static activites
#  
# Dynamic activities, e.g. investment, influence state transition
# Static activities, e.g. price setting, do not not influence state transition
#  
# We assume that payoffs are additively separable in static in dynamic activities
#                  g(a) = g.dyn(a.dyn) + g.stat(a.stat)

# This means that also cheating payoffs are additively separable:
#                  c(a) = c.dyn(a.dyn) + c.stat(a.stat)
#                  
# An action profile can be implemented given continuation payoffs U if and only if
#                  (1-delta)G(a)+ delta*EU(a) >= C(a)
#                  (1-delta)[G.dyn(.)+G.stat(.)] + delta*EU(a) >= C.dyn(.)+(1-delta)*C.stat
#                  (1-delta)[C.stat-G.stat] <= (1-delta)[G.dyn]+delta*EU(a) - C.dyn 
#                  (1-delta)L.stat(.) <= U-V
#                  L.stat <= U_V / (1-delta)
                 
                  

make.static.admiss.and.change = function(m, U_V.ax, delta=m$delta, ret.e=NULL, ret.i = NULL, tol=m$tol) {
  
  restore.point("make.static.admiss.and.change")

  stopifnot(!is.null(m$stat.ag))
  
  # equilibrium path
  static.G = rep(0,NROW(m$admiss))
  
  ag = m$stat.ag
  for (x in 1:m$nx) {
    rows = rows.to.v.ind(m$ind.admiss,x)
    static.G[rows] = ag$opt[[x]]$fun.e(U_V.ax[rows] / (1-delta))
    rows = rows + 1
  }
  infeasible = which(!is.finite(static.G))

  # punishment
  static.ci = matrix(NA,length(m$admiss),m$n)
  for (i in 1:m$n) {        
    ci = rep(0,length(m$admiss))
    if (is.null(m$stat.ag)) {
      ci[U_V.ax<0] = Inf
      static.ci[,i] = ci
      break
    }  
    
    ag = m$stat.ag
    for (x in 1:m$nx) {
      rows = rows.to.v.ind(m$ind.admiss,x)
      ci[rows] = ag$opt[[x]]$fun.i[[i]](U_V.ax[rows] / (1-delta))
      rows = rows + 1
    }  
    static.ci[,i] = ci
  }

  # Check whether payoff for previously optimal profiles changed
  if (!is.null(ret.e)) {
    m$change.e = !approxeq(m$static.G[ret.e$ax.admiss],
                           static.G[ret.e$ax.admiss],tol=tol)                           
  }
  
  if (!is.null(ret.i)) {
    m$change.i = rep(FALSE,m$n)
    for (i in 1:m$n) {
      m$change.i[i] = !approxeq(m$static.ci[ret.i[[i]]$ax.admiss],
                                static.ci[ret.i[[i]]$ax.admiss],tol=tol)
    }                           
  }
  
  # Update admissible action profiles
  
  if (length(infeasible) >0) {
    m$admiss = m$admiss[-infeasible]
    m$ind.admiss = remove.elements(m$ind.admiss,infeasible)
    m$static.G = static.G[-infeasible]
    m$static.ci = static.ci[-infeasible,]
    m$dyn.ci = m$dyn.ci[-infeasible,]
  } else {
    m$static.G = static.G
    m$static.ci = static.ci
  }
  
  names(m$static.G) = rownames(m$static.ci) = names(m$admiss) = m$ax.lab[m$admiss]

}

solve.game = function(m,delta,tol = m$tol, tol.feasible = tol, augment.sol.mat = FALSE) {
  
  restore.point("solve.game")
  ms = new.dynsol(m)
  m$delta = delta
  ms$delta = delta
  ms$name = paste("delta=",delta, m$name)
    
  ret.i = list()
  Ue.x = rep(NA,m$nx) 
  infeas.k = rep(TRUE,m$n+1)
  
 # admissible ax action profiles
  m$admiss = 1:m$nax
  m$ind.admiss = VectorListInd(m$na)

  m$dyn.ci = matrix(NA,NROW(m$admiss),m$n)

  U_V.ax = rep(Inf,m$nax)
  make.static.admiss.and.change(m,U_V.ax)
  m$change.e = TRUE
  m$change.i = rep(TRUE,m$n)

  no.change.count = 0
  
  iter = 0
  while(TRUE) {
    iter = iter+1
    print("")
    
    # Calculate optimal equilibrium state actions
    # if some of the previous action profiles became infeasible
    if (m$change.e) {
      print("get.highest.joint.payoff")
      ret.e = get.highest.joint.payoff(m)
      Ue.x = ret.e$Ue
    }
    
    V.x = rep(0,m$nx)

    for (i in 1:m$n) {
      # Calculate optimal punishment profiles for player i
      # if some of the previous optimal punishment profiles were infeasible
      if (m$change.i[i]) {
        print(paste("get.harshest.punishment",i))
        ret.i[[i]] = get.harshest.punishment(m,i)
        
        #Cheating payoffs for all admisisble ax given the just
        #calculated punishment payoffs for all states x
        m$dyn.ci[,i] = get.cheating.payoffs(m,i,delta=delta,v = ret.i[[i]]$vi)
      }
      V.x = V.x + ret.i[[i]]$vi  # Punishment payoff per state
    }
  
        
    Ud.ax = as.numeric((1-delta)*m$G[m$admiss] + delta * (m$tau[m$admiss,] %*% Ue.x))
    Vd.ax  = rowSums(m$dyn.ci)
    U_V.ax = Ud.ax-Vd.ax
    
    # Need this line to avoid rounding errors
    # Marginally infeasible action profiles
    # will still be considered feasible
    U_V.ax[U_V.ax<0 & U_V.ax > -tol.feasible] = 0

    old.change.e = m$change.e; old.change.i = m$change.i;
    make.static.admiss.and.change(m,U_V.ax,ret.e=ret.e, ret.i=ret.i)
     
    # Stop if in two consequtive rounds there was no change
    if (!m$change.e & (!any(m$change.i))) {             
      if (no.change.count == 0) {
        no.change.count = 1;
        m$change.e = !old.change.e;
        m$change.i = !old.change.i;
      } else {
       break;
      }
    } else {
      no.change.count = 0;
    }
    
    plot(1:m$nx,Ue.x,ylim=range(c(Ue.x,V.x)),main=paste(iter,ms$name),col="blue",pch=3)
    points(1:m$nx,V.x,col="red")
  }
  
  
  
  # Check whether an equilibrium exists
  if (any(Ue.x==-Inf) | any(V.x==Inf))
    stop("There does not exist a subgame perfect equilibrium")

  
  plot(1:m$nx,Ue.x,ylim=range(c(Ue.x,V.x)),main=paste("*",iter,ms$name),col="blue",pch=3)
  points(1:m$nx,V.x,col="red")
  
  ax.i = vi = list()
  for (i in 1:m$n) {
    ax.i[[i]] = ret.i[[i]]$ax
    vi[[i]] = ret.i[[i]]$vi
  }
  
  sol = list(ax.e = ret.e$ax, ax.i = ax.i, Ue = ret.e$Ue, vi = vi)
  mat = matrix(unlist(sol),m$nx,(m$n+1)*2)
  colnames(mat)=c("ax.e",paste("ax.",1:m$n,sep=""),"U",paste("v",1:m$n,sep=""))
  rownames(mat) = m$x.lab
  if (augment.sol.mat) {
    x = 1:m$nx
    delta = rep(ms$delta,m$nx)
    mat = cbind(mat,x,delta)
  }
    
  ms$sol.mat = mat
       
  
  return(ms)
}    


get.highest.joint.payoff = function(m) {
  
  restore.point("get.highest.joint.payoff")
  
  # Call with all action profiles
  # to improve speed it would be better to call only with profiles  
  # that are not infeasible (static.G > -Inf)
  T = m$tau[m$admiss,]
  R = m$G[m$admiss] + m$static.G
  
  ret = policy_iteration(T=T,R=R,delta=m$delta,na=m$ind.admiss)
  names(ret$p) = names(ret$V) =  m$x.lab
  #print(rbind(label.a(m,opt.a), ret$V))

  return(list(Ue=ret$V, ax = m$admiss[ret$p], ax.admiss = ret$p))
}

get.harshest.punishment = function(m,i, tol=1e-8) {
  
  restore.point("get.harshest.punishment")
  delta = m$delta
  static.ci = m$static.ci[,i]
  
  # Start with action profiles that minimize player i's static cheating payoff
  cheat.pay = get.cheating.payoffs(m,i,delta = 0) + static.ci
  ax.admiss = which.RowMins(m$ind.admiss, cheat.pay)
  
  names(cheat.pay) = m$ax.lab[m$admiss]
  names(ax.admiss) = m$ax.lab[m$admiss[ax.admiss]]
  
  # Get corresponding cheating payoffs
  # We solve a MDP for player i
  v = get.full.dyn.vi(m,i,ax.admiss)$vi
  
  iter = 0
  
  tol.counter = 0
  while (TRUE) {
    iter = iter +1
    # Get action profiles that are optimal given the 
    # actual values of player i's continuation payoffs v
    cheat.pay = get.cheating.payoffs(m,i,delta = delta, v=v) + (1-delta)*static.ci
    old.v = v
    old.ax.admiss = ax.admiss
    
    mat = as.matrix(m$ind.admiss,cheat.pay)
    rm1 = RowMins(m$ind.admiss, cheat.pay)
    rm2 = rowMins(mat)
    rma = RowMaxs(m$ind.admiss, cheat.pay)
    
    rownames(mat) = names(rm1) = names(rm2) = m$x.lab
    rbind(rm1,rm2,rma)
    
    
    ax.admiss = which.RowMins(m$ind.admiss, cheat.pay)
    
    names(ax.admiss) = m$ax.lab[m$admiss[ax.admiss]]

    #print(label.a(m,act.a))
    v = get.full.dyn.vi(m,i,ax.admiss)$vi
    #plot(act.a)
    change = max(old.v-v)
    if (change <= 0) break;
    if (change <= tol) { 
      tol.counter = tol.counter +1
    } else {
      tol.counter = 0
    }
    if (tol.counter > 10) {
      warning("get.harshest.punishment stoped since change was 10 times in a row below tol. Yet it was still positive")
      break
    }
  }
  #if (change > -tol) {
  #  warning(paste("get.harshest.punishment: v increased in last step by ", change))
  #}
  ax.admiss = old.ax.admiss
  v = old.v
  return(list(vi=v,ax.admiss = ax.admiss, ax = m$admiss[ax.admiss]))
}    
 

# Returns for every ax the dynamic cheating payoff of player i if continuation payoffs are given by v
# Not yet effectively implemented: don't check which action profiles are admissible
get.cheating.payoffs = function(m, i, delta = 0,v=rep(0,m$nx),admiss=m$admiss) {
  
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
  if (!is.null(admiss)) {
    return(cheat.payoff[m$ind.ax.to.ax_i[[i]][admiss]])
  } else {
    return(cheat.payoff[m$ind.ax.to.ax_i[[i]]])
  }   
}    

# ax.admiss is a nx x 1 vector of action profiles
# one profile for every state x
# static cheating payoffs are already added
get.full.dyn.vi = function(m,i,ax.admiss) {
  
  restore.point("get.full.dyn.vi")
  
  ax_i = m$ind.ax.to.ax_i[[i]][m$admiss[ax.admiss]]
  replies = m$replies[[i]][rows.to.v.ind(m$ind.replies.by.ax_i[[i]],ax_i)]

  ind.replies = make.child.VectorListInd(m$ind.replies.by.ax_i[[i]],rows = ax_i)  
  static.ci = const.rows(ind.replies,m$static.ci[ax.admiss,i])

  names(replies) = names(static.ci) = m$ax.lab[replies]
  names(ax.admiss) = m$ax.lab[m$admiss[ax.admiss]]
  
  # Reward function
  R = m$g[replies,i] + static.ci  
  # Transition function between states
  T = m$tau[replies,,drop=FALSE]

  ret = policy_iteration(T,R,delta=m$delta,na=ind.replies)
  names(ret$p) = names(ret$V) =  m$x.lab
  opt.ax = replies[ret$p]
  return(list(vi=ret$V, ax=opt.ax))
}  






