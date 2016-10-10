# the payoffs we would get if ae were played
# in every state but no payments would be allowed
compute.no.payments.u = function(sol, regime="e") {
  restore.point("compute.no.payments.u")
  
  sol.mat = sol$sol.mat
  m = sol$get.m()
  
  player = 1
  # equilibrium regime action policy
  p = sol.mat[,paste0("a",regime)]
  
  # transmission probabilities
  Tp = m$tau[p,,drop=FALSE]
  nx = m$nx
  delta = m$delta
  npu = lapply(1:2, function(player) {
    Rp = m$g[p,player, drop=FALSE]
    u = solve(diag(nx)-delta*Tp, (1-delta)*Rp)
  })
  npu

}

plot.payoff.set = function(state = 1, sol, rsg.sol = NULL, all.iter=!is.null(rsg.sol$ipoints), tight.approx=FALSE, main=paste0("State ", state), show.np.e = !is.null(rsg.sol)) {
  restore.point("plot.payoff.set")
  
  np.e = compute.no.payments.u(sol,"e")

  G.e = sol$get.m()$g[sol$sol.mat[state,"ae"],]
  # plot for state 
  tpoints = with(as.data.frame(sol$sol.mat[state,,drop=FALSE]),rbind(
    c(v1,v2),c(v1,U-v1),c(U-v2,v2),c(v1,v2)
  ))
  
  if (all.iter){
    ipoints = rsg.sol$ipoints[[state]]
    xlim = range(c(tpoints[,1],ipoints[,1]))
    ylim = range(c(tpoints[,2],ipoints[,2]))
  } else {
    xlim = range(c(tpoints[,1]))
    ylim = range(c(tpoints[,2]))
  }
  plot(xlim,ylim, col="white", type="s", xlim=xlim,ylim=ylim,xlab="u1",ylab="u2")
  if (all.iter) {
    points(ipoints,col="grey",lty=2)
    rsg.sol.outside(state, sol, rsg.sol)
  }
  points(tpoints, col="red", type="o", lwd=2)
  if (isTRUE(tight.approx)) tight.approx = tighter.notransfer.payoff.set(sol)
    
  if (is.list(tight.approx)) {
    pset = tight.approx[[state]]
    points(pset, col="orange", lwd=2, type="o")
  }
  if (!is.null(rsg.sol)) {
    lines(rsg.sol$points[[state]], col="blue", type="o",lty=2)
  }
  
  if (show.np.e) {
    points(np.e[[1]][state],np.e[[2]][state], col="red", lwd=2,cex=1.8)
  }
}

rsg.sol.outside = function(state,sol, rsg.sol, cat=TRUE) {
  restore.point("rsg.sol.outside")
  
  tpoints = with(as.data.frame(sol$sol.mat[state,,drop=FALSE]),rbind(
  c(v1,v2),c(v1,U-v1),c(U-v2,v2),c(v1,v2)
))
  ipoints = rsg.sol$ipoints[[state]]

  tol = 1e-12
  inside = 
  ipoints[,1] >= min(tpoints[,1])-tol &
  ipoints[,1] <= max(tpoints[,1])+tol &
  ipoints[,2] >= min(tpoints[,2])-tol &
  ipoints[,2] <= max(tpoints[,2])+tol &
  ipoints[,1]+ipoints[,2] <= max(tpoints[,1]+tpoints[,2])+tol 
  
  outside = !inside
  rev.outside = unique(rsg.sol$revolution[outside])
  
  #rows = rsg.sol$revolution==63 & outside
  #points(ipoints[rows,],col="blue")
  
  share.outside = sum(outside)/nrow(ipoints)
  share.rev.outside = length(rev.outside)/(max(rsg.sol$revolution)+1)
  if (cat) {
    cat("\nFraction of pivots outside equilibrium payoff set with transfers: ", round(share.outside,3))  
    cat("\nFraction of revolutions with pivots outside equilibrium payoff set with transfers: ", round(share.rev.outside,3))
    
  }
  
  invisible(nlist(share.outside,share.rev.outside))
}

dyngame.to.json.rsg = function(m, file) {
  rsg = dyngame.to.rsg(m)
  json = jsonlite::toJSON(rsg, pretty=TRUE)
  writeLines(json,file)
  rsg
}

dyngame.to.rsg = function(m) {

  states = lapply(seq.int(m$nx), function(x) {
    ai_lab = m$act.fun(m$xv.mat[x,])$lab
    numActions = sapply(ai_lab, length)
    
    # transpose action profile order
    tax = as.vector(t(matrix(which(m$ax[,"x"]==x), numActions[1], numActions[2])))
    
    list(
      numActions = numActions,
      payoffs = m$g[tax,,drop=FALSE],
      transition = m$tau[tax,,drop=FALSE]
    )
  })
  li = list(
    delta = m$delta,
    numPlayers = m$n,
    numStates = m$nx,
    states = states
  )
  return(li)
}


tighter.notransfer.payoff.set = function(sol) {
  restore.point("tighter.notransfer.payoff.set")
  
  m = sol$get.m()
  Umin = get.lowest.joint.payoff.no.transfer(m)
  u1 = get.highest.ui.notransfer(m,1)
  u2 = get.highest.ui.notransfer(m,2)
  ma = sol$sol.mat
  v1 = ma[,"v1"]
  v2 = ma[,"v2"]
  U = ma[,"U"]
  
  u1 = pmin(u1, U-v2)
  u2 = pmin(u2, U-v1)
  Umin = pmax(Umin, v1+v2)
  #data.frame(v1,v2,U,Umin,u1,u2)
  
  # loop through states
  pset.li = lapply(seq.int(NROW(ma)), function(state) {
    pset = matrix(ncol=2, byrow=TRUE,c(
      v1[state], v2[state],
      v1[state], u2[state],
      U[state]-u2[state], u2[state],
      u1[state], U[state]-u1[state],
      u1[state],v2[state],
      v1[state], v2[state]
    ))
  })  
  pset.li
}

#' Calculates the lowest joint payoff for each state
#' can be used to tighten bound on payoff set
#' for games without transfers
get.lowest.joint.payoff.no.transfer = function(m) {
  restore.point("get.lowest.joint.payoff.no.transfer")
  
  T = m$tau[m$admiss,,drop=FALSE]
  R = -m$G[m$admiss]
  
  ret = policy_iteration(T=T,R=R,delta=m$delta,na=m$admiss.ind)
  names(ret$p) = names(ret$V) =  m$x.lab
  return(ret$V)
}

#' Calculates the lowest joint payoff for each state
#' can be used to tighten bound on payoff set
#' for games without transfers
get.highest.ui.notransfer = function(m,i) {
  
  restore.point("get.highest.joint.payoff")
  
  T = m$tau[m$admiss,,drop=FALSE]
  R = m$g[m$admiss,i]
  
  ret = policy_iteration(T=T,R=R,delta=m$delta,na=m$admiss.ind)
  names(ret$p) = names(ret$V) =  m$x.lab
  return(ret$V)
}

# Stops time to solve a game
stop.game.solve.time = function(game, times=1, just.gk=FALSE, use.cache=TRUE, seed=NULL, run.id = 1:times) {
  library(dyngame)
  library(dplyr)
  library(RSGSolve)
  library(R.cache)
  restore.point("stop.game.solve.time")
  
  set.storing(FALSE)
  
  m = init.game(game)
  if (!is.null(seed))
    set.seed(seed)

  m.org = as.environment(as.list(m))
  res = lapply(run.id, function(run) {
    # solve with transfers
    cachedEval(add.expr.to.key = FALSE, overwrite=!use.cache,key=list(game=game, run=run, seed=seed, just.gk=just.gk),
    {
      start = proc.time()
      m = as.environment(as.list(m.org))
      log = capture.output(sol <- suppressMessages(solve.game(m,plots = FALSE, verbose = FALSE)))
      gk = proc.time()- start
      
      if (just.gk) {
        return(data_frame(delta=m$delta,gk=gk["elapsed"], gk.has.eq = sol$sol.exists,numStates = m$nx, numA=m$nax))
      }
      
      # solve without transfers
      rsg = dyngame.to.rsg(m)
      start = proc.time()
      rsg.sol = solveSG(rsg=rsg,noreturn=TRUE)
      abs = proc.time()- start
    
      data_frame(delta=m$delta,gk=gk["elapsed"],abs=abs["elapsed"], gk.has.eq = sol$sol.exists, abs.solved = rsg.sol$solved==1, numStates = m$nx, numA=m$nax)
    })
  })
  set.storing(TRUE)
  bind_rows(res)
  
}

cachedEval <- function(expr, key=NULL, envir=parent.frame(), overwrite=FALSE, add.expr.to.key=TRUE,character.key = TRUE,...) {
  expr <- substitute(expr);

  library(R.cache)
  if (add.expr.to.key)
    key <- c(list(expr=expr), key);

  if (character.key) {
    key = as.list(as.character(key))
  }

  # Look for cached results
  resList <- loadCache(key=key,...);
  restore.point("cachedEval")
  if (!overwrite && !is.null(resList)) {
    # Return the results of the memoized evaluation
    return(resList$result);
  }

  env <- new.env(parent=envir);
  res <- eval(expr, envir=env);
  resList <- list(results=res);
  saveCache(resList, key=key, ...);
  res;
} # cachedEval()


