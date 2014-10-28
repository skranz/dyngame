# Tools for analysing multiple models

get.sol.from.model.list = function(mlist) {
  sol = sublist(mlist,"sol.mat")
  sol = rbind.list(sol)
}

# Cannot yet change delta!
solve.multi.model = function(m,delta=m$delta,par = NULL,g.fun=m$g.fun,tau.fun=m$tau.fun) {
 
  restore.point("solve.multi.model")
  type=""
  if (length(delta)>1) {
    nm = length(delta)
    par.seq = delta
    par.name=type="delta"
  }
  if (!is.null(par)) {
    if (type=="delta" | length(par)>1) {
      stop("So far we can only change one parameter.")
    }
    par.name = names(par)[1]
    par.seq = par[[1]]
    type = "par"
  }
  
  nm = length(par.seq)
  ms = list()
  # Parameter that changes payoffs
  if (type!="delta") {
    for (b in 1:nm) {
      par = par.seq[b]
      # Assign the new parameter value to the environment
      # that is used by g.fun
      if (!is.null(g.fun))
        assign(par.name, par, envir = environment(g.fun))
      if (!is.null(tau.fun))
        assign(par.name, par, envir = environment(tau.fun))

      # Update stage game payoffs
      set.g.and.tau(m,g.fun=g.fun,tau.fun=tau.fun) 
      m$name = paste(par.name,round(par,3))
      #Solve game
      ms[[b]] = solve.game(m,delta=delta)
    }   
  } # Not yet implemented for changes in delta
    
  
  # Get solution matrix
  sol = sublist(ms,"sol.mat")
  sol = rbind.list(sol)
  sol = cbind(rep(1:nm,each=m$nx),rep(par.seq,each=m$nx), rep(1:m$nx,times=nm),sol)
  colnames(sol)[1:3] = c("model","par","x")
  mm = ListEnv(ms=ms,sol=sol,par.name = par.name, par.seq = par.seq, delta = delta)
  class(mm)=c("multimodel","ListEnv")
  return(mm)
}

update.multimod.sol = function(mm) {
  nx = mm$ms[[1]]$nx
  nm = NROW(mm$par.seq)
  # Get solution matrix
  sol = sublist(mm$ms,"sol.mat")
  sol = rbind.list(sol)
  sol = cbind(rep(1:nm,each=nx),rep(mm$par.seq,each=nx), 
              rep(1:nx,times=nm),sol)
  colnames(sol)[1:3] = c("model","par","x")
  mm$sol = sol
}  


get.sol.from.model.list = function(mlist) {
  sol = sublist(mlist,"sol.mat")
  sol = rbind.list(sol)
}

add.transitions = function(m) {
  
  restore.point("add.transitions")
  m$tau.e = m$tau[m$sol.mat[,"ae"],]
  
  vec = as.vector(t(m$tau.e))
  gi = GridInd(dim=dim(m$tau.e))
  rows = which(vec != 0)
  
  if (length(rows>0)) {
    mat = cbind(v.ind.to.x.ind(gi,rows),vec[rows])
    mat = mat[mat[,1] != mat[,2],,drop = FALSE]
    colnames(mat) = c("x.from","x.to","prob")
    rownames(mat) = paste(m$x.lab[mat[,1]],m$x.lab[mat[,2]])
    
    m$trans.e = mat
  } else {
    m$trans.e = matrix(NA,0,3)
    colnames(m$trans.e) = c("x.from","x.to","prob")
  }
}
  
get.x.long.run = function(m) {
  tau.e = m$tau[m$sol.mat[,"ax.e"],]
  prob = solve(diag(m$nx) - tau.e)
  return(prob)
}
