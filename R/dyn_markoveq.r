#' Try to find a Markov Perfect Equilibrium
#' using a Cournot adjustment process
#'
#' Idea:
#' 1. Start with some policy.
#' 2. Calculate each player's continuation payoff u
#' 3. Get for player i a locally best reply holding fix the other players' strategies and assume continuation payoffs are u
#' 4. repeat steps 2 and 3 for different players
#' Variant if there are static actions
#' 0. Solve for a Nash equilibrium in every state
#' 1.-4
#' @export
markov.pm.algo = function(m,delta=m$delta,start.a = rep(1,m$nx), eps = 10^-8, sym2=TRUE) {
  
  restore.point("markov.pm.algo")
  
  #u = matrix(0,m$nx,m$n)
  a = start.a
  ax = a.to.ax(m,a)

  x.for.2  = xv.to.x(m,m$xv.mat[,c(2,1)])
  #ax.for.2 = a.to.ax(m,m$ax[,"a"],x = x.for.2[m$ax[,"x"]])
  
  counter = 0
  i = 1
  u.old = matrix(-Inf,m$nx,m$n)
  while (TRUE) {
	
		ax.tau = as.matrix(m$tau[ax,])
		
    # Payoffs for every player in every state corresponding to the actual policy ax
    u = solve(diag(m$nx)-delta*ax.tau, (1-delta)*m$g[ax,])
    
		plot(u[,1], main = counter,col="blue",xlab="x")
    points(u.old[x.for.2,2], pch=3,col="grey")
    points(u.old[,1],col="grey")
    points(u[x.for.2,2], pch=3,col="red")
    
    if (approxeq(u,u.old,tol=eps)) {
      break;
    }		
    # Get best-reply actions and values for every player
    ax = get.best.replies(m,i,delta,u[,i],ax)
    
    counter = counter +1
    i = (i %% m$n)+1
    u.old = u
    
    #rbind(m$x.lab,m$x.lab[x.for.2],u[,1],u[x.for.2,2])
  }
  
  return(list(ax=ax,u=u))
}    


get.best.replies = function(m, i, delta = 0,u,ax) {
  
  restore.point("get.best.replies")
  
	# Convert ax to ax_i
	ax_i = m$ind.ax.to.ax_i[[i]][ax]
	
	# Get replies and corresponding index for the rows ax_i
	ind.all = m$ind.replies.by.ax_i[[i]]
	# ax indices of the replies
	replies = m$replies[[i]][rows.to.v.ind(ind.all,ax_i)]
	ind.replies = make.child.VectorListInd(ind.all,ax_i)
	
	# Get payoffs for each reply
  g.replies = m$g[replies,i]
	replies.tau = as.matrix(m$tau[replies,])
	u.replies = replies.tau %*% u 
  val = (1-delta)*g.replies + delta*u.replies
	
	# Select best reply for each ax_i
	ax.opt = replies[which.RowMaxs(ind.replies,val)]
	
  return(ax.opt)
}    



markov.pm.algo.sym = function(m,delta=m$delta,start.a = rep(1,m$nx), eps = 10^-8, sym2=TRUE) {
  
  restore.point("markov.pm.algo")
  
  #u = matrix(0,m$nx,m$n)
  a = start.a
  ax = a.to.ax(m,a)

  x.for.2  = xv.to.x(m,m$xv.mat[,c(2,1)])
  #ax.for.2 = a.to.ax(m,m$ax[,"a"],x = x.for.2[m$ax[,"x"]])
  
  counter = 0
  i = 1
  u.old = matrix(-Inf,m$nx,m$n)
  while (TRUE) {
	
		ax.tau = as.matrix(m$tau[ax,])
		
    # Payoffs for every player in every state corresponding to the actual policy ax
    u = solve(diag(m$nx)-delta*ax.tau, (1-delta)*m$g[ax,])
    
		plot(u[,1], main = counter,col="blue",xlab="x")
    points(u.old[x.for.2,2], pch=3,col="grey")
    points(u.old[,1],col="grey")
    points(u[x.for.2,2], pch=3,col="red")
    
    if (approxeq(u,u.old,tol=eps)) {
      break;
    }
    
		
		
    # Get best-reply actions and values for every player
    ax = get.best.replies(m,i,delta,u[,i],ax)
    
    counter = counter +1
    i = (i %% m$n)+1
    u.old = u
    
    #rbind(m$x.lab,m$x.lab[x.for.2],u[,1],u[x.for.2,2])
  }
  
  return(list(ax=ax,u=u))
}