# Quickly solving dynamic games with perfect monitoring in which state transitions are completely exogenous
# We can also (partially) solve macro models 

# Idea: Calculate for every state x the highest joint payoffs U and punishment payoffs v for every level of liquidity
#  
                
init.exo = function(m) {
  
  m = m
	restore.point("init.static.action.group")
 
  # Go through all state groups
  # and init replies then go through all states
  # in this group and solve the corresponding repeated game
  init.x.groups(m)
	init.action.group(m)
	
	m$opt = list()
	m$opt.n = matrix(0,m$nx,m$n+1)

  for (ig in 1:m$nxg) {
    xg = m$xg[[ig]]
    xg$replies = list()
    for (i in 1:m$n) {
      xg$replies[[i]] = make.replies.gi(xg$av.gi,k.i = which(xg$iav == i))
    }
    # Go through all states
    for (x in xg$x) {
      init.exo.x(m,xg,x)
    }    
    # If we have a memory problem, we could now remove xg
    #m$xg[[i]] = NULL
  }
}

init.exo.x = function(m,xg,x) {
	restore.point("init.exo.x")

  # Calculate static payoffs and cheating payoffs for state x
  # and the static action group
  xv = t(m$xv.mat[x,])
  
  av.mat = make.grid.matrix(xg$av.val)
  
	# Stage game payoffs
	g = m$g.fun(av.mat,xv)
  G = rowSums(g)

	# Stage game cheating payoffs
	c.mat = matrix(NA,xg$na,m$n)  
  for (i in 1:m$n) {
    gi.rep = xg$replies[[i]]
    g.replies = matrix(g[gi.rep$reply.mat,i],nrow=NROW(gi.rep$reply.mat))
    g.max = rowMaxs(g.replies)
    c.mat[,i] = g.max[gi.rep$cp.ind]
  }
	
	# Liquidity for every action profile that can be chosen in the stage game
  L = rowSums(c.mat) - G

  m$opt[[x]] = list()
  
	# Equilibrium phase: Compute optimal profiles for every L
	opt.ind = sk.pareto.frontier(-L,G)
	opt.mat = cbind(L[opt.ind],G[opt.ind],opt.ind)
	colnames(opt.mat) = c("L","G","a")
	#rownames(opt.mat) = xg$a.lab[opt.ind]
  m$opt[[x]]$mat.e = opt.mat
  m$opt.n[x,1] = NROW(opt.mat)  #plot(opt.mat[,"L"],opt.mat[,"G"],ylab="G",xlab="L",main=paste(m$x.lab[x], " static G(L)"),col="blue")
  
	
	# Calculate optimal punishment payoffs for every player i and every L
	
  m$opt[[x]]$mat.i = list()
	for (i in 1:m$n) {
  	opt.ind = sk.pareto.frontier(-L,-c.mat[,i])

  	opt.mat = cbind(L[opt.ind],c.mat[opt.ind,i],opt.ind)
  	colnames(opt.mat) = c("L","ci","a")
  	#rownames(opt.mat) = xg$a.lab[opt.ind]
  
    m$opt[[x]]$mat.i[[i]] = opt.mat
		m$opt.n[x,i+1] = NROW(opt.mat)  
	}
	
  plot(m$opt[[x]]$mat.e[,"L"],m$opt[[x]]$mat.e[,"G"],ylab="U",xlab="L",main=paste(m$x.lab[x], " "),col="blue",
       xlim= range(c(m$opt[[x]]$mat.e[,"L"],m$opt[[x]]$mat.i[[m$n]][,"L"])),
       ylim= range(c(m$opt[[x]]$mat.e[,"G"],m$opt[[x]]$mat.i[[m$n]][,"ci"])))
  points(opt.mat[,"L"],opt.mat[,"ci"],col="red",pch=3)
}  

# Get G,L,a for every state and phase for given positions in the list of
# optimal action profiles
get.from.pos = function(m,pos.k) {
	G = C = ae = Le = rep(NA,nx)
	ai = Li = matrix(NA,n,nx)
	
	pos.e = pos.k[,1]
	pos.i = pos.k[,-1]
	for (x in 1:nx) {	
		
		G[x] <- m$opt[[x]]$mat.e[pos.e[x],"G"]			
		ae[x] <- m$opt[[x]]$mat.e[pos.e[x],"ae"]			
		Le[x] <- m$opt[[x]]$mat.e[pos.e[x],"L"]	

		C_temp = 0		
		for (i in 1:n) {
			C_temp = C_temp+m$opt[[x]]$mat.i[[i]][pos.i[x,i],"ci"]
			ai[x,i] = m$opt[[x]]$mat.i[[i]][pos.i[x,i],"ai"]
			Li[x,i] = m$opt[[x]]$mat.i[[i]][pos.i[x,i],"L"]		
		}
		C[x] = C_temp
		
	}
	return(list(G=G,C=C,ae=ae,Le=Le,ai=ai,Li=Li,ak=cbind(ae,ai),Lk=cbind(Le,Li)))
	U
}


solve.exo.dyngame = function(m, delta = m$delta) {
	nx = m$nx; n= m$n;
	
	pos.k = matrix(1,nx,n+1)
	success=FALSE
	if (m$tau.type=="exo")
		tau = m$tau
		
	while (TRUE) {
		p = get.from.pos(m,pos.k)
		
		if (m$tau.type=="macro")
			tau = get.tau(m,act$ae)
		
		# Find states and phases for which ak(x) can be implemented
		LHS = (diag(nx)-delta*tau)%*%p$Lk
		RHS = delta*(p$G-p$C)
		RHS = matrix(RHS,ncol=NCOL(LHS))
		notok = LHS > RHS
		
		# All optimal action profiles can be implemented
		if (sum(notok)==0) {
			success=TRUE
			break
		}
		
		# There is a state and phase in which an action profile can't be 
		# implemented
		pos.k = pos.k + notok
		if (any(pos.k<m$opt.n)) {
			success = FALSE;
			break
		}
	}

	if (success) {
		return(p)
	} else {
		print("No equilibrium found")
		return(NULL)
	}
}
