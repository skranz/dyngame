# Init a dynamic game


# payoffs g and the transition matrix tau have a 3 dimensional structure
# I do not store them in a 3d array however, since R does not have such nice and quick support for arrays than for matrices.

# Therefore, I store these values in a semi flat matrix:
# g:   n   cols for players, a, x
# tau: xn  cols for destination prob, a, x

FUN_NAMES = c("stat.act.fun","stat.g.fun","tau.fun","act.fun","g.fun",
              "stat.x.group", "x.group","x.group", "x.stage.fun","extra.sol.fun")

new.dyngame = function(...) {
	m = ListEnv(...)	  
  class(m)=c("dyngame","ListEnv","environment")
  m
}

new.dynsol = function(m) {
	ms = new.ListEnv(parent=m)	  
  class(ms)=c("dynsol","dyngame","ListEnv","environment")
  ms$get.m = function() m
  ms$sol.exists = TRUE
  ms
}


# There are three different ways how payoffs can be initialized
init.game = function(my.game=NULL,n=2, name="",
                     xv.val=NULL, nxv=NULL, xv.names = names(xv.val), xv.lab = xv.val, ixv = NULL,
                     symmetric = FALSE, a.lab.sep = "|", x.lab.sep = ";", tol = 10e-8, tau.type="endo",tau.sparse=FALSE,para=list()) {
 	
	# Copy all variables of my.game into local environment
	# Have to do this before 
	# In order to keep assignment of parameter default values
	copy.into.env(source=my.game, exclude="my.game")
	print(xv.val)
	print(xv.lab)
	
	
	
	restore.point("init.game")

  m = new.dyngame()
  
	m$no.labels = NO.LABELS
	
  	# Transform values
  if (!is.list(xv.val)) xv.val=list(xv.val)
  if (!is.list(xv.lab)) xv.lab=list(xv.lab)

  #Copy all local variables into m
	copy.into.env(dest=m, exclude="my.game")
	
  # For symmetric games certain operations only
  # have to be conducted for player 1
  
	m$para = para
    
  m$nxv = length(xv.val)
  m$xv.dim = sapply(xv.val,NROW)
  m$nx = prod(m$xv.dim)
  
  m$xv.mat = make.grid.matrix(x=xv.val)
	if (is.null(m$xv.names))
	  m$xv.names = paste("xv",1:m$nxv,sep="")
		
	if (!NO.LABELS) {	
		colnames(m$xv.mat) = m$xv.names
	}
	
  gm = make.grid.matrix(m$xv.lab)
  m$x.lab = rep(gm[,1])
  if (m$nxv > 1) {
    for (z in 2:m$nxv)
      m$x.lab = paste(m$x.lab,gm[,z],sep=x.lab.sep)
  }
	
	if (!NO.LABELS) {
		rownames(m$xv.mat) = m$x.lab
		colnames(m$xv.mat) = names(m$xv.val) = names(m$xv.lab) = m$xv.names
  }

	if (m$symmetric) {
    m$n.sym = 1
		m$ixv = ixv
		if (is.null(m$ixv)) {
			if (m$nxv==m$n) {
				mywarning("You should specify ixv if symmetric. I use the default ixv=1:m$n")
				m$ixv = 1:m$n
			} else {
				stop("You have to specify ixv or set symmetric=FALSE")
			}
		}
  } else {
    m$n.sym = m$n
  }
  
	m$tau.sparse = tau.sparse
	m$tau.type = tau.type
	
	if (tau.sparse) {
		m$tau.sub = function(tau.sub) as.matrix(tau.sub)
	} else {
		m$tau.sub = function(tau.sub) tau.sub
	}
	


	if (tau.type!="endo") {	
		stop('works in the moment only for tau.type = "endo"')
	}

	init.x.groups(m)
	init.ax(m)
	init.replies(m)
	set.g.and.tau(m)
  make.ax.labels(m,"dyn")
	adapt.for.symmetry(m) 	
	
	if (FALSE) {
	# State transitions are endogenous and all actions can change transition probabilities
	# if (tau.type!="endo") {
		# init.action.group(m) 
	# # Action profiles can be split into two dimensions:
	# # actions that change transition probabilities, e.g. capacity investments
	# # actions that don't affect state transitions, e.g. prices 
	# } else if (tau.type =="dynstat") {
		
		# init.action.group(m) # Also inits transition probabilities
		# init.static.action.group(m)
  
	# # State transitions are completely independent from actions		
	# } else if (tau.type =="exo") {
		# init.exo(m)
	
	# # Many small relationships in a large economy. Players consider state
	# # transitions as exogenous in their decisions, but tau depends on
	# # aggregated actions
	# } else if (tau.type =="macro") {
		# init.exo(m)	
	# }
	}

  # In repeated multi-stage games we need to map every state to a stage
  # in order to adapt payoffs such that we won't have discounting between stages
  m$x.stage = NULL
  if (!is.empty(m,"x.stage.fun")) {
    m$x.stage = m$x.stage.fun(m$xv.mat,m)
		if (!NO.LABELS)
			names(m$x.stage) = m$x.lab
    m$n.stage = max(m$x.stage)
  } else {
    m$n.stage = 1
  }
  m$ax = ind.row.matrix(m$ind.ax.by.x,colnames=c("ax","x"))

	m$i.max = m$n
  # Some consistency checks
  check.m(m)
  return(m)
}

make.ax.labels = function(ag,type) {
  
 	restore.point("make.ag.labels")
  m = ag
  if (NO.LABELS) return();
  if (type=="dyn") { 
    for (ig in 1:ag$nxg) {
      xg = ag$xg[[ig]]
      
      gm = make.grid.matrix(xg$av.lab)
      xg$a.lab = paste.matrix.cols(gm)
      
      if (is.empty(xg,"ax")) {
        ax.rows = as.numeric(rows.to.v.ind(ag$ind.ax.by.x,xg$x))    
        ax = cbind(rep(xg$x,each=xg$na),rep(1:xg$na,times=xg$nx))
        colnames(ax) = c("x","a")
      } else {
        ax.rows = xg$ax.rows
        ax = xg$ax
      }
          
      # Labels
      ag$ax.lab[ax.rows] = paste(m$x.lab[ax[,"x"]],xg$a.lab[ax[,"a"]],sep=" ")
			ag$ax.lab.short[ax.rows] = xg$a.lab[ax[,"a"]]

    }
    rownames(m$tau) = rownames(ag$g) = names(ag$G) = ag$ax.lab
    colnames(m$tau) = m$x.lab
  } else {
    # Labels for static action profiles
  }
}


# Initializes or Resets Payoffs and Transmission Matrix
# This can be done without changing the rest of the model structure
set.g.and.tau = function(ag,g.fun=ag$g.fun,tau.fun=ag$tau.fun,
                        recalc.g=!is.null(g.fun),recalc.tau = !is.null(tau.fun),para=NULL) {
 	restore.point("set.g.and.tau")

	
	if (!recalc.g & (!recalc.tau)) {
		return
	}
  m = ag
  if (!is.null(g.fun)) {
    m$g.fun = g.fun
  }
  if (!is.null(tau.fun)) {
    m$tau.fun = tau.fun
  }

	# For what ever reason it goes MUCH quicker if we locally define tau and g
  # instead of writing into m$tau and m$g	
	if (ag$tau.sparse) {
		tau = Matrix(0,m$nax,m$nx,sparse=TRUE)
	} else {
		tau = matrix(0,m$nax,m$nx)
	}
	g = matrix(0,m$nax,m$n)
	
	# Do it to have better profiling
	call.g.fun = function(avm,xvm,m) {
		m$g.fun(avm,xvm,m)
	}	
	# Do it to have better profiling
	call.tau.fun = function(avm,xvm,m) {
		m$tau.fun(avm,xvm,m)
	}	
	
	# Copy new parameters into environments of g and tau
	if (!is.null(para)) {
	 warning("Adapt environment, but do not perform initial calculations in make.my.game")
	 env = environment(m$tau.fun)
	 copy.into.env(source=para, dest=env)
	 
	 envg = environment(m$g.fun)
	 if (!identical(env,envg))
	 copy.into.env(source=para, dest=envg)
	}
		
  for (ig in 1:ag$nxg) {
    xg = ag$xg[[ig]]    
		if (is.empty(xg,"ax")) {
			ax.rows = as.numeric(rows.to.v.ind(ag$ind.ax.by.x,xg$x))    
			ax = cbind(rep(xg$x,each=xg$na),rep(1:xg$na,times=xg$nx))
			colnames(ax) = c("x","a")
		} else {
			ax.rows = xg$ax.rows
			ax = xg$ax
		}
   
    # Payoffs
    av.mat = make.grid.matrix(xg$av.val)
    avm = av.mat[ax[,"a"],,drop=FALSE]
    xvm = m$xv.mat[ax[,"x"],,drop=FALSE]
  
    if (recalc.g) {
			if (ag$nxg==1) {
				g = call.g.fun(avm,xvm,m)
			} else  {
				g[ax.rows,] = call.g.fun(avm,xvm,m)
			}
      
    }
    # Transition Probabilities
    # Go through all dynamic state groups 
    if (recalc.tau & m$tau.type!="exo") {
			if (ag$nxg==1) {
				tau = call.tau.fun(avm,xvm,m)
			} else  {
				X = call.tau.fun(avm,xvm,m)
				tau[ax.rows,] <- X
			}
	  } 
	}
	
	if (recalc.tau & m$tau.type=="exo") {
  	m$tau = m$tau.fun(avm=NULL,m$xv.mat,m)
	}
  if (recalc.tau & m$tau.type=="macro") {
		m$get.tau = function(m,ae) {
			m$tau[ae,]
		}
	}
	

	if (recalc.g) {
		ag$g = g
		ag$G = rowSums(ag$g)
	}
	
	if (recalc.tau) {
		ag$tau = tau
	}
}


# So far only store permutations of states
# This will allow to solve for punishment values only once
# We should also store permutations of ax. This allows to solve fewer action profiles on the equilibrium path
adapt.for.symmetry = function(m) {
 	restore.point("adapt.for.symmetry")
	
	# OLD STUFF BELOW
	if (!m$symmetric) {
		m$ax.sym = 1:m$nax
		m$x.sym = 1:m$nx
		m$nx.sym = m$nx
		m$nax.sym = m$nax
		return;
	}
		
	
	# Generate permutation matrix for states x
	x.perm = matrix(0,m$nx,m$n)
	x.perm[,1] = 1:m$nx
	# We just store the permutation in which we simply swap player 1 and player i
	# Punishment payoffs of player i in the state of row r and column i will be the same
	# as those of player 1 in the state r 		
	for (i in 2:m$n) {
		cols = 1:m$n
		cols[1]=i; cols[i]=1;
		xv.mat = m$xv.mat[,cols]
		x.perm[,i] = xv.to.x(m,xv.mat)
	}
	m$x.perm = x.perm

	m$x.sym = which(x.perm[,1]<=rowMins(x.perm[,-1,drop=FALSE]))
	m$nx.sym = length(m$x.sym)
	
	# Now create a permutation matrix for ax
	ax.perm = matrix(0,m$nax,m$n)
	ax.perm[,1] = 1:m$nax
	
  for (ig in 1:m$nxg) {
		xg = m$xg[[ig]]
		x = xg$x[1]
		
		iav = xg$iav
		stopifnot(length(iav)==m$n) # Can't yet handle multidimensional actions		
		xg$av.val
		
		ax1 = unlist(lapply(xg$x,x.to.ax,m=m))
		for (i in 2:m$n) {
			perm.x  = m$x.perm[x,i]
			perm.xg = m$xg[[m$ind.xg.of.x[perm.x]]]
						
			iav.perm = perm.xg$iav
			stopifnot(length(iav)==m$n) # Can't yet handle multidimensional actions
			
			iav.perm[1] = i
			iav.perm[i] = 1
			
			perm.rows = grid.matrix.permutation(xg$av.val,iav.perm)
			
			# These are the original ax values of the permuted states
			axi = unlist(lapply(x.perm[xg$x,i],x.to.ax,m=m))

			# Add the indices given by perm.rows and remove the original ordering 1:na
			axi = axi + rep(perm.rows,length.out = NROW(axi)) - rep(1:NROW(perm.rows),length.out = NROW(axi))
			
			ax.perm[ax1,i] = axi
		}
	}
	m$ax.perm = ax.perm
	
	m$ax.sym = which(1:m$nax<=rowMins(m$ax.perm[,-1,drop=FALSE]))
	m$nax.sym = length(m$ax.sym)
		
	#matrix(m$ax.lab[as.numeric(ax.perm)],NROW(ax.perm),NCOL(ax.perm))
}

init.x.groups = function(ag) {
	m = ag
	x.group=ag$x.group
	act.fun=ag$act.fun
	
	# Groups of states that may have different action profiles
  if (is.function(x.group)) {
    x.group = x.group(m$xv.mat,m)
  } else if (length(x.group)==1) {
    x.group = rep(1,m$nx)
  }
	m$x.group = x.group
  m$nxg = length(unique(x.group))

  if (!isTRUE(all.equal(1:ag$nxg,sort(unique(x.group))))) {
    stop("x.group must be of form 1,...,nxg (Starting with 1)")
  }
  ag$xg = list()
  ag$ind.xg.of.x = x.group
  ag$na = rep(NA,m$nx)
  
  for (ig in 1:ag$nxg) {
    x = which(x.group==ig)[1]
    # Activities corresponding to the actual group of states
    xg = new.ListEnv()
    
    ret = act.fun(m$xv.mat[x,],m)
    
    xg$av.val = ret$val
    xg$av.gi = GridInd(xg$av.val)
    xg$nav = length(xg$av.val)
    xg$av.dim = sapply(xg$av.val,NROW)
    xg$na = prod(xg$av.dim)
    xg$x = which(x.group==ig)
    xg$nx = length(xg$x)
    if (is.null(ret$i)) {
      stop("Error: Return value of act.fun does not include the argument i that maps players to action vectors. Please include this value in your function that specifies the game.")
    }
    
    xg$iav = ret$i
    
    if (is.null(ret$lab)){
      xg$av.lab = ret$val
    } else {
      xg$av.lab = ret$lab
    }
    
    
    ag$xg[[ig]] = xg
    
    ag$na[which(x.group == ig)] = xg$na
  }
}

init.ax = function(m) {
  
	restore.point("init.ax")
	
	if (is.null(m$act.fun)) {
		print("No act.fun exists. No action can influence state transitions")
		return(NULL)
	}
	ag = m
	
	act.fun = m$act.fun
	x.group = m$x.group

	init.x.groups(m)
  # Go through all state groups
  # and init actions, labels, payoffs and replies
  # for all states in the group

  ag$ind.ax.by.x = VectorListInd(ncols=ag$na)
  ag$start.ax.by.x = c(1,cumsum(ag$na)+1)#[-(m$nx+1)]
	
	
  ag$nax = length(ag$ind.ax.by.x)
  ag$ax.lab = rep(NA,ag$nax)
	ag$ax.lab.short = rep(NA,ag$nax)
  ag$g = matrix(NA,ag$nax,ag$n)
  
  
  for (ig in 1:ag$nxg) {
    xg = ag$xg[[ig]]
    ax.rows = as.numeric(rows.to.v.ind(ag$ind.ax.by.x,xg$x))    
    ax = cbind(rep(xg$x,each=xg$na),rep(1:xg$na,times=xg$nx))
    colnames(ax) = c("x","a")
    xg$ax.rows = ax.rows
    xg$ax = ax
	}

} 


init.replies = function(m) {
  
	restore.point("init.replies")
	
	if (is.null(m$act.fun)) {
		print("No act.fun exists. No action can influence state transitions")
		return(NULL)
	}
	ag = m
	#m$dyn.ag = ag  # Commented out to avoid recursion, because 
  
  ag$replies = list()
  ag$ind.ax.to.ax_i = list()
  ag$ind.ax_i.by.x = list()
  ag$ind.replies.by.ax_i = list()
  ag$num.replies.by.x = list()
  ag$num.ax_i.by.x = list()

  for (i in 1:m$n) {
    #ag$replies[[i]] = VectorList(nrow=m$nx)
    ag$replies[[i]] = list()
    ag$ind.ax.to.ax_i[[i]] = list()
    ag$ind.ax_i.by.x[[i]] = rep(NA,m$nx)
    ag$ind.replies.by.ax_i[[i]] = list()
    ag$num.replies.by.x[[i]] = list()
    ag$num.ax_i.by.x[[i]] = rep(NA,m$nx)
  }
  
  for (ig in 1:ag$nxg) {
    xg = ag$xg[[ig]]    
    # Labels
    # ag$ax.lab[ax.rows] = paste(m$x.lab[ax[,"x"]],xg$a.lab[ax[,"a"]],sep=" ")
    x.offset = c(0,cumsum(ag$na)[-length(ag$na)])  
    # Generate ax_i indices of player i's replies
    for (i in 1:m$n) {
      repl = make.replies.gi(xg$av.gi,k.i = which(xg$iav == i))
			repl.vec = as.vector(t(repl$reply.mat))
			nx = NROW(xg$x)
			num.replies = rep(NCOL(repl$reply.mat),NROW(repl$reply.mat))
			ag$num.replies.by.x[[i]][xg$x] = replicate(nx,num.replies,simplify=FALSE)
			ag$ind.ax.to.ax_i[[i]][xg$x] = replicate(nx,repl$cp.ind,simplify=FALSE)
			ag$ind.ax_i.by.x[[i]][xg$x] = NROW(repl$reply.mat)	
			ag$replies[[i]][xg$x] = lapply(xg$x, function(x) { repl.vec + x.offset[x] } )
    }
  }
    
  for (i in 1:m$n) {
    ax_i.offset = c(0,cumsum(ag$ind.ax_i.by.x[[i]][-m$nx]))
		ag$ind.ax.to.ax_i[[i]] = lapply(1:m$nx, 
		              function(x) { ag$ind.ax.to.ax_i[[i]][[x]] + ax_i.offset[x] } )
    ag$ind.ax.to.ax_i[[i]] = unlist(ag$ind.ax.to.ax_i[[i]]) # m$nax rows 
    ag$ind.ax_i.by.x[[i]] = VectorListInd(ag$ind.ax_i.by.x[[i]])  
    ag$replies[[i]] = do.call("c",ag$replies[[i]]) # Store replies in one big vector
    #ncols = sapply(ag$replies[[i]],NROW)
    ag$ind.replies.by.ax_i[[i]] = VectorListInd(do.call("c",ag$num.replies.by.x[[i]]))
  }
  ag$num.replies.by.x = NULL
  	
  #adapt.for.symmetry(ag)
} 


init.static.action.group = function(m) {
  

	if (is.empty(m,"stat.act.fun")) {
		print("No stat.act.fun exists. Purely dynamic game")
		return(NULL)
	}
	ag = new.ListEnv(parent=m)
	m$stat.ag = ag
	act.fun = m$stat.act.fun
	g.fun = m$stat.g.fun
	ag$g.fun = g.fun

	x.group = m$stat.x.group
	ag$x.group = x.group
	ag$act.fun = act.fun
      
  class(ag)=c("ActionGroup","dyngame","ListEnv","environment")  
	init.x.groups(ag,x.group)

	# Go through all state groups
  # and init replies then go through all states
  # in this group and solve the corresponding repeated game
  ag$opt = list()
  for (ig in 1:ag$nxg) {
    xg = ag$xg[[ig]]
    xg$replies = list()
    for (i in 1:m$n) {
      xg$replies[[i]] = make.replies.gi(xg$av.gi,k.i = which(xg$iav == i))
    }
    # Go through all states
    for (x in xg$x) {
      init.static.frontiers(ag,xg,x)
    }    
    # If we have a memory problem, we could now remove xg
    #ag$xg[[i]] = NULL
  }

	make.ag.labels(ag,type)  
}

init.static.frontiers = function(ag,xg,x) {
  
	restore.point("init.static.frontiers")

  # Calculate static payoffs and cheating payoffs for state x
  # and the static action group
  xv = t(ag$xv.mat[x,])
  
  av.mat = make.grid.matrix(xg$av.val)
  g = ag$g.fun(av.mat,xv)
  G = rowSums(g)
  c.mat = matrix(NA,xg$na,ag$n)  
  for (i in 1:ag$n) {
    gi.rep = xg$replies[[i]]
    g.replies = matrix(g[gi.rep$reply.mat,i],nrow=NROW(gi.rep$reply.mat))
    g.max = rowMaxs(g.replies)
    c.mat[,i] = g.max[gi.rep$cp.ind]
  }  
  L = rowSums(c.mat) - G

  ag$opt[[x]] = list()
  
	# Equilibrium state
	#opt.ind = rev(sk.pareto.frontier(G,-L))
	opt.ind = sk.pareto.frontier(-L,G)
	opt.mat = cbind(L[opt.ind],G[opt.ind],opt.ind)
	colnames(opt.mat) = c("L","G","a")
	#rownames(opt.mat) = xg$a.lab[opt.ind]
	opt.fun = approxfun(x=opt.mat[,1], y = opt.mat[,2], method="constant",
          yleft=-Inf, yright=opt.mat[NROW(opt.mat),2], f = 1, ties = "ordered")

  ag$opt[[x]]$mat.e = opt.mat
  ag$opt[[x]]$fun.e = opt.fun
  
  #plot(opt.mat[,"L"],opt.mat[,"G"],ylab="G",xlab="L",main=paste(ag$x.lab[x], " static G(L)"),col="blue")
  
  ag$opt[[x]]$mat.i = list()
  ag$opt[[x]]$fun.i = list()

	# Punishment states
	for (i in 1:ag$n) {
  	#opt.ind = rev(sk.pareto.frontier(-c.mat[,i],-L))
  	opt.ind = sk.pareto.frontier(-L,-c.mat[,i])

  	opt.mat = cbind(L[opt.ind],c.mat[opt.ind,i],opt.ind)
  	colnames(opt.mat) = c("L","ci","a")
  	#rownames(opt.mat) = xg$a.lab[opt.ind]
  	opt.fun = approxfun(x=opt.mat[,1], y = opt.mat[,2], method="constant",
            yleft=Inf, yright=opt.mat[NROW(opt.mat),2], f = 1, ties = "ordered")
  
    ag$opt[[x]]$mat.i[[i]] = opt.mat
    ag$opt[[x]]$fun.i[[i]] = opt.fun
  }
  plot(ag$opt[[x]]$mat.e[,"L"],ag$opt[[x]]$mat.e[,"G"],ylab="G",xlab="L",main=paste(ag$x.lab[x], " static G(L)"),col="blue",
       xlim= range(c(ag$opt[[x]]$mat.e[,"L"],ag$opt[[x]]$mat.i[[ag$n]][,"L"])),
       ylim= range(c(ag$opt[[x]]$mat.e[,"G"],ag$opt[[x]]$mat.i[[ag$n]][,"ci"])))
  points(opt.mat[,"L"],opt.mat[,"ci"],col="red",pch=3)
}  


# Assume imperfect monitoring of static signals
set.beta.I = function(m,beta.I) {

  if (is.empty(m,"beta.I")) {
        old.beta.I = 1
  } else { 
    old.beta.I = m$beta.I
  }
  ag = m$stat.ag
  m$beta.I = beta.I
  for (x in 1:m$nx) {
    ag$opt[[x]]$mat.e[,"L"] = ag$opt[[x]]$mat.e[,"L"] *
                            (old.beta.I / m$beta.I)
    opt.mat = ag$opt[[x]]$mat.e
    ag$opt[[x]]$fun.e = approxfun(x=opt.mat[,1], y = opt.mat[,2], method="constant",
          yleft=-Inf, yright=opt.mat[NROW(opt.mat),2], f = 1, ties = "ordered")
    
    for (i in 1:m$n) {
      ag$opt[[x]]$mat.i[[i]][,"L"] = ag$opt[[x]]$mat.i[[i]][,"L"] *
                        (old.beta.I / m$beta.I)
      opt.mat = ag$opt[[x]]$mat.i[[i]]
      ag$opt[[x]]$fun.i[[i]] = approxfun(x=opt.mat[,1], y = opt.mat[,2], method="constant",
            yleft=Inf, yright=opt.mat[NROW(opt.mat),2], f = 1, ties = "ordered")
    }
  }
}                

# Some checks whether model is consistent
check.m = function(m, stop.on.error = TRUE) {
  
  restore.point("check.m")
  
  if (!all.eq(rowSums(m$tau),rep(1,m$nax))) {
    warning("check.m: m$tau is not correct. Not all rows add to 1! Check tau.fun!")
  }
	
  if (any(m$tau<0)) {
    warning("check.m: m$tau is not correct. Some transition probabilities are negative! Check tau.fun!")
  }
  if (!is.null(m$x.stage)) {
    if (length(m$x.stage) != m$nx) {
      warning("check.m: m$x.stage does not have length nx. Check x.stage.fun!")
    }
    for (stage in 1:m$n.stage) {
      rows = which(m$x.stage == stage)
      if (NROW(rows) == 0) {
        warning(paste("check.m: stage",stage,"not specified. Check x.stage.fun!"))
      } else {
        rows = m$ind.ax.by.x[rows]
        
        next.stage = stage +1 
        if (stage == m$n.stage) next.stage = 1
        next.states = which(colSums(m$tau[rows,,drop=FALSE])>0)
        if (any(m$x.stage[next.states] != next.stage))
          warning(paste("check.m: Wrong state transitions. Not all states in stage",
                        stage, "transit for sure to states in stage",next.stage,
                        "Check tau.fun or x.stage.fun!"))
      }                 
    }
  }
}
    