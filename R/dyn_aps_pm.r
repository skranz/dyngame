
# Computes for every ax the joint payoff that can be decomposed
# Uo is an nx*1 vector that describes an outer approximation for the
# maximum equilibrium values of U(x)
aps.pm.decompose.U = function(m,Uo=m$Uo) {
	Ud = (1-m$delta)*m$G+ m$delta * m$tau %*% Uo
	return(as.numeric(Ud))
}

# Computes for every ax the punishment payoff that can be decomposed
aps.pm.decompose.vi = function(m,voi,i=1, g.replies=m$g[m$replies[[i]],i]) {
	delta = m$delta
  V.ax = as.numeric(m$tau %*% voi)
  V.replies = V.ax[m$replies[[i]]]

  # cheat.payoffs is a n.ax_i x 1 vector 	
  cheat.payoff = RowMaxs(m$ind.replies.by.ax_i[[i]],
                    (1-delta)*g.replies + delta*V.replies)
	
	#	the ax_i will be mapped back to ax 	
	return(cheat.payoff[m$ind.ax.to.ax_i[[i]]])
}

aps.pm.solve.outer = function(m,eps=10^(-6),max.iter = 1000,Uo=NULL,vo=NULL,
                              feasible = rep(TRUE,m$nax), return.for.imp=FALSE) {
	("aps.pm.solve.outer")
	#re("aps.pm.solve.outer")
																
	print("")
	print("")
	print("Perfect Monitoring Outer approximation")
	
	# If not given, initialize Um as an upper bound on U...
	if (is.null(Uo)) {
		Uo = rep(max(m$G),m$nx)
	}
	#... and vm as a lower bound on v
	if (is.null(vo)) {
		min.g = colMins(m$g)
		vo = matrix(min.g,m$nx,m$n,byrow=TRUE)
	}
	
	counter = 0
	decompose.U = decompose.v = TRUE
	n.sym = m$n
	vd = matrix(Inf,m$nax,m$n)
	ai.opt = matrix(NA,m$nx,m$n)
	while(TRUE) {
		counter = counter+1
	
		# Find for every ax the maximum U that can be decomposed given Uo
		if (decompose.U)
			Ud = aps.pm.decompose.U(m,Uo)
			
		
		# Find for every ax the mimimal vi that can be decomposed given vo
		if (decompose.v) {
			for (i in 1:m$n) {
				vd[,i] = aps.pm.decompose.vi(m,vo[,i],i)
			}
			Vd = rowSums(vd)
		}
		
		# Determine which ax are feasible, i.e. can be implemented given the previous
		# Uo and vo
		feasible = feasible & (Ud >= Vd)
		
		Ud[!feasible]  = -Inf
		vd[!feasible,] = Inf
		
		# Calculate new Uo and new vo and optimal action profiles
		Uo.old = Uo; vo.old = vo;
		
		ae.opt = which.RowMaxs(m$ind.ax.by.x,Ud)
		Uo = Ud[ae.opt]


		for (i in 1:m$n) {
			ai.opt[,i] = which.RowMins(m$ind.ax.by.x,vd[,i])
			vo[,i] = vd[ai.opt[,i],i]
		}
    #ai.opt.lab = matrix(m$ax.lab[ai.opt],ncol=2)
		
		
		change.Uo = max(Uo.old-Uo)
		change.vo = max(vo-vo.old)
		
		info = c(counter,round(mean(feasible)*100,1),mean(Uo),mean(vo)*m$n,change.Uo,change.vo)
		info = round(info,5)
		names(info)=c("r","feas%","Uo","Vo","diff Uo","diff vo")
		print(info)
		
		if (any(!is.finite(Uo))) {
			mywarning(paste("Error: Found infinite value for Uo in states",m$x.lab[which(!is.finite(Uo))]))
			break;
		}
		if (max(change.Uo,change.vo*m$n)<eps) {
			print(paste("Outer approximation with precision eps=",eps))
			break;
		}
		if (counter>max.iter){
			mywarning("max.iter reached in outer approximation")
			break;
		}
	}
	
	# Check whether the resulting optimal action plan can be solved
	inner = pm.solve.action.plan(m,ae.opt,ai.opt)

	Vo = rowSums(vo)	
	mat = cbind(ae.opt,ai.opt,inner$U,Uo,rowSums(inner$vi),Vo,inner$vi,vo,Uo.old-Uo,vo-vo.old,inner$feas)
	colnames(mat) = c("ae",paste("a",1:m$n,sep=""),
	                 "U","Uo","V","Vo",
	                  paste("vi",1:m$n,sep=""),paste("vo",1:m$n,sep=""),"change.Uo", paste("change.vo",1:m$n,sep=""),
										"feas.e",paste("feas.",1:m$n,sep=""))
	rownames(mat) = paste(m$ax.lab[ae.opt],m$ax.lab.short[ai.opt[,1]],m$ax.lab.short[ai.opt[,2]])
	return(list(feasible = feasible,opt.mat=mat,min.slack=min(inner$feas)))
}

# Calculates joint equilibrium and punishment payoffs for a given action plan and determines whether the solution is feasible
pm.solve.action.plan = function(m,ae,ai) {
	("pm.solve.action.plan")
	#re("pm.solve.action.plan")
	

	U = pm.get.U.from.ae(m,ae)
	vi = matrix(NA,m$nx,m$n)
	for (i in 1:NCOL(ai)) {
		vi[,i] = pm.get.vi.from.ai(m,i,ai[,i])
	}
	
	# Check feasibility. Not yet very efficient
	# would be better if I can assign a subset of action profiles
	ax = c(ae,as.numeric(ai))
	
	Ud = aps.pm.decompose.U(m,U)
	vd = matrix(Inf,m$nax,m$n)
	for (i in 1:m$n) {
		vd[,i] = aps.pm.decompose.vi(m,vi[,i],i)
	}
	Vd = rowSums(vd)
	feas = Ud[ax]-Vd[ax]
	feas = matrix(feas,ncol=1+m$n)
	return(list(feas=feas,U=U,vi=vi))
}

# ai is an nx*1 vector of player i's punishment actions indexed as ax
pm.get.vi.from.ai = function(m,i,ai) {
	return(get.full.dyn.vi(m,i,m$ind.ax.to.ax_i[[i]][ai])$vi)
}

# ae is an nx*1 vector of equilibrium phase actions indexed as ax
pm.get.U.from.ae = function(m,ae) {
	tau.e = as.matrix(m$tau[ae,])
	U = solve(diag(m$nx)-m$delta*tau.e, (1-delta)*m$G[ae])
	return(U)
}  

imp.solve.outer = function(m,eps=10^(-4),max.iter = 1000000,Uo=NULL,vo=NULL,
							 feasible = rep(TRUE,m$nax),opt.mat = NULL,Uin=-Inf,vin=rep(Inf,m$n),plot.freq=20) {
	("imp.solve.outer")
  #re("imp.solve.outer")		
	
	remove.LP_PROBS()
	# If no outer bounds are given, we first solve the problem with perfect monitoring
	if (is.null(Uo) | is.null(vo)) {
		if (is.empty(m,"sol.mat")) {
			m.pm = solve.game(m)
		} else {
			m.pm = m
		}
		Uo = m.pm$sol.mat[,"U"]
		vo = m.pm$sol.mat[,paste("v",1:m$n,sep="")]
		feasible = feasible & (1:m$nax %in% m.pm$admiss)
		print(paste(sum(feasible), " ax feasible (", round(sum(feasible)*100 / m$nax,2),"%)",sep=""))
	}
	
	print("")
	print("")
	print("Imperfect monitoring, outer approximation")

	delta = m$delta;
	
	Ud = rep(-Inf,m$nax)
	vd = matrix(Inf,m$nax,m$n)
	U.lp = vector(mode = "list", length = m$nax)
	vi.lp = list()
	n.sym = m$n.sym
	
	for (i in 1:n.sym) {
		vi.lp[[i]] = vector(mode = "list", length = m$nax)
	}
	
	U.max = -Inf
	vi.min = rep(Inf,m$n)
	
	ae.opt = rep(NA,m$nx)
	ai.opt = matrix(NA,m$nx,m$n)
	
	Uo.change = rep(Inf,m$nx)
	vo.change = matrix(Inf,m$nx,m$n)

	
	Ud.imp = vd.imp = NULL
	
	# As first bounds on what can be decomposed, we take the values under perfect monitoring
	Ud.bound = aps.pm.decompose.U(m,Uo)
	# Find for every ax the mimimal vi that can be decomposed given vo
	vd.bound = matrix(Inf,m$nax,m$n)
	for (i in 1:n.sym) 
		vd.bound[,i] = aps.pm.decompose.vi(m,vo[,i],i)

		
	ylim = range(c(Uo,vo))
	counter = 0
	
	Vo = rowSums(vo)
	plot(Uo,ylim=ylim,main=paste(counter,"."),col="blue")
	for (i in 1:n.sym)
		points(vo[,i],col="red",pch=paste(i,sep=""),cex=0.7)


	# Table from which we select the next pair of phase k and state x that we update
	max.diff = diff(range(c(m$G,as.numeric(m$g))))*m$n
	xk.tab = NULL
	for (i in 1:n.sym) 
		xk.tab = rbind(xk.tab,cbind(i,1:m$nx,vo[,i],0,max.diff+1)) # Punishment phase
		
	xk.tab = rbind(cbind(0,m$x.sym,Uo[m$x.sym],0,max.diff),xk.tab) # Equilibrium phase
	colnames(xk.tab) = c("k","x","val","rounds","change")
	
	xk.ind = 0
	outer.counter = 0
	total.skipped = 0
	total.looked.at = 0
	num.change.wrong = 0
	while(TRUE) {
		("imp.solve.outer2")
		#re("imp.solve.outer2")		
	
		outer.counter = outer.counter +1
		# We may consider a more sophisticated scheme to assign the next x,k that will be decomposed
		xk.ind = xk.ind + 1
		if (xk.ind > NROW(xk.tab))
			xk.ind = 1

		x = xk.tab[xk.ind,"x"]
		k = xk.tab[xk.ind,"k"]

		
		##############################################################
		# Decompose U
		##############################################################
		if (k == 0) {
			

			ax = rows.to.v.ind(m$ind.ax.by.x,x)
		
			feas.ax = feasible[ax]
			
			# sort ax for the given state x
			ord.ax = ax[order(-feas.ax,-Ud.bound[ax])[1:sum(feas.ax)]]
			
			# Can ignore ax action profiles that are symmetric
			ord.ax = intersect(ord.ax,m$ax.sym)
			
			
			num.success = 0
			best.U = -Inf
			counter = 0
			skipped = 0
			for (ax in ord.ax) {
				#ax = ord.ax[1]
				# Stop the loop if even with the upper bound on payoffs does not do better
				# than the currently best decomposed value
				if (Ud.bound[ax] < best.U) {
					skipped = NROW(ord.ax)-counter
					break;	
				}
				
				print (paste(counter,". ax=",ax, sep=""))
				counter = counter+1
				ret = imp.decompose.U.ax(m=m,ax=ax,Uo=Uo,vo=vo,lp=U.lp[[ax]],x=x,update.lp=TRUE, warn.if.larger=TRUE)
				U.lp[[ax]] = ret$lp
				
				# ax can be implemented 
				if (ret$feasible) {
					num.success = num.success+1
					Ud[ax] = Ud.bound[ax] = ret$Ud
					
					# Best decomposed value so far for state x
					if (Ud[ax] > best.U) {
						best.U = Ud[ax]
						ae.opt[x] <- ax
					}
					points(x,Ud[ax],col="green")
					
					if (m$symmetric) {
						xs  = m$x.perm[x,-1]
						axs = m$ax.perm[ax,-1]
						Ud[axs] = Ud.bound[axs] =Ud[ax]
						ae.opt[xs] <- axs
					}	
				# ax cannot be implemented
				} else {
					feasible[ax] <- FALSE
					if (m$symmetric) {
						axs = m$ax.perm[ax,-1]
						feasible[axs] <- FALSE
					}	
				}
			}
			
			Uo.old = Uo[x]
			Uo[x] = best.U
			
			if (m$symmetric) {
				xs  = m$x.perm[x,-1]
				axs = m$ax.perm[ae.opt[x],-1]
				Uo[xs] <- Uo[x]
				ae.opt[xs] <- axs
			}
			
			Uo.change[x] <- Uo.old-Uo[x]	
			change = Uo.change[x]

			print(paste("Uo x=",m$x.lab[x]," ", skipped, "/", NROW(ord.ax), "=", 
					         round(skipped*100/NROW(ord.ax),1), "% action profiles skipped.",sep=""))
			print("...")

		##############################################################
		# Decompose vi
		##############################################################		
		} else {
			i = k
		

		
			ax = rows.to.v.ind(m$ind.ax.by.x,x)
			feas.ax = feasible[ax]
			# sort ax for the given state x
			ord.ax = ax[order(-feas.ax,vd.bound[ax,i])[1:sum(feas.ax)]]
		
			num.success = 0
			# Action profiles will be checked in descending order from
			# the calculated values or bounds on what can be decomposed
			#ord.ax = order(-feasible,vd.bound[[i]])[1:sum(feasible)]
			counter = 0
			skipped = 0
			best.vi = Inf
			for (ax in ord.ax) {
				# Stop the loop if even with the upper bound on payoffs does not do better
				# than the currently best decomposed value
				if (vd.bound[ax,i] > best.vi) {
					skipped = NROW(ord.ax)-counter
					break;	
				}
				counter = counter +1
				
				# Decompose vi
				ret = imp.decompose.vi.ax(m=m,ax=ax,i=1,Uo=Uo,vo=vo,lp=vi.lp[[i]][[ax]],x=x,update.lp=TRUE, warn.if.larger=TRUE)
				vi.lp[[i]][[ax]] = ret$lp
				# ax can be implemented
				if (ret$feasible) {
					vd[ax,i] = vd.bound[ax,i] = ret$vd
					num.success = num.success+1
					# The best decomposed value so far for this state
					if (vd[ax,i] < best.vi) {
						best.vi <- vd[ax,i]
						ai.opt[x,i] <- ax
					}
					points(x,vo[x,i],col="orange",pch=paste(i,sep=""),cex=0.7)
				# ax cannot be implemented
				} else {
					feasible[ax] <- FALSE
					if (m$symmetric) {
						axs = m$ax.perm[ax,-1]
						feasible[axs] <- FALSE
					}	

				}
			}
			
		
			print(paste(" v",i," x=",m$x.lab[x]," ",skipped, "/", NROW(ord.ax), "=", 
					         round(skipped*100/NROW(ord.ax),1), "% action profiles skipped.",sep=""))
			print("...")

			if (num.success == 0) {
				mywarning(paste("Could not implement any action profile in state ",m$x.lab[x]))
				mywarning("No PPE in pure strategies exists")
				return(NULL)
			}

			vo.old = vo[x,i]
			vo[x,i] = best.vi
			if (m$symmetric) {
				xs  = m$x.perm[x,-1]
				axs = m$ax.perm[ai.opt[x,i],-1]
				vo[xs,-1] <- vo[x,1]
				ai.opt[xs,-1] <- axs
			}	
			vo.change[x,i] = change = vo[x,i]-vo.old 	
		}
		
		num.change.wrong = num.change.wrong+ret$larger.set
		total.skipped = total.skipped + skipped
		total.looked.at = total.looked.at + NROW(ord.ax)
		

			
		if (outer.counter %% plot.freq == 0) {
			Vo = rowSums(vo)
			plot(Uo,ylim=ylim,main=paste(outer.counter,". change=", change),col="blue")
			for (i in 1:m$n)
				points(vo[,i],col="red",pch=paste(i,sep=""),cex=0.7)
		}
		
		info = c(outer.counter,k,x,round(mean(feasible)*100,1),round(mean(Uo),4),round(mean(vo)*m$n,4),round(change,5))
		names(info)=c("r","k","x","feas%","Uo","Vo","change")
		print(info)

		if (any(!is.finite(Uo))) {
			mywarning("Error: Infinite Uo found.")
			break;
		}
		
		if (xk.ind == NROW(xk.tab)) {
			max.change = max(c(Uo.change[m$x.sym],as.numeric(vo.change[,1:n.sym]*n.sym)))
			print("##################################################################")
			print(paste("Max.change",max.change))
			print("##################################################################")			
			if (max.change<eps) {
				print(paste("Outer approximation with precision eps=",eps))
				break;
			}
		}
		if (outer.counter>max.iter){
			mywarning("max.iter reached in outer approximation")
			break;
		}
	}
	
	print(paste("Average skipped:", round(total.skipped *100 / total.looked.at,1),"%"))
	print(paste(num.change.wrong, " times decomposed a larger set"))
	remove.LP_PROBS()
	Vo = rowSums(vo)	
	mat = cbind(ae.opt,ai.opt,Uo,Vo,vo,Uo.old-Uo,vo-vo.old)
	colnames(mat) = c("ae",paste("a",1:m$n,sep=""),
	                 "Uo","Vo",
	                  paste("vo",1:m$n,sep=""),"change.Uo", paste("change.vo",1:m$n,sep=""))
	rownames(mat) = paste(m$ax.lab[ae.opt],m$ax.lab.short[ai.opt[,1]],m$ax.lab.short[ai.opt[,2]])
	
	remove.LP_PROBS()
	return(list(feasible = feasible,opt.mat=mat))
}

# We only specify payments for states that can potentially be reached
# This means if fewer states can be reached from ax, we have fewer variables and fewer constraints
imp.make.rhs.dir.ax = function(m,ax,U,v,x=ax.to.x(m,ax),rx = get.reachable.states.from.ax.replies(m,ax)) {
	("imp.make.rhs.dir.ax ")
	#re("imp.make.rhs.dir.ax ")
	
	n = m$n; nx = m$nx; nax = m$nax; delta = m$delta
	# The different lists that will describe the constraints
		
	# states that can be directly reached with positive probability
	nrx = NROW(rx)
	# Pick continuation payoffs u in which player 2:n just get their punishment payoffs
	# and player 1 gets all excess liquidity
	u = v
	u[,1] = U-rowSums(v[,-1,drop=FALSE])
	u = u[rx,,drop=FALSE]
	v = v[rx,,drop=FALSE]

	
	################################################################
  # Action Constraints
	################################################################
	repl = mapply(get.replies,i=1:m$n,MoreArgs = list(m=m,ax=ax,keep.ax=FALSE),SIMPLIFY=FALSE)
	nrepl = sapply(repl,length)	
	nac = sum(nrepl) # Number of possible replies by all players

	ac.rhs = list()
	start.row = 0
	
	# State transitions given ax
	tau.ax = as.numeric(m$tau[ax,rx])
	# Loop through all players
	for (i in 1:m$n) {
		if (nrepl[i]>0) {
			# Matrix of state transitions for each reply of player i
			tau.repl  = as.matrix(m$tau[repl[[i]],rx,drop=FALSE])
			tau.ax.mat = matrix(tau.ax,nrepl[i],nrx,byrow=TRUE)
			# The right hand side
			
			ac.rhs[[i]] = (1-delta)*(m$g[repl[[i]],i]-m$g[ax,i])+
										delta*(tau.repl-tau.ax.mat) %*% u[,i]
		} else {
			ac.rhs[[i]]  = numeric(0)
		}
	}
	
	
  ################################################################
  # Payment Constraints
	################################################################
	
	# We have one payment constraint per player and rx
	# Each constraint has just a single column
	npc = nrx*n
	pc.rhs = c(u[,1]-v[,1],rep(0,nrx*(n-1)))
	nbc = nrx
	bc.rhs       = rep(0,nbc)
	
	rhs = as.numeric(c(unlist(ac.rhs),pc.rhs,bc.rhs))
	dir = c(rep(">=",nac),rep("<=",npc),rep(">=",nbc))

	return(list(rhs=rhs,dir=dir))
}

# We only specify payments for states that can potentially be reached
# This means if fewer states can be reached from ax, we have fewer variables and fewer constraints
imp.make.constr.ax = function(m,ax,U,v,x=ax.to.x(m,ax),rx = get.reachable.states.from.ax.replies(m,ax)) {
	("imp.make.constr.ax")
	#re("imp.make.constr.ax")
	
	n = m$n; nx = m$nx; nax = m$nax; delta = m$delta
	# The different lists that will describe the constraints
	
	# Bp >= k
	
	# Triplet matrix of coefficients 
	coef.row = list()
	coef.col = list()
	coef.val = list()
	
	# direction of comparison and rhs
	rhs = list()
	dir = list()
	
	
	# states that can be directly reached with positive probability
	nrx = NROW(rx)

	# Pick continuation payoffs u in which player 2:n just get their punishment payoffs
	# and player 1 gets all excess liquidity
	u = v
	u[,1] = U-rowSums(v[,-1,drop=FALSE])
	u = u[rx,,drop=FALSE]
	v = v[rx,,drop=FALSE]

	
	################################################################
  # Action Constraints
	################################################################
	repl = mapply(get.replies,i=1:m$n,MoreArgs = list(m=m,ax=ax,keep.ax=FALSE),SIMPLIFY=FALSE)
	nrepl = sapply(repl,length)
	
	
	nac = sum(nrepl) # Number of possible replies by all players

	ac.coeff.val=ac.coeff.col = ac.coeff.row = ac.rhs = ac.lab = list()
	start.row = 0
	
	# State transitions given ax
	tau.ax = as.numeric(m$tau[ax,rx])

	# Loop through all players
	for (i in 1:m$n) {
		if (nrepl[i]>0) {
			# Matrix of state transitions for each reply of player i
			tau.repl  = as.matrix(m$tau[repl[[i]],rx,drop=FALSE])
			tau.ax.mat = matrix(tau.ax,nrepl[i],nrx,byrow=TRUE)
		
			# Payments are negative and weighted by the probability that they take place
			ac.coeff.val[[i]] = delta*(1-delta)*as.numeric(t((tau.repl-tau.ax.mat)))
			
			ac.coeff.col[[i]] = rep(1:nrx+((i-1)*nrx),times=length(repl[[i]]))
			ac.coeff.row[[i]] = rep(1:length(repl[[i]]),each=nrx) + start.row
			
			start.row = start.row + nrepl[[i]]
			
			# The right hand side
			ac.rhs[[i]] = (1-delta)*(m$g[repl[[i]],i]-m$g[ax,i])+
										delta*(tau.repl-tau.ax.mat) %*% u[,i]
										
			if (!m$no.labels) {
				ac.lab[[i]] = paste("AC.",i," ",m$ax.lab[ax],"->",m$ax.lab.short[repl[[i]]],sep="")
			}
		} else {
				# Payments are negative and weighted by the probability that they take place
			ac.coeff.val[[i]] = ac.coeff.col[[i]] = ac.coeff.row[[i]] = ac.rhs[[i]]  = numeric(0)
		}
	}
	
	
  ################################################################
  # Payment Constraints
	################################################################
	
	# We have one payment constraint per player and rx
	# Each constraint has just a single column
	npc = nrx*n

	pc.coeff.row = 1:npc + start.row
	#pc.coeff.col = as.numeric(sapply(1:m$n, function(i) 1:nrx + (i-1)*nrx)) # perhaps simply 1:npc
	pc.coeff.col = 1:npc
	pc.coeff.val = rep(1-delta,npc) # 
	pc.rhs = c(u[,1]-v[,1],rep(0,nrx*(n-1)))
	
	if (!m$no.labels) {
		pc.lab = paste("PC.",rep(1:n,each=nrx)," ",rep(m$x.lab[rx],times=n),sep="")
	}


	start.row = start.row + npc
  ################################################################
  # Budget Constraints
	################################################################

	nbc = nrx
	bc.coeff.row = rep(1:nbc,each = n)+start.row
	bc.coeff.col = as.numeric(sapply(1:nrx, function(x) x + ((1:n)-1)*nrx))
	bc.coeff.val = rep(1,NROW(bc.coeff.col))
	bc.rhs       = rep(0,nbc)
	
	if (!m$no.labels) {
		bc.lab = paste("BC ",m$x.lab[rx],sep="")
	}

	
	start.row = start.row + nbc
	
	coeff.row = c(unlist(ac.coeff.row),pc.coeff.row,bc.coeff.row)
	coeff.col = c(unlist(ac.coeff.col),pc.coeff.col,bc.coeff.col)
	coeff.val = c(unlist(ac.coeff.val),pc.coeff.val,bc.coeff.val)
	rhs = c(unlist(ac.rhs),pc.rhs,bc.rhs)
	dir = c(rep(">=",nac),rep("<=",npc),rep(">=",nbc))
	
	lab = NULL
	if (!m$no.labels) {
		lab = c(unlist(ac.lab),pc.lab,bc.lab)
	}
	
	
	return(list(coeff=list(row=coeff.row,col=coeff.col,val=coeff.val),rhs=rhs,dir=dir,lab=lab))
}

imp.make.objective.U.ax = function(m,ax,rx=get.reachable.states.from.ax.replies(m,ax)) {
	("imp.make.objective.U.ax")
	#re("imp.make.objective.U.ax")
	
	# Minimize the total expected amount of money burning
	# This means the coefficients are given by the probabilities, we reach a particular state
	tau.ax = as.numeric(m$tau[ax,rx])
	ret = rep(tau.ax,times=m$n)
	return(ret)
}

imp.make.objective.vi.ax = function(m,ax,i=1,rx=get.reachable.states.from.ax.replies(m,ax)) {
	("imp.make.objective.vi.ax")
	#re("imp.make.objective.vi.ax")
	
	nrx = NROW(rx)
	# Minimize the total expected continuation payoff for player i
	# This is equivalent to maximizing expected payments of player i
	tau.ax = as.numeric(m$tau[ax,rx])
	ret = rep(0,nrx*m$n)
	ret[(1:nrx)+((i-1)*nrx)]=tau.ax
	return(ret)
}

imp.make.U.lp = function(m,ax,Uo,vo) {
	("imp.make.U.lp")
	#re("imp.make.U.lp")
		
	n = m$n; nx = m$nx;

	rx=get.reachable.states.from.ax.replies(m,ax)
	nrx = NROW(rx)
	
	obj    = imp.make.objective.U.ax(m,ax,rx=rx)
	constr = imp.make.constr.ax(m,ax,U=Uo,v=vo,rx = get.reachable.states.from.ax.replies(m,ax)) 
	
	var.lab = paste("p",rep(1:m$n,each=nrx),"(",m$x.lab[rx],")",sep="")
	prob.name = paste("LP-U ax=",m$ax.lab[ax],sep="")
	lp = make.glpk.lp(obj, constr$coeff, constr$dir, constr$rhs, obj.max = FALSE, 
		             row.lab=constr$lab, col.lab=var.lab, prob.name=prob.name, prob = NULL, lp=NULL,
		             set.names = !m$no.labels)
	#info.glpk.lp(lp)
	lp= imp.update.U.lp(m,lp,ax,Uo,vo)
	return(lp)
}

imp.make.vi.lp = function(m,ax,i=1,Uo,vo) {
	("imp.make.vi.lp")
	#re("imp.make.vi.lp")
			
	n = m$n; nx = m$nx;

	rx  = get.reachable.states.from.ax.replies(m,ax)
	nrx = NROW(rx)
	
	obj    = imp.make.objective.vi.ax(m=m,ax=ax,i=i,rx=rx)
	constr = imp.make.constr.ax(m,ax,U=Uo,v=vo,rx = get.reachable.states.from.ax.replies(m,ax)) 
	
	var.lab = paste("p",rep(1:m$n,each=nrx),"(",m$x.lab[rx],")",sep="")
	prob.name = paste("LP-v",i," ",m$ax.lab[ax],sep="")
	lp = make.glpk.lp(obj, constr$coeff, constr$dir, constr$rhs, obj.max = TRUE, 
		             row.lab=constr$lab, col.lab=var.lab, prob.name=prob.name, prob = NULL, lp=NULL,
		             set.names = !m$no.labels)

	#info.glpk.lp(lp)
	return(lp)
}

imp.update.U.lp = function(m,lp,ax,Uo,vo) {
	#remove.lp(lp)
	#return(imp.make.U.lp(m,ax,Uo,vo))

	ret = imp.make.rhs.dir.ax(m=m,ax=ax,U=Uo,v=vo)
	lp.set.dir.and.rhs(lp=lp,rhs=ret$rhs,dir=ret$dir)
	return(lp)
}

imp.update.vi.lp = function(m,lp,ax,i=1,Uo,vo) {
	("imp.update.vi.lp")
	#re("imp.update.vi.lp")
	
	#remove.lp(lp)
	#return(imp.make.vi.lp(m,ax,i=i,Uo,vo))
	
	ret = imp.make.rhs.dir.ax(m=m,ax=ax,U=Uo,v=vo)
	lp.set.dir.and.rhs(lp=lp,rhs=ret$rhs,dir=ret$dir)
	return(lp)
}


imp.decompose.U.ax = function(m,ax,Uo,vo,lp=NULL,x=ax.to.x(m,ax),update.lp=TRUE, warn.if.larger=TRUE,tol.decomp = 10^(-14)) {
	("imp.decompose.U.ax")
	#re("imp.decompose.U.ax")
	
	# Generate or update linear program
	if (is.null(lp)) {
		lp = imp.make.U.lp(m,ax,Uo,vo)
	} else {
		if (update.lp)
			lp = imp.update.U.lp(m,lp,ax,Uo,vo)
	}
	# Solve the linear program
	sol = solve.glpk.lp(lp, retry.with.standard.basis = !TRUE, warning.if.no.solution = !TRUE, should.always.solve = FALSE)
	#info.glpk.lp(lp, fname="lp_sol.txt")
	
	# ax can be implemented 
	if (sol$status == 0) {
		# Expected amount of money burning
		B = sol$objval
		print(paste("B=",B))
		
		# Joint continuation payoff that is decomposed
		tau.ax = as.numeric(m$tau[ax,])
		Ud <- (1-m$delta)*m$G[ax] + m$delta*(sum(tau.ax*Uo)-(1-m$delta)*B)
		larger.set = Ud>Uo[x]+tol.decomp
		if (larger.set & warn.if.larger) {
				warning(paste("Decomposed U x=",m$x.lab[x], " to ", Ud, " but Uo[x]=",Uo[x], " is smaller! That should not be possible at this point."))
				#Ud.pm = aps.pm.decompose.U(m,Uo=Uo)[ax]
				#stop()
		}
		return(list(lp=lp,feasible = TRUE,larger.set = larger.set, Ud = Ud,B=B))
	} else {
		return(list(lp=lp,feasible = FALSE,larger.set = FALSE, Ud = NULL,B=NULL))
	}
}

imp.decompose.vi.ax = function(m,ax,i=1,Uo,vo,lp=NULL,x=ax.to.x(m,ax),update.lp=TRUE, warn.if.larger=TRUE,tol.decomp = 10^(-14)) {
	("imp.decompose.vi.ax")
	#re("imp.decompose.vi.ax")
	lp = NULL
	if (is.null(lp)) {
		lp = imp.make.vi.lp(m,ax,i=i,Uo,vo)
	} else {
		if (update.lp)
			lp = imp.update.vi.lp(m,lp,ax,i=i,Uo,vo)
	}
	# Solve the linear program
	sol = solve.glpk.lp(lp, retry.with.standard.basis = !TRUE, warning.if.no.solution = !TRUE, should.always.solve = FALSE)
			
	#info.glpk.lp(lp)
	# ax can be implemented 
	if (sol$status == 0) {
		# Expected payment for player i
		pi = sol$objval
		# Joint continuation payoff that is decomposed
		tau.ax = as.numeric(m$tau[ax,])
		
		# Excess liquidity of all states is assigned to player 1
		uo = vo
		uo[,1] = Uo-rowSums(vo[,-1,drop=FALSE])

		vd <- (1-m$delta)*m$g[ax,i] + m$delta*(sum(tau.ax*uo[,i])-(1-m$delta)*pi)
		larger.set = vd<vo[x,i] - tol.decomp
		if (larger.set & warn.if.larger) {
				warning(paste("Decomposed v",i," x=",m$x.lab[x], " to ", vd, " but vo[x,i]=",vo[x,i], " is larger! That should not be possible at this point."))
				
				if (FALSE) {
				# Diagnostics
				vd.pm = aps.pm.decompose.vi(m,voi=vo[,i],i=i)[ax]
				rx = get.reachable.states.from.ax.replies(m,ax)
				nrx = NROW(rx)
				
				p1.sol = rep(0,m$nx)
				names(p1.sol) = m$x.lab
				p1.sol[rx] = sol$solution[1:nrx]
				
				
				repl = get.replies(m=m,ax=ax,keep.ax=TRUE,i=i)
				tau.repl  = as.matrix(m$tau[repl,,drop=FALSE])

				vd.sol.all <- sapply(repl,function(ax) {
					tau.ax = as.numeric(m$tau[ax,])
					print(tau.ax)
					(1-m$delta)*m$g[ax,i] + m$delta*(sum(tau.ax*uo[,i])) 
						
					-m$delta*(1-m$delta)*sum(tau.ax*p1.sol)
				})
				names(vd.sol.all) = m$ax.lab[repl]
				vd.sol.all
				br.sol = repl[which.max(vd.sol.all)]
				m$ax.lab[br.sol]
				vd.sol = max(vd.sol.all)
				
				vd.up = max(vd.man.all)
				p1.up = (Uo - rowSums(vo)) / (1-m$delta)
				p1.up[setdiff(1:m$nx,rx)] = 0
				
				vd.up.all <- sapply(repl,function(ax) {
					tau.ax = as.numeric(m$tau[ax,])
					print(tau.ax)
					(1-m$delta)*m$g[ax,i] + m$delta*(sum(tau.ax*uo[,i])) 
					-m$delta*(1-m$delta)*sum(tau.ax*p1.up)
				})	
				br.up = repl[which.max(vd.up.all)]
				m$ax.lab[br.up]
				vd.up = max(vd.up.all)
				
				c(vd.sol,vd.up,vd,vd.pm,vo[x,i])
				cbind(vd.sol.all,vd.up.all,vd.up.all-vd.sol.all)
				stop()
				}
		}
		return(list(lp=lp,feasible = TRUE,larger.set = larger.set, vd = vd))
	} else {
		return(list(lp=lp,feasible = FALSE,larger.set = FALSE, vd = NULL))
	}
}