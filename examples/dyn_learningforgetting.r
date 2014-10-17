# Try to implement the learning by doing model of Besanko et. al. Econometrica (2010) 
# In each period there is some probability phi 
# that firms are able to refill their storage at cost 0
# With probability eta the embargo ends and oil markets are perfectly competitive


# Quicker translation function
xv.to.x = function(m,xvm) {
  return( (xvm[,1]-1)*m$xv.dim[[2]]+xvm[,2] )
}


new.learning.sim = function(...) {
  sim = new.dyngame.sim(...)
  sim$parent.result.fun = sim$result 
  
  # Adapt average discounted extra sol values
  # Discount to average values during crisis
  sim$result = function(sim,para) {
    store.objects("sim.result",force=TRUE)
    #restore.objects("sim.result")
    
    delta = para$delta
    eta = para$eta
    ms= sim$ms
    if (!ms$sol.exists) {
      return(NULL)
    }
    #adprob = get.average.discounted.prob(m=ms,delta=(1-eta)*delta)
    #ms$extra.sol.ad = adprob %*% ms$extra.sol.cur
    #colnames(ms$extra.sol.ad) = paste(colnames(ms$extra.sol.cur),"ad",sep="")
    sim$parent.result.fun(sim,para)
  }
  return(sim) 
}

# Learning by doing model
make.my.game = function(para=list(x.max=10, x.learn=10, kappa=10, rho=0.85,gamma=0.1,sigma=1,delta=0.5,p.min=0,p.max=12,p.step=1,p.cheat.step=0.1,tau.sparse=TRUE,pred.pricing=TRUE, F = 1, p0=10),symmetric=TRUE) {
  if (!is.null(para)) {
    copy.into.env(source=para)
  }
  # Cost for producing one unit
	eta = log(rho,2)

  store.local.objects()

	xv.val = list(1:x.max,1:x.max)
	
  # xv will be one state vector
  act.fun = function(xv,m=NULL) {
    val = seq(p.min,p.max,by=p.step)
		val = list(val,val)
    list(val=val,lab=val, i=c(1,2))
  }
  x.group = function(xvm,m=NULL) {
		x.group = rep(1,NROW(xvm))
		#1:NROW(xvm)
	}
    
  
  # Some additional statistics of the solution that we might be interested in
  # Here we want to know about price, total output, joint profits, consumer surplus,
  # and total welfare 
  extra.sol.fun = function(avm,xvm,m=NULL) {
    p   = avm
		p.diff = p[,1]-p[,2]
		
		eta = log(rho,2)
		cost = kappa*xvm^eta
    mat = cbind(xvm,p,p.diff,p.avg=(p[,1]+p[,2])/2,cost)
		
    colnames(mat)=c("xv1","xv2","p1","p2","p.diff","p.avg","MC1","MC2")
    return(mat)
  }
    
  # Now we want to calculate average expected values

  # We do not consider a multistage game
  x.stage.fun = NULL
  
	g.fun = function(avm,xvm,m=NULL) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  	
		p = avm
		#D = 1 / (1+exp((p[,1:2]-p[,2:1])/sigma)) # Expected probability to sell one unit
		
		# Demand with an outside competitor
		D = exp(-p/sigma) / 
			  (exp(-p[,1]/sigma)+exp(-p[,2]/sigma)+exp(-p0/sigma))
		
		
		# Cost for producing one unit
		eta = log(rho,2)
		cost = kappa*xvm^eta
		cost[xvm>x.learn] = kappa*x.learn^eta
		
		pi = D*(p-cost)
		#pi = pi- F*(xvm<=x.max)
		if (!pred.pricing) {
		  pi[p<cost] = -10000 # Fine for pricing below marginal cost
		}
		
		return(pi)
  }

	
  # State transitions
  tau.fun = function(avm,xvm,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")

		# Extremely important for speeding up things...
		rownames(xvm)=rownames(avm)=NULL
		
		p = avm
		#D = 1 / (1+exp((p[,1:2]-p[,2:1])/sigma)) # Probability to sell one unit
		D = exp(-p/sigma) / 
			  (exp(-p[,1]/sigma)+exp(-p[,2]/sigma)+exp(-p0/sigma))

		Gamma = 1-(1-gamma)^xvm # Probability to forget one unit		
    # Create transition matrix
		# Sparse transition matrix
		offset = m$xv.dim[[2]]

		# Column indices of states that can be reached
		jCC = ((xvm[,1]-1))*offset+(xvm[,2])
		jCM = jCC-1;                         
		jCP = jCC+1
		jMC = jCC-offset;
		jMP = jMC+1
		jPC = jCC+offset;
		jPM = jPC-1;

		# Probablities for one player gaining, the other loosing experience
		pPM = D[,1]*(1-Gamma[,1])*Gamma[,2]
		pMP = D[,2]*(1-Gamma[,2])*Gamma[,1]

		# Probablities for one player gaining, the other player keeping experience
		pPC = D[,1]*(1-Gamma[,1])*(1-Gamma[,2])
		pCP = D[,2]*(1-Gamma[,2])*(1-Gamma[,1])

		# Probablities for one player keeping the other loosing experience
		pCM = D[,1]*(Gamma[,1])*Gamma[,2]
		pMC = D[,2]*(Gamma[,2])*Gamma[,1]
		
		# Probability for no change in probabilities
		pCC = 1-pPM-pMP-pPC-pCP-pCM-pMC
		
		# Adjust probabilities for floor and ceiling
		floor = xvm==1
		ceiling = xvm==x.max

		rows = floor[,1]
		pCP[rows]=pCP[rows]+pMP[rows]
		pCC[rows]=pCC[rows]+pMC[rows]

		rows = floor[,2]
		pPC[rows]=pPC[rows]+pPM[rows]
		pCC[rows]=pCC[rows]+pCM[rows]
		
		rows = ceiling[,1]
		pCM[rows]=pCM[rows]+pPM[rows]
		pCC[rows]=pCC[rows]+pPC[rows]

		rows = ceiling[,2]
		pMC[rows]=pMC[rows]+pMP[rows]
		pCC[rows]=pCC[rows]+pCP[rows]

			
		X = pCC+pCM+pMC+pPC+pCP+pPM+pMP
		range(X)
		
		# States to which we do not transist
		del  = which(c(floor[,1]  , floor[,1]  | ceiling[,2],
													floor[,2],rep(FALSE,NROW(xvm)),    ceiling[,2],
						 ceiling[,1]| floor[,2],ceiling[,1]))
						 
		if (length(del)>0) {	 
			i    = rep(1:NROW(xvm),times=7) [-del]
			j    = c(jMC,jMP, jCM,jCC,jCP, jPM,jPC)[-del]
			prob = c(pMC,pMP, pCM,pCC,pCP, pPM,pPC)[-del]
		} else {
			i    = rep(1:NROW(xvm),times=7) 
			j    = c(jMC,jMP, jCM,jCC,jCP, jPM,jPC)
			prob = c(pMC,pMP, pCM,pCC,pCP, pPM,pPC)
		}
	
		
		if (tau.sparse) {
			tau = sparseMatrix(i=i,j=j,x=prob)
		} else {
			tau = as.matrix(sparseMatrix(i=i,j=j,x=prob))
		}
	
		#tau = Matrix(tau,sparse=TRUE)
   	if (min(tau)<0)
			stop("Negative tau!")
		tau
  }
	
  integrated = FALSE  
  list(n=2,
       delta=delta,
       integrated = integrated,
       xv.val = xv.val,
       act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.stage.fun = x.stage.fun,
       x.group=x.group,
       extra.sol.fun=extra.sol.fun,
			 tau.sparse = tau.sparse,
			 ixv = 1:2,
			 symmetric = symmetric,
			 para = para)
}


# Parameters and their default values
STORE_OBJECTS = TRUE
NO.LABELS = !TRUE

para=list(x.max=10, x.learn=10, kappa=5, rho=0.5,gamma=0.1,sigma=0.5,delta=0.6,p.min=0,p.max=8,p0=8,F=1,p.step=1,tau.sparse = TRUE, pred.pricing = TRUE)
copy.into.env(source=para)

my.game = make.my.game(para=para,symmetric=TRUE)
m = init.game(my.game = my.game)

print(object.size(m$tau),units="Mb")

# Collusive solution
STORE_OBJECTS = TRUE
m$delta = 0.6
#mc = solve.game(m)

STORE_OBJECTS = !TRUE
ret.imp = imp.solve.outer(m,eps=10^(-3))

ret.pm = aps.pm.solve.outer(m,eps=10^(-5))
print.sol(mc)

#mm = markov.pm.algo(m,delta=m$delta,start.a = rep(1,m$nx), eps = 10^-8, sym2=TRUE)

par(mar=c(5,5,2,2))
state.levelplot(mc,z=mc$extra.sol.cur[,"p.avg"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Avg. prices",
                xlab="State 1", ylab="State 2",reverse.colors=TRUE)

m$integrated = TRUE
set.g.and.tau(m,recalc.g=TRUE,para=para)
mm = solve.game(m)
par(mar=c(5,5,2,2))
state.levelplot(mm,z=mm$extra.sol.cur[,"p.avg"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Avg. prices",
                xlab="State 1", ylab="State 2",reverse.colors=TRUE)
print.sol(mm)
								
para=list(x.max=30, x.learn=20, kappa=15, rho=0.15,gamma=0.1,sigma=0.1,delta=0.5,p.min=0,          p.max=20,p.step=1,tau.sparse = TRUE, pred.pricing = !TRUE)
set.g.and.tau(m,recalc.g=TRUE,para=para)

par(mar=c(5,5,2,2))
state.levelplot(mc,z=mc$extra.sol.cur[,"p2"],arrows=FALSE,cuts=10,
                xrange=list(c(1,5),c(1,5)), main="p2",
                xlab="State 1", ylab="State 2",reverse.colors=TRUE)

								
								
# Monopoly solution
# mon = clone(m)
# mon$integrated = TRUE
# mon = solve.game(mon,delta=delta)


								
STORE_OBJECTS = !TRUE
NO.LABELS = TRUE
options(error=traceback)


para=list(x.max=30, x.learn=15, kappa=20, rho=0.5,gamma=0.1,sigma=0.1,delta=0.5,p.min=0,p.max=30,p.step=1,tau.sparse = TRUE, pred.pricing = !TRUE)

no.g.para   = c("gamma","delta")
no.tau.para = c("x.learn","kappa","rho","sigma","delta","pred.pricing")

block = list(def = list(), 
             val = list(gamma = c(0.1,0.2,0.4),
												kappa = c(10,20),
												rho = c(1/4,1/2,3/4),
												sigma = c(0.05,0.1,1),
												pred.pricing = c(FALSE,TRUE),
                        delta = c(1/2,3/4,0.9)
                    )) 

sim = new.learning.sim(para.def=para, blocks=list(block), make.my.game=make.my.game,
                       no.g.para=no.g.para, no.tau.para = no.tau.para)
											 

run.sim(sim,file.name="learnsim.csv",append=TRUE)


sim = load.sim(file.name="learnsim.csv")

#save.sim(sim)

x = sim$sol[,"x"]
x.max = sim$sol[,"x.max"]
xv1 = ceiling(x / x.max)
xv2 = ((x-1) %% x.max)+1
sim$sol = cbind(sim$sol,xv1,xv2)


para=list(x.max=30, x.learn=15, kappa=20, rho=0.5,gamma=0.1,sigma=0.1,delta=0.5,p.min=0,p.max=25,p.step=1,tau.sparse = TRUE, pred.pricing = !TRUE)
para = c(names(para),"xv1","xv2")
draw.sim.GUI(sim$sol,para=para)


STORE_OBJECTS = TRUE
rows.pred = which(sim$mod.mat[,"pred.pricing"]==1)
rows.nopred = which(sim$mod.mat[,"pred.pricing"]==0)

for (i in 1:NROW(rows.pred)) {
	
	split.screen(c(1,2), erase = TRUE)
	row = rows.pred[i]
	lab = row

	screen(1)
  sol.levelplot(sim$sol,sim$para.mat[row,],xyz=c("xv1","xv2","p.avg"),
          xlab="x1", ylab="x2",
          main=paste("Pred. Prices",lab), reverse.colors=TRUE)

	row = rows.nopred[i]
	lab = row
	screen(2)
	sol.levelplot(sim$sol,sim$para.mat[row,],xyz=c("xv1","xv2","p.avg"),
          xlab="x1", ylab="x2",
          main=paste("No.pred Prices",lab), reverse.colors=TRUE)

}

dt = as.data.table(sim$sol)
split(dt,sim$sol["delta"])

para=as.list(sim$para.mat[2,])
copy.into.env(source=para)

my.game = make.my.game(para=para)
m = init.game(my.game = my.game)
m = solve.game(m)
print.sol(m)














funs = make.funs(K=K,x.max,eta=eta)
#funs = make.funs() 
xv.val = list(-1:x.cap,0:x.cap)
                                        
Rprof(tmp <- "profile.out", memory.profiling=FALSE)
m = init.game(n=2, name="Embargo", symmetric = FALSE,nxv=2,
              xv.val=xv.val, functions=funs)
m$delta = delta
Rprof()
summaryRprof(tmp)
STORE.OBJECTS=!TRUE

                
# First best solution
mfb = clone(m)
mfb$integrated = TRUE
mfb$goal = "W"
set.g.and.tau(mfb,recalc.g=TRUE)
mfb = solve.game(mfb,delta=delta)
print.sol(mfb)


# Monopoly solution
mon = clone(m)
mon$integrated = TRUE
mon = solve.game(mon,delta=delta)

                
state.levelplot(mon,z=get.Pae,arrows=!FALSE,cuts=10,xrange=c(0,20),main="Monopoly",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=TRUE)

			
# Collusive solution
mc = clone(m)
mc = solve.game(mc,delta=delta)  
print.sol(mc)

get.average.discounted.prob(mc)


STORE_OBJECTS = !TRUE

par(mar=c(5,5,2,2))

state.levelplot(mc,z=get.Pae,arrows=FALSE,cuts=10,xrange=c(0,20),main="Collusion",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=TRUE)
                
state.levelplot(mon,z=get.Pae,arrows=FALSE,cuts=10,xrange=c(0,20),main="Monopoly",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=TRUE)
                
state.levelplot(mfb,z=get.Pae,arrows=!FALSE,cuts=10,xrange=c(0,20),main="First Best",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=TRUE)

nm = 21
eta.seq = seq(0.1,0.9,length=nm)
#gamma.seq = c(0.1,0.3)
delta = 0.8
mmc = solve.multi.model(m=mc,delta = delta,par = list(eta=eta.seq))
mmon = solve.multi.model(m=mon,delta = delta,par = list(eta=eta.seq))
mmfb = solve.multi.model(m=mfb,delta = delta,par = list(eta=eta.seq))

mmlist = list(c=mmc,m=mmon,fb=mmfb)

# Also discount with eta




for (imm in 1:3) {
  for (im in 1:nm) {
    eta = eta.seq[im]
    m = mmlist[[imm]]$ms[[im]]
    adprob = get.average.discounted.prob(m=m,delta=(1-eta)*delta)
    m$extra.sol = adprob %*% m$extra.sol.current
    #m$extra.sol = m$extra.sol.current *10
    cols = 1:(NCOL(m$sol.mat)-NCOL(m$extra.sol))
    m$sol.mat[,-cols] = m$extra.sol
  }
  update.multimod.sol(mmlist[[imm]])
}

cols = -(1:3)
colnames(mmc$sol)[cols]  = paste(colnames(mmc$sol[,cols]),".c",sep="")
colnames(mmon$sol)[cols]= paste(colnames(mmon$sol[,cols]),".m",sep="")
colnames(mmfb$sol)[cols] = paste(colnames(mmfb$sol[,cols]),".fb",sep="")


mm = clone(mmc)
sol = cbind(mmc$sol,mmon$sol,mmfb$sol)
mm$sol = sol

par(mar=c(5,2,1,1))
xv = c(10,10)
mat = sol[sol[,"x"]==xv.to.x(mm$ms[[1]],xv),]
yvar="W"
yvars = paste(yvar,c("fb","c","m"),sep=".")
plot.multi.lines(mat,xvar="par",yvar=yvars,ylab="",xlab="Probability that crisis ends",
 yname=c("Social Planner","Collusion","Monopoly"),col=c("green","blue","red"),
 legend.pos = "topleft",legend.title="Welfare",lwd=2)
 




print.sol(mm,order="x")

STORE_OBJECTS=TRUE
plot.multi.model(mm,y="U",type="p",xv=c(10,10))


mlist = list()
Rprof(tmp <- "profile.out", memory.profiling=FALSE)
for (b in 1:nm) {
  gamma = gamma.seq[b] # Global environment will be searched for gamma
  set.g.and.tau(m,g.fun=funs$g.fun) # Update stage game payoffs
  m$name = paste("Gamma",round(gamma,3))
  mlist[[b]] = solve.game(m,delta=delta)
}



	
	 # C_g_fun = '
		# NumericVector p(av_);
		# IntegerVector xv(xv_);

		# double sigma = par["sigma"];
		# double rho = par["rho"];
		# double kappa = par["kappa"];
		# double x.learn = par["x.learn"];

		# double D0 = 1 / (1+exp((p[0]-p[1])/sigma);
		# double D1 = 1 / (1+exp((p[1]-p[0])/sigma);
						
		# // Cost for producing one unit
		# double eta = log(rho,2);
		# double cost0 = kappa*min(xv[0],x.learn)^eta;
		# double cost1 = kappa*min(xv[1],x.learn)^eta;

		# double pi0 = D0*(p[0]-cost0);
		# double pi1 = D1*(p[1]-cost1);
				
		# return(NumericVector::create(pi0,pi1))
	# '

	 # C_sparse_tau_fun = '
		# NumericVector p(av_);
		# IntegerVector xv(xv_);

		
		
		# double sigma = par["sigma"];
		# double rho = par["rho"];
		# double kappa = par["kappa"];
		# double x.learn = par["x.learn"];

		# double D0 = 1 / (1+exp((p[0]-p[1])/sigma);
		# double D1 = 1 / (1+exp((p[1]-p[0])/sigma);
						
		# // Cost for producing one unit
		# double eta = log(rho,2);
		# double cost0 = kappa*min(xv[0],x.learn)^eta;
		# double cost1 = kappa*min(xv[1],x.learn)^eta;

		# double pi0 = D0*(p[0]-cost0);
		# double pi1 = D1*(p[1]-cost1);
				
		# return(NumericVector::create(pi0,pi1))
	# '


# Public Goods Game with long and short term investments
make.my.game = function(para=list(x.max=10, x.learn=10, kappa=10, rho=0.85,gamma=0.1,sigma=1,delta=0.5,p.min=0,p.max=12,p.step=1,p.cheat.step=0.1,tau.sparse=TRUE,pred.pricing=TRUE),symmetric=TRUE) {
  if (!is.null(para)) {
    copy.into.env(source=para)
  }
  # Cost for producing one unit
	eta = log(rho,2)

  store.local.objects()

	xv.val = list(1:x.max,1:x.max)
	
  # xv will be one state vector
  act.fun = function(xv,m=NULL) {
    val = seq(p.min,p.max,by=p.step)
		val = list(val,val)
    list(val=val,lab=val, i=c(1,2))
  }
  x.group = function(xvm,m=NULL) {
		x.group = rep(1,NROW(xvm))
		#1:NROW(xvm)
	}
    
  
  # Some additional statistics of the solution that we might be interested in
  # Here we want to know about price, total output, joint profits, consumer surplus,
  # and total welfare 
  extra.sol.fun = function(avm,xvm,m=NULL) {
    p   = avm
		p.diff = p[,1]-p[,2]
		
		eta = log(rho,2)
		cost = kappa*xvm^eta
    mat = cbind(xvm,p,p.diff,p.avg=(p[,1]+p[,2])/2,cost)
		
    colnames(mat)=c("xv1","xv2","p1","p2","p.diff","p.avg","MC1","MC2")
    return(mat)
  }
    
  # Now we want to calculate average expected values

  # We do not consider a multistage game
  x.stage.fun = NULL
  
	g.fun = function(avm,xvm,m=NULL) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  	
		p = avm
		D = 1 / (1+exp((p[,1:2]-p[,2:1])/sigma)) # Expected probability to sell one unit
		
		# Cost for producing one unit
		eta = log(rho,2)
		cost = kappa*xvm^eta
		cost[xvm>x.learn] = kappa*x.learn^eta
		
		pi = D*(p-cost)
		if (!pred.pricing) {
		  pi[p<cost] = -10000 # Fine for pricing below marginal cost
		}
		
		return(pi)
  }

	
  # State transitions
  tau.fun = function(avm,xvm,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")

		# Extremely important for speeding up things...
		rownames(xvm)=rownames(avm)=NULL
		
		p = avm
		D = 1 / (1+exp((p[,1:2]-p[,2:1])/sigma)) # Probability to sell one unit

		Gamma = 1-(1-gamma)^xvm # Probability to forget one unit
		
    # Create transition matrix
		# Sparse transition matrix
		offset = m$xv.dim[[2]]

		# Column indices of states that can be reached
		jCC = ((xvm[,1]-1))*offset+(xvm[,2])
		jCM = jCC-1;                         
		jCP = jCC+1
		jMC = jCC-offset;
		jMP = jMC+1
		jPC = jCC+offset;
		jPM = jPC-1;

		# Probablities for one player gaining, the other loosing experience
		pPM = D[,1]*(1-Gamma[,1])*Gamma[,2]
		pMP = D[,2]*(1-Gamma[,2])*Gamma[,1]

		# Probablities for one player gaining, the other player keeping experience
		pPC = D[,1]*(1-Gamma[,1])*(1-Gamma[,2])
		pCP = D[,2]*(1-Gamma[,2])*(1-Gamma[,2])

		# Probablities for one player keeping the other loosing experience
		pCM = D[,1]*(Gamma[,1])*Gamma[,2]
		pMC = D[,2]*(Gamma[,2])*Gamma[,1]
		
		# Probability for no change in probabilities
		pCC = 1-pPM-pMP-pPC-pCP-pCM-pMC
		
	
		
		# Adjust probabilities for floor and ceiling
		floor = xvm==1
		ceiling = xvm==x.max

		rows = floor[,1]
		pCP[rows]=pCP[rows]+pMP[rows]
		pCC[rows]=pCC[rows]+pMC[rows]

		rows = floor[,2]
		pPC[rows]=pPC[rows]+pPM[rows]
		pCC[rows]=pCC[rows]+pCM[rows]
		
		rows = ceiling[,1]
		pCM[rows]=pCM[rows]+pPM[rows]
		pCC[rows]=pCC[rows]+pPC[rows]

		rows = ceiling[,2]
		pMC[rows]=pMC[rows]+pMP[rows]
		pCC[rows]=pCC[rows]+pCP[rows]

			
		X = pCC+pCM+pMC+pPC+pCP+pPM+pMP
		range(X)
		
		# States to which we do not transist
		del  = which(c(floor[,1]  , floor[,1]  | ceiling[,2],
													floor[,2],rep(FALSE,NROW(xvm)),    ceiling[,2],
						 ceiling[,1]| floor[,2],ceiling[,1]))
						 
		if (length(del)>0) {	 
			i    = rep(1:NROW(xvm),times=7) [-del]
			j    = c(jMC,jMP, jCM,jCC,jCP, jPM,jPC)[-del]
			prob = c(pMC,pMP, pCM,pCC,pCP, pPM,pPC)[-del]
		} else {
			i    = rep(1:NROW(xvm),times=7) 
			j    = c(jMC,jMP, jCM,jCC,jCP, jPM,jPC)
			prob = c(pMC,pMP, pCM,pCC,pCP, pPM,pPC)
		}
	
		
		if (tau.sparse) {
			tau = sparseMatrix(i=i,j=j,x=prob)
		} else {
			tau = as.matrix(sparseMatrix(i=i,j=j,x=prob))
		}
	
		#tau = Matrix(tau,sparse=TRUE)
   	tau
  }


	
	 
  # State transitions
	# Try to write a function that works quick if xvm is a single vector
	# while avm can have multiple entries
	if (1==0) {
  # tau.fun = function(avm,xvm,m=NULL) {
    # store.local.objects("tau.fun")
  	# #restore.local.objects("tau.fun")

		# # Extremely important for speeding up things...
		# rownames(xvm)=rownames(avm)=NULL
		
		# p = avm
		# D1 = 1 / (1+exp((p[,1]-p[,2])/sigma))
		# D2 = 1-D1

		# Gamma1 = 1-(1-gamma)^xvm[,1] # Probability to forget one unit
		# Gamma2 = 1-(1-gamma)^xvm[,2] # Probability to forget one unit

		# # Adjust probabilities for floor and ceiling
		# floor   = xvm==1
		# ceiling = xvm==x.max

		# pr1.minus = ((1-D1)*Gamma1)
		# pr1.plus  = (D1*(1-Gamma1))
		# pr1.const = (1-pr1.minus-pr1.plus)
		
		# pr2.minus = (1-D2)*Gamma2
		# pr2.plus  = D2*(1-Gamma2)
		# pr2.const = 1-pr2.minus-pr2.plus



		# pr1.const  = pr1.const+pr1.plus*ceiling[,1]+pr1.minus*floor[,1]
		# pr1.minus  = pr1.minus*(1-floor[,1])
		# pr1.plus   = pr1.plus*(1-ceiling[,1])
		
		# pr2.const  = pr2.const+pr2.plus*ceiling[,2]+pr2.minus*floor[,2]
		# pr2.minus  = pr2.minus*(1-floor[,2])
		# pr2.plus   = pr2.plus*(1-ceiling[,2])
		
		# pr1.const  = pr1.const+pr1.plus*ceiling[,1]+pr1.minus*floor[,1]

		# # Column indices of states that can be reached
		# jCC = ((xvm[,1]-1))*offset+(xvm[,2])
		# jCM = jCC-1;                         
		# jCP = jCC+1;
		# jMC = jCC-offset;
		# jMM = jMC-1;
		# jMP = jMC+1
		# jPC = jCC+offset;
		# jPM = jPC-1;
		# jPP = jPC+1
			
		# # Probablities for each column that can be reached
		# pCC = pr1.const*pr2.const
		# pCM = pr1.const*pr2.minus
		# pCP = pr1.const*pr2.plus
		# pMC = pr1.minus*pr2.const
		# pMM = pr1.minus*pr2.minus
		# pMP = pr1.minus*pr2.plus
		# pPC = pr1.plus *pr2.const
		# pPM = pr1.plus *pr2.minus
		# pPP = pr1.plus *pr2.plus


		# # States to which we do not transist
		# del  = c(floor[,1]  | floor[,2],floor[,1]  ,floor[,1]  | ceiling[,2],
			       # floor[,2],rep(FALSE,NROW(xvm)),    ceiling[,2],
						 # ceiling[,1]| floor[,2],ceiling[,1],ceiling[,1]| ceiling[,2])
							 
		# if (NROW(xvm)>1) {					 
			# i    = rep(1:NROW(xvm),times=9) [!del]
			# j    = c(jMM,jMC,jMP, jCM,jCC,jCP, jPM,jPC,jPP)[!del]
			# prob = c(pMM,pMC,pMP, pCM,pCC,pCP, pPM,pPC,pPP)[!del]
		# } else {
			# i    = rep(1:NROW(avm),each=sum(!del))
			# j    = rep(c(jMM,jMC,jMP, jCM,jCC,jCP, jPM,jPC,jPP)[!del],times=NROW(avm))
			# prob = c(pMM,pMC,pMP, pCM,pCC,pCP, pPM,pPC,pPP)[!rep(del,times=]
		# }
		# if (tau.sparse) {
			# tau = sparseMatrix(i=i,j=j,x=prob)
		# } else {
			# tau = as.matrix(sparseMatrix(i=i,j=j,x=prob))
		# }
	
		# #tau = Matrix(tau,sparse=TRUE)
   	# tau
  # }
  }

	
  integrated = FALSE
  
  list(n=2,
       delta=delta,
       integrated = integrated,
       xv.val = xv.val,
       act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.stage.fun = x.stage.fun,
       x.group=x.group,
       extra.sol.fun=extra.sol.fun,
			 tau.sparse = tau.sparse,
			 ixv = 1:2,
			 symmetric = symmetric,
			 para = para)
}

