# Story: A market is hit by an oil embargo
#
# There are two firms that have oil storage capacity of x.cap
#
# In each period there is some probability phi 
# that firms are able to refill their storage at cost 0
# With probability eta the embargo ends and oil markets are perfectly competitive

new.oil.sim = function(...) {
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
    adprob = get.average.discounted.prob(m=ms,delta=(1-eta)*delta)
    ms$extra.sol.ad = adprob %*% ms$extra.sol.cur
    colnames(ms$extra.sol.ad) = paste(colnames(ms$extra.sol.cur),"ad",sep="")
    sim$parent.result.fun(sim,para)
  }
  return(sim) 
}

para.const = list(inc.type=c("det","poisson"),
                  sol.type=c("coll","mon","fb.w","fb.cs"))
string.const.to.int.const(para.const,prefix="") # Generates integer const in uppercase
FB.W


# Public Goods Game with long and short term investments
make.my.game = function(x.cap=10, K=10, x.inc=2, eta=0.05,delta=0.5,
                        sol.type=COLL,inc.type=DET,para=NULL) {
  if (!is.null(para)) {
    copy.into.env(source=para)
  }
    
  store.local.objects()                     
  # xv will be one state vector
  act.fun = function(xv,m=NULL) {
    store.local.objects("act.fun")
  	#restore.local.objects("act.fun");restore.local.objects("make.funs")
  	if (xv[1]==-1) {
  	  val=list(0,0)
  	  lab = list("-","-")
	  } else {
      val = list(0:xv[1],0:xv[2])
  	  lab = val
	  }
    list(val=val,lab=lab, i=c(1,2))
  }
  x.group = function(xvm,m=NULL) {
    #x.group = rep(1,NROW(xvm))
    x.group = 1:NROW(xvm)
    x.group
  } 
  
  C.g.fun = '
		NumericVector q;
		double P,Q,CS,Pi;
		double K = par["K"];
		int sol_type = par["sol_type"];
	
		q = av;
		P = K-q[0]-q[1]
  	if (P<0) P = 0;
  	
  	pi = P*q 
  	
  	# For first best solution, we may want to study alternative goals
  	if (sol_type == FB_W || sol_type == FB.CS) {
    	Q  = q[0]+q[1]
    	CS = (K-P)*min(Q,K)*(1/2)
    	Pi = pi[,0]+pi[,1]
    	if (sol.type==FB.CS) {
      	pi[,0] = CS
      	pi[,1] = 0
    	} else if (sol.type==FB.W) {
      	pi[,0] = CS
      	pi[,1] = Pi
    	}
  	}
      	 	
  	# If the embargo ends, profits are zero
  	pi[xvm[,1]==-1,] = 0

  	pi

	'
	
  g.fun = function(avm,xvm,m=NULL) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  	
  	q = avm
  	P = K-q[,1]-q[,2]
  	P[P<0]=0
  	
  	pi = P*q 
  	
  	# For first best solution, we may want to study alternative goals
  	if (sol.type %in% c(FB.W,FB.CS)) {
    	Q  = q[,1]+q[,2]
		
    	CS = (K-P)*pmin(Q,K)*(1/2)
    	Pi = pi[,1]+pi[,2]
    	if (sol.type==FB.CS) {
      	pi[,1] = CS
      	pi[,2] = 0
    	} else if (sol.type==FB.W) {
      	pi[,1] = CS
      	pi[,2] = Pi
    	}
  	}
      	 	
  	# If the embargo ends, profits are zero
  	pi[xvm[,1]==-1,] = 0

  	pi
  }
  
  # Some additional statistics of the solution that we might be interested in
  # Here we want to know about price, total output, joint profits, consumer surplus,
  # and total welfare 
  extra.sol.fun = function(avm,xvm,m=NULL) {
    q   = avm
    P   = K-q[,1]-q[,2]
    P[P<0]=0
    pi = P*q 
    Q  = q[,1]+q[,2]
    CS = (K-P)*pmin(Q,K)*(1/2)
    Pi = pi[,1]+pi[,2]
    W  = CS+Pi
    
      
    mat = cbind(q,P,Q,CS,Pi,W)
    colnames(mat)[1:2]=c("q1","q2")
    return(mat)
  }
    
  # Now we want to calculate average expected values

  # We do not consider a multistage game
  x.stage.fun = NULL
    
  # State transitions
  tau.fun = function(avm,xvm,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")

  	tau = matrix(0,NROW(avm),m$nx)  	  	
    # New state after investment
    # For the moment assume, firms allways get one unit of supply restored
    stopped = xvm[,1] == -1
    
    xvm = with.floor.and.ceiling(xvm - avm + x.inc ,floor=0,ceiling=x.cap)

    xvm[stopped,1] = -1
    # Translate state values to state indices
    ind.mat = cbind(1:NROW(xvm),xv.to.x(m,xvm))
  	tau[ind.mat]=1-eta

  	# Embargo Stopped
  	xvm[,1] = -1
  	xvm[,2] = 0
  	ind.mat = cbind(1:NROW(xvm),xv.to.x(m,xvm))
  	tau[ind.mat]=tau[ind.mat]+eta
  	
   	tau
  }
  
  integrated = FALSE
  if (!is.null(para$sol.type)) {
    integrated = para$sol.type != COLL
  }
  
  list(n=2,
       delta=delta,
       integrated = integrated,
       xv.val = list(-1:x.cap,0:x.cap),
       act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.stage.fun = x.stage.fun,
       x.group=x.group,
       extra.sol.fun=extra.sol.fun)
}




# Parameters and their default values
STORE_OBJECTS = !TRUE
NO.LABELS = TRUE

Rprof(tmp <- "profile.out", memory.profiling=FALSE)
para     = list(delta=0.8,eta=0.2,K=20,x.cap=20,x.inc=2,
                inc.type=DET,sol.type=COLL)
copy.into.env(source=para)


my.game = make.my.game(para=para)
mc = init.game(my.game = my.game)

Rprof()
summaryRprof(tmp)

# Collusive solution
STORE_OBJECTS = !TRUE
mc = solve.game(mc)  
#print.sol(mc)


m = clone(mc)

# Monopoly solution
mon = clone(m)
mon$integrated = TRUE
mon = solve.game(mon,delta=delta)

par(mar=c(5,5,2,2))
state.levelplot(mc,z=mc$extra.sol.cur[,"P"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Prices under Collusion",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",zlim=c(8,20),
                reverse.colors=TRUE)

state.levelplot(mc,z=mc$extra.sol.cur[,"P"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Prices under Collusion",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",zlim=c(8,20),
                reverse.colors=TRUE,col.scheme = "grey")

V = mc$sol.mat[,"v1"]+mc$sol.mat[,"v2"]
state.levelplot(mc,z=V,arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Sum of punishment payoffs",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=!TRUE,col.scheme = "grey")
                
V = mc$sol.mat[,"v1"]+mc$sol.mat[,"v2"]
state.levelplot(mc,z=mc$sol.mat[,"v2"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Punishment payoffs player 1",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=!TRUE,col.scheme = "grey")
                
                
state.levelplot(mc,z=mc$extra.sol.cur[,"P"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Preise unter Kollusion",
                xlab="Ölreserven Firma 1", ylab="Ölreserven Firma 2",zlim=c(8,20),
                reverse.colors=TRUE)
                
                
par(mar=c(5,5,2,0))
state.levelplot(mon,z=mon$extra.sol.cur[,"P"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Prices under Monopoly",zlim=c(8,20),
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=TRUE)
                
                
V = mc$sol.mat[,"v1"]+mc$sol.mat[,"v2"]


state.levelplot(mc,z=V,arrows=FALSE,cuts=10,xrange=c(0,20),main="Punishment Collusion",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=TRUE)

par(mar=c(5,5,2,2))
                
state.levelplot(mon,z="P",arrows=FALSE,cuts=10,xrange=c(0,20),main="Monopoly",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=TRUE)
                


								
STORE_OBJECTS = !TRUE
NO.LABELS = TRUE
options(error=traceback)


para     = list(delta=0.9,eta=0.2,K=20,x.cap=10,x.inc=4,
                inc.type=DET,sol.type=COLL)
para.change = list(delta="",eta="",K="",x.cap="ax",x.inc="",
                   inc.type=DET,sol.type=COLL)
block = list(def = list(), 
             val = list(sol.type = c(COLL,MON),
                        x.inc = c(4)
                    )) 
STORE_OBJECTS = !TRUE
sim = new.oil.sim(para.def=para, blocks=list(block), make.my.game=make.my.game)
run.sim(sim,file.name="sim20.csv",append=TRUE)


load.sim(sim)

#save.sim(sim)


for (i in 1:NROW(sim$mod.mat)) {
  lab = paste("+",sim$mod.mat[i,"x.inc"],"eta=",sim$mod.mat[i,"eta"],
              para.const$sol.type[sim$mod.mat[i,"sol.type"]],sep=" ")
  sol.levelplot(sim$sol,sim$para.mat[i,],xyz=c("xv1","xv2","P"),
          xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
          main=paste("Prices",lab),xrange=c(0,Inf), reverse.colors=TRUE)
}

delta = 0.9
eta = 0.2
K = 40
x.cap = 20


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

