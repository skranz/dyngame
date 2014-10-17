# Story: Small firms have a relational contract
#        Total surplus created can have an impact on macro level

new.relmac.sim = function(...) {
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


gamma = 0.00001
a = seq(0,1.4,length=200)
k = (a^(1+gamma))/(1+gamma)
Pi = a-k
Pid = diff(Pi)
plot(a,Pi,type="l",col="blue")
lines(a[-1],Pid,col="red")
abline(h=0,col="grey")

# Parameters of effort choice game
# Cost function has the form a^gamma
make.my.game = function(para = list(num.a=10, gamma=0.01,delta=0.5,eps.tau = 0.1)) {
  if (!is.null(para)) {
    copy.into.env(source=para)
  }
    
  store.local.objects()

  xv.val = list(0:1,seq(0,1,length=num.a))
	
  # xv will be one state vector
  act.fun = function(xv,m=NULL) {
    store.local.objects("act.fun")
  	#restore.local.objects("act.fun");restore.local.objects("make.funs")
    val = list(seq(0,1,length=num.a),0)
    list(val=val,lab=val, i=c(1,2))
  }
  x.group = function(xvm,m=NULL) {
    rep(1,NROW(xvm))
  } 
  g.fun = function(avm,xvm,m=NULL) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  	
		effort = avm[,1]
		k  = (1/(1+gamma))*effort^(1+gamma)
		pi = cbind(- k,xvm[,1]*effort) 
  	pi
  }
  
  # Some additional statistics of the solution that we might be interested in
  # Here we want to know about price, total output, joint profits, consumer surplus,
  # and total welfare
  extra.sol.fun = function(avm,xvm,m=NULL) {
		effort = avm[,1]
		k  = (1/(1+gamma))*effort^(1+gamma)
		pi = cbind(xvm[,1]*effort, - k) 
		Pi = pi[,1]+pi[,2]
    mat = cbind(effort,k,Pi)
    return(mat)
  }
    

  # State transitions are aggregated on a macro level
  tau.fun = function(avm,xvm,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")

  	tau = matrix(0,NROW(avm),m$nx)  	  	

		# Deterministic global state
		xglob = avm[,1]


		# Combination of geting an order and global state
		xvm = cbind(1,xglob)
  	tau[xv.to.x(m,xvm)]=xglob
		
		# Combination of geting NO order and global state
		xvm = cbind(0,xglob)
  	tau[xv.to.x(m,xvm),]=1-xglob

		# Make noise tau matrix
		noise.tau = matrix(0,NROW(avm),m$nx)
		xglob = xv.vak[[2]] # Set of possible global states
		xvm = cbind(1,xglob)
  	noise.tau[xv.to.x(m,xvm)]=xglob
		xvm = cbind(0,xglob)
  	noise.tau[xv.to.x(m,xvm)]=1-xglob
		
		tau = (1-eps.tau)*tau + eps.tau * noise.tau 
   	tau
  }
  
  integrated = FALSE
  if (!is.null(para$sol.type)) {
    integrated = para$sol.type != COLL
  }
  
  list(n=2,
			 tau.type = "macro",
       delta=delta,
       integrated = integrated,
       xv.val = xv.val,
       act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.group=x.group,
       extra.sol.fun=extra.sol.fun)
}


# Parameters and their default values
STORE_OBJECTS = TRUE
NO.LABELS = !TRUE

Rprof(tmp <- "profile.out", memory.profiling=FALSE)
para = list(num.a=10, gamma=0.01,delta=0.5,eps.tau = 0.1)
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

