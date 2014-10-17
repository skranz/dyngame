# Industry Dynamics á la Besanko & Doraszelski (RAND, 2004)

# Two firms compete á la Cournot and can invest into capacity

# States:  production capacity of each player, between 0 and x.max
# Actions: investment size 0 or 1
#          output between 0 and A and at most at capacity qi <=x1

# Demand: P = (100-Q)/beta

            
##################################################################
# Section 5: Cournot Games
##################################################################

library(repgames)
# Section 5.1 ##############################################
init.cournot = function(n=2,q.seq, A, B, MC) {
  P.fun = function(Q) {A-B*Q}
  cost.fun = function(q) {MC*q}  

  g.fun = function(qm) {
    Qm = rowSums(qm)
    Pm = P.fun(Qm)
    Pm[Pm<0]=0
    g = matrix(NA,NROW(qm),n)
    for (i in 1:n) {
      g[,i] = Pm*qm[,i]-cost.fun(qm[,i])
    }
    g
  }
  name=paste(n,"Cournot Game (g.fun)")
  m = repgames::init.game(n=2,g.fun=g.fun,action.val = q.seq,
                name=name, lab.ai=round(q.seq,2))
  return(m)
}
# Parameters of the model
n=2; A=100; B=1; MC=10;
q.seq = seq(0,100,by=1); # Grid of possible quantities
m = init.cournot(n=n,q.seq=q.seq,A=A,B=B,MC=MC)
m = repgames::solve.game(m)
m.rep = m
plot(m,legend.pos = "right")


make.funs = function(beta,MC.inv,MC.q,q.step,
                    x.min,x.max,x.step, prob.dep, prob.inv) {
  store.local.objects()
                      
  # output decisions
  stat.act.fun = function(xv,m=NULL) {
    store.local.objects("stat.act.fun")
  	#restore.local.objects("stat.act.fun");restore.local.objects("make.funs")

    val = list(seq(0,xv[1],by=q.step),
               seq(0,xv[2],by=q.step))
    i   = 1:2
    names(val) = c("q1","q2")
    
    list(val=val,i=i)
  }
  
  stat.g.fun = function(q,xvm,m=NULL) {
  	P = pmax((100-q[,1]-q[,2])/beta,0)
    (P-MC.q)*q
  }  
  
  # investment decisions
  act.fun = function(xv,m=NULL) {
    val = replicate(2,(-1):1,simplify=FALSE)
    names(val) = c("inv1","inv2")
    i = 1:2    
    list(val=val,i=i)
  }
  
  g.fun = function(avm,xvm,m=NULL) {
  	-abs(avm) * MC.inv
  }

  stat.x.group = function(xvm,m) 1:m$nx
  x.group = 0
  
  # state transitions  
  tau.fun = function(avm,xv,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun");restore.local.objects("make.funs")

  	#tau = Matrix(0,NROW(avm),m$nx) # Generates a Sparse Matrix
  	tau = matrix(0,NROW(avm),m$nx) # Generates a traditional matrix
  	
  	ind.mat = cbind(1:NROW(xv),NA)
  	
    # New state after investment
    xv = with.floor.and.ceiling(xv + avm*x.step,floor=0,ceiling=x.max)

    # Translate state values to state indices
    ind.mat[,2] = xv.to.x(m,xv)
  	tau[ind.mat]=1-prob.dep
  	
  	# New state if deprecation takes place
    xv = with.floor(xv - x.step,floor=x.min)
    
    # Translate state values to state indices
    ind.mat[,2] = xv.to.x(m,xv)
  	tau[ind.mat]=tau[ind.mat]+ prob.dep
  	
    tau
  }
  
  list(stat.act.fun=stat.act.fun,
       stat.g.fun=stat.g.fun,
       act.fun = act.fun,
       g.fun = g.fun,
       tau.fun=tau.fun,
       stat.x.group=stat.x.group,
       x.group=x.group)
}


STORE_OBJECTS = TRUE
reset.call.counts()

MC.inv = 100
MC.q = 10
beta = 1
q.step = 1

x.min=0;x.max=100;x.step = 10
prob.dep=0; prob.inv=1

xv.val = seq(x.min,x.max,by=x.step)
funs = make.funs(beta,MC.inv,MC.q,q.step,
                    x.min,x.max,x.step, prob.dep, prob.inv) 

                    
                    
Rprof(tmp <- "profile.out", memory.profiling=FALSE)
m = init.game(n=2, name="B&D", symmetric = TRUE,nxv=2, xv.val=xv.val, functions=funs)
Rprof()
summaryRprof(tmp)
#library(proftools)
#printProfileCallGraph(readProfileData(tmp))

m.init = clone(m)
m = clone(m.init)

#tracemem(m)      
ms = solve.game(m,delta=0.1)        
print.sol(ms)

# Cournot Equilibrium
P.c = (100 + 2*MC.q) / 3
q.c = (100 - MC.q) / 3
Pi.c = 2*(P.c-MC.q)*q.c 
P.m = (100+MC.q) / 2
Q.m = (100-MC.q) / 2
Pi.m = (P.m-MC.q)*Q.m
c(Pi.c,Pi.m)


#delta.seq = c(0.5,0.6,0.62,0.65,0.67,0.7,0.8,0.9)
delta.seq = seq(0,0.4,by=0.01)
mlist = list()
Rprof(tmp <- "profile.out", memory.profiling=FALSE)
for (b in 1:NROW(delta.seq)) {
  mlist[[b]] = solve.game(m,delta=delta.seq[b], augment.sol.mat = TRUE)  
}
sol = sublist(mlist,"sol.mat")
sol = rbind.list(sol)

Rprof()
#summaryRprof(tmp)

plot(m.rep,yvar="Ue",xlim=c(0,0.5))
solx = sol[sol[,"x"]==121,]
lines(solx[,"delta"],solx[,"U"],type="l")
solx = sol[m$x.lab[sol[,"x"]]=="50;0",]
lines(solx[,"delta"],solx[,"U"],type="l",col="green")



# Expected profits under 0 investment

U = sapply(mlist, function(m) m$sol.mat[1,"U"])
plot(delta.seq,U)

for (b in 1:length(mlist)) {
  m = mlist[[b]]
  
  plot(rowSums(m$xv.mat),m$sol.mat[,"U"],main=m$name,xlab="Total Capacity",ylab="Total Profits")
}



for (b in 1:length(mlist)) {
  m = mlist[[b]]
  plot(rowSums(m$xv.mat),m$sol.mat[,"U"],main=m$name,xlab="Total Capacity",ylab="Total Profits")
}


sum.av = function(m) {
  ave.mat = ax.to.av(m,m$sol.mat[,"ax.e"])
  rownames(ave.mat) = m$x.lab
  AVe = rowSums(ave.mat)
  AVe
}  
sum.av(m)

state.plot = function(m,f,name="") {
  sk.levelplot(grid.xyz = cbind(m$xv.mat,f(m)),xlab="x1",ylab="x2",main=paste(name,m$name))
}  


# Show total investment
for (b in 1:length(mlist)) state.plot(mlist[[b]], f= sum.av, name="Tot. Inv")

# Show total profits
for (b in 1:length(mlist))
   state.plot(mlist[[b]], f= function(m) m$sol.mat[,"U"], name="U")



av1.mat = a.to.av(m,m$sol.mat[,"a1"])
sk.levelplot(grid.xyz = cbind(m$xv.mat,av1.mat[,2]),xlab="x1",ylab="x2",main="a1.2")
sk.levelplot(grid.xyz = cbind(m$xv.mat,m$sol.mat[,"U"]),xlab="x1",ylab="x2",main="U")


eq.dist = get.eq.dist(m,100)
sol.mat.e(m)


library(proftools)
printProfileCallGraph(readProfileData(tmp))


unlink(tmp)


# Dynamic action group: Investments
ag = ListEnv()
inv.min = 0; inv.max = 1;
ag$av.val = inv.min:inv.max;
ag$g.fun  = g.funs[[1]]
ag$type = "dyn"
ag.dyn = ag

# Static action group: Output
q.min = 0; q.max = 100; q.step = 10; q.seq = seq(q.min,q.max,by = q.step)
ag = ListEnv() # Need to generate a new action group
ag$av.val = q.seq;
ag$g.fun  = g.funs[[2]]
ag$type = "static"
stat.ag = ag


# States: capacity of each firm
x.min = 0; x.max = 100; x.step = 10; x.seq = seq(x.min,x.max,by=x.step);
xv.val = list(x.seq,x.seq);
tau.fun = make.tau.fun(x.min  = x.min,x.max=x.max,x.step=x.step, p.dep = p.dep)










# Prisoners' Dilemma Game with investment that increases externalities of cooperation


# Example: Prisoner's Dilemma Game with payoffs
# An activity profile shall have the following structure:
# z[,1],z[,2]: Actions: Cooperate or Defect of player 1 and 2
# z[,3],z[,4]: Investments 0 = no investment, 1 = investment


make.g.fun = function(s=1,d=1,x.gain = 1, inv.cost = 1) {
  g.fun = function(avm,xvm) {
    store.local.objects("g.fun"); store.local.objects("make.g.fun",parent.num=-2);
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  
    g = matrix(NA,NROW(avm),2)
  
    C = 1; D = 2;
    g1.mat = matrix(c(
      1, -s,
      1+d, 0),2,2,byrow=TRUE)
   
    g[,1] = g1.mat[cbind(avm[,1],avm[,2])]
    g[,2] = g1.mat[cbind(avm[,2],avm[,1])]
    
    g[,1] = g[,1] + (avm[,2]==C & xvm[,1] == 2)*x.gain
    g[,2] = g[,2] + (avm[,1]==C & xvm[,2] == 2)*x.gain
    
    g[,1] = g[,1] - avm[,3]*inv.cost
    g[,2] = g[,2] - avm[,4]*inv.cost
    g
  }
  return(g.fun)
}

s = 1; d=1;
x.gain = 1;
inv.cost = 1;
g.fun = make.g.fun(s=s,d=d,x.gain=x.gain,inv.cost=inv.cost)


dep.prob = 0.5
inv.success = 0.9
# Let av and xv have 
# returns a vector of transition probabilities for each row
tau.fun = function(avt,xv,m=NULL) {
  store.local.objects("tau.fun")
	#restore.local.objects("tau.fun")

  L = 1; H = 2
  i = 1
  pH = rep(NA,NROW(avt))

  for (i in 1:2) {  
    pH[xv[,i] == H & avt[,i] == 0] = 1-dep.prob
    pH[xv[,i] == H & avt[,i] == 1] = (1-dep.prob) + dep.prob * inv.success
    pH[xv[,i] == L & avt[,i] == 0] = 0
    pH[xv[,i] == L & avt[,i] == 1] = inv.success
    if (i == 1) pH1 = pH
    if (i == 2) pH2 = pH
  }
  tau = matrix(1,NROW(avt),m$nx)
  
  xH1 = m$xv.mat[,1] == H
  xH2 = m$xv.mat[,2] == H
  tau[,xH1 & xH2] = pH1*pH2
  tau[,xH1 & !xH2] = pH1*(1-pH2)
  tau[,!xH1 & xH2] = (1-pH1)*pH2
  tau[,!xH1 & !xH2] = (1-pH1)*(1-pH2)
  
  tau
}

av.names = c("act1","act2","in1","in2")
av.val = list( c(1,2), c(1,2),c(0,1),c(0,1))
av.lab = list( c("C","D"), c("C","D"), c(0,1),c(0,1))

avt = c(3,4)

xv.lab = list(c("L","H"),c("L","H"))
xv.val = list( c(1,2), c(1,2) )
xv.names = NULL
 
iav = c(1,2,1,2)



m = init.game(n=2,g.fun=g.fun, tau.fun = tau.fun, iav = iav, avt=avt,
              av.val=av.val,av.names = av.names, av.lab = av.lab,
              xv.val=xv.val,xv.names = xv.names, xv.lab = xv.lab, 
              name="Investment PD", symmetric = TRUE)
              
m$delta = 0.4

solve.game(m)
