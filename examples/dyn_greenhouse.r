# A green house gas model

# Two countries excerting negative pollution externalities on their neighbours

# They can invest into more efficient industry
#

# Need: optimal solution if players maximize only their own utility
#       optimal solution under fully cooperative behavior

make.funs = function(k.R,k.I,x.max) {
  # Copy global variables into local environment
  # copy.into.env(source = globalenv(), names = ls.vars(globalenv()))  
  store.local.objects() 
                      
  # xv will be one state vector
  act.fun = function(xv,m=NULL) {
    store.local.objects("stat.act.fun")
  	#restore.local.objects("stat.act.fun");restore.local.objects("make.funs")
  	val = list(c(-1,0,1),c(-1,0,1))
  	lab = list(c("-","0","+"),c("-","0","+"))  	   
    list(val=val,lab=lab, i=1:2)
  }
  x.group = function(xvm,m=NULL) {
    x.group = rep(1,NROW(xvm))
    x.group
  } 
  
  # Total polution level is given by
  # ((2*x.max - x1 -x2)/(2*x.max)) * 100 
  g.fun = function(avm,xvm,m=NULL) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 

  	pol = ((2*x.max - xvm[,1] -xvm[,2])/(2*x.max)) * 100  
    g = matrix(100-pol,NROW(avm),2)
    
  	pol1 = ((x.max - xvm[,1])/(x.max)) * 100 
  	pol2 = ((x.max - xvm[,2])/(x.max)) * 100 
    g = matrix(c(100-pol2,100-pol1),NROW(avm),2)
    
    # Cost of reducing technological state
    g = g - k.R* (avm==-1)    
    # Cost of reducing technological state
    g = g - k.I * (avm==1)

    g    
  }
  
  # We do not consider a multistage game
  x.stage.fun = NULL
    
  # State transitions
  tau.fun = function(avm,xvm,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")

  	tau = matrix(0,NROW(avm),m$nx)

  	#For the moment, we abstract from depreciation
  	xvm.dest = with.floor.and.ceiling(xvm + avm, floor = 0, ceiling  = x.max)
  	
  	x.dest = xv.to.x(m,xvm.dest)
  	
  	# From any offered price we move to state with xv = -3 if accepted
  	rows = 1:NROW(xvm)
  	tau[cbind(rows,x.dest)] = 1

   	tau
  }
  
  list(act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.stage.fun = x.stage.fun,
       x.group=x.group)
}


STORE_OBJECTS = TRUE

k.R =  110  # costs to reduce   technological state
k.I =  50   # costs to increase technological state

x.max = 10
funs = make.funs(k.R=k.R,k.I=k.I,x.max = x.max) 

xv.val = list(0:x.max,0:x.max)
m = init.game(n=2, name="Greenhouse", symmetric = FALSE,nxv=2,
              xv.val=xv.val, functions=funs)
m$g

              
delta = 0.85
ms = solve.game(m,delta=delta)  
print.sol(ms)

delta.seq = c(0.5,0.6,0.62,0.65,0.67,0.7,0.8,0.9)
delta.seq = seq(0.5,0.9,by=0.02)
mlist = list()
Rprof(tmp <- "profile.out", memory.profiling=FALSE)
for (b in 1:NROW(delta.seq)) {
  mlist[[b]] = solve.game(m,delta=delta.seq[b])
}

Rprof()
#summaryRprof(tmp)


for (b in 1:length(mlist)) {
  m = mlist[[b]]
  plot(rowSums(m$xv.mat),m$sol.mat[,"U"],main=m$name)
}


sum.av = function(m) {
  ave.mat = a.to.av(m,m$sol.mat[,"ae"])
  rownames(ave.mat) = m$x.lab
  AVe = rowSums(ave.mat)
}  
state.plot = function(m,f,name="") {
  sk.levelplot(grid.xyz = cbind(m$xv.mat,f(m)),xlab="x1",ylab="x2",main=paste(name,m$name))
}  

# Show total investment
lapply(mlist,state.plot,f=sum.av,name="Tot. Inv")

lapply(mlist,state.plot,f=function(m) m$sol.mat[,"U"],name="U")



av1.mat = a.to.av(m,m$sol.mat[,"a1"])
sk.levelplot(grid.xyz = cbind(m$xv.mat,av1.mat[,2]),xlab="x1",ylab="x2",main="a1.2")
sk.levelplot(grid.xyz = cbind(m$xv.mat,m$sol.mat[,"U"]),xlab="x1",ylab="x2",main="U")


eq.dist = get.eq.dist(m,100)
sol.mat.e(m)


library(proftools)
printProfileCallGraph(readProfileData(tmp))


unlink(tmp)











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



# Copy global variables into local environments
x = 5
make.f =  function() {
  copy.into.env(source = globalenv(), names = ls.vars(globalenv()))
  f = function() x;
  return(f)
}
f = make.f()
f()
x = 10
f()


