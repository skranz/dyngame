# Investment into abatement reduction

# Two players cause negative externalities to each other, e.g. air pollution
# A player can invest into an abatement technology that reduces the negative
# externality on the other player. 
# Investments will have a permanent (long-lasting) effect

# States:  technology level of each player between 0 and x.max
# Actions: investment size between 0 and x.max: by how much shall 
#          the abatement technology be improved


# Perhaps associate to every action vector a set of players and a set of states
# That should yield a 

with.floor = function(mat,floor=0) {
  mat[mat<floor] = floor
  mat
}

with.ceiling = function(mat,ceiling) {
  mat[mat>ceiling] = ceiling
  mat
}


with.floor.and.ceiling = function(mat,floor=0,ceiling = 1) {
  mat[mat<floor] = floor
  mat[mat>ceiling] = ceiling
  mat
}



make.g.fun = function(kb,kd) {
  store.local.objects("make.g.fun");
  g.fun = function(avm,xvm) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  
    g = matrix(NA,NROW(avm),2) 
    g[,1] = -kb * (avm[,1] == 1) - kd * (avm[,1] == -1)
    g[,2] = xvm[,1]
    g
  }
  return(g.fun)
}

make.tau.fun = function(pi) {
  tau.fun = function(avm,xv,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")
    xv = xv[,1]
    avm = avm[,1]
  	
  	tau = matrix(0,NROW(avm),m$nx)
  	# Investment in state x=0: Move to state x=1
  	tau[avm == 1 & xv == 0,2] = 1
  	# No investment in state x=0: stay in to state x=0
  	tau[avm != 1 & xv == 0,1] = 1
  	# Remove investment in state x=1: Move to state x=0
  	tau[avm == -1 & xv == 1,1] = 1
  	# No disinvestment in state x=1: depreciation can take place
  	# Move to state x=0 with prob. pi and remain in state x=1 with prob 1-pi
  	tau[avm != -1 & xv == 1,1] = pi
  	tau[avm != -1 & xv == 1,2] = 1-pi

  	tau
  }
  return(tau.fun)
}

activity.fun = function(

STORE_OBJECTS = TRUE
reset.call.counts()

pi = 0.1
delta = 0.5

kb.0 = delta / ((1-delta)+(pi*delta))
kb.1 = delta / (1+(pi*delta))
kb.0
kb.1
U0 = c(0,((1-delta)) / ((1-delta+pi*delta)))
U1 = (1/(1+pi*d))*(
(d-(1-d)k_{b}-pdk_{b})
1-pdk_{b}
)
U1 = (1/(1+pi*delta)) * c(delta,1)
U0
U1

kb = 0.5
kd = 0.1

g.fun = make.g.fun(kb,kd)
tau.fun = make.tau.fun(pi=pi)

activities = list(
  val =   

av.val = list((-1):1, 0);
xv.val = list(0:1);
av.lab = av.val; xv.lab = xv.val;

m = init.game(n=2,g.fun=g.fun, tau.fun = tau.fun,
              av.val=av.val, av.lab = av.lab,
              xv.val=xv.val,xv.lab = xv.lab, 
              name="Coase", symmetric = TRUE)


ms = solve.game(m,delta=delta)  
ms

      
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