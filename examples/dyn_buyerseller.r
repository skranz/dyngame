# A seller buyer model

# A period consists of a sequential game
#
# Player 1 (seller) can offer a service at a price p
# Player 2 (buyer) can accept or reject the offer
#          If rejected the period ends and players get a payoff of 0 in the actual period

# If accepted, the seller can choose low or high effort
#          costs for player 1 are kL and kH, respectively
#          gains for player 2 are vL and vH, respectively

# We assume kH>kL and vH>vL and kH-kL<vH-vL

# We assume p in [kL,vH]

# Transformation into a dynamic game
# States: x=-1 Beginning of a period, i.e. seller is going to make an offer
#         x>=0 Buyer faces an offer of size x
#         x=-3 Seller chooses effort

# Payoff adaption: p = p / delta  cost and v are divided by delta^2 

make.funs = function(p.grid,kL,kH,wL,wH) {
  # Copy global variables into local environment
  # copy.into.env(source = globalenv(), names = ls.vars(globalenv()))  
  store.local.objects() 
                      
  # xv will be one state vector
  act.fun = function(xv,m=NULL) {
    store.local.objects("stat.act.fun")
  	#restore.local.objects("stat.act.fun");restore.local.objects("make.funs")
  	
  	# Seller makes offer
  	if (xv == -1) {
    	val = list(p.grid,0)
    	lab = list(p.grid,"")
    # Buyer faces offer of p and can accept or reject
    } else if (xv >= 0) {
      val = list(0,c(0,1))
      lab = list("",c("A","R"))
    # Seller can decide on low or high effort
  	} else if (xv == -3) {
    	val = list(c(0,1),0)
    	lab = list(c("L","H"),"")
    # Stage 3 after offer has been rejected: nobody can do anything
    } else if (xv == -4) {
    	val = list(0,0)
    	lab = list("","")
  	} else {
    	stop(paste("unknown state xv = ",xv))
  	}
   
    list(val=val,lab=lab, i=1:2)
  }
  x.group = function(xvm,m=NULL) {
    x.group = rep(2,NROW(xvm))
    x.group[xvm==-1] = 1
    x.group[xvm==-3] = 3
    x.group[xvm==-4] = 4
    x.group
  }
  
  g.fun = function(avm,xvm,m=NULL) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  
    g = matrix(0,NROW(avm),2)
    
    # Price offers that are accepted
    rows = xvm[,1] > 0 & avm[,2] == 0
    g[rows,1] = xvm[rows,1] # Seller gets price
    g[rows,2] = -xvm[rows,1] # Buyer pays price
    
    # Seller that picks high quality
    rows = (xvm[,1] == -3) & avm[,1] == 1
    g[rows,1] = -kH 
    g[rows,2] = wH  

    # Seller that picks low quality
    rows = (xvm[,1] == -3) & avm[,1] == 0
    g[rows,1] = -kL 
    g[rows,2] = wL 

    g    
  }
  
  # How many stages have been before the actual state
  # for which there is no discounting between the stages
  x.stage.fun = function(xvm,m=NULL) {
    no.discount = rep(1,NROW(xvm))
    no.discount[xvm[,1]>=0]  = 2
    no.discount[xvm[,1]==-3] = 3 # Stage 3 if a contract is accepted
    no.discount[xvm[,1]==-4] = 3 # Stage 3 if a contract is rejected
    no.discount
  }
    
  # State transitions
  tau.fun = function(avm,xv,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")

  	tau = matrix(0,NROW(avm),m$nx)
  	
  	# From state -1 we move to the state that corresponds to the offered price
  	# Note that avm and xvm have the same values in this case (the offered price)
  	rows = which(xv[,1] == -1)
  	if (length(rows)>0) {
   	  x.dest = xv.to.x(m,avm[rows,1])
  	  tau[cbind(rows,x.dest)] = 1
	  }
	  
  	# From any offered price we move to state with xv = -3 if accepted
  	x.dest = xv.to.x(m,-3)
  	rows = (xv[,1] >= 0) & avm[,2] == 0
  	tau[rows,x.dest] = 1
  	
  	# From any offered price we move to state xv=-4 if rejected
  	x.dest = xv.to.x(m,-4)
  	rows = (xv[,1] >= 0) & avm[,2] == 1
  	tau[rows,x.dest] = 1
  	
  	# From state xv=-3 or xv=-4 we always move to xv=-1
  	x.dest = xv.to.x(m,-1)
  	rows = (xv[,1] == -3 | xv[,1] == -4)
  	tau[rows,x.dest] = 1

   	tau
  }
  
  list(act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.stage.fun = x.stage.fun,
       x.group=x.group)
}


STORE_OBJECTS = TRUE

kL = 0 # a sensible normalization
kH = 1

wL = 0.5
wH = wL + 2

p.grid = seq(kL,wH,length=11) # Offered prices can lie between kL and wH
funs = make.funs(p.grid=p.grid,kL=kL,kH=kH,wL=wL,wH=wH) 

xv.val = c(-1,-3,-4,p.grid)
xv.lab = c("St1","St3A","St3R",paste("p",round(p.grid,2),sep=""))

m = init.game(n=2, name="Buyer/Seller", symmetric = FALSE,nxv=1,
              xv.val=xv.val, xv.lab = xv.lab, functions=funs)
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


