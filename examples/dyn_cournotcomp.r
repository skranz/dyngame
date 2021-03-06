# Cournot example with stochastic water reservers from the paper

# Some constants used in the game specification

# type of reserve increment
DET = 1  # deterministic reserve increment
UNIF = 2 # uniformely distributed reserve increment

# type of solution
COLL = 1 # collusive solution (optimal dynamic game equilibrium)
MON = 2  # integrated monopoly solution


# a quicker translation function from states vectors to state indices
xv.to.x = function(m,xvm) {
  return( (xvm[,1])*m$xv.dim[[2]]+xvm[,2]+1)
}

# Cournot duopoly with stochastic water reserves. Example from paper
cournot.reserves.game = function(
           x.cap=20, K=20, x.inc=2, delta=2/3,inc.min=0,inc.max=2,
           sol.type=COLL,inc.type=DET,para=NULL) 
{
  
  # If parameters are given as a list just copy them all
  # into the local environment
  if (!is.null(para)) {
    copy.into.env(source=para)
  }
    
  # act.fun specifies the action set for a given state xv
  # val = numerical values of actions
  # lab = labels of actions
  act.fun = function(xv,m=NULL) {
    restore.point("act.fun")
		val = list(0:xv[1],0:xv[2])
		lab = val
		
    list(val=val,lab=lab, i=c(1,2))
  }
  
  # States can be grouped into sets of states with same
  # action sets.
  # In our game each state has a different action set
  x.group = function(xvm,m=NULL) {
    #x.group = rep(1,NROW(xvm))
    x.group = 1:NROW(xvm)
    x.group
  } 
	
  # Stage game payoff function
  # Specifies the stage game payoffs as a function
  # of the action profile av and state xv
  # The function must be vectorized, i.e. avm and xvm
  # are matrices of action profiles and states
  g.fun = function(avm,xvm,m=NULL) {
    restore.point("g.fun");
  	
		rownames(avm)=rownames(xvm)=NULL
  	
    # compute outputs and prices
		q = avm
  	P = K-q[,1]-q[,2]
  	P[P<0]=0
    # return profits
  	pi = P*q 
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
    
  # We do not consider a multistage game
  x.stage.fun = NULL
    
  # State transitions
  # For a matrix of action profiles and states specifies
  # the matrix of state transitions
  tau.fun = function(avm,xvm,m=NULL) {
  	#restore.point("tau.fun")

		rownames(avm)=rownames(xvm)=NULL
  	tau = matrix(0,NROW(avm),m$nx)  	  	
 
		# Firms reserves are restored by a fixed amount x.inc
		if (inc.type==DET) {
      
      # matrix of resulting state values
			new.xvm = with.floor.and.ceiling(xvm - avm + x.inc ,floor=0,ceiling=x.cap)

			# Translate state values to state indices
      new.x = xv.to.x(m,new.xvm)
      
      # Set transition probabilities to 1 for resulting states 
			ind.mat = cbind(1:NROW(xvm),new.x)
			tau[ind.mat]=1
			
		# Alternatively, firms reserves are restored by a uniformely 
    # distributed integer amount between 0 and x.inc
		} else if (inc.type==UNIF) {
      
      # loop through all combinations of independently
      # distributed restore amount draws
			for (x.add1 in inc.min:inc.max) {
				for (x.add2 in inc.min:inc.max) {
					x.add = c(x.add1,x.add2)				
					xvm1 = with.floor.and.ceiling(xvm[,1] - avm[,1] + x.add1 ,floor=0,ceiling=x.cap)
					xvm2 = with.floor.and.ceiling(xvm[,2] - avm[,2] + x.add2 ,floor=0,ceiling=x.cap)
				
					# Translate state values to state indices
          new.x = xv.to.x(m,cbind(xvm1,xvm2))
          
          # Add probability of restore draw to transition matrix
					ind.mat = cbind(1:NROW(xvm),new.x)
					tau[ind.mat]=tau[ind.mat] +1/((inc.max-inc.min+1)^2)
				}
			}	
		}
  	
   	tau
  }
  
  integrated = sol.type == MON
  
  # return required information for the dynamic game
  list(n=2,
       delta=delta,
       integrated = integrated,
       xv.val = list(0:x.cap,0:x.cap),
       act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.stage.fun = x.stage.fun,
       x.group=x.group,
       extra.sol.fun=extra.sol.fun)
}


# Run the inner code of this function manually
examples.cournot.reserves.game = function() {
  # Turn off restore points to speed up computation
  set.storing(FALSE)

  # Specify parameters 
  para = list(
    #delta=2/3,
    #K=20,
    #x.cap=20,
    #x.inc=3,
    #inc.min=3,
    #inc.max=4,
    delta=0.85,
    K=4,
    x.cap=4,
    x.inc=1,
    inc.min=0,
    inc.max=2,
    inc.type=UNIF,
    sol.type=COLL
  )
  
  # init game
  my.game = cournot.reserves.game(para=para)
  m = init.game(my.game = my.game)
  # solve game
  sol = solve.game(m)
  head(print.sol(sol))

  # solve game without transfers
  library(RSGSolve)
  setwd("D:/libraries/dyngame")
  rsg = dyngame.to.json.rsg(m, file="cournot.json")
  #rsg = loadJsonSG("cournot.json")
  rsg.sol = NULL
  rsg.sol = solveSG(delta=rsg$delta, states=rsg$states, all.iter = TRUE)

    
  set.storing(TRUE)
  plot.payoff.set(1, sol, rsg.sol, tight.approx = TRUE)
  plot.payoff.set(2, sol, rsg.sol, tight.approx = TRUE)
  plot.payoff.set(3, sol, rsg.sol, tight.approx = TRUE)
  plot.payoff.set(4, sol, rsg.sol, tight.approx = TRUE)
  

  # Analyse solution graphically
  par(mar=c(5,5,2,2))
  
  # Equilibrium prices
  state.levelplot(mc,z=mc$extra.sol.cur[,"P"],
                  arrows=FALSE,cuts=10,xrange=c(0,20),zlim=c(7,20),
                  main="Prices under Collusion",
                  xlab="Reserves firm 1", ylab="Reserves firm 2",
                  reverse.colors=TRUE)
  
  state.levelplot(mc,z=mc$extra.sol.cur[,"P"],
                  arrows=FALSE,cuts=10,xrange=c(0,20),zlim=c(8,20),
                  main="Prices under Collusion",
                  xlab="Reserves firm 1", ylab="Reserves firm 2",
                  reverse.colors=TRUE,col.scheme = "grey")
  
  # Punishment payoffs
  
  V = mc$sol.mat[,"v1"]+mc$sol.mat[,"v2"]
  state.levelplot(mc,z=V,arrows=FALSE,cuts=10,xrange=c(0,20),
                  main="Sum of punishment payoffs",
                  xlab="Reserves firm 1", ylab="Reserves firm 2",
                  reverse.colors=!TRUE,col.scheme = "grey")
                  
  V = mc$sol.mat[,"v1"]+mc$sol.mat[,"v2"]
  state.levelplot(mc,z=mc$sol.mat[,"v1"],arrows=FALSE,cuts=10,xrange=c(0,20),
                  main="Punishment payoffs player 1",
                  xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                  reverse.colors=!TRUE,col.scheme = "grey")
                                    
  # Monopoly solution
  mon = clone(mc)
  mon$integrated = TRUE
  mon = solve.game(mon,delta=delta)
  
  par(mar=c(5,5,2,0))
  state.levelplot(mon,z=mon$extra.sol.cur[,"P"],arrows=FALSE,cuts=10,xrange=c(0,20),
                  main="Prices under Monopoly",zlim=c(8,20),
                  xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                  reverse.colors=TRUE)
                  
   
 }
  
