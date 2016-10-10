# Two state Prisoners' Dilemma game
# from Abreu, Brooks and Sannikov (2016)
abs.pd.game = function(delta=2/3) 
{

  # act.fun specifies the action set for a given state xv
  # val = numerical values of actions
  # lab = labels of actions
  act.fun = function(xv,m=NULL) {
    restore.point("act.fun")
		val = list(1:2,1:2)
		lab = list(c("C","D"),c("C","D"))
    list(val=val,lab=lab, i=c(1,2))
  }
  
  # Stage game payoff function
  # Specifies the stage game payoffs as a function
  # of the action profile av and state xv
  # The function must be vectorized, i.e. avm and xvm
  # are matrices of action profiles and states
  g.fun = function(avm,xvm,m=NULL) {
    restore.point("g.fun");
  	
		# Payoff matrices state 1
		g1 = matrix(c(2, 3, -1,0),2,2) # player 1
		g2 = t(g1) # player 2
		
		g = cbind(g1[avm], g2[avm])
		
		# In state 2 payoffs are just by 2 larger
		g[xvm==2] = g[xvm==2] + 2
		return(g)
  }

  # State transitions
  # For a matrix of action profiles and states specifies
  # the matrix of state transitions
  tau.fun = function(avm,xvm,m=NULL) {
  	#restore.point("tau.fun")

		rownames(avm)=rownames(xvm)=NULL
		# set default of 1/2
		
		tau = matrix(1/2,NROW(avm),m$nx)
  	## CC in state 1
    row = avm[,1] == 1 & avm[,2] == 1 & xvm == 1
  	tau[row,] = c(1/3, 2/3)
  	## DD in state 1
    row = avm[,1] == 2 & avm[,2] == 2 & xvm == 1
  	tau[row,] = c(1/3, 2/3)

  	## CC in state 2
    row = avm[,1] == 1 & avm[,2] == 1 & xvm == 2
  	tau[row,] = c(2/3, 1/3)
  	## DD in state 2
    row = avm[,1] == 2 & avm[,2] == 2 & xvm == 2
  	tau[row,] = c(2/3, 1/3)
   	tau
  }
  
  # return required information for the dynamic game
  list(n=2,
       delta=delta,
       xv.val = list(1:2), # states
       act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun
  )
}


# Run the inner code of this function manually
examples.abs.pd.game = function() {
  # Turn off restore points to speed up computation
  set.storing(TRUE)

  # init game
  my.game = abs.pd.game(delta=0.357)
  m = init.game(my.game = my.game)
  # solve game with transfers
  sol = solve.game(m)
  print.sol(sol)

  tighter.notransfer.payoff.set(sol)
    # solve game without transfers
  library(RSGSolve)
  setwd("D:/libraries/dyngame")
  rsg = dyngame.to.json.rsg(m, "pd_abs.json")
  rsg.sol = solveSG(rsg=rsg)

  plot.payoff.set(state=1, sol, rsg.sol,tight.approx = TRUE)
  plot.payoff.set(state=2, sol, rsg.sol,tight.approx = TRUE)
}
