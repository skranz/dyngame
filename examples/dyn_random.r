# Cournot example with stochastic water reservers from the paper

# Cournot duopoly with stochastic water reserves. Example from papa
random.game = function(delta=2/3, numStates=3, numActions=2) 
{

  # act.fun specifies the action set for a given state xv
  # val = numerical values of actions
  # lab = labels of actions
  act.fun = function(xv,m=NULL) {
    restore.point("act.fun")
		val = list(1:numActions,1:numActions)
    list(val=val,lab=val, i=c(1,2))
  }
  
  # Stage game payoff function
  # Specifies the stage game payoffs as a function
  # of the action profile av and state xv
  # The function must be vectorized, i.e. avm and xvm
  # are matrices of action profiles and states
  g.fun = function(avm,xvm,m=NULL) {
    restore.point("g.fun");
		rownames(avm)=rownames(xvm)=NULL
		nc = 2; nr = NROW(avm)
		g = matrix(runif(nc*nr),nr,nc)
		g
  }

  # State transitions
  # For a matrix of action profiles and states specifies
  # the matrix of state transitions
  tau.fun = function(avm,xvm,m=NULL) {
  	#restore.point("tau.fun")

		rownames(avm)=rownames(xvm)=NULL
		nr = NROW(avm)
		nc = numStates
		tau = matrix(runif(nr*nc),nr,nc)
		tau = tau / rowSums(tau)
   	tau
  }
  
  # return required information for the dynamic game
  list(
    n=2,
    delta=delta,
    #integrated = FALSE,
    xv.val = 1:numStates,
    act.fun=act.fun,
    g.fun=g.fun,
    tau.fun=tau.fun
  )
}


# Run the inner code of this function manually
examples.random.game = function() {
  # Turn off restore points to speed up computation
  set.storing(!TRUE)

  numStates = 5
  numActions = 5
  # init game
  my.game = random.game(delta=0.7,numActions=numActions, numStates=numStates)
  m = init.game(my.game = my.game)
  # solve with transfers
  m.sol = solve.game(m)
  head(print.sol(m.sol))

  # solve game without transfers
  library(RSGSolve)
  setwd("D:/libraries/dyngame")
  rsg = dyngame.to.json.rsg(m, "random.json")
  rsg.sol = solveSG(rsg=rsg)
  set.storing(TRUE)
  plot.payoff.set(1, m.sol, rsg.sol, tight.approx = TRUE)
  plot.payoff.set(2, m.sol, rsg.sol, tight.approx = TRUE)
  plot.payoff.set(3, m.sol, rsg.sol, tight.approx = TRUE)
}
