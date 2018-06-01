# dyngame: an R package to solve discounted  dynamic games with transfers

The R package `dyngame` can be used to solve discounted stochastic games with perfect monitoring in which monetary transfers can be used. The package computes the payoff set of all pure strategy subgame perfect equilibria and characterizes Pareto-optimal SPE. It implements the  policy elimination algorithm described in the paper "Discounted Stochastic Games with Voluntary Transfers" by Susanne Goldluecke and Sebastian Kranz.


The code below shows how to install the required R packages and contains the dynamic Cournot game example from the paper. Look in the example section in the packages github folder:

https://github.com/skranz/dyngame/tree/master/examples

for the source code files of the example.

Unfortunatley, the package documentation is still not well developed.



# Installation


To install the dyngame package run the following code in your R console:
```r
install.packages("dyngame",repos = c("https://skranz-repo.github.io/drat/",getOption("repos")))
```
Note that `dyngame` is not hosted on CRAN but on my own Github based repository. You therefore have to specify the `repos` argument as in the call above. 

If you use Windows, binaries are available in my repository for some R version (e.g. 3.3.x and 3.5.x), but not for all. If there is no binary version available and installation fails, you also have to install RTools for Windows to compile the package from source. See

[https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)


It is also possible to install the packages from Github and the dependencies manually from CRAN by pasting the following code (RTools will then always be necessary on Windows):

```r
if (!require(devtools))
  install.packages("devtools")

devtools::install_github("skranz/skUtils")
devtools::install_github("skranz/rgmpl")
devtools::install_github("skranz/dyngame")
```


# Cournot example with stochastic water reservers from the paper


## Game definition



```r
library(dyngame)

# Some constants used in the game specification

# type of reserve increment
DET = 1  # deterministic reserve increment
UNIF = 2  # uniformely distributed reserve increment

# type of solution
COLL = 1  # collusive solution (optimal dynamic game equilibrium)
MON = 2  # integrated monopoly solution


# a quicker translation function from states vectors to state indices
xv.to.x = function(m, xvm) {
    return((xvm[, 1]) * m$xv.dim[[2]] + xvm[, 2] + 1)
}

# Cournot duopoly with stochastic water reserves. Example from paper
cournot.reserves.game = function(x.cap = 20, K = 20, x.inc = 2, delta = 2/3, 
    inc.min = 0, inc.max = 2, sol.type = COLL, inc.type = DET, para = NULL) {
    
    # If parameters are given as a list just copy them all into the local
    # environment
    if (!is.null(para)) {
        copy.into.env(source = para)
    }
    
    # act.fun specifies the action set for a given state xv val = numerical
    # values of actions lab = labels of actions
    act.fun = function(xv, m = NULL) {
        restore.point("act.fun")
        val = list(0:xv[1], 0:xv[2])
        lab = val
        
        list(val = val, lab = lab, i = c(1, 2))
    }
    
    # States can be grouped into sets of states with same action sets. This will
    # speed up computations.  In our game each state has a different action set
    x.group = function(xvm, m = NULL) {
        # return a different group index for each state
        x.group = 1:NROW(xvm)
        x.group
    }
    
    # Stage game payoff function Specifies the stage game payoffs as a function
    # of the action profile av and state xv The function must be vectorized,
    # i.e. avm and xvm are matrices of action profiles and states
    g.fun = function(avm, xvm, m = NULL) {
        restore.point("g.fun")
        
        rownames(avm) = rownames(xvm) = NULL
        
        # compute outputs and prices
        q = avm
        P = K - q[, 1] - q[, 2]
        P[P < 0] = 0
        # return profits
        pi = P * q
        pi
    }
    
    # Some additional statistics of the solution that we might be interested in.
    # Here we want to know about price, total output, joint profits, consumer
    # surplus, and total welfare
    extra.sol.fun = function(avm, xvm, m = NULL) {
        q = avm
        P = K - q[, 1] - q[, 2]
        P[P < 0] = 0
        pi = P * q
        Q = q[, 1] + q[, 2]
        CS = (K - P) * pmin(Q, K) * (1/2)
        Pi = pi[, 1] + pi[, 2]
        W = CS + Pi
        
        
        mat = cbind(q, P, Q, CS, Pi, W)
        colnames(mat)[1:2] = c("q1", "q2")
        return(mat)
    }
    
    # We do not consider a multistage game
    x.stage.fun = NULL
    
    # State transitions For a matrix of action profiles and states specifies the
    # matrix of state transitions
    tau.fun = function(avm, xvm, m = NULL) {
        # restore.point('tau.fun')
        
        rownames(avm) = rownames(xvm) = NULL
        tau = matrix(0, NROW(avm), m$nx)
        
        # Firms reserves are restored by a fixed amount x.inc
        if (inc.type == DET) {
            
            # matrix of resulting state values
            new.xvm = with.floor.and.ceiling(xvm - avm + x.inc, floor = 0, ceiling = x.cap)
            
            # Translate state values to state indices
            new.x = xv.to.x(m, new.xvm)
            
            # Set transition probabilities to 1 for resulting states
            ind.mat = cbind(1:NROW(xvm), new.x)
            tau[ind.mat] = 1
            
            # Alternatively, firms reserves are restored by a uniformely distributed
            # integer amount between 0 and x.inc
        } else if (inc.type == UNIF) {
            
            # loop through all combinations of independently distributed restore amount
            # draws
            for (x.add1 in inc.min:inc.max) {
                for (x.add2 in inc.min:inc.max) {
                  x.add = c(x.add1, x.add2)
                  xvm1 = with.floor.and.ceiling(xvm[, 1] - avm[, 1] + x.add1, 
                    floor = 0, ceiling = x.cap)
                  xvm2 = with.floor.and.ceiling(xvm[, 2] - avm[, 2] + x.add2, 
                    floor = 0, ceiling = x.cap)
                  
                  # Translate state values to state indices
                  new.x = xv.to.x(m, cbind(xvm1, xvm2))
                  
                  # Add probability of restore draw to transition matrix
                  ind.mat = cbind(1:NROW(xvm), new.x)
                  tau[ind.mat] = tau[ind.mat] + 1/((inc.max - inc.min + 1)^2)
                }
            }
        }
        
        tau
    }
    
    integrated = sol.type == MON
    
    # return required information for the dynamic game
    list(n = 2, delta = delta, integrated = integrated, xv.val = list(0:x.cap, 
        0:x.cap), act.fun = act.fun, g.fun = g.fun, tau.fun = tau.fun, x.stage.fun = x.stage.fun, 
        x.group = x.group, extra.sol.fun = extra.sol.fun)
}
```


## Initializing and solving game

We now solve the game using the parameterizations in our paper.

```r
NO.LABELS = FALSE

# Turn off restore points to speed up computation
set.storing(FALSE)

# Specify parameters
para = list(delta = 2/3, K = 20, x.cap = 20, x.inc = 3, inc.min = 3, inc.max = 4, 
    inc.type = UNIF, sol.type = COLL)

# init game
my.game = cournot.reserves.game(para = para)
mc = init.game(my.game = my.game)

# solve for collusive solution
mc = solve.game(mc)
```

```
## 
## Game successfully solved for delta =  0.666666666666667
```

```r
# Solution table with one row per state
sol = print.sol(mc)
head(sol)
```

```
##   x      ae      a1      a2     U    v1    v2 q1 q2  P Q   CS Pi    W
## 1 1 0;0 0 0 0;0 0 0 0;0 0 0 60.33 30.17 30.17  0  0 20 0  0.0  0  0.0
## 2 2 0;1 0 1 0;1 0 1 0;1 0 1 66.67 30.17 36.50  0  1 19 1  0.5 19 19.5
## 3 3 0;2 0 2 0;2 0 2 0;2 0 2 72.33 30.17 42.17  0  2 18 2  2.0 36 38.0
## 4 4 0;3 0 3 0;3 0 3 0;3 0 3 77.33 30.17 47.17  0  3 17 3  4.5 51 55.5
## 5 5 0;4 0 4 0;4 0 4 0;4 0 4 81.67 30.17 51.50  0  4 16 4  8.0 64 72.0
## 6 6 0;5 0 5 0;5 0 5 0;5 0 5 85.33 30.17 55.17  0  5 15 5 12.5 75 87.5
##     U-V G.ae
## 1 -0.01    0
## 2  0.00   19
## 3 -0.01   36
## 4 -0.01   51
## 5  0.00   64
## 6 -0.01   75
```



## Analyse solution graphically

The code below performs the graphical analysis

```r
# Analyse solution graphically
par(mar=c(5,5,2,2))

# Equilibrium prices
state.levelplot(mc,z=mc$extra.sol.cur[,"P"],
                arrows=FALSE,cuts=10,xrange=c(0,20),zlim=c(7,20),
                main="Prices under Collusion",
                xlab="Reserves firm 1", ylab="Reserves firm 2",
                reverse.colors=TRUE)
```

```
## Loading required package: lattice
```

![](dyn_cournotreserves_files/figure-html/unnamed-chunk-4-1.png) 

```r
state.levelplot(mc,z=mc$extra.sol.cur[,"P"],
                arrows=FALSE,cuts=10,xrange=c(0,20),zlim=c(8,20),
                main="Prices under Collusion",
                xlab="Reserves firm 1", ylab="Reserves firm 2",
                reverse.colors=TRUE,col.scheme = "grey")
```

![](dyn_cournotreserves_files/figure-html/unnamed-chunk-4-2.png) 

```r
# Punishment payoffs

V = mc$sol.mat[,"v1"]+mc$sol.mat[,"v2"]
state.levelplot(mc,z=V,arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Sum of punishment payoffs",
                xlab="Reserves firm 1", ylab="Reserves firm 2",
                reverse.colors=!TRUE,col.scheme = "grey")
```

![](dyn_cournotreserves_files/figure-html/unnamed-chunk-4-3.png) 

```r
V = mc$sol.mat[,"v1"]+mc$sol.mat[,"v2"]
state.levelplot(mc,z=mc$sol.mat[,"v1"],arrows=FALSE,cuts=10,xrange=c(0,20),
                main="Punishment payoffs player 1",
                xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                reverse.colors=!TRUE,col.scheme = "grey")
```

![](dyn_cournotreserves_files/figure-html/unnamed-chunk-4-4.png) 

## Monopoly solution


```r
  # Monopoly solution
  mon = clone(mc)
  mon$integrated = TRUE
  mon = solve.game(mon,delta=delta)
  
  par(mar=c(5,5,2,0))
  state.levelplot(mon,z=mon$extra.sol.cur[,"P"],arrows=FALSE,cuts=10,xrange=c(0,20),
                  main="Prices under Monopoly",zlim=c(8,20),
                  xlab="Oil reserves firm 1", ylab="Oil reserves firm 2",
                  reverse.colors=TRUE)
```

                  
   
 }
  
