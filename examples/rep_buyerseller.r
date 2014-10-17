# Repeated Game Version of the seller buyer model

# The following simultaneous move game will be repeated in every period: period consists of a sequential game
#
# Player 1 (seller) can offer a service at a price p
#          and chooses low or high effort L or H
# Player 2 (buyer) specifies a minimum price pB


# Payoffs are: If p<=pB a contract is created that has the following features
#          buyer pays p to seller 
#          costs for player 1 are kL or kH (depending on effort)
#          gains for player 2 are vL or vH (depending on effort)

# We assume kH>kL and vH>vL and kH-kL<vH-vL
# We assume p in [kL,vH]

# Transformation into a dynamic game
# States: x=-1 Beginning of a period, i.e. seller is going to make an offer
#         x>=0 Buyer faces an offer of size x
#         x=-3 Seller chooses effort

# Payoff adaption: p = p / delta  cost and v are divided by delta^2 


##################################################################
# Section 6: Hotelling Games
##################################################################

library(repgames)
kL = 0 # a sensible normalization
kH = 1
wL = 1.5
wH = wL + 2
p.grid = seq(kL,wH,length=11) # Offered prices can lie between kL and wH

np = NROW(p.grid)
na1   = 2*NROW(p.grid)

pS.mat  = matrix(p.grid,np,np,byrow=!TRUE)
pB.mat  = matrix(p.grid,np,np,byrow=TRUE)
doTrade = pS.mat <= pB.mat 


gS  = rbind( (pS.mat - kL) * doTrade,
             (pS.mat - kH) * doTrade)
                 
gB  = rbind( (-pS.mat + wL) * doTrade,
             (-pS.mat + wH) * doTrade)
                 
colnames(gS)=colnames(gB) = round(p.grid,2)
rn = paste( c(rep("L",np),rep("H",np)), c(round(p.grid,2),round(p.grid,2)))
rownames(gS)=rownames(gB) = rn

mr = repgames::init.game(g1=gS,g2=gB,name="Seller Buyer")

# Solve the model 
mr = repgames::solve.game(mr, keep.only.opt.rows=TRUE)
mr

# Note that minmax payoffs (0,0) are also Nash equilibrium payoffs of the stage game!


# Let us specify the repeated game using the dyngame package


make.funs = function(p.grid,kL,kH,wL,wH) {
  # Copy global variables into local environment
  # copy.into.env(source = globalenv(), names = ls.vars(globalenv()))  
  store.local.objects() 
                      
  # xv will be one state vector
  act.fun = function(xv,m=NULL) {
    store.local.objects("stat.act.fun")
  	#restore.local.objects("stat.act.fun");restore.local.objects("make.funs")
  	val = list(p.grid,p.grid,c(0,1))
    lab = list(p.grid,p.grid,c("L","H"))   
    list(val=val,lab=lab, i=c(1,2,1))
  }
  x.group = 0
  
  g.fun = function(avm,xvm,m=NULL) {
    store.local.objects("g.fun");
  	#restore.local.objects("g.fun"); restore.local.objects("make.g.fun"); 
  
    g = matrix(0,NROW(avm),2)
    
    # Accepted contracts with low effort
    rows = which((avm[,1] <= avm[,2]) & avm[,3] == 0)   
    g[rows,1] =  avm[rows,1] -kL # Seller
    g[rows,2] = -avm[rows,1] +wL # Buyer 

        
    # Accepted contracts with high effort
    rows = which((avm[,1] <= avm[,2]) & avm[,3] == 1)   
    g[rows,1] =  avm[rows,1] -kH # Seller
    g[rows,2] = -avm[rows,1] +wH # Buyer 

    g    
  }
    
  # State transitions
  tau.fun = function(avm,xv,m=NULL) {
    store.local.objects("tau.fun")
  	#restore.local.objects("tau.fun")

  	tau = matrix(1,NROW(avm),m$nx)
    return(tau)
  }
  
  list(act.fun=act.fun,
       g.fun=g.fun,
       tau.fun=tau.fun,
       x.group=x.group)
}


STORE_OBJECTS = TRUE
funs = make.funs(p.grid=p.grid,kL=kL,kH=kH,wL=wL,wH=wH) 

xv.val = c(1)
xv.lab = c("x")

m = init.game(n=2, name="Buyer/Seller", symmetric = FALSE,nxv=1,
              xv.val=xv.val, xv.lab = xv.lab, functions=funs)
m$g
         
delta = 0.26
ms = solve.game(m,delta=delta)  
print.sol(ms)
