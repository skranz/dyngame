#' We want to store action profiles in a parsimonious way
#' The underlying structure will be activity vectors
#' To every activity vector, we assign a set of players and states
#' Every assigned player can choose an activity from the vector in the assigned states

#' Every activity vector has a label and if it is assigned to multiple players
#' the player number will be added
 
make.Activities = function(av, x = list(), i = list(), av.names=NULL, av.lab = NULL) {
  ai = list()
  ai$av = av
  ai$x = x
  ai$i = i
  if (is.null(av.names))
    av.names = names(av)
  if (is.null(av.names))
    av.names = paste("av",1:NROW(av))
    
  names(ai$av) = av.names
  ai$av.names = av.names
  ai$av.lab = av.lab
  
      
  class(gi) = c("GridInd","list")

	return(gi)
}

make.action.mat = function() {
  
}

# x value to an index
x.to.x.ind = function(gi,x) {
  if (is.matrix(x)) {
    x.ind = matrix(n,nrow=NROW(x),ncol=NCOL(x))
    for (k in 1:gi$m)
      x.ind[,k] = match(x[,k],gi$val[[k]])
    return(x.ind)
  }
  x.ind = numeric(length(x))
  for (k in 1:gi$m)
    x.ind[k] = match(x[k],gi$val[[k]])
  return(x.ind)
}



# matrix of x.indeces to x.values
x.ind.to.x = function(gi,x.ind,k=1:gi$m) {
  if (is.matrix(x.ind)) {
    x = matrix(NA,nrow=NROW(x),ncol=NCOL(x))
    for (k.act in k)
      x[,k.act] = gi$val[[k.act]][x.ind]
  } else if (length(k)==1) {
    x = gi$val[[k]][x.ind]
  } else {
    x = numeric(length(x))
    for (i in 1:length(x)) {
      k.act = k[i]
      x[i] = gi$val[[k.act]][x.ind[i]]
    }
  }
  return(x.ind)
}


# x value to an index
x.to.x.ind = function(gi,x,k=1:gi$m) {
  if (is.matrix(x)) {
    x.ind = matrix(NA,nrow=NROW(x),ncol=NCOL(x))
    for (k.act in k)
      x.ind[,k.act] = match(x[,k.act],gi$val[[k.act]])
  } else if (length(k)==1) {
    x.ind = match(x,gi$val[[k]])
  } else {
    x.ind = numeric(length(x))
    for (i in 1:length(x)) {
      k.act = k[i]
      x.ind[i] = match(x[i],gi$val[[k.act]])
    }
  }
  return(x.ind)
}

x.to.v.ind = function(gi,x) {
  x.ind = x.to.x.ind(gi,x)
  x.ind.to.v.ind(gi,x.ind)
}

x.ind.to.v.ind = function(gi,x.ind) {
  if (is.matrix(x.ind)) {
  	v = rep(1,NROW(x.ind))
  	for (k in 1:gi$m) {
  		v = v + (x.ind[,k]-1)*gi$shift[k]
  	}
  	return(v)
  } else {
    return(sum((x.ind-1)*gi$shift)+1)
  }
}

v.ind.to.x.ind = function(gi,v.ind=1:gi$n, k=1:gi$m,as.matrix = length(k)>1) {
  
  restore.point("v.ind.to.x.ind")
  n.k = length(k)
  if (as.matrix) {
    x.ind = matrix(0,NROW(v.ind),n.k)
  	for (i in 1:n.k) {
  		x.ind[,i] = ceiling((((v.ind-1) %% gi$modulus[k[i]])+1) / gi$shift[k[i]])
  	}
  	return(x.ind)
  } else {
    if (n.k > 1) {
      x.ind = numeric(n.k)
  	  for (i in 1:n.k) {
    		x.ind[i] = ceiling((((v.ind-1) %% gi$modulus[k[i]])+1) / gi$shift[k[i]])
    	}
    	return(x.ind)
  	} else {
      x.ind = ceiling((((v.ind-1) %% gi$modulus[k])+1) / gi$shift[k])
    	return(x.ind)
    }      	
	}
}

v.ind.to.x = function(gi,v.ind=1:gi$n,k=1:gi$m, as.matrix = TRUE) {
  x.ind = v.ind.to.x.ind(gi,v.ind,k,as.matrix)
  x.ind.to.x(gi,x.ind,k)    
}

gi.x.ind.k = function(gi,k) {
  if (length(k)>1) 
    stop("only works for a single index k")
  x = 1:gi$dim[k]
  mydim = c(1,gi$dim,1,1)
  rep(rep(x,each=prod(mydim[(k+2):(gi$m+2)])), times = prod(mydim[1:k]))
}  


# Assume that only a subset of the m vectors will vary
make.replies.gi = function(gi, k.i) {
  k_i = (1:gi$m)[-k.i]
  cgi = make.child.gi(gi,k_i)
  cgi$k.i = k.i
  cgi$reply.mat = matrix(order(cgi$cp.ind),nrow=cgi$n,byrow=TRUE)
  return(cgi)
}


  

# Assume that only a subset of the m vectors will vary
make.child.gi = function(gi, k) {
  
  restore.point("make.child.gi")

  cgi = GridInd(gi$val[k])
  # Generate index with cgi$parent$n rows. Each element gives the corresponding index number in cgi
  # This means vec[row] will be the index in cgi where the columns k are the same as in gi
  cgi$cp.ind = rep(1,gi$n)
  for (i in 1:length(k)) {
    x.ind = gi.x.ind.k(gi,k[i])
    cgi$cp.ind = cgi$cp.ind + (x.ind-1) * cgi$shift[i]
  }
  cgi$parent = gi

  return(cgi)
}



# Some simple tests
# grid.matrix  = function(gi) {
#   make.grid.matrix(x=gi$val)
# }

# gi = GridInd(1:2,m=3)
# gi
# grid.matrix(gi)
# v.ind.to.x.ind(gi)

# gi.reply.mat(gi,k.var=c(1))
# cgi = make.child.gi(gi,k=c(2,3))
# grid.matrix(cgi)
# cbind(cgi$cp.ind,grid.matrix(gi))

# rep.gi = make.replies.gi(gi,k.i=c(1,2))
# rep.gi
# cbind(rep.gi$cp.ind,grid.matrix(gi))
