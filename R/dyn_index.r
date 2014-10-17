# ##############################################################################################################
# Search functions to access certain action profiles
# ##############################################################################################################

# Finding optimal joint payoffs

# Notation: x physical state; z activity

# g will be an array:
# 1st dimension action profile a
# 2nd dimension player i
# 3rd dimension state x

# We denote by z an activity
# activities are


x.to.ax = function(m,x) {
	if (NROW(x)==1) {
		return(m$start.ax.by.x[x]:(m$start.ax.by.x[x+1]-1))
	} else {
		stop("x.to.ax so far implemented only for a single state x")
	}
}


ax.to.x = function(m,ax) {
	findInterval(ax,m$start.ax.by.x)
}

a.to.ax = function(m,a,x=1:m$nx) {
  m$start.ax.by.x + a -1
}
  
# Transform action indices to action labels
label.a = function(m,mat) {
  mat[] = m$a.lab[mat]
  mat
}

# Transform ax action indices to action labels    
label.ax = function(m,mat) {
  mat[] = m$ax.lab[mat]
  mat
}

# A unique key generated from a set of matrix columns
make.mult.col.key = function(mat,col=1:NCOL(mat),sep = "#") {
  key = mat[,col[1]]
  if (NROW(col) > 1) {
    for (i in 2:NROW(col))
      key = paste(key,mat[,col[i]],sep=sep)
  }
  return(key)
  #as.character(interaction(as.data.frame(mat[,col])))
}
   
# Match matrix m1 in matrix m2 using the specified columns as keys 
match.mat = function(m1,m2,col1=1:NCOL(m1),col2=1:NCOL(m2),key.char = "#", nomatch=NA) {
  
  restore.point("match.mat")

  if (!is.matrix(m1))
    m1 = t(m1)
  
  key1 = make.mult.col.key(m1,col1,key.char)
  key2 = make.mult.col.key(m2,col2,key.char)
  return(match(key1,key2,nomatch=nomatch))
}


# Translates a matrix of xv vectors into an index
xv.to.x = function(m,xvm) {
  if (!is.matrix(xvm)) {
    if (NCOL(m$xv.mat)>1) {
      xvm = matrix(xvm,1,NCOL(m$xv.mat))
    } else {
      xvm = matrix(xvm,NROW(xvm),1)
    }
  }
  
  j = NCOL(xvm)
  x = match(xvm[,j],m$xv.val[[j]])
  
  mult = m$xv.dim[j]
  if (j >= 2) { 
  	for (j in (NCOL(xvm)-1):1) {
    	xv.ind = match(xvm[,j],m$xv.val[[j]])
    	x = x+(xv.ind-1)*mult
    	mult = mult*m$xv.dim[j]
    }
  }
  return(x)
}

get.av.mat = function(m,ax=NULL,x=NULL,xv=NULL) {
	if (!is.null(ax)) {
		return(ax.to.av(m,ax,as.matrix=TRUE))
	} else if (!is.null(xv)) {
		stopifnot(is.matrix(xv))
		return(m$act.fun(xv)$val)
	} else {
		stopifnot(length(x)==1)
		return(m$act.fun(m$xv.mat[x,])$val)
	}
}

# Not very quick
ax.to.av = function(m,ax, as.matrix = TRUE) {
  
  restore.point("ax.to.av")
  
  av.list = list()
  xa = v.ind.to.rowcol(m$ind.ax.by.x,ax)
  old.av.val = NULL
  for (i in 1:NROW(xa)) {
    av.val = m$act.fun(m$xv.mat[xa[i,1],])$val
    if (!identical(av.val, old.av.val))
      av.mat = make.grid.matrix(av.val)
    old.av.val = av.val
    av.list[[i]] = av.mat[xa[i,2],]
  }
  vl = VectorList(av.list)
  if (as.matrix) {
    return(as.matrix(vl))
  } else {
    return(vl)
  }    
}


a.to.av = function(m,a) {
  return(m$av.mat[a,])
}


nice.ax = function(m,vec,a=NULL,x=NULL) {
  if (!is.null(a)) {
    if (!is.numeric(a)) {
      a  = match(a,m$a.lab)
    }
  }
  if (!is.null(x)) {
    if (!is.numeric(x)) {
      x  = match(x,m$x.lab)
    }
  }
  if (is.vector(vec) & is.null(a) & is.null(x)) {
    mat = matrix(vec,m$nx,m$na)
    rownames(mat) = m$x.lab
    colnames(mat) = m$a.lab
    return(mat)
  }
  
  if (is.vector(vec)) {
    names(vec) = m$ax.lab
    use = rep(TRUE,NROW(vec))
    if (!is.null(a))
      use = m$ax[,"a"] %in% a
    if (!is.null(x))
      use = use & (m$ax[,"x"] %in% x)
    return(vec[use])
  }  else if (is.matrix(vec)) {
    rownames(vec) = m$ax.lab
    use = rep(TRUE,NROW(vec))
    if (!is.null(a))
      use = m$ax[,"a"] %in% a
    if (!is.null(x))
      use = use & (m$ax[,"x"] %in% x)
    return(vec[use,])
  }    
}



  
sol.mat.e = function(m,sol=m$sol) {
  mat = matrix(unlist(sol),m$nx,(m$n+1)*2)
  colnames(mat)=c("ae",paste("a",1:m$n,sep=""),"U",paste("v",1:m$n,sep=""))
  
  mav = a.to.av(m,mat[,"ae"])
  mav = cbind(mav,mat[,"U"])
  rownames(mav) = m$x.lab
  colnames(mav) = c(m$av.lab,"U")
  mav
}

# Gives the indices of possible replies of player i holding fixed the actions of the other players
# Very slow implementation in the moment
get.replies = function(m,i,ax=1:m$nax,keep.ax = TRUE) {
	ax_i = m$ind.ax.to.ax_i[[i]][ax]
	act.repl = m$replies[[i]][rows.to.v.ind(m$ind.replies.by.ax_i[[i]],ax_i)]
	if (keep.ax) {
		return(act.repl)
	} else {
		stopifnot(NROW(ax)==1)
		return(setdiff(act.repl,ax))
	}
}  

get.union.of.replies = function(m,ax) {
	store.objects("get.union.of.replies")
	#restore.objects("get.union.of.replies")
	
  repl = NULL
	for (i in 1:m$n) {
		ax_i = m$ind.ax.to.ax_i[[i]][ax]
		act.repl = m$replies[[i]][rows.to.v.ind(m$ind.replies.by.ax_i[[i]],ax_i)]
		repl = c(repl,act.repl)
	}
	return(unique(repl))
}

# states that can be possibly reached from  state x given any ax
get.reachable.states.from.ax.replies = function(m,ax) {
	store.objects("get.reachable.states.from.ax.replies")
	#restore.objects("get.reachable.states.from.ax.replies")
	
	all.ax = get.union.of.replies(m,ax)
	tau.x  = as.matrix(m$tau[all.ax,,drop=FALSE])
	max.tau = colMaxs(tau.x)
	
	return(which(max.tau>0))
}


