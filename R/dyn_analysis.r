
add.transitions = function(m) {
  
  restore.point("add.transitions")
  m$tau.e = m$tau[m$sol.mat[,"ae"],]
  
  vec = as.vector(t(m$tau.e))
  gi = GridInd(dim=dim(m$tau.e))
  rows = which(vec != 0)
  
  if (length(rows>0)) {
    mat = cbind(v.ind.to.x.ind(gi,rows),vec[rows])
    mat = mat[mat[,1] != mat[,2],,drop = FALSE]
    colnames(mat) = c("x.from","x.to","prob")
    rownames(mat) = paste(m$x.lab[mat[,1]],m$x.lab[mat[,2]])
    
    m$trans.e = mat
  } else {
    m$trans.e = matrix(NA,0,3)
    colnames(m$trans.e) = c("x.from","x.to","prob")
  }
}
  
get.x.long.run = function(m) {
  tau.e = m$tau[m$sol.mat[,"ax.e"],]
  prob = solve(diag(m$nx) - tau.e)
  return(prob)
}

# Simulated probabilities
get.average.discounted.prob = function(m=NULL,delta = m$delta,ax=m$sol.mat[,"ae"], tau = m$tau) {
	store.objects("get.average.discounted.prob")
	#restore.objects("get.average.discounted.prob")
																											
  tau.e = m$tau[ax,]  
  prob = (1-delta)*solve(diag(NROW(tau.e)) - delta*tau.e)
  return(prob)
}  


replace.ax.with.ax.lab = function(m,mat) {
  
  restore.point("replace.ax.with.ax.lab")
  
  if (!is.data.frame(mat)) {
    df = data.frame(mat)
  } else {
    df = mat
  }
  cols = intersect(colnames(df),c("ae","a1","a2","a3","a4","a5"))
  for (k in cols) { 
    df[,k] = m$ax.lab[df[,k]]
  }
  df
}

print.sol = function(m,stack.sym.x = TRUE,order=NULL, round=2, x=NULL, digits=2,extra.sol=TRUE) {
  
  restore.point("print.sol")
  
  
  if ("multimodel" %in% class(m)) {
    mat = m$sol
    
    if (!is.null(digits)) mat = round(mat,digits)
    if (!is.null(x)) {
      mat = mat[mat[,"x"] %in% x,]
    }
    if (!is.null(order)) {
      if (identical(order,"x")) {
        order = order(mat[,"x"],mat[,"model"])
      }
      mat = mat[order,]
    }
    df = replace.ax.with.ax.lab(mm$ms[[1]],mat)
    return(df)
  }
    
  mat = cbind(m$sol.mat, m$extra.sol.cur)
  if (!is.null(digits)) mat = round(mat,digits)
  
  if ("a1" %in% colnames(mat)) {
    v.col = paste("v",1:m$n,sep="")
    add = cbind(mat[,"U"]-rowSums(mat[,v.col]),
                m$G[mat[,"ae"]])
    colnames(add) = c("U-V","G.ae")
    df = as.data.frame(round(cbind(mat,add),2))
    rownames(df) = NULL #rownames(m$sol.mat)
  } else {
    df = as.data.frame(mat)
  }
  df = replace.ax.with.ax.lab(m,df)  
 
  if (!is.null(order)) {
    return(df[order,])
  }
  if (stack.sym.x & m$symmetric) {
    order = rev(get.sym.order(m$x.sym.ind))
    return(df[order,])
  } 
  return(df)
}

