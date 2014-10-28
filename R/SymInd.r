# Store information about symmetry
# We have a vector and assign symmetry structures

SymInd = function(x.mat=NULL, perm=NULL,x.val = NULL, ix = NULL) {
  si = list()
  if (!is.null(perm)) {
    si$np = NCOL(perm)
    si$nv = NROW(perm)
    si$perm = perm
  } else {  
    if (!is.null(x.val)) {
      x.mat = make.grid.matrix(x.val)
    }
    if (!is.null(ix)) {
      si$np = max(ix)
      x.mat.old = x.mat
      x.mat = matrix(NA,NROW(x.mat),si$np)
      for (i in 1:si$np) {
        x.mat[,i] = paste.matrix.cols(x.mat.old,which(ix==i),sep="§")
      }
    }
    stopifnot(!is.null(x.mat))        
    
    si$np = NCOL(x.mat)
    si$nv = NROW(x.mat)
    
    si$perm  = matrix(NA,si$nv,si$np)
    si$perm[,1] = 1:si$nv
    key1 = paste(x.mat[,1],x.mat[,2],sep="\r")
    key2 = paste(x.mat[,2],x.mat[,1],sep="\r")
    si$perm[,2] = match(key2,key1)
  }
  stopifnot(si$np==2)

  si$strong.sym = which(si$perm[,1] == si$perm[,2])
  # Those indices that do first appear in perm[,1]
  si$distinct = which(si$perm[,1]<=si$perm[,2])
  class(si) = c("SymInd","list")
  si
}


as.matrix.SymInd = function(si,v, empty = NA, rows=NULL, distinct=is.null(rows)) {
  stopifnot(si$np==2)
  if (is.null(rows)) {
    stopifnot(NROW(v) == si$nv)
    if (distinct) {          
      return(cbind(v[si$distinct],v[si$perm[si$distinct,2]]))
    } else {
      return(cbind(v,v[si$perm[,2]]))
    }
  } else {
    stopifnot(!distinct)
    return(cbind(v,v[si$perm[rows,2]]))
  }
}    

get.sym.v = function(si,v,p.source=1,p.dest=2, rows=NULL) {
  stopifnot(si$np==2, is.null(rows))
  if (is.null(rows)) {
    return(v[si$perm[si$perm[,p.source],p.dest]])
  }
}

get.sym.ind = function(si,rows=1:si$nv,p.source=1,p.dest=2) {
  stopifnot(si$np==2)
  return(si$perm[si$perm[rows,p.source],p.dest])
}

get.sym.order = function(si) {
  restore.point("get.sym.order")
  stopifnot(si$np==2)
  dns = setdiff(si$distinct,si$strong.sym)
  ord = as.vector(rbind(si$perm[dns,1],si$perm[dns,2]))
  return(c(si$strong.sym,ord))
}
# 
# 
# mat = make.grid.matrix(list(0:1,0:1,0:1,0:1))
# si = SymInd(x.mat=mat,ix=c(1,1,2,2))
# as.matrix(si,1:16)
# get.sym.ind(si,1:16,1,2)
