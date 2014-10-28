#' A VectorList shall be a vector representation of a list of vectors
#' The class allows subsetting as for a matrix 

# Need to rewrite the functions from the base package
NROW <- function(x) {
  d = dim(x)
  if (!is.null(d)) d[1] else length(x)
}
NCOL <- function(x) {
  c(dim(x),1,1)[2]
}

NCOLS <- function(vli) {
  vli[["ncols"]]
}


RowMaxs <- function(...) {
  UseMethod("RowMaxs")
}
RowMins <- function(...) {
  UseMethod("RowMins")
}

which.RowMaxs <- function(...) {
  UseMethod("which.RowMaxs")
}

which.RowMins <- function(...) {
  UseMethod("which.RowMins")
}

RowSums <- function(...) {
  UseMethod("RowSums")
}
ColSums <- function(...) {
  UseMethod("ColSums")
}
RowSums.default <- function(...) {
  rowSums(...)
}
ColSums <- function(...) {
  colSums(...)
}

VectorList = function(list=NULL,ncols=NULL,v=NULL,empty=NA,nrow=NULL,skel=NULL) {
  ("VectorList")
  #re("VectorList")
	
  vl = list()
  if (!is.null(list)) {
    vl$nrow = length(list)
    vl$ncols = sapply(list,length)
    vl$li = list
  } else if (!is.null(ncols)) {
    vl$nrow = NROW(ncols)
    vl$ncols = ncols
		if (is.null(skel)) {
			skel = lapply(1:vl$nrow, function(x) {1:vl$ncols[x]})
			skel[ncols==0] = replicate(sum(ncols==0),NULL)
		}
    if (is.null(v)) {
      vl$li = skel   
    } else {
      vl$li = relist(v,skel)
    }
  } else if (!is.null(nrow)) {
    vl$nrow = nrow
    vl$ncol = rep(0,nrow)
    vl$li = list()
  } else {
    stop("Cannot construct with given parameters")
  }
  vl$empty = empty
  class(vl) = c("VectorList","list")
  return(vl)
}

#' Makes a VectorListInd for the selected rows
make.child.VectorList <- function(vl,rows) {
  vlc = list()
	vlc$li = vl$li[rows]
	vlc$nrow = NROW(rows)
	vlc$ncols = vl$ncols[rows]
  vl$empty = vl$empty
  class(vlc) = class(vl)
  return(vlc)
}


as.vector.VectorList = function(vl) {
  unlist(vl$li)
}

unlist.VectorList = function(vl) {
  unlist(vl$li)
}

as.list.VectorList = function(vl) {
  vl$li
}


as.matrix.VectorList <- function(vl, empty=NA, rows = 1:vli$nrow) {
  mat = matrix(empty,NROW(vl),max(vl$ncols))
  for (i in 1:NROW(mat)) {    
    if (vl$ncols[i]>0) {
      mat[i,1:vl$ncols[i]] = vl$li[[i]]
    }
  }
  mat
}

set.vector = function(vl,vec) {
  vl$li = relist(vl$list,vec)
  vl
}

length.VectorList = function(x) sum(x[["ncols"]])
dim.VectorList <- function(x) c(x$nrow, max(x[["ncols"]]))



row.to.v.ind = function(vl,i) {
  return(vl$rowstart[i]:(vl$rowstart[i]+vl$ncols[i]))
}


"[.VectorList" <- function (vl,i,j,cols,empty,drop) {
  
  ijcommas = length(sys.call()) -3 - (!missing(cols)) - (!missing(empty)) - (!missing(drop))
  mycall = sys.call()
  #print(as.character(mycall))
  #print(ijcommas)
  
  if (missing(i) & missing(j) & missing(cols))  {
    # Called [,]
    if (ijcommas>0)
      return(as.matrix.VectorList(vl))
    # Called []
    return(unlist(vl$li))
  }  
  if (missing(drop)) drop = TRUE    
  if (missing(empty)) empty = vl$empty

  if (missing(j) & missing(cols)) {
    if (ijcommas == 0 | (NROW(i) == 1 & drop==TRUE)) {
      return(unlist(vl$li[i]))
    } else {
      return(as.matrix.VectorList(vl,rows=i))
    }
  }
  if (NROW(i) == 1 & drop == TRUE) {
    if (!missing(j)) {
      return(vl$li[[i]][j])
    } else {
      return(vl$li[[i]][cols])
    }
  }
  if (missing(i)) {
    i = 1:vl$nrow
  }  
  # I need to speed up the indexing performed below
  if (!missing(j)) {
    if (NROW(j)>1) {
      cols = matrix(j,NROW(i),NROW(j),byrow=TRUE)
    } else {
      cols = rep(j,NROW(i))
    }
  }
  
  # Have a vector of i and a vector or matrix of cols
  # Need to speed up the procedure below at best using c code
  vec = unlist(vl$li)
  rowstart = cumsum(vl$ncols)
  if (!is.matrix(cols)) {
    ret = rep(empty,NROW(cols))
    ind = cols + rowstart[i]
    ok  = which(cols<=vl$ncols[i])
    ret[ok] = vl$v[i[ok]]
  } else {
    ret = matrix(empty,NROW(cols),NCOL(cols))
    for (k in 1:NCOL(cols)) {
      ind = cols[,k] + rowstart[i]
      ok  = which(cols[,k]<=vl$ncols[i])
      ret[ok,k] = vec[i[ok],k]
    }
  }            
  return(ret)
}


"[<-.VectorList" <- function (vl,i,j,value,cols,empty,drop) {
  
  ijcommas = length(sys.call()) -3 - (!missing(cols)) - (!missing(empty)) - (!missing(drop))
  mycall = sys.call()
  
  if (missing(j)) {
    # Need to make this piece of code quicker...
    if (is.matrix(value)) {
      for (k in 1:length(i)) {
        vl$li[[i[k]]] <- value[k,]
      }
      vl$ncols[i] = NCOL(value)
    } else if (is.list(value)) {
      vl$li[i] = value
      vl$ncols[i] = sapply(value,length)
    } else if (is.vector(value)) {
      vl$li[i] = replicate(NROW(i),value,simplify=FALSE)
      vl$ncols[i] = length(value)
    }
    #vl$max.cols = max(vl$ncols)  
    
    return(vl)
  } else {
    stop("Assignment of VectorList not yet fully implemented")
 }
}

examples.VectorList = function() {
  x = VectorList(ncols=c(1,2,3))
  x
  x[1:2] = 1:5
  x
  x[2:3] = matrix(1:8,2,4)
  x
}


# Still quite slow...
as.matrix.VectorList = function(vl,empty=vl$empty,rows = 1:vl$nrow) {
  mat = matrix(empty,NROW(rows),max(vl$ncols))
  for (i in 1:length(rows)) {
    mat[i,1:vl$ncols[rows[i]]] = vl$li[[rows[i]]]
  }
  mat
} 


RowSums.VectorList <- function(vl,vec = vl[]) {
	C_RowSums_VectorList(vec,vl$ncols)
}
RowMaxs.VectorList <- function(vl,vec = vl[]) {
	C_RowMaxs_VectorList(vec,vl$ncols)
}



VectorListInd = function(ncols) {
  vli = list()
  
  vli$ncols = ncols
  vli$nrow = length(vli$ncols)
  class(vli) = c("VectorListInd","VectorList","list")  
  vli
}

length.VectorListInd <- function(vli) {
  return(sum(vli$ncols))
}
dim.VectorListInd <- function(x) c(x$nrow, max(x$ncols))

# Need to speed up this function
rows.to.v.ind = function(vli,rows) {
  restore.point("rows.to.v.ind")
  v.ind.start = c(1,cumsum(vli$ncols)[-vli$nrow]+1)
  unlist(sapply(rows, function(i) (1:vli$ncols[i])+v.ind.start[i]-1, simplify=FALSE))
}


# Need to speed up this function
rowcol.to.v.ind = function(vli,rows) {
  restore.point("rows.to.v.ind")
  v.ind.start = c(1,cumsum(vli$ncols)[-vli$nrow]+1)
  unlist(sapply(rows, function(i) (1:vli$ncols[i])+v.ind.start[i]-1, simplify=FALSE))
}

v.ind.to.rowcol = function(vli, ind) {
  v.ind.start = c(1,cumsum(vli$ncols)+1)
  row = findInterval(ind,v.ind.start)
  col = ind - v.ind.start[row]+1
  return(cbind(row,col))
}


v.ind.to.row = function(vli, ind) {
  v.ind.start = c(1,cumsum(vli$ncols)+1)
  row = findInterval(ind,v.ind.start)
  return(row)
}


vli = VectorListInd(c(2,5,3,5))
v.ind.to.rowcol(vli,c(2,4,9))

v.ind.to.cols = function(vli, ind, x = 1:NROW(ind)) {
  stopifnot(TRUE)
}

#' Makes a VectorListInd for the selected rows
make.child.VectorListInd <- function(vli,rows) {
  VectorListInd(vli[["ncols"]][rows])
}

#' Fills every row with a constant
const.rows = function(vli,row.val) {
  
  restore.point("const.rows")

  unlist(sapply(1:NROW(vli),function (i) rep(row.val[i],vli[["ncols"]][i]),simplify = FALSE))
}

#vli = VectorListInd(c(2,4,5))
#constant.rows(vli,c(3,2,1))

remove.elements = function(vli, v.ind, remove.empty.rows = FALSE) {
  
  restore.point("remove.elements")  
  
  if (length(v.ind)==0) 
    return(vli)
  
  stopifnot(!remove.empty.rows, min(v.ind)>=1, max(v.ind)<= length(vli))

  
  
  v.ind.start = c(1,cumsum(vli$ncols)[-vli$nrow]+1)
  rows.remove = findInterval(unique(v.ind),v.ind.start)
  vli$ncols = vli$ncols - tabulate(rows.remove,nbin = vli$nrow)
  return(vli)
}
 
# Gives a length(vli) * 2 matrix
# The second column gives the row that corresponds to the index in vli
ind.row.matrix <- function(vli,colnames=NULL) {
  
  restore.point("ind.row.matrix")
  
  f = function(k) (rep(k,vli$ncols[k]))
  mat = cbind(1:length(vli),unlist(sapply(1:NROW(vli),f,simplify=FALSE)))
  if (!is.null(colnames))
    colnames(mat)=colnames
  mat
}

as.matrix.VectorListInd <- function(vli, vec, empty=NA, rows = 1:vli$nrow) {
  ncol = max(vli$ncols)
  mat = matrix(empty,NROW(rows),max(vli$ncols))
  v.ind.start = c(1,cumsum(vli$ncols)+1)
  for (i in 1:length(rows)) {
    ri = rows[i]
    if (vli$ncols[ri]>0) {
      mat[i,1:vli$ncols[ri]] = vec[v.ind.start[ri]:(v.ind.start[ri+1]-1)]
    }
  }
  mat
}


RowSums.VectorListInd <- function(vl,vec = vl[]) {
  C_RowSums_VectorList(vec,vl$ncols)
}

RowMaxs.VectorListInd <- function(vl,vec = vl[]) {
  C_RowMaxs_VectorList(vec,vl$ncols)
}

RowMins.VectorListInd <- function(vl,vec = vl[]) {
  -C_RowMaxs_VectorList(-vec,vl$ncols)
}


which.RowMaxs.VectorListInd <- function(vl,vec = vl[],return.v.ind=TRUE) {
  stopifnot(return.v.ind)
  C_which_RowMaxs_VectorList(vec,vl$ncols)
}
which.RowMins.VectorListInd <- function(vl,vec = vl[],return.v.ind=TRUE) {
  stopifnot(return.v.ind)
  C_which_RowMaxs_VectorList(-vec,vl$ncols)
}



"[.VectorListInd" <- function (vl,i,j,empty,drop) {
  
  ijcommas = length(sys.call()) -3 - (!missing(empty)) - (!missing(drop))
  mycall = sys.call()
  
  #print(as.character(mycall))
  #print(ijcommas)
  if (missing(empty)) empty = NA
  if (missing(drop)) drop = TRUE    
   
  if (missing(i) & missing(j))  {
    # Called [,]
    if (ijcommas>0)
      return(as.matrix.VectorListInd(vl,vec=1:length(vl),empty=empty))
    # Called []
    return(1:length(vl))
  }  

  if (missing(j)) {
    if (ijcommas == 0 | (NROW(i) == 1 & drop==TRUE)) {
      return(rows.to.v.ind(vl,rows=i))
    } else {
      return(as.matrix.VectorListInd(vl,vec = rows.to.v.ind(vl,rows=i), rows=i,empty=empty))
    }
  }
  if (!missing(i) & !missing(j)) {
    if (NROW(i) == 1) {
      return(rows.to.v.ind(vl,rows=i)[j])
    }
  }
  stop("Column indexing not yet fully implemented for VectorListInd")
}

examples.VectorListInd = function() {
  vl = VectorListInd(c(1,2,3,4))
  vl[2]
  vl[2,1]
  
  x = sample(1:10,10)
  x
  as.matrix(vl,x)
  
  RowMins(vl,x)
  which.RowMins(vl,x)
  RowMins.VectorListInd(vl,x)
  which.RowMins.VectorListInd(vl,x)
    
  RowSums(vl,c(1:10))
}
