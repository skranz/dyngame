
shorten.arrows = function(mat,x0=mat[,1],y0=mat[,2],x1=mat[,3],y1=mat[,4],
                          clip.x,clip.y) {
  # Length of arrow
  len = ((x1-x0)^2+(y1-y0)^2)^(1/2)
  
  #if (length(clip.x)==1) clip.x = c(clip.x,clip.x)
  #if (length(clip.y)==1) clip.y = c(clip.y,clip.y)
    
  # Relative clipping
  rel.clip.x = clip.x / len
  rel.clip.y = clip.y / len
  cx0 = x0+(x1-x0)*rel.clip.x
  cy0 = y0+(y1-y0)*rel.clip.y
  cx1 = x1-(x1-x0)*rel.clip.x
  cy1 = y1-(y1-y0)*rel.clip.y
  mat=cbind(cx0,cy0,cx1,cy1)
  colnames(mat)=c("x0","y0","x1","y1")
  return(mat)
}



state.levelplot = function(m,z=m$sol.mat[,"U"],name="", arrows = TRUE,
                           xlab="x1", ylab="x2",main=paste(name,m$name),xrange=NULL,
                           clip.x=NULL,clip.y=NULL,...) {
  
  restore.point("state.levelplot")

                             
  if (is.function(z)) 
    z = z(m)
    
  xv.mat = m$xv.mat
  if (!is.null(xrange)) {
		if (is.list(xrange)) {
			x.shown = which(xv.mat[,1]>=xrange[[1]][1] & xv.mat[,1]<=xrange[[1]][2] &
			                xv.mat[,2]>=xrange[[2]][1] & xv.mat[,2]<=xrange[[2]][2])
    } else {
			x.shown = which(xv.mat[,1]>=xrange[1] & xv.mat[,1]<=xrange[2])
		}
		xv.mat = xv.mat[x.shown,]
    z = z[x.shown]
  }
  
  sk.levelplot(grid.xyz = cbind(xv.mat,z),xlab=xlab,ylab=ylab,main=main,...)

  #sk.levelplot(grid.xyz = cbind(xv.mat,z),xlab=xlab,ylab=ylab,main=main)
  
  # Need to adapt for xrange
  if (arrows) {
    if (is.empty(m,"trans.e"))
      add.transitions(m)
      
    #colnames(m$trans.e) = c("x.from","x.to","prob")      
    mat = m$trans.e
    
    if (!is.null(xrange)) {
      # Remove transitions outside shown area
      arrow.shown = mat[,"x.from"] %in% x.shown & mat[,"x.to"] %in% x.shown
      mat = mat[arrow.shown,]
    }  
    
    xv.mat = m$xv.mat   
    if (NROW(mat)>0) {
      trellis.focus("panel", 1, 1)
      
 
      # How much will be cut from an arrow?
      if (is.null(clip.x)) {
        clip.x = abs((m$xv.val[[1]][2]-m$xv.val[[1]][1]) / 5)
      }
      if (is.null(clip.y)) {
        clip.y = abs((m$xv.val[[2]][2]-m$xv.val[[2]][1]) / 5)
      } 
      
      arr.mat = cbind(xv.mat[mat[,1],1],xv.mat[mat[,1],2],
                      xv.mat[mat[,2],1],xv.mat[mat[,2],2])
      # Shorten arrows      
      arr.mat = shorten.arrows(mat = arr.mat,clip.x=clip.x,clip.y=clip.y)
      
      larrows(x0=arr.mat[,1],y0=arr.mat[,2],
             x1=arr.mat[,3],y1=arr.mat[,4],
             lwd = mat[,3]*1,length=0.05 )
    }
  }
}

#xv=mm$ms[[1]]$xv.mat[1,];type="l";xlab=mm$par.name
plot.multi.model = function(mm,y,xv=mm$ms[[1]]$xv.mat[1,],type="l",xlab=mm$par.name,add=TRUE,...) {
  
  restore.point("plot.multi.model")
  
  m = mm$ms[[1]]
  if (is.character(y)) {
    y = mm$sol[,y]
  }
  x = xv.to.x(m,xv)
  rows = which(mm$sol[,"x"] == x)
  
  if (!add) {
    plot(mm$par.seq,y[rows],type=type,xlab=xlab,...)
  } else if (add) {
    lines(mm$par.seq,y[rows],type=type,...)
  }
}
