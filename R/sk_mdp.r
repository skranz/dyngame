# Solving Markov Decision Processes using Policy Iteration with Value Determination
#
# Main function: policy_iteration(T,R,delta,use_val_iter=FALSE,oldp=NULL,V=NULL)
#
# T(s,a,s') = prob(s' | s, a)
# R(s,a)


# POLICY_ITERATION
# [new_policy, V, Q, niters] = policy_iteration(T, R, delta, use_val_iter, old_policy)
#
# If use_val_iter is not specified, we use value determination instead of value iteration.
# If the old_policy is not specified, we use an arbitrary initial policy.
#
# T(s,a,s') = prob(s' | s, a)
# R(s,a)
policy_iteration = function(T, R, delta, na, oldp=NULL, V = NULL, tol = 1e-8) {
  
  restore.point("policy_iteration")
  #stop()
  
  # na is vector with as many rows as states. na[x] is the number of
  # action profiles in state x
  if (! is(na,"VectorListInd"))
    na = VectorListInd(na)

  n = length(na)  
  nx = NROW(na)  
  Q = rep(-Inf,n);

  # Default policy is the one that maximizes immediate rewards
  if (is.null(oldp))
    oldp = which.RowMaxs(na,R)  # na describes the row indices
  
  p = oldp
  iter = 0;
  done = 0;
	
	sparse = !is.matrix(T)
	
		
  
  if (any(!is.finite(R[p]))) {
    stop("Even under optimal policy value is -Inf")
    return(NULL)
  }
  
  while (TRUE) {
    iter = iter + 1;
    
    # Determine Value
		if (!sparse) {
		  Tp = T[p,,drop=FALSE]
		} else {
			Tp = as.matrix(T[p,,drop=FALSE])
		}

    Tp = T[p,,drop=FALSE]
    Rp = R[p]
    
    # V = (1-delta)*R + delta*T*V  => (I-delta*T)V = (1-delta)*R  =>
    # V = inv(I-gT)*(1-delta)*R    
    V = solve(diag(nx)-delta*Tp, (1-delta)*Rp)
    # V = solve(Diagonal(nx)-delta*Tp, (1-delta)*Rp) #For sparse matrices

    # Get value for every (a,x) pair  
    oldQ = Q
		
		if (!sparse) {
			Q = (1-delta)*R + delta * T %*% V;
    } else {
		  Q = (1-delta)*R + delta * as.vector(T %*% V);
		}
    # Get optimal policy
    p = which.RowMaxs(na,Q, return.v.ind=TRUE)
    
    #print(V)
    #names(V) = names(p)= m$x.lab;
    #rbind(V,label.ax(m,p),oldV,label.ax(m,oldp))
    
    if (identical(p, oldp) | approxeq(Q, oldQ, tol)) {
      # if we just compare p and oldp, it might oscillate due to ties
      # However, it may converge faster than Q
      break();
    }
    oldp = p;
    oldQ = Q;
  }
  #print(paste(iter-1, " policy iterations"))
  return(list(p=p,V=as.numeric(V)))
}
    