clone = function (...) UseMethod("clone")
clone.default = function(x) x


names.ListEnv = function(le,rev=!TRUE) {
  if (rev) {
    rev(ls(le))
  } else {
    ls(le)
  }
}

print.ListEnv = function(le) {
  NextMethod()
  print(as.list(le))
}

unlist.ListEnv = function(le,recursive = TRUE, use.names = TRUE) {
  unlist(as.list(le))
}

new.ListEnv = function(hash = TRUE, parent = emptyenv(), size = 29L) {
  le = new.env(hash,parent,size)
  class(le) = c("ListEnv","environment")
  le
}

#' Check if a particuar item does not exist or is NULL
is.empty = function (...) UseMethod("is.empty")

is.empty.ListEnv = function(le,name, inherits = TRUE) {
  if (!exists(name, envir=le, inherits = inherits))
    return(TRUE)
  if (is.null(le[[name]]))
    return(TRUE)
  return(FALSE)
}

"$.ListEnv" <- function (x,i) {
  return(get(i,envir=x))
}


"[[.ListEnv" <- function (x,i) {
  if (!is.numeric(i)) {
    return(get(i,envir=x))
  } else {
    return(get(names(x)[i],envir=x))    
  }
}

"[.ListEnv" <- function (x,i) {
  if (is.numeric(i)) {
    i = names(x)[i]
  }
  as.list(x)[i]  
}


#'
summary.ListEnv = function(le) {
  ls.str(le)
}

#' Create a new ListEnv
#' A ListEnv is an environment, in particular it is copied by reference
#' Stil it behaves in some aspects like a list, e.g. there is a function names
#' or we can use unlist
ListEnv = function(...) {
  le = new.ListEnv(hash = TRUE, parent = parent.frame(), size = 29L)
  li = list(...)
  names = names(li)
  for (n in names) {
      le[[n]] <- li[[n]]
  }
  le
}


as.ListEnv = function(x) {
  if (is.list(x)) {
    le = new.ListEnv(hash = TRUE, parent = parent.frame(), size = 29L)
    names = names(x)
    for (n in names) {
        le[[n]] <- x[[n]]
    }
    return(le)
  } else if (is.environment(x)) {
    class(x) = c("ListEnv","environment")
    return(x)
  }
  stop("Can only convert lists or environment to ListEnv!")
}

examples.ListEnv = function() {
  env = new.env()
  env$a = "Hi"
  env$b = 1:10
  env[["a"]]
  env$a
  env["a"]
  names(env)
  
  le = as.ListEnv(env)

  le$a = "Hi"
  le$b = 1:10
  le[["a"]]
  le$a
  le["a"]
  le[1:2]
  le[[2]]
  names(env)
}