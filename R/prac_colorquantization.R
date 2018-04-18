dyn.load("./fortran/bar.so")
.Fortran("bar", n=as.integer(5), x=as.double(rnorm(5)))
