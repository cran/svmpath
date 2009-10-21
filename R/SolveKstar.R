SolveKstar<-function(Kstar,eps=1e-10){
  onestar<-rep(1,ncol(Kstar))
  onestar[1]<-0
 # solve(Kstar+diag(length(onestar))*eps,onestar) #goes wrong
solve(Kstar,onestar)
}
