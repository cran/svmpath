".First.lib" <-
  function (lib, pkg) 
  library.dynam("svmpath", pkg, lib)
Balanced.Initialization<-function(K,y,nplus,nminus=nplus,eps=1e-12){
  n<-2*nplus
  f<-K%*%y
  Iplus<-seq(y)[y>0]
  Iminus<-seq(y)[y<0]
  fmax<-max(f[Iplus])
  fmin<-min(f[Iminus])
  iplus<-Iplus[match(fmax,f[Iplus],0)]
  iminus<-Iminus[match(fmin,f[Iminus],0)]
###this seems to take into account ties
  lambda<- (fmax-fmin)/2
  beta0<-1-fmax/lambda
  alpha0<-beta0*lambda
### package parameters for the left of the start
  alpha00<-c(slope=beta0,intercept=0)
  
  list(Elbow=c(iplus,iminus),lambda=lambda,alpha0=alpha0,alpha00=alpha00,alpha=rep(1,n))
}
  
coef.svmpath<-function(object,lambda,...){
  if(missing(lambda)){
    alpha<-object$alpha
    lambda<-object$lambda
    alpha0<-object$alpha0
    }
  else{
    alphs<-predict(object,lambda=lambda,type="alpha")
    alpha<-alphs$alpha
    alpha0<-alphs$alpha0
    }
  alpha<-alpha*object$y
  beta<-scale(t(object$x)%*%alpha,FALSE,lambda)
  beta0<-alpha0/lambda
  list(beta=beta,beta0=beta0,lambda=lambda)
}
  
DowndateKstar<-function(Kstar,index){
  index<-index+1
  Kstar[-index,-index,drop=FALSE]
}
"enlist" <-
function(...)
{
	result <- list(...)
	if((nargs() == 1.) & is.character(n <- result[[1.]])) {
		result <- as.list(seq(n))
		names(result) <- n
		for(i in n)
			result[[i]] <- get(i)
	}
	else {
		junk <- sys.call()
		n <- NULL
		for(i in junk[-1.])
			n <- c(n, deparse(i))
		if(!is.null(n2 <- names(result))) {
			which <- n2 != ""
			n[which] <- n2[which]
		}
		names(result) <- n
	}
	result
}
get.svmstep<-function(lambda,lambda.old){
      lambda.range<-range(c(lambda,lambda.old))
      lambda.cuts<-c(lambda.range[2]+1,lambda.old,lambda.range[1]-1)
      as.numeric(cut(-lambda,-lambda.cuts))-1
    }
      
"modulus" <-
  function(x, n)
  n * (x/n - trunc(x/n))
OptInit.alpha<-function(mbig,msmall,Rmat,c.objective){
  a<-rep(1,mbig)
  bl<-c(rep(0,mbig),msmall)
  bu<-c(rep(1,mbig),msmall)
  nclin<-1
  istate <- rep(0, mbig + nclin)
  storage.mode(a)<-"double"
  storage.mode(bl)<-"double"
  storage.mode(bu)<-"double"
  storage.mode(c.objective)<-"double"
  storage.mode(istate)<-"integer"
  storage.mode(Rmat)<-"double"
  lenw<-2*mbig*mbig+10*mbig+6
  alpha<-rep(1,mbig)*msmall/mbig
  storage.mode(alpha)<-"double"
  fit<-.Fortran("lssol",
              as.integer(mbig),
              as.integer(mbig),
              as.integer(nclin),
              as.integer(nclin),
              as.integer(mbig),
              a,
              bl,
              bu,
              c.objective,
              istate=istate,
              integer(mbig),
              alpha=alpha,
              Rmat,
              double(mbig),
              inform=integer(1),
              iter=integer(1),
              obj=double(1),
              clambda=double(mbig+1),
              iw=integer(mbig),
              as.integer(mbig),
              double(lenw),
              as.integer(lenw),
              PACKAGE="svmpath"
              )
  list(alpha=fit$alpha,obj=fit$obj)
}
"poly.kernel" <-
  function(x, y=x, param.kernel = 1)
{
  if(is.null(param.kernel))
    param.kernel <- 1
  if(param.kernel == 1)
    x %*% t(y)
  else (x %*% t(y) + 1)^param.kernel
}
"predict.svmpath" <-
  function(object,newx,lambda,type=c("function","class","alpha"),...){
    type<-match.arg(type)
    oalpha<-object$alpha
    oalpha0<-object$alpha0
    if(missing(lambda)){
      lambda<-object$lambda
      alpha<-oalpha
      alpha0<-oalpha0
    }
    else{
      olambda<-object$lambda
      nalpha<-length(object$y)
      minl<-min(olambda);maxl<-max(olambda)
      anysmaller<-seq(along=lambda)[lambda<minl]
      if(length(anysmaller))lambda[anysmaller]<-minl
      lmax<-max(lambda)
      if(lmax>maxl){# we need to modify our alphas
        alpha00<-object$alpha00
        maxl<-lmax
        oalpha<-cbind(oalpha[,1],oalpha)
        oalpha0<-c(alpha00["slope"]*lmax+alpha00["intercept"],oalpha0)
        olambda<-c(maxl,olambda)
      }
      lfrac<-(lambda-minl)/(maxl-minl)
      olfrac<-(olambda-minl)/(maxl-minl)
      coord<-approx(olfrac,seq(olambda),lfrac)$y
      left<-floor(coord);right<-ceiling(coord)
      alpha <- outer(rep(1,nalpha),olfrac[right] - lfrac) * oalpha[,left , drop = FALSE] + 
                outer(rep(1,nalpha),lfrac - olfrac[left]) * oalpha[,right , drop = FALSE]
      alpha<-scale(alpha,  FALSE, olfrac[right] - olfrac[left])
      alpha[,left == right] <- oalpha[,left[left == right] ]
      alpha0 <- ((olfrac[right] - lfrac) * oalpha0[left] + 
                 (lfrac - olfrac[left]) * oalpha0[right])/(olfrac[right] - 
                                                           olfrac[left])
      alpha0[left == right] <- oalpha0[left[left == right] ]
      }
    if(type=="alpha"){
      attr(alpha,"scaled:scale")<-NULL
      object<-list(alpha0=alpha0,alpha=drop(alpha),lambda=lambda)
    }
    else{
      if(missing(newx))newx<-object$x
      K<-object$kernel(newx,object$x,object$param.kernel)
      fit<-K%*%(alpha*object$y)
      fit<- scale(fit,-alpha0,lambda)
      attr(fit,"scaled:center")<-NULL
      attr(fit,"scaled:scale")<-NULL
    }
    switch(type,
           "function"=fit,
           "class"=sign(fit),
           alpha=object
           )
  }
PrintPath<-function(trace,step,obs,moveto,movefrom,lambda,digits,stats){
  if(trace){
    moveto<-rep(moveto,length=length(obs))
    movefrom<-rep(movefrom,length=length(obs))
  for(i in seq(along=obs)){
    cat(step,":\tObs ",obs[i],"\t",movefrom[i],"->",moveto[i],"  lambda = ",format(lambda, digits=digits,nsmall=digits),"  Margin = ",format(round(stats$margin,2))," Elbow = ",stats$selbow," Error = ",stats$error,"\n",sep="")
  }
}
  invisible()
}
print.svmpath<-function(x,digits=6,...){
  for( i in seq(x$Step)){
    step<-x$Step[i]
    stats<-list(selbow=x$Size.Elbow[step],error=x$Error[step],margin=x$SumEps[step])
PrintPath(TRUE,x$Step[i],x$Obs.step[i],x$Moveto[i],x$Movefrom[i],x$lambda[step],digits,stats)
  }
  invisible()
}
  
"radial.kernel" <-
  function(x, y=x, param.kernel = 1)
{

  ###Note param.kernel is now gamma 
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
  exp( - a* param.kernel)
}
"SnapPath" <-
  function(step,x, y,f, alpha, alpha0, lambda, Elbow, kernel.function,param.kernel, linear.plot,Size=60,movie=FALSE,movie.root="./",...)
{
  if(movie)jpeg(file=paste(movie.root,step,".jpg",sep=""),quality=90,width=540,height=540,bg="wheat",...)
  stats<-StatPath(y,f,Elbow)
###only for 2 dim inputs
  x <- x[, 1:2]
  n<-length(y)
    ss<-abs(diff(range(x[,2])))
    soffset<-x*0
    soffset[,2]<-ss/50
    x <- x[, 1:2]
if(!linear.plot){
        xr <- apply(x, 2, range)
      xg <- apply(xr, 2, function(x, Size)
              seq(from = x[1], to = x[2], length = Size), Size)
      xG <- data.matrix(expand.grid(as.list(data.frame(xg))))
      Knew <- kernel.function(x, xG, param.kernel=param.kernel,...)
      fhat <- ((alpha * y) %*% Knew + alpha0)/lambda
      fhat <- matrix(fhat, Size, Size)
      }

    plot(x[,  ], type = "n",xlab="X1",ylab="X2")
   title(paste("              Step: ",format(step,digits=3),"   Error:",format(round(stats$error),digits=3), "   Elbow Size:",format(round(stats$selbow),digits=2),"   Margin:",format(round(stats$margin,2),digits=7)),adj=0)
    pointid <- seq(y)
    points(x[y == 1,  ], col = 2,pch="*",cex=2)
    points(x[y == -1,  ], col = 4,pch="*",cex=2)

    if(n<15)text((x-soffset)[y == 1,  ], labels = paste(pointid[y == 1]), col = 2)
    if(n<15)text((x-soffset)[y == -1,  ], labels = paste(pointid[y == -1]), col = 4)
  
    if(n<15&&length(Elbow))text((x-soffset)[Elbow,  ], labels = paste(pointid[Elbow]), col = 3)
  if(linear.plot){
    beta <- (alpha * y) %*% x
    abline( - alpha0/beta[2],  - beta[1]/beta[2], col = 3, lwd = 3)
    abline(lambda/beta[2] - alpha0/beta[2],  - beta[1]/beta[2], col = 3, 
           lwd = 2, lty = 3)
    abline( - lambda/beta[2] - alpha0/beta[2],  - beta[1]/beta[2], col = 3, 
           lwd = 2, lty = 3)
  }
  else{
      contour(xg[, 1], xg[, 2], fhat, levels = 0, add = TRUE, labels = "",col=3,lwd=3)
      contour(xg[, 1], xg[, 2], fhat, levels = c(-1,1), add = TRUE, labels = c("",""),col=3,lwd=1)
    }
  if(movie)dev.off()
    invisible()
    
}
SolveKstar<-function(Kstar,eps=.0001){
  onestar<-rep(1,ncol(Kstar))
  onestar[1]<-0
#  solve(Kstar+diag(length(onestar))*eps,onestar)
solve(Kstar,onestar)
}
StatPath<-function(f,y,Elbow){
    yhat<-y*f
  error<-sum(yhat<0)
  margin<-sum((1-yhat)[yhat<1])
  selbow<-length(Elbow)
list(error=error,margin=margin,selbow=selbow)
  }
summary.svmpath<-function(object,nsteps=5,digits=6,...){
  cat("\nCall:\n")
  cat(paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
  N<-length(object$Error)
  nsupp<-apply(object$alpha,2,function(x)sum(x>0))
  m<-sort(unique(c(pretty(seq(N),nsteps),1,N)))
  m<-seq(N)[match(m,seq(N),0)]
  cat("Number of steps:",N,"\n")
  cat("Selected steps:\n")
  data.frame(Lambda=round(object$lambda[m],digits),Error=object$Error[m],Size.Elbow=object$Size.Elbow[m],Support=nsupp[m],SumEps=object$SumEps[m],row.names=m)
  

}
"svmpath" <-
  function(x, y, K = kernel.function(x, x,param.kernel=param.kernel), kernel.function=poly.kernel,param.kernel=1, trace = FALSE, plot.it = FALSE, linear.plot=missing(kernel.function),eps = 
           1e-10, Nmoves = 3 * n, digits=6,lambda.min=1e-4,...)
{
### Function to compute the entire SVM path of solutions for a two-class
### classification problem
### Copyright Trevor Hastie, May 2003
###
### Initializations
### y must be -1 and 1
###
  this.call<-match.call()
  if(plot.it&&(ncol(x) > 2))stop("Plotting only for 2-dim X")
  n <- length(y)
  yvals <- table(y)
  if(length(yvals) != 2)
    stop("SvmPath works with binary problems only")
  nplus <- yvals[2]
  nminus <- yvals[1]	### Initialize the sets of observations
  Right <- Elbow <- NULL
  Left <- seq(n)
  Kscript <- K * outer(y, y)	
### We start with a maximum of 2*n moves, but these can be increased
###
### Initializations of counters
  alpha <- matrix(1, n, Nmoves)
  alpha0 <- double(Nmoves)
  SumEps <- double(Nmoves)
  Elbow.list<-as.list(seq(Nmoves))
  Size.Elbow<-integer(Nmoves)
  Error<-integer(Nmoves)
  Step<-integer(2*Nmoves)
  Obs.step<-integer(2*Nmoves)
  Movefrom<-character(2*Nmoves)
  Moveto<-character(2*Nmoves)
  lambda <- double(Nmoves)
  Kstar<-matrix(0,1,1)
###Initialization of path
### Two cases nplus=nminus or else not
  if(nplus == nminus){
    init<-Balanced.Initialization(K, y,nplus,nminus)
    Elbow <- init$Elbow
    Left <- setdiff(Left, Elbow)
  }
  else{
    init<-Unbalanced.Initialization(K, y, nplus,nminus)
    Elbow<-init$Elbow
    Right<-init$Right
    Left<-init$Left
  }
  Elbow.list[[1]]<-Elbow
  lambda0 <- init$lambda
  alpha0[1] <- init$alpha0
  alpha[,1] <-init$alpha
  alpha00<-init$alpha00 #safekeeping
  Kstar<-UpdateKstar(Kstar,Kscript[Elbow,Elbow],NULL,y[Elbow])
  lambda[1] <- lambda0
  fl <- (K %*% (alpha[, 1] * y) +init$alpha0)/lambda0
  stats<-StatPath(y,fl,Elbow)
  SumEps[1]<-stats$margin
  Error[1]<-stats$error
  Size.Elbow[1]<-stats$selbow
  nobs<-seq(along=Elbow)
  Step[nobs]<-1
  Obs.step[nobs]<-Elbow
  Movefrom[nobs]<-" ";Moveto[nobs]<-"E"
  move.counter<-length(nobs)
  PrintPath(trace,1,Elbow,"E"," ",lambda0,digits,stats)
  k <- 1
  if(plot.it)
    SnapPath(k,x, y, fl, alpha[, k], alpha0[k], lambda[k], Elbow,kernel.function,param.kernel,linear.plot,...)
  while((k < Nmoves)&& (lambda[k]>lambda.min)) {
### Now we implement the updates in Section 4.0
    if(length(Elbow)==0){
### The elbow has become empty; need to resort to an initial condition
      if(sum(y[Left])!=0)stop("Unbalanced data in interior empty elbow situation")
      init<-Balanced.Initialization(K[Left,Left],y[Left],length(Left)/2)
      lambda0 <- init$lambda
      alpha0[k+1] <- init$alpha0
      Elbow <- Left[init$Elbow]
      Left <- setdiff(Left, Elbow)
      lambda[k+1] <- lambda0
      alpha[,k+1]<-alpha[,k]
      Kstar<-UpdateKstar(Kstar,Kscript[Elbow,Elbow],NULL,y[Elbow])
      fl<- (lambda[k]/lambda[k + 1]) * (fl + (alpha0[k+1]-alpha0[k])/lambda[k])
      movefrom<-" ";moveto<-"E";obs<-Elbow
    }
    else{
      bstar<-SolveKstar(Kstar)
      b0 <- bstar[1]
      b <- bstar[-1]
### Now find the first event
### Check for immobile margin
      gl <- K[, Elbow] %*% (y[Elbow] * b) + b0
      dl<-fl-gl
      immobile<-sum(abs(dl))/n < eps
### now check for exits from Elbow
      temp <-  - alpha[Elbow, k] + lambda[k] * b
      lambda.left<-(1+temp)/b
      lambda.left[abs(b)<eps]<- -1 #anything negative
      lambda.right<-temp/b
      lambda.right[abs(b)<eps]<- -1
      lambda01 <- c(lambda.right,lambda.left)
      lambda.exit <- max(lambda01[lambda01 < lambda[k] - eps])
### Check to see if we leave the margin when it is immobile
      if(immobile&(lambda.exit < eps))break
### Now check for entries
      if(!immobile){
      lambdai <- (lambda[k] * (dl))/(y - gl)
      lambdai[abs(y-gl)<eps]<- -Inf
      lambda.entry <- max(lambdai[lambdai < lambda[k] - eps])
    }
    else lambda.entry<- -1 #any negative will do
      lambda.max <- max(lambda.entry, lambda.exit)	
### update lambda, alphas and fit
      lambda[k + 1] <- lambda.max
      alpha[, k + 1] <- alpha[, k]
      alpha[Elbow, k + 1] <- alpha[Elbow, k] - (lambda[k] - 
                                                lambda.max) * b
      alpha0[k + 1] <- alpha0[k] - (lambda[k] - lambda[k + 1]) * b0
      fl <- (lambda[k]/lambda[k + 1]) * (dl) + gl	
### update active sets
      if(lambda.entry > lambda.exit) {
        
###point joins the elbow
        i <- match(lambda.entry, lambdai, 0)[1]
        obs<-i
###assumes for now there is only 1
        moveto<-"E"
        if(match(i, Left, FALSE)) {
          Left <- setdiff(Left, i)
          movefrom<-"L"
        }
        else {
          Right <- setdiff(Right, i)
          movefrom<-"R"
        }
        Kstar<-UpdateKstar(Kstar,Kscript[i, i],drop(Kscript[i, Elbow]),y[i])
        Elbow <- c(Elbow, i)
      }
      else {
###point(s) leaves the elbow; can be more than one
        movefrom<-"E"
        moveto<-NULL
        idrop<-Leaveright<-NULL
        i<-Elbow[abs(lambda.right-lambda.exit)< eps]
        if(length(i)>0){
          Leaveright<-rep(TRUE,length(i))
          idrop<-i
        }
        i<-Elbow[abs(lambda.left-lambda.exit)< eps]
        if(length(i)>0){
          Leaveright<-c(Leaveright,rep(FALSE,length(i)))
          idrop<-c(idrop,i)
        }
        obs<-idrop
        for(j in seq(along=idrop)){
          if(Leaveright[j]) {
            moveto<-c(moveto,"R")
            Right <- c(Right, idrop[j])
          }
          else {
            moveto<-c(moveto,"L")
            Left <- c(Left, idrop[j])
          }
          mi<-match(idrop[j],Elbow)
          Kstar<-DowndateKstar(Kstar,mi)
          Elbow <- Elbow[-mi]
        }
      }
    }
    k <- k + 1
    stats<-StatPath(y,fl,Elbow)
    SumEps[k]<-stats$margin
    Error[k]<-stats$error
    Size.Elbow[k]<-stats$selbow
    nobs<-seq(along=obs)
    Moveto[move.counter+nobs]<-moveto
    Movefrom[move.counter+nobs]<-movefrom
    Step[move.counter+nobs]<-k
    Obs.step[move.counter+nobs]<-obs
    move.counter<-move.counter+length(nobs)
    Elbow.list[[k]]<-Elbow
    PrintPath(trace,k,obs,moveto,movefrom,lambda[k],digits,stats)
    if(plot.it)
    SnapPath(k,x, y,fl, alpha[, k], alpha0[k], lambda[k], Elbow, kernel.function,param.kernel, linear.plot,...)
  }
  obj<-list(alpha=alpha[,seq(k)],alpha0=alpha0[seq(k)],lambda=lambda[seq(k)],alpha00=alpha00,Error=Error[seq(k)],SumEps=SumEps[seq(k)],
            Size.Elbow=Size.Elbow[seq(k)],Elbow=Elbow.list[seq(k)],Moveto=Moveto[seq(move.counter)],Movefrom=Movefrom[seq(move.counter)],Obs.step=Obs.step[seq(move.counter)],Step=Step[seq(move.counter)],kernel=kernel.function,param.kernel=param.kernel,x=x,y=y,call=this.call)
class(obj)<-"svmpath"
  obj
}

Unbalanced.Initialization<-function(K,y,nplus,nminus,eps=1e-8){
pos.init<-function(K,y,nplus,nminus,eps){
  
### Build up initial Left and dominant class indices, 
  Iplus<-seq(y)[y>0]
  Iminus<-seq(y)[y<0]

  Left<-Iminus
  Kscript <- K * outer(y, y)	
  Rmat<-Kscript[Iplus,Iplus]
  c.objective<-Kscript[Iplus,Iminus]%*%rep(1,nminus)
  
  alpha.opt<-OptInit.alpha(nplus,nminus,Rmat,c.objective)$alpha
  alpha<-rep(1,length(y))
  alpha[Iplus]<-alpha.opt
  
  zeros<-alpha.opt<eps
  Right<-if(any(zeros))Iplus[zeros] else NULL
  elbows<-(!zeros)&(alpha.opt<1-eps)
  Elbow<-if(any(elbows))Iplus[elbows] else NULL
  Iplus.Left<-setdiff(Iplus,c(Right,Elbow))
  Left<-c(Left,Iplus.Left)
  
  f<-K%*%(y*alpha)

  fmin<-min(f[Iminus])
    
  if(length(Elbow)>0)
    fmax<-max(f[Elbow])
  else
    {
      fmax<- max(f[Iplus.Left])
      Elbow <-Iplus.Left[abs(f[Iplus.Left]-fmax)<eps]
      Left<-setdiff(Left,Elbow)
    }
  iminus<-Iminus[abs(f[Iminus]-fmin)<eps]
  Elbow<-c(Elbow,iminus)
  Left<-setdiff(Left,iminus)
  lambda<- (fmax-fmin)/2
  beta0<-1-fmax/lambda
  list(Elbow=Elbow,Right=Right,Left=Left,lambda=lambda,alpha=alpha,beta0=beta0)
}

###Find out which class is dominant
more.y <- +1
if(nminus>nplus){
  more.y<- -1
  ### We like to work with + class dominant
  init<-pos.init(K,-y,nminus,nplus,eps)
  ###Fix up some of the entries
  init$beta0<- -init$beta0
}
else init<- pos.init(K,y,nplus,nminus,eps)
alpha0<-init$beta0*init$lambda
 
### package parameters for the left of the start
### alpha0 is linear in lambda, so just solve a system of equations
alpha00<-c(slope = more.y, intercept=init$lambda*(init$beta0-more.y))
c(init,list(alpha0=alpha0,alpha00=alpha00))
}

      
UpdateKstar<-function(Kstar,Kell,Krow,y){
  if(length(y)==1){
    advec<-c(y,Krow)
    Kstar<-rbind(Kstar,advec)
    cbind(Kstar,c(advec,Kell))
  }
  else{
    adrect<-cbind(y,Krow)
    Kstar<-rbind(Kstar,adrect)
    cbind(Kstar,rbind(t(adrect),Kell))
  }
}
  
