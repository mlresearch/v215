##Packages

library(quadprog)

###Functions

##Function to draw Figure 1
ex_I_vec<-function(x,e){
  erg<-c()
  for (i in 1:length(e)) {
    erg<-c(erg, length(x)*mean(x-e[i])^2/(t(x-e[i])%*%(x-e[i])/length(x)))
  }
  return(erg)
}

##Function to find optimal e and objective function 
#in the case of Omega=S => depends on quadprog

Q_n_S <- function(h,e_min,e_max){
  Dmat<-2*solve(var(h)) 
  # times 2 because quadprog uses 1/2*x'Dx-d'x as objective function 
  Amat<-cbind(diag(-1,ncol(h)),diag(1,ncol(h)))
  #Amat'%*%x >= bvec is used by quadprog => -1 is used for upper bounds mean(h)-e_min
  bvec<-c(e_min-apply(h,FUN=mean,MARGIN=2),apply(h,FUN=mean,MARGIN=2)-e_max)
  #transformed feasible region I_ex to mean(h)-I_ex, see Section 2.2.
  return(solve.QP(Dmat,rep(0,ncol(h)),Amat,bvec))
  #d is set to zero as our objective function is of the form 1/2*x'Dx
}

##Function to calculate the objective function in the case Omega=Sigma
#based on the objective function value in the case Omega=S

Q_n_Sigma<-function(x){
  x/((n-1)/n+x/n) #See (proof of) Theorem 4
}


### Drawing Figure 1

set.seed(1783)
y<-rnorm(50,3,1)
curve(ex_I_vec(y,x),xlim=c(mean(y)-6,mean(y)+6),
      ylab="value of objective function",
      xlab="external value")


### Simulation study

##Scenarios: 

#n = 30 , 50
#moments = EY, VARY , EX
#=> 3 single moments , 3 moment duos , 1 moment triple
#=> 7 combinations
#in each scenario the moment(s) could be correctly or incorrectly specified 
#(e_0 in I_ex or not)
#3 versions of the test 
#=> 2x2x7x3 = 84 scenarios

mu<-4 # expected value of X
sigma<-2 #standard deviation of X
beta<-c(16,1) #regression betas
err_sd<-sqrt(60) #standard deviation of regression errors

ey<-c(beta%*%c(1,mu)) #expected value of Y
vary<-beta[2]^2*sigma^2+err_sd^2 #true variance of Y
ex<-mu

nsim=10000 #number of Simulation-loops per scenario

results<-list()

set.seed(12345) #seed to recalculate results in the paper

for(n in c(30,50)){ #sample size loop
  true_m<-c(ey,vary,ex)
  for(j in list(1,2,3,c(1,2),c(2,3),c(1,3),1:3)){ #moment(s) loop
    true_m_used<-true_m[j]
    df=length(true_m_used)
    #correct specified I_ex:
    e_min=0.95*true_m_used
    e_max=1.05*true_m_used
    #incorrectly specified I_ex:
    e_min_mis=1.2*true_m_used
    e_max_mis=1.3*true_m_used
    erg<-c()
    erg_mis<-c()
    for (i in 1:nsim) {#Simulation loop per scenario
      #data generation
      X<-cbind(1,rnorm(n,mu,sigma))
      y<-X%*%beta+rnorm(n,0,err_sd)
      #calculate h^:
      sample_mom<-cbind(y,n/(n-1)*(y-mean(y))^2,X[,2])
      #calculate and collect the two chi^2_min in the correctly specified case
      QS<-Q_n_S(as.matrix(sample_mom[,j]),e_min,e_max)$value*n
      QSig<-Q_n_Sigma(QS)
      erg<-rbind(erg,c(QS,QSig))
      #calculate and collect the two chi^2_min in the incorrectly specified case
      QS<-Q_n_S(as.matrix(sample_mom[,j]),e_min_mis,e_max_mis)$value*n
      QSig<-Q_n_Sigma(QS)
      erg_mis<-rbind(erg_mis,c(QS,QSig))
    }
    results<-c(results,n,list(j),
      list(cbind(correct<-c( 
        #compare chi^2_min to critical values 
        #and calculate type I error rates (correctly specified case)
        sum(erg[,1]>qchisq(0.95,df=df))/nsim, #SH(S)
        sum(erg[,2]>qchisq(0.95,df=df))/nsim, #SH(Sigma)
        sum(erg[,1]*(n-df)/(df*(n-1))>qf(0.95,df1=df,df2=n-df))/nsim #IDC
      ),
        misspeci<-c( 
        #compare chi^2_min to critical values 
        #and calculate power (incorrectly specified case)
        sum(erg_mis[,1]>qchisq(0.95,df=df))/nsim, #SH(S) 
        sum(erg_mis[,2]>qchisq(0.95,df=df))/nsim, #SH(Sigma)
        sum(erg_mis[,1]*(n-df)/(df*(n-1))>qf(0.95,df1=df,df2=n-df))/nsim #IDC
    ))))
  }
}

results 
#prints list of results 

#sorted by n (30 and 50) > moments (1 = EY, 2= Vary , 3= EX) 

#in each matrix, row = test and column = correctly or incorrectly specified
#left column = type I errors (correctly specified case) => compare to alpha=0.05
#right column = power (incorrectly specified case)
#first row = SH(S) ; second row = SH(Sigma) ; third row = IDC
