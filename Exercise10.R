rm(list=ls()) #remove global environment
setwd("/Users/Ashley/Documents/Biocomputing_2018/Biocomp-Fall2018-181109-Exercise10")

data=read.csv("data.txt", header = TRUE)
data

# plot our observations
library(ggplot2) 
ggplot(data,aes(x=x,y=y))+geom_point()+theme_classic()

### Custom likelihood functions
nllike<-function(p,x,y){
  B0=p[1] 
  B1=p[2] 
  sigma=exp(p[3])
  expected=B0+B1*x
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll) 
}

nllike2<-function(p,x,y){
  B0=p[1] 
  B1=p[2] 
  B2=p[3]
  sigma=exp(p[4])
  expected=B0+B1*x +B2*x^2
  nll2=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll2) 
}

### estimate parameters by minimizing the negative log likelihood 
initialGuess=c(1,1,1,1) 
fit = optim(par=initialGuess,fn=nllike,x=data$x,y=data$y)
fit2 = optim(par=initialGuess,fn=nllike2,x=data$x,y=data$y)
# fit is a variable that contains a list describing the result of minimization
# $par is a vector of most likely parameter values in the order specified by
#   your custom likelihood function
# $value is the minimum negative log likelihood found by optim
# $counts gives us the number of parameter combinations tried by optim
# $convergence is a code with information about the optim run; 0 means success
print(fit)
print(fit2)









#2
tumorSim<-function(t,y,p){
  N=y[1]
  T=y[2]
  
  RN=p[1]
  a11=p[2]
  RT=p[3]
  a22=p[4]
  a21=p[5]
  a12=p[6]
  
  dNdt=RN*(1-N*a11 - T*a12) *N
  dTdt=RT*(1-T*a22 - N*a21) *T
  
  return(list(c(dNdt,dTdt)))
}

# case 1
times=1:100
y0=c(0.2,0.8)
params2=c(0.5,.015,0.5,0.015,.01,0.01)
sim2=ode(y=y0,times=times,func=tumorSim,parms=params2)
out2=data.frame(time=sim2[,1],normal=sim2[,2],tumor=sim2[,3])
ggplot(out2,aes(x=time,y=normal))+geom_line()+geom_line(data=out2,mapping=aes(x=time,y=tumor),col='red')+theme_classic()

# case 2 (a21>a22)
times=1:100
y0=c(0.2,0.8)
params3=c(0.5,.015,0.5,0.015,.02,0.01)
sim3=ode(y=y0,times=times,func=tumorSim,parms=params3)
out3=data.frame(time=sim3[,1],normal=sim3[,2],tumor=sim3[,3])
ggplot(out3,aes(x=time,y=normal))+geom_line()+geom_line(data=out3,mapping=aes(x=time,y=tumor),col='red')+theme_classic()

# case 3 (a12>a11)
times=1:100
y0=c(0.2,0.8)
params4=c(0.5,.015,0.5,0.015,.01,0.02)
sim4=ode(y=y0,times=times,func=tumorSim,parms=params4)
out4=data.frame(time=sim4[,1],normal=sim4[,2],tumor=sim4[,3])
ggplot(out4,aes(x=time,y=normal))+geom_line()+geom_line(data=out4,mapping=aes(x=time,y=tumor),col='red')+theme_classic()
