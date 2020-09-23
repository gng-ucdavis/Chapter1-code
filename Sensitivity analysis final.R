####For the sensitivity analysis, I will rerun the analysis for all 4 scenarios but swap out the other parameters each time by halving or doubling the value individually. But no interaction of value so here's the total number of analyses
# a* high and low
# b* high and low
# dh* high and low
# dp* high and low
# f* high and low 
##### so five different parameters, each with high and low so 10 different parameter combinations over 4 scearios = 40 different graphs

###adding beep noises
# install.packages('beepr') 
library(beepr) 
setwd('~/Dropbox/Projects/Sensory stress basal/Paper/Final scripts/Swapped m images 8_18_19/Sensitivity analysis')

library('deSolve')

low=1
high=2
m.domain= c(0.2, 2)#c(0.2,2)#c(0.1, 5)#c(0.5,5)
n.domain=c(0.002, 25)
####Pulling out the parameters I am tweaking to make it consistent from run to run
##The values below are the original ones used for the model and annotated out to the right are the different values used for this sensitivity analysis
a=2.17			#2.17/2	#2*2.17 #		#2.17 #done low a, low m; low a, high m; high a, low m; high a, high m
b=0.1			#0.1/2	#0.1*2			#0.1 #done low b, low m; low b, high m; high b, low m; high b, high m
f=1 				#2#1/2					#1 #done low f, low m; low f, high m; high f, low m; high f, high m
dh=0.056			#0.056/2 0.056*2		#0.056 #done low dh, low m; low dh, high m; high dh, low m; high dh, high m
dp=0.01			#0.01/2	 0.01*2			#0.01	#done low dp, low m; low dp, high m; high dp low m; high dp, high m				
m=m.domain[high]							#Still need to change each round
n.time=1000								#Used to speed up models or slow down if models are not running

name.low=('low dp low n high m.csv') ##name low ##swap out the name as appropriate depending on which parameter is being tested
name.high=('low dp high n high m.csv') ##name high


final.n=NULL

B.final=NULL
H.final=NULL
P.final=NULL

a.percap.final=NULL
b.percap.final=NULL
a.tot.final=NULL
b.tot.final=NULL
r.b.final=NULL
a.percap.std.final=NULL
b.percap.std.final=NULL
b.percap.a.final=NULL
a.percap.B.final=NULL

B=0.8
H=0.5
P=0.2
# # # # # # s.domain=seq(from =0, to =0.45, by =0.001) #n=high, pred more sensitive, right this is why you wanted to split out s and n
s.domain=seq(from =0, to =6000, by =2) #n=low, prey more sensitive #changing s to be more coarse just to speed up time #original length is to 2000 but used 6000 for some long runs
for(i in 1: length(s.domain))
{
n=n.domain[low]#0.002#25
s=s.domain[i]


parameters = c(a, f, n , m , s, b, dh, dp)
state=c(B=B[i]+0.01,H=H[i]+0.01,P=P[i]+0.01)
# state=c(B=0.8,H=0.5,P=0.2)

Lorenz=function(t, state, parameters){
	with(as.list(c(state,parameters)), {
		dB=B*(1-B)-a*H*B/(1+f*m*P/(1+s))
		dH=a*B*H/(1+f*m*P/(1+s))-b*P*H/((1+f*P/(1+s))*(1+n*s))-dh*H
		dP=b*P*H/((1+f*P/(1+s))*(1+n*s))-dp*P
	
		list(c(dB,dH,dP))
	})
}


times = seq(0, 20000, by =n.time) #dropping the by from 1000 to 1 can make the model run in certain cases


out=lsoda(y=state, times=times, func=Lorenz, parms = parameters, rtol=1e-8, atol=1e-8)
plot(out,main=c(paste(paste('s =',s), paste(','), paste('i =',i), paste(','), paste('f=',f), paste(','), paste('m=',m))))

out=as.data.frame(out)
B.equil=out$B[length(out$B)]
H.equil=out$H[length(out$H)]
P.equil=out$P[length(out$P)]
r.b=B.equil*(1-B.equil)
a.percap=a*B.equil/(1+f*m*P.equil/(1+s))
b.percap=b*H.equil/((1+f*P.equil/(1+s))*(1+n*s))
a.tot=a*B.equil*H.equil/(1+f*m*P.equil/(1+s))
b.tot=b*P.equil*H.equil/((1+f*P.equil/(1+s))*(1+n*s))
a.percap.std=a/(1+f*m*P.equil/(1+s))
b.percap.std=b/((1+f*P.equil/(1+s))*(1+n*s))
b.percap.a=b*P.equil/((1+f*P.equil/(1+s))*(1+n*s))
a.percap.B=a*H.equil/(1+f*m*P.equil/(1+s))

B=c(B,B.equil)
H=c(H, H.equil)
P=c(P, round(P.equil, 50))

B.final=c(B.final, B.equil)
H.final=c(H.final, H.equil)
P.final=c(P.final, P.equil)
r.b.final=c(r.b.final,r.b)
a.percap.final=c(a.percap.final,a.percap)
b.percap.final=c(b.percap.final, b.percap)
a.tot.final=c(a.tot.final, a.tot)
b.tot.final=c(b.tot.final, b.tot)
a.percap.std.final=c(a.percap.std.final,a.percap.std)
b.percap.std.final=c(b.percap.std.final,b.percap.std)
b.percap.a.final=c(b.percap.a.final, b.percap.a)
a.percap.B.final=c(a.percap.B.final, a.percap.B)
}

final.s=(cbind(B.final, H.final, P.final, a.percap.final, b.percap.final, a.tot.final, b.tot.final, r.b.final, a.percap.std.final, b.percap.std.final,b.percap.a.final,a.percap.B.final))
final.n=cbind(final.n,final.s)


final=NULL
final=as.data.frame(final.n)

####################Plotting results
s.domain=seq(from =0, to =6000, by =2) 
par(mfrow=c(2,2))
en.1=12

B.domain=seq(from =1, to =ncol(final), by =en.1)
plot(final[,B.domain[1]]~s.domain, type='n', xlab=("Effective stress"), ylab=('Basal Resource'))
for( i in 1:length(B.domain))
{
lines(final[,B.domain[i]]~s.domain, col=(3))	
}


H.domain=seq(from =2, to =ncol(final), by =en.1)
plot(final[,H.domain[1]]~s.domain, type='n', xlab=("Effective stress"), ylab=('Prey'))
for( i in 1:length(H.domain))
{
lines(final[,H.domain[i]]~s.domain, col=(4))	
}


P.domain=seq(from =3, to =ncol(final), by =en.1)
plot(final[,P.domain[1]]~s.domain, type='n', xlab=("Effective stress"), ylab=('Predator'))
for( i in 1:length(P.domain))
{
lines(final[,P.domain[i]]~s.domain, col=(2))	
}

plot(final[,B.domain]~s.domain,type='n', ylim=c(0,1), xlab=("Sensory stress"), ylab=('Standardized population abundance'))
for( i in 1:length(B.domain))
{
	Br.new=final[,B.domain[i]]/max(final[,B.domain[i]])
lines(Br.new~s.domain, col=(3), lwd=4)	
}
for( i in 1:length(H.domain))
{
	H.new=final[,H.domain[i]]/max(final[,H.domain[i]])
lines(H.new~s.domain, col=(4), lwd=4)	
}
for( i in 1:length(P.domain))
{
	P.new=final[,P.domain[i]]/max(final[,P.domain[i]])
lines(P.new~s.domain, col=(2), lwd=4)	
}

write.table(final, name.low , sep=',')
final1=final #just a placeholder for plotting 
# final=final1

#############High n
low=1
high=2
m.domain= c(0.2, 2)#c(0.2,2)#c(0.1, 5)#c(0.5,5)
n.domain=c(0.002, 25)
# f.domain=c(0.1,10)

final.n=NULL


B.final=NULL
H.final=NULL
P.final=NULL

a.percap.final=NULL
b.percap.final=NULL
a.tot.final=NULL
b.tot.final=NULL
r.b.final=NULL
a.percap.std.final=NULL
b.percap.std.final=NULL
b.percap.a.final=NULL
a.percap.B.final=NULL

B=0.8
H=0.5
P=0.2
s.domain=seq(from =0, to =0.7, by =0.001) #n=high, pred more sensitive, right this is why you wanted to split out s and n #original end is 0.45 but used 0.7 for some long runs
for(i in 1: length(s.domain))
{
n=n.domain[high]#0.002#25
s=s.domain[i]


parameters = c(a, f, n , m , s, b, dh, dp)
state=c(B=B[i]+0.01,H=H[i]+0.01,P=P[i]+0.01)
# state=c(B=0.8,H=0.5,P=0.2)

Lorenz=function(t, state, parameters){
	with(as.list(c(state,parameters)), {
		dB=B*(1-B)-a*H*B/(1+f*m*P/(1+s))
		dH=a*B*H/(1+f*m*P/(1+s))-b*P*H/((1+f*P/(1+s))*(1+n*s))-dh*H
		dP=b*P*H/((1+f*P/(1+s))*(1+n*s))-dp*P
	
		list(c(dB,dH,dP))
	})
}


times = seq(0, 20000, by =n.time)


out=lsoda(y=state, times=times, func=Lorenz, parms = parameters, rtol=1e-8, atol=1e-8)
plot(out,main=c(paste(paste('s =',s), paste(','), paste('i =',i), paste(','), paste('f=',f), paste(','), paste('m=',m))))

out=as.data.frame(out)
B.equil=out$B[length(out$B)]
H.equil=out$H[length(out$H)]
P.equil=out$P[length(out$P)]
r.b=B.equil*(1-B.equil)
a.percap=a*B.equil/(1+f*m*P.equil/(1+s))
b.percap=b*H.equil/((1+f*P.equil/(1+s))*(1+n*s))
a.tot=a*B.equil*H.equil/(1+f*m*P.equil/(1+s))
b.tot=b*P.equil*H.equil/((1+f*P.equil/(1+s))*(1+n*s))
a.percap.std=a/(1+f*m*P.equil/(1+s))
b.percap.std=b/((1+f*P.equil/(1+s))*(1+n*s))
b.percap.a=b*P.equil/((1+f*P.equil/(1+s))*(1+n*s))
a.percap.B=a*H.equil/(1+f*m*P.equil/(1+s))

B=c(B,B.equil)
H=c(H, H.equil)
P=c(P, round(P.equil, 50))

B.final=c(B.final, B.equil)
H.final=c(H.final, H.equil)
P.final=c(P.final, P.equil)
r.b.final=c(r.b.final,r.b)
a.percap.final=c(a.percap.final,a.percap)
b.percap.final=c(b.percap.final, b.percap)
a.tot.final=c(a.tot.final, a.tot)
b.tot.final=c(b.tot.final, b.tot)
a.percap.std.final=c(a.percap.std.final,a.percap.std)
b.percap.std.final=c(b.percap.std.final,b.percap.std)
b.percap.a.final=c(b.percap.a.final, b.percap.a)
a.percap.B.final=c(a.percap.B.final, a.percap.B)
}

final.s=(cbind(B.final, H.final, P.final, a.percap.final, b.percap.final, a.tot.final, b.tot.final, r.b.final, a.percap.std.final, b.percap.std.final,b.percap.a.final,a.percap.B.final))
final.n=cbind(final.n,final.s)


final=NULL
final=as.data.frame(final.n)


####################Plotting results
par(mfrow=c(2,2))
en.1=12

B.domain=seq(from =1, to =ncol(final), by =en.1)
plot(final[,B.domain[1]]~s.domain, type='n', xlab=("Effective stress"), ylab=('Basal Resource'))
for( i in 1:length(B.domain))
{
lines(final[,B.domain[i]]~s.domain, col=(3))	
}


H.domain=seq(from =2, to =ncol(final), by =en.1)
plot(final[,H.domain[1]]~s.domain, type='n', xlab=("Effective stress"), ylab=('Prey'))
for( i in 1:length(H.domain))
{
lines(final[,H.domain[i]]~s.domain, col=(4))	
}


P.domain=seq(from =3, to =ncol(final), by =en.1)
plot(final[,P.domain[1]]~s.domain, type='n', xlab=("Effective stress"), ylab=('Predator'))
for( i in 1:length(P.domain))
{
lines(final[,P.domain[i]]~s.domain, col=(2))	
}

plot(final[,B.domain]~s.domain,type='n', ylim=c(0,1), xlab=("Sensory stress"), ylab=('Standardized population abundance'))
for( i in 1:length(B.domain))
{
	Br.new=final[,B.domain[i]]/max(final[,B.domain[i]])
lines(Br.new~s.domain, col=(3), lwd=4)	
}
for( i in 1:length(H.domain))
{
	H.new=final[,H.domain[i]]/max(final[,H.domain[i]])
lines(H.new~s.domain, col=(4), lwd=4)	
}
for( i in 1:length(P.domain))
{
	P.new=final[,P.domain[i]]/max(final[,P.domain[i]])
lines(P.new~s.domain, col=(2), lwd=4)	
}

write.table(final, name.high ,sep=',')
beep('facebook')