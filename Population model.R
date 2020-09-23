###Using values by McCann and Yodzis
#8_18_19 Now I am going to swap m out between reduced in foraging and predation to make the interpretation of the equation easier
#Changelog 8_18_19 in line 59 and 60, m was moved from 1+fmp under predation rate to 1+fmp under foraging rate
#Changelog 8_20_10 I'm breaking the code up into 2 pieces: one for high n and one for low n to make it easier to run. Also, definitely putting the swtiched m into manuscript so need to pretty the figures up. Also changed the vital rates with m
#Notes for 8_20_19 I'm using 0.2 and 2 for m in the final results
##Here pred is more sensitive
setwd('~/Dropbox/Projects/Sensory stress basal/Paper/Final scripts/Swapped m images 8_18_19/Model output for figures for swapped m')

library('deSolve')

low=1
high=2
m.domain= c(0.2,2)#c(0, 5)#c(0.5,5)
n.domain=c(0.002, 25)
# f.domain=c(0.1,10)

final.n=NULL


f.domain=c(1)
for(j in 1:length(f.domain))
{
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
s.domain=seq(from =0, to =0.3, by =0.001) #n=high, pred more sensitive
# # # # # # s.domain=seq(from =0, to =2000, by =1) #n=low, prey more sensitive
for(i in 1: length(s.domain))
{
a=2.17#0.7 #It works!!!!#140#278#5
b=0.1#3.7#7.8#53.69
# f= 0.1#10 #0.1
f=f.domain[j]#0.1, 10
m=m.domain[high]#0.1 #5
n=n.domain[high]#0.002#25
s=s.domain[i]
dh=0.056#0.018 
dp=0.01#0.2#16.77

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

# times = seq(0, 1000000, by =10000)
times = seq(0, 20000, by =1000)


out=lsoda(y=state, times=times, func=Lorenz, parms = parameters, rtol=1e-8, atol=1e-8)
plot(out,main=c(paste(paste('s =',s), paste(','), paste('i =',i), paste(','), paste('f=',f))))

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
}

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




# B.domain=seq(from =1, to =ncol(final), by =en.1)
# plot(final[,B.domain[1]]~s.domain, type='n', ylim=c(min(c(final[,B.domain[1]],final[,B.domain[2]])),max(c(final[,B.domain[1]],final[,B.domain[2]] ))), xlab=("Effective stress"), ylab=('Basal Resource'))
# legend('topright', c( "Low fear", "High fear"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(B.domain))
# {
# lines(final[,B.domain[i]]~s.domain, col=(i+2))	
# }


# H.domain=seq(from =2, to =ncol(final), by =en.1)
# plot(final[,H.domain[2]]~s.domain, type='n',ylim=c(min(c(final[,H.domain[1]],final[,H.domain[2]])),max(c(final[,H.domain[1]],final[,H.domain[2]] ))), xlab=("Effective stress"), ylab=('Prey'))
# legend('topright', c( "Low fear", "High fear"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(H.domain))
# {
# lines(final[,H.domain[i]]~s.domain, col=(i+2))	
# }


# P.domain=seq(from =3, to =ncol(final), by =en.1)
# plot(final[,P.domain[1]]~s.domain, type='n', ylim=c(min(c(final[,P.domain[1]],final[,P.domain[2]])),max(c(final[,P.domain[1]],final[,P.domain[2]] ))), xlab=("Effective stress"), ylab=('Predator'))
# legend('topright', c( "Low fear", "High fear"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(P.domain))
# {
# lines(final[,P.domain[i]]~s.domain, col=(i+2))	
# }

#Vital rates
par(mfrow=c(2,4))
en.1=12
a.percap.std.domain=seq(from =9, to =ncol(final), by =en.1)
plot(final[,a.percap.std.domain[1]]~s.domain,ylim=c(min(final[,a.percap.std.domain[1]]),max(final[,a.percap.std.domain[1]])), type='n', xlab=("Effective stress"), ylab=('Proportion of resource consumed per capita'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(a.percap.std.domain))
{
lines(final[,a.percap.std.domain[i]]~s.domain, col=(i+2))	
}

a.percap.domain=seq(from =4, to =ncol(final), by =en.1)
plot(final[,a.percap.domain[1]]~s.domain, ylim=c(min(final[,a.percap.domain[1]]),max(final[,a.percap.domain[1]])), type='n', xlab=("Effective stress"), ylab=('Number of resource consumed per capita'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(a.percap.domain))
{
lines(final[,a.percap.domain[i]]~s.domain, col=(i+2))	
}

a.percap.B.domain=seq(from =12, to =ncol(final), by =en.1)
plot(final[,a.percap.B.domain[1]]~s.domain, type='n',ylim=c(min(final[,a.percap.B.domain[1]]),max(final[,a.percap.B.domain[1]])),  xlab=("Effective stress"), ylab=('Probability of a resource being consumed')) #This may equal total proportion of resource consumed
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(a.percap.B.domain))
{
lines(final[,a.percap.B.domain[i]]~s.domain, col=(i+2))	
}

a.tot.domain=seq(from =6, to =ncol(final), by =en.1)
plot(final[,a.tot.domain[1]]~s.domain, type='n',ylim=c(min(final[,a.tot.domain[1]]),max(final[,a.tot.domain[1]])), xlab=("Effective stress"), ylab=('Number of resource consumed total'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(a.tot.domain))
{
lines(final[,a.tot.domain[i]]~s.domain, col=(i+2))	
}


b.percap.std.domain=seq(from =10, to =ncol(final), by =en.1)
plot(final[,b.percap.std.domain[1]]~s.domain, type='n',ylim=c(min(final[,b.percap.std.domain[1]]),max(final[,b.percap.std.domain[1]])) ,xlab=("Effective stress"), ylab=('Proportion of prey consumed per capita'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(b.percap.std.domain))
{
lines(final[,b.percap.std.domain[i]]~s.domain, col=(i+2))	
}

b.percap.domain=seq(from =5, to =ncol(final), by =en.1)
plot(final[,b.percap.domain[1]]~s.domain, type='n',ylim=c(min(final[,b.percap.domain[1]]),max(final[,b.percap.domain[1]])), xlab=("Effective stress"), ylab=('Number of prey consumed per capita'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(b.percap.domain))
{
lines(final[,b.percap.domain[i]]~s.domain, col=(i+2))	
}


b.percap.a.domain=seq(from = 11, to=ncol(final), by =en.1)
plot(final[,b.percap.a.domain[1]]~s.domain,ylim=c(min(final[,b.percap.a.domain[1]]),max(final[,b.percap.a.domain[1]])), type='n', xlab=("Effective stress"), ylab=('Probability of a prey being consumed')) #May be the same as total proportion of prey being consumed
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(b.percap.a.domain))
{
lines(final[,b.percap.a.domain[i]]~s.domain, col=(i+2))	
}


b.tot.domain=seq(from =7, to =ncol(final), by =en.1)
plot(final[,b.tot.domain[1]]~s.domain, type='n', xlab=("Effective stress"),ylim=c(min(final[,b.tot.domain[1]]),max(final[,b.tot.domain[1]])), ylab=('Total number of prey consumed'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
for( i in 1:length(b.tot.domain))
{
lines(final[,b.tot.domain[i]]~s.domain, col=(i+2))	
}






# r.b.domain=seq(from =8, to =ncol(final), by =en.1)
# plot(0, type='n', ylim=c(0.0,0.27), xlim=c(0, en), xlab=("Effective stress"), ylab=('Growth rate of B pop'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(r.b.domain))
# {
# lines(final[,r.b.domain[i]]~s.domain, col=(i+2))	
# }

#writing data
write.table(final,'high n high m.csv', sep=',')