#####This script solves for a tri-trophic population model under varying sensory environment where both predator and prey have difficulty detecting each other
# setwd('~/Dropbox/Projects/Sensory stress basal/Paper/Scripts for Yodzis and Innes/Images new style/Raw data for f=1')
#Change log 1/19/21: m is swapped to be in the prey foraging rates to be consistent with the manuscript

library('deSolve')

##This model examines a variety of conditions: when either predator or prey are more vulnerable to the sensory stress and when the fear response to prey is either tied or decoupled to its foraging capability
#Low and high denote the two different conditions for both sensory vulnerability and nature of the fear response in prey
low=1
high=2
m.domain=c(0.5,5) # m denotes whether the fear response in prey is tied to prey foraging (low m) or decoupled (high m)
n.domain=c(0.002, 25) # n denotes whether predator (high n) or prey (low n) is more vulnerable to sensory stress

final.n=NULL
f.domain=1 #The default fear intensity is set to 1 for this simulation but in previous runs, f has been allowed to vary

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

#Initial starting conditions for the three populations; B is the basal resource, H is the prey, and P is the predator
B=0.8
H=0.5
P=0.2

#The sensory domain in which the tri-trophic system exists varies depending on the value of n. Therefore, the correct s.domain line must be run depending on whether n has a high or low value, see line 46
s.domain=seq(from =0, to =0.3, by =0.001) #n=high 
# s.domain=seq(from =0, to =2000, by =1) #n=low
for(i in 1: length(s.domain))
{
a=2.17 #Fixed parameter values for attack rates, birth rates
b=0.1
f=f.domain[j]#0.1, 10
m=m.domain[low]#0.5 #5
n=n.domain[high]#0.002#25
s=s.domain[i]
dh=0.056
dp=0.01

parameters = c(a, f, n , m , s, b, dh, dp)
state=c(B=B[i]+0.01,H=H[i]+0.01,P=P[i]+0.01)

#The Lorenz function encapsulates the model equations, which is a modified Lotka-Volterra predatory-prey model. Prey are able to respond to the presence of predators but at a cost of reduced foraging on the basal resource. Sensory stress (s) is introduced into the environment where as s increases, predation rate decreases but concurrently, prey dampen their response to the presence of predators. The overarching question of this model is to determine how do trophic cascades change when sensory stress is present
Lorenz=function(t, state, parameters){
	with(as.list(c(state,parameters)), {
		dB=B*(1-B)-a*H*B/(1+f*m*P/(1+s))
		dH=a*B*H/(1+f*m*P/(1+s))-b*P*H/((1+f*P/(1+s))*(1+n*s))-dh*H
		dP=b*P*H/((1+f*P/(1+s))*(1+n*s))-dp*P
	
		list(c(dB,dH,dP))
	})
}

times = seq(0, 20000, by =100)

#The lsoda function numerically solves for the population abundance after 20000 time steps
out=lsoda(y=state, times=times, func=Lorenz, parms = parameters, rtol=1e-8, atol=1e-8)
plot(out,main=c(paste(paste('s =',s), paste(','), paste('i =',i), paste(','), paste('f=',f)))) #This line of code plots the results of the three populations over time and is a way of visualizing the progress of the function

#Once the population equilibria have been established, the vital rates of the populations can be calculated
out=as.data.frame(out)
B.equil=out$B[length(out$B)] #The equilibrium basal resource population
H.equil=out$H[length(out$H)]#The equilibrium prey population
P.equil=out$P[length(out$P)]#The equilibrium predator population
r.b=B.equil*(1-B.equil) #The growth rate of the basal resource population
a.percap=a*B.equil/(1+f*P.equil/(1+s)) #The per capita foraging rate of prey
b.percap=b*H.equil/((1+f*m*P.equil/(1+s))*(1+n*s)) #The per capita predation rate of predators
a.tot=a*B.equil*H.equil/(1+f*P.equil/(1+s)) #The total foraging rate of prey population i.e. the amount of the basal resource the prey population is consuming
b.tot=b*P.equil*H.equil/((1+f*m*P.equil/(1+s))*(1+n*s)) #The total predation rate of predator population
a.percap.std=a/(1+f*P.equil/(1+s)) #The per capita foraging rate of prey but standardized by the basal resource population (The proportion of the basal resource population each prey is consuming)
b.percap.std=b/((1+f*m*P.equil/(1+s))*(1+n*s)) #The per capita predation rate of predators but standardized by the prey population
a.percap.B=a*H.equil/(1+f*P.equil/(1+s)) #The proportion of the basal resource that is being consumed by the prey population
b.percap.a=b*P.equil/((1+f*m*P.equil/(1+s))*(1+n*s)) #The proportion of the prey that is being consumed by the predator population

B=c(B,B.equil)
H=c(H, H.equil)
P=c(P, round(P.equil, 50))

#This creates a data frame of the equilibrium population and vital rates levels as the for loop loops through increasing sensory stress levels
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

####################Preliminary plots of results
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


### Plotting of vital rates
# par(mfrow=c(2,4))
# en.1=12
# a.percap.std.domain=seq(from =9, to =ncol(final), by =en.1)
# plot(final[,a.percap.std.domain[1]]~s.domain,ylim=c(min(final[,a.percap.std.domain[2]]),max(final[,a.percap.std.domain[1]])), type='n', xlab=("Effective stress"), ylab=('Proportion of resource consumed per capita'))
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(a.percap.std.domain))
# {
# lines(final[,a.percap.std.domain[i]]~s.domain, col=(i+2))	
# }

# a.percap.domain=seq(from =4, to =ncol(final), by =en.1)
# plot(final[,a.percap.domain[1]]~s.domain, ylim=c(min(final[,a.percap.domain[2]]),max(final[,a.percap.domain[1]])), type='n', xlab=("Effective stress"), ylab=('Number of resource consumed per capita'))
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(a.percap.domain))
# {
# lines(final[,a.percap.domain[i]]~s.domain, col=(i+2))	
# }

# a.percap.B.domain=seq(from =12, to =ncol(final), by =en.1)
# plot(final[,a.percap.B.domain[1]]~s.domain, type='n',ylim=c(min(final[,a.percap.B.domain[1]]),max(final[,a.percap.B.domain[1]])),  xlab=("Effective stress"), ylab=('Probability of a resource being consumed')) #This may equal total proportion of resource consumed
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(a.percap.B.domain))
# {
# lines(final[,a.percap.B.domain[i]]~s.domain, col=(i+2))	
# }

# a.tot.domain=seq(from =6, to =ncol(final), by =en.1)
# plot(final[,a.tot.domain[2]]~s.domain, type='n',ylim=c(min(final[,a.tot.domain[2]]),max(final[,a.tot.domain[1]])), xlab=("Effective stress"), ylab=('Number of resource consumed total'))
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(a.tot.domain))
# {
# lines(final[,a.tot.domain[i]]~s.domain, col=(i+2))	
# }


# b.percap.std.domain=seq(from =10, to =ncol(final), by =en.1)
# plot(final[,b.percap.std.domain[2]]~s.domain, type='n',ylim=c(min(final[,b.percap.std.domain[2]]),max(final[,b.percap.std.domain[1]])) ,xlab=("Effective stress"), ylab=('Proportion of prey consumed per capita'))
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(b.percap.std.domain))
# {
# lines(final[,b.percap.std.domain[i]]~s.domain, col=(i+2))	
# }

# b.percap.domain=seq(from =5, to =ncol(final), by =en.1)
# plot(final[,b.percap.domain[2]]~s.domain, type='n',ylim=c(min(final[,b.percap.domain[2]]),max(final[,b.percap.domain[2]])), xlab=("Effective stress"), ylab=('Number of prey consumed per capita'))
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(b.percap.domain))
# {
# lines(final[,b.percap.domain[i]]~s.domain, col=(i+2))	
# }


# b.percap.a.domain=seq(from = 11, to=ncol(final), by =en.1)
# plot(final[,b.percap.a.domain[2]]~s.domain,ylim=c(min(final[,b.percap.a.domain[2]]),max(final[,b.percap.a.domain[1]])), type='n', xlab=("Effective stress"), ylab=('Probability of a prey being consumed')) #May be the same as total proportion of prey being consumed
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(b.percap.a.domain))
# {
# lines(final[,b.percap.a.domain[i]]~s.domain, col=(i+2))	
# }


# b.tot.domain=seq(from =7, to =ncol(final), by =en.1)
# plot(final[,b.tot.domain[2]]~s.domain, type='n', xlab=("Effective stress"),ylim=c(min(final[,b.tot.domain[2]]),max(final[,b.tot.domain[1]])), ylab=('Total number of prey consumed'))
# # legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(b.tot.domain))
# {
# lines(final[,b.tot.domain[i]]~s.domain, col=(i+2))	
# }






# r.b.domain=seq(from =8, to =ncol(final), by =en.1)
# plot(0, type='n', ylim=c(0.0,0.27), xlim=c(0, en), xlab=("Effective stress"), ylab=('Growth rate of B pop'))
# legend('topright', c( "Predator does better", "Prey does better"), lty=c(1,1), col=c(3,4))
# for( i in 1:length(r.b.domain))
# {
# lines(final[,r.b.domain[i]]~s.domain, col=(i+2))	
# }

#writing data
# write.table(final,'low n high m.csv', sep=',')