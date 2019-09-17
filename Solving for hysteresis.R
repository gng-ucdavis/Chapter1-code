library('rootSolve')

##Annotating the crap out of this.
#This here is to set your parameters for m,n, and f
low=1
high=2
magnitude=low
m.domain=c(0.5,5)
n.domain=c(0.002, 25)
f.domain=c(0.1,10)

#Pick the right s.domain that corresponds to n level
par(mfrow=c(2,2))
# s.domain=seq(from =0, to =0.3, by =0.001) #high
s.domain=seq(from =0, to =2000, by =1) #low
H.star=matrix(nrow=length(s.domain), ncol=4) #creating a blank matrix to fill in root solutions (assuming that you don't get more than 4 equil values)

####Actual parameter values
for(j in 1: length(s.domain))  #we are going to run this over the range of ss values (no surprises here)
{
#Parameters
a=2.17
b=0.1
n=n.domain[low]#0.002#25 #we've done low, high and high high, low low, [High,low has true alternative equilibrium values! Those turn out to be unstable but the P=0 is stable. Time to try high fear]
m=m.domain[low]#0.5 #5
f=f.domain[high] ##Writing down the results of high fear here. High low also has true alternative equil values and P=0 is stable. # High high has no ass# low high has none either, nor low low. Completed stability analysis
s=s.domain[j]
dh=0.056
dp=0.01

PforH.equil=function(P) #A function to find at what P level will the population rate of H be 0. The math for this equation is in iPad. Solve everything in terms of P.
{
	H.equil=a*((-a*((dp*(n*s + 1) * (f*m*P + s + 1))/(b*(s+1)))*(s + 1) + f*P + s + 1)/(f*P+s+1))*((dp*(n*s + 1) * (f*m*P + s + 1))/(b* (s + 1)))/(1+(f*P)/(1+s))-b*((dp*(n*s + 1) * (f*m*P + s + 1))/(b* (s + 1)))*P/((1+(f*m*P)/(1+s))*(1+n*s))-dh*((dp*(n*s + 1) * (f*m*P + s + 1))/(b* (s + 1)))
	return(H.equil)
}
P.test=seq(-200,200, by =10) ##This parameter here is just for plotting purposes
H.equil=PforH.equil(P.test) #We take P.test and run it through the function to see what kind of dh/dt values we get. Obvi, want the ones that goes through 0
plot(H.equil~P.test, main=c(paste(paste('s=',s), paste(','), paste('i =',i), paste(','), paste('f=',f)))) #Plotting the analysis
sol=uniroot.all(PforH.equil, c(-20,30)) #This is the important function. Uniroot.all looks for multiple roots over the range of -20 to 30. This range is the one you will have to tweak if you aren't getting solutions. Numeric(0) means you get nothing; not 0 is a solution
H.star[j,1:length(sol)]=sol #Filling in the matrix with the solutions
}

H.star=cbind(s.domain,H.star) #Adding a row for sensory stress
H.starr=as.data.frame(H.star) #Making the matrix a dataframe but we still have the issues of having negative solutions
head(H.starr) 
#####You're probably safe running the script up to this point

# H.starr.temp=H.starr[1,][H.starr[1,]>=0] #Alright, this works but I need it to work for the entire matrix

#This part of the script here is to remove negative solutions
H.starr.final=matrix(nrow=length(s.domain), ncol=5) ##Create new matrix of blanks
for(i in 1: nrow(H.starr))
{
	H.starr.temp=H.starr[i,][H.starr[i,]>=0] ###This is a cool bit of code that filters out all negative values from H.starr for a particular row. The syntas is weird af
	H.starr.final[i,1:length(H.starr.temp)]=H.starr.temp #Take the row of non negative values and add it to H.starr.final
}
H.starr.final=as.data.frame(H.starr.final)
head(H.starr.final)
H.starr.final=H.starr.final[complete.cases(H.starr.final$V2),] #What we have here is the P value needed for H to be in equilibrium at a given S. So what's H? when it is in equil? Plug this into the P equation #Complete.cases is a function that removes NA that is different than is.na (better?)

####Run the following script only if you get mutliple equil values #Look up the script "checking unstable equilibrium" on what to do. 
#This is to remove other equilibrium points if you happen to get more than one set
# test=c(H.starr.final$V2[1:181], H.starr.final$V3[182:225]) #Obvi, the numbers of H.starr.final will change
# H.starr.final$V2=test


##Solving for prey quilibrium value using predator equation
# 0=b*H*P/((1+(f*m*P)/(1+s))*(1+n*s))-d*P
##Now that we have the P values at equilibiru, let's solve for H

#Create blank matrix
P.star=matrix(nrow=length(s.domain), ncol=4)

par(mfrow=c(2,2))
for(j in 1: nrow(H.starr.final))
{
	s=H.starr.final$V1[j] #Call up relevant sensory stress values in H.starr.finals
	# s=H.starr.final$s.domain[j]
	P=H.starr.final$V2[j] #Call up relevant P.equil values
HforP.equil=function(H) #Here's the function to solve for H.equil when P population isn't changing
{
	# P.equil=b*H*P/((1+(f*m*P)/(1+s))*(1+n*s))-dp*P
	P.equil=b*H/((1+(f*m*P)/(1+s))*(1+n*s))-dp
	return(P.equil)
}
H.test=seq(-200,200, by =1) #For plotting purposes
P.equil=HforP.equil(H.test)
plot(P.equil~H.test, main=c(paste(paste('s=',s), paste(','), paste('i =',i), paste(','), paste('f=',f))))
sol=uniroot.all(HforP.equil, c(-20,200)) #Again, you might need to change the range
P.star[j,1:length(sol)]=sol #Filling out the matrix
}

P.star=as.data.frame(P.star)
P.star.final=P.star[complete.cases(P.star$V1),]

###Now that we have 2 data frames with solutions, time to combine them into final data set
final=NULL
final$s=H.starr.final$V1 #for sensory stress
final$P=H.starr.final$V2 #For p.equil
final$H=P.star.final$V1 #For h.equil
final=as.data.frame(final)


###So now that I have P and H that causes P and H to be in equil, I need to know what the B is.

B.star=matrix(nrow=length(s.domain), ncol=4)

par(mfrow=c(2,2))
for(j in 1: nrow(final))
{
	s=final$s[j]
	P=final$P[j]
	H=final$H[j]
BforB.equil=function(B)
{
	B.equil=B*(1-B)-a*B*H/(1+(f*P)/(1+s))	 #Function to solve B given H and P. Note that with the way this is set up, B = 0 is always a solution and should be ignored.
	return(B.equil)
}
B.test=seq(-200,200, by =0.1)
B.equil=BforB.equil(B.test)
plot(B.equil~B.test, main=c(paste(paste('s=',s), paste(','), paste('i =',i), paste(','), paste('f=',f))), ylim=c(-0.5,0.5), xlim =c(0,1))
sol=uniroot.all(BforB.equil, c(0,1)) #Don't think you have to tweak this since B is constrained by 0 and 1
B.star[j,1:length(sol)]=sol
}
B.star=as.data.frame(B.star)
B.star.final=B.star[complete.cases(B.star$V2),] #Remove NAs
final$B=B.star.final$V2 #Add to final data frame

##I think we're done here, but let's check how well this works
B.test=NULL
H.test=NULL
P.test=NULL
for(k in 1:nrow(final))
{
s=final[k,1]
P=final[k,2]
H=final[k,3]
B=final[k,4] #Getting the respective B, H, and B values to test if they come out to 0 when plugged in

B.test[k]=B*(1-B)-a*B*H/(1+(f*P)/(1+s)) #db
H.test[k]=a*(B)*(H)/(1+(f*P)/(1+s))-b*(H)*P/((1+(f*m*P)/(1+s))*(1+n*s))-dh*(H) #dh
P.test[k]=b*H*P/((1+(f*m*P)/(1+s))*(1+n*s))-dp*P #dp
}
round(B.test,4)
round(H.test,4) ##Round to get pretty 0 values
round(P.test,4)

##Just check stability of 'correct' equilibrium value
##We've basically runge kutta our equations and now we check for their stability although they should all be stable.
e.values=NULL
for(i in 1: nrow(final))
{
	s=final[i,1]
	B=final[i,4] #Take your respective s,B,H, and P equil values... (I named the columns in final but used numbered columns here for notations. Inconsistent on my part)
	H=final[i,3]
	P=final[i,2]
	
# aa=a*H*(1+s)/(1+s+f*P)-2*B+1 #Bad derivation
aa=-(a*H*(s + 1))/(f*P+s+1)-2*B+1

# bb=a*B*(1+s)/(1+s+f*P) #Bad
 bb=-(a*B*(s+1))/(f*P+s+1)

# cc=a*f*B*H*(1+s)/(1+s+f*P)^2
 cc=(a*B*f*H*(s+1))/(f*P+s+1)^2

# dd=a*H*(1+s)/(1+s+f*P)
dd=(a*H*(s+1))/(f*P+s+1)

# ee=a*B/(1+f*P/(1+s))-b*P/((1+f*m*P/(1+s))*(1+n*s))-dh
ee=(a*B)/((f*P)/(s+1)+1)-(b*P)/((n*s+1)*((f*m*P)/(s+1)+1))-dh

# ff=-a*f*B*H/((1+s)*(1+f*P/(1+s))^2)-b*H/((1+n*s)*(1+(f*m*P)/(1+s)))+b*f*m*H*P/((1+s)*(1+n*s)*(1+(f*m*P)/(1+s))^2)
ff=-(a*B*f*H)/((s+1)*((f*P)/(s+1)+1)^2)-(b*H)/((n*s+1)*((f*m*P)/(s+1)+1))+(b*f*H*m*P)/((s+1)*(n*s+1)*((f*m*P)/(s+1)+1)^2)

gg=0

# hh=b*P*(1+s)/((1+s+f*m*P)*(1+n*s))
hh=(b*P*(s+1))/((n*s+1)*(f*m*P+s+1))

# ii=b*H/((1+f*m*P/(1+s))*(1+n*s))-b*f*m*H*P/((1+s)*(1+n*s)*(1+f*m*P/(1+s))^2)-dp
ii=(b*H)/((n*s+1)*((f*m*P)/(s+1)+1))-(b*f*H*m*P)/((s+1)*(n*s+1)*((f*m*P)/(s+1)+1)^2)-dp

jmat=(matrix(c(aa,bb,cc,dd,ee,ff,gg,hh,ii), ncol=3, byrow=T))

e.values=rbind(e.values, eigen(jmat)$values)  ##whole bunch of derivations and eigen values up top
}
e.values.fin=cbind(final$s,e.values)
e.values.fin=as.data.frame(e.values.fin)
e.values.fin #all values should have negative real parts here

####Part 2
#####I want to test for when P is 0
#The simple solution is to set P to 0, find the equil values and construct a matrix of repeating values. However, i want to make sure that this works so I am setting P to 0 and then actively solve for it over ss range values (just a way to double check my coding. )

B.equil.P0=matrix(nrow=length(s.domain), ncol=4) #Empty matrix

alt.equil=NULL #Empty matrix?

#Parameters (should already be done up top) since we are looking at the same scenarios
# a=2.17
# b=0.1
# n=n.domain[low]#0.002#25
# m=m.domain[high]#0.5 #5
# f=f.domain[low]
# s=s.domain[j]
# dh=0.056
# dp=0.01
P=0 #The one addition of a new constant


for(j in 1: nrow(H.starr.final))
{
	s=final$s[j]
	# s=H.starr.final$s.domain[j]
HforP.equil=function(B)
{
	# P.equil=b*H*P/((1+(f*m*P)/(1+s))*(1+n*s))-dp*P
	H.equil=a*B/(1+f*P/(1+s))-b*P/((1+f*m*P/(1+s))*(1+n*s))-dh #solve for B when dh/dt is 0 and when P is set to 0
	return(H.equil)
}
B.test=seq(-200,200, by =1)
H.equil=HforP.equil(B.test)
plot(H.equil~B.test, main=c(paste(paste('s=',s), paste(','), paste('i =',i), paste(','), paste('f=',f))))
sol=uniroot.all(HforP.equil, c(-20,200)) ##I guess I should set this range to be 0 to 1
B.equil.P0[j,1:length(sol)]=sol
}

B.equil.P0=as.data.frame(B.equil.P0)
alt.equil=cbind(final$s, B.equil.P0[complete.cases(B.equil.P0$V1),1],rep(0, length(H.starr.final$V2))) ##Creates a data frame of ss, B.equil, and column of 0s to represent P
# alt.equil=final$s
# alt.equil$B=B.equil.P0[complete.cases(B.equil.P0$V1),1]
# alt.equil$P=rep(0, length(H.starr.final$V2))

alt.equil=as.data.frame(alt.equil)
colnames(alt.equil)=c('s','B','P') ##This is a cool function to name columns on dataframe.
head(alt.equil)

##Now we solve for prey values. Again, this isn't hard nor do we need to do this (It's also the same across all scenarios)
H.equil.B=matrix(nrow=length(s.domain), ncol=4)
##Found basal resource equilibrium. Now what is H?
for(j in 1: nrow(H.starr.final))
{
	s=final$s[j]
	B=alt.equil$B[j]
HforB.equil=function(H)
{
	# P.equil=b*H*P/((1+(f*m*P)/(1+s))*(1+n*s))-dp*P
	B.equil=(1-B)-a*H/(1+f*P/(1+s))
	return(B.equil)
}
H.test=seq(-200,200, by =1)
B.equil=HforB.equil(H.test)
plot(B.equil~H.test, main=c(paste(paste('s=',s), paste(','), paste('i =',i), paste(','), paste('f=',f))))
sol=uniroot.all(HforB.equil, c(-20,2))
H.equil.B[j,1:length(sol)]=sol
}
H.equil.B

H.equil.B=as.data.frame(H.equil.B)
alt.equil$H=H.equil.B[complete.cases(H.equil.B$V1),1] ##adding the prey values to alt.equil dataframe
# colnames(H.starr.final)=c('s','P','B','H')


alt.final=cbind(alt.equil$s,alt.equil$B,alt.equil$H,alt.equil$P) ##Let's rearrange them and name them
colnames(alt.final)=c('s','B','H','P')
alt.final=as.data.frame(alt.final)

##Let's find eigenvalues of the second equilbrium
e.values=NULL
for(i in 1: nrow(alt.final))
{
	s=alt.final[i,1] ##again, names should work here but I'm calling by numbers. (shakes head)
	B=alt.final[i,2]
	H=alt.final[i,3]
	P=alt.final[i,4]
	
# aa=a*H*(1+s)/(1+s+f*P)-2*B+1 #Bad derivation
aa=-(a*H*(s + 1))/(f*P+s+1)-2*B+1

# bb=a*B*(1+s)/(1+s+f*P) #Bad
 bb=-(a*B*(s+1))/(f*P+s+1)

# cc=a*f*B*H*(1+s)/(1+s+f*P)^2
 cc=(a*B*f*H*(s+1))/(f*P+s+1)^2

# dd=a*H*(1+s)/(1+s+f*P)
dd=(a*H*(s+1))/(f*P+s+1)

# ee=a*B/(1+f*P/(1+s))-b*P/((1+f*m*P/(1+s))*(1+n*s))-dh
ee=(a*B)/((f*P)/(s+1)+1)-(b*P)/((n*s+1)*((f*m*P)/(s+1)+1))-dh

# ff=-a*f*B*H/((1+s)*(1+f*P/(1+s))^2)-b*H/((1+n*s)*(1+(f*m*P)/(1+s)))+b*f*m*H*P/((1+s)*(1+n*s)*(1+(f*m*P)/(1+s))^2)
ff=-(a*B*f*H)/((s+1)*((f*P)/(s+1)+1)^2)-(b*H)/((n*s+1)*((f*m*P)/(s+1)+1))+(b*f*H*m*P)/((s+1)*(n*s+1)*((f*m*P)/(s+1)+1)^2)

gg=0

# hh=b*P*(1+s)/((1+s+f*m*P)*(1+n*s))
hh=(b*P*(s+1))/((n*s+1)*(f*m*P+s+1))

# ii=b*H/((1+f*m*P/(1+s))*(1+n*s))-b*f*m*H*P/((1+s)*(1+n*s)*(1+f*m*P/(1+s))^2)-dp
ii=(b*H)/((n*s+1)*((f*m*P)/(s+1)+1))-(b*f*H*m*P)/((s+1)*(n*s+1)*((f*m*P)/(s+1)+1)^2)-dp

jmat=(matrix(c(aa,bb,cc,dd,ee,ff,gg,hh,ii), ncol=3, byrow=T))

e.values=rbind(e.values, eigen(jmat)$values)
}
e.values.fin=cbind(alt.final$s,e.values)
e.values.fin=as.data.frame(e.values.fin)
e.values.fin

#So what do we expect? None of these should be stable (some positive eigen values)
