# install.packages('rootSolve')
library('rootSolve')

###The purpose of this script is to check for the stability of the solutions to the tri-trophic model and whether alternative stable states exist. To do so, we first find the possible solutions to the tri-trophic model.

#This here is to set your parameters for m and n. I ran this looking at all four combinations of high/low for m and n corresponding to the four scenarios in the manuscript
low=1
high=2
m.domain=c(0.2,2)
n.domain=c(0.002, 25)

#Pick the right s.domain that corresponds to n level
s.domain=seq(from =0, to =0.3, by =0.001) #n=high, pred more sensitive
# s.domain=seq(from =0, to =2000, by =1) #n=low, prey more sensitive
H.star=matrix(nrow=length(s.domain), ncol=4) #creating a blank matrix to fill in root solutions (assuming that you don't get more than 4 equil values)

##The following forloop uses the uniroot.all function to find the the value of P such that the population growth rate of H is 0.
for(j in 1: length(s.domain))  
{
#Parameters
a=2.17
b=0.1
n=n.domain[high] 
m=m.domain[low]
f=1 
s=s.domain[j]
dh=0.056
dp=0.01

PforH.equil=function(P)
{
	H.equil=a*((-a*((dp*(n*s + 1) * (f*P + s + 1))/(b*(s+1)))*(s + 1) + f*m*P + s + 1)/(f*m*P+s+1))*((dp*(n*s + 1) * (f*P + s + 1))/(b* (s + 1)))/(1+(f*m*P)/(1+s))-b*((dp*(n*s + 1) * (f*P + s + 1))/(b* (s + 1)))*P/((1+(f*P)/(1+s))*(1+n*s))-dh*((dp*(n*s + 1) * (f*P + s + 1))/(b* (s + 1))) 
	return(H.equil)
} 
H.equil=PforH.equil(P.test) #We take P.test and run it through the function to see what kind of dh/dt values we get. We are identifying the predator numbers that cause H.equil (or dH/dt) to become 0.
sol=uniroot.all(PforH.equil, c(0,30)) 
# sol=max(uniroot.all(PforH.equil, c(0,30))) When both m and n are high, multiple solutions are found so this function allows us to test for the stability of those solutions, which turns out to be unstable
H.star[j,1:length(sol)]=sol 
}

H.star=cbind(s.domain,H.star) #Adding a row for sensory stress
H.starr=as.data.frame(H.star) #Making the matrix a dataframe 
head(H.starr) 


#This part of the script here is to remove negative solutions
H.starr.final=matrix(nrow=length(s.domain), ncol=5) ##Create new matrix of blanks
for(i in 1: nrow(H.starr))
{
	H.starr.temp=H.starr[i,][H.starr[i,]>=0] ### This will automatically filter out the sensory stress values where the predator population collapses and we don't have a non-zero solution
	H.starr.final[i,1:length(H.starr.temp)]=H.starr.temp 
}
H.starr.final=as.data.frame(H.starr.final)
head(H.starr.final)
H.starr.final=H.starr.final[complete.cases(H.starr.final$V2),] ##This should result in a data frame with predator numbers that result in the prey population being at equilibrium
##Now that we have the P values at equilibirum, let's solve for the corresponding H when dP/dt = 0

#Create blank matrix
P.star=matrix(nrow=length(s.domain), ncol=4)

for(j in 1: nrow(H.starr.final))
{
	s=H.starr.final$V1[j] #Call up relevant sensory stress values in H.starr.finals
	P=H.starr.final$V2[j] #Call up relevant P.equil values
HforP.equil=function(H) #Here's the function to solve for H.equil when P population isn't changing
{
	P.equil=b*H/((1+(f*P)/(1+s))*(1+n*s))-dp 
	return(P.equil)
}
P.equil=HforP.equil(H.test)
sol=uniroot.all(HforP.equil, c(0,200)) 
P.star[j,1:length(sol)]=sol 
}

P.star=as.data.frame(P.star)
P.star.final=P.star[complete.cases(P.star$V1),]

###Now that we have 2 data frames with solutions, combine them into final data set
final=NULL
final$s=H.starr.final$V1 #for sensory stress
final$P=H.starr.final$V2 #For p.equil
final$H=P.star.final$V1 #For h.equil
final=as.data.frame(final)

###With the solutions for P and H, we can calculate the abundance of the basal resource at equilbrium
B.star=matrix(nrow=length(s.domain), ncol=4)

for(j in 1: nrow(final))
{
	s=final$s[j]
	P=final$P[j]
	H=final$H[j]
BforB.equil=function(B)
{
	B.equil=B*(1-B)-a*B*H/(1+(f*m*P)/(1+s))	 #Function to solve B given H and P. Note that with the way this is set up, B = 0 is always a solution and should be ignored. 
	return(B.equil)
}
B.equil=BforB.equil(B.test)
sol=uniroot.all(BforB.equil, c(0,1)) 
B.star[j,1:length(sol)]=sol
}
B.star=as.data.frame(B.star)
B.star.final=B.star[complete.cases(B.star$V2),] #Remove NAs
final$B=B.star.final$V2 #Add to final data frame

##With the solutions for the three populations, we can ensure that they are correct by plugging them into the original equation. If they are the solution, the results should equal 0
B.test=NULL
H.test=NULL
P.test=NULL
for(k in 1:nrow(final))
{
s=final[k,1]
P=final[k,2]
H=final[k,3]
B=final[k,4] #Getting the respective B, H, and B values to test if they come out to 0 when plugged in

B.test[k]=B*(1-B)-a*H*B/(1+f*m*P/(1+s)) 
H.test[k]=a*B*H/(1+f*m*P/(1+s))-b*P*H/((1+f*P/(1+s))*(1+n*s))-dh*H
P.test[k]=b*P*H/((1+f*P/(1+s))*(1+n*s))-dp*P
}
round(B.test,4)
round(H.test,4) ##Round to get 0 values
round(P.test,4)


###Now, I'm creating the Jacobian matrix of the tri-trophic model and filling it with the solutions calculated 
e.values=NULL
for(i in 1: nrow(final))
{
	s=final[i,1]
	B=final[i,4] 
	H=final[i,3]
	P=final[i,2]

###The following are the partial derivatives of the tri-trophic model 
aa=-a*H*(1+s)/(1+s+f*m*P)-2*B+1 ###This is the partial derivative of the basal resource equation wrt the basal resource.

bb=-a*B*(1+s)/(1+s+f*m*P)  #This is the partial derivative of the basal resource wrt the prey population

cc=a*f*m*B*H*(1+s)/(1+s+f*m*P)^2 #This is the partial derivative of the basal resource wrt the predator population

dd=a*H*(1+s)/(1+s+f*m*P) #This is the partial derivative of the prey wrt the basal resource

ee=a*B/(1+((f*m*P)/(1+s)))-b*P/((1+f*P/(1+s))*(1+n*s))-dh #This is the partial derivative of the prey wrt the prey

ff=-a*f*m*B*H/((1+s)*(1+f*m*P/(1+s))^2)-b*H/((1+n*s)*(1+(f*P)/(1+s)))+b*f*H*P/((1+s)*(1+n*s)*(1+(f*P)/(1+s))^2) #This is the partial derivative of the prey wrt the predator

gg=0 #This is the partial derivative of the predator equation wrt the basal resource

hh=b*P*(1+s)/((1+s+f*P)*(1+n*s)) #This is the partial derivative of the predator equation wrt the prey

ii=b*H/((1+f*P/(1+s))*(1+n*s))-b*f*H*P/((1+s)*(1+n*s)*(1+f*P/(1+s))^2)-dp #This is the partial derivative of the predator equation wrt the predator

jmat=(matrix(c(aa,bb,cc,dd,ee,ff,gg,hh,ii), ncol=3, byrow=T))

e.values=rbind(e.values, eigen(jmat)$values)  #This function calculates the eigenvalues of the matrix
}
e.values.fin=cbind(final$s,e.values)
e.values.fin=as.data.frame(e.values.fin)

#Eigenvalues should all have negative real numbers for the solutions to be stable
head(e.values.fin)
Re(e.values.fin$V2)[Re(e.values.fin$V2)>0]
Re(e.values.fin$V3)[Re(e.values.fin$V3)>0]
Re(e.values.fin$V4)[Re(e.values.fin$V4)>0]

