####9/17/19 This script is to look for eigenvalues in your Jacobian matrix. Solving for this will let me know if my solutions are stable or not (doesn't say if it is the only stable one)
###This is the first part of the stability analysis. The script 'solving for hyseresis' comes later
##The following are the partial derivatives of my jacobian matrix (corrected for the flipped m). 
###And that's mostly it. This has the up to date parital derivative (9/17/19) that will be used for stability analyses
aa=-a*H*(1+s)/(1+s+f*m*P)-2*B+1 ###This is the partial derivative of the basal resource equation wrt the basal resource. The reason why this looks weird is because I flipped the (1+s) to the numerator to keep the equations neater. *shrug
bb=-a*B*(1+s)/(1+s+f*m*P) #This is the partial derivative of the basal resource wrt the prey population
cc=a*f*m*B*H*(1+s)/(1+s+f*m*P)^2 #This is the partial derivative of the basal resource wrt the predator population
dd=a*H*(1+s)/(1+s+f*m*P) #This is the partial derivative of the prey wrt the basal resource
ee=a*B/(1+((f*m*P)/(1+s)))-b*P/((1+f*P/(1+s))*(1+n*s))-dh #This is the partial derivative of the prey wrt the prey
ff=-a*f*m*B*H/((1+s)*(1+f*m*P/(1+s))^2)-b*H/((1+n*s)*(1+(f*P)/(1+s)))+b*f*H*P/((1+s)*(1+n*s)*(1+(f*P)/(1+s))^2) #This is the partial derivative of the prey wrt the predator
gg=0 #This is the partial derivative of the predator equation wrt the basal resource
hh=b*P*(1+s)/((1+s+f*P)*(1+n*s)) #This is the partial derivative of the predator equation wrt the prey
ii=b*H/((1+f*P/(1+s))*(1+n*s))-b*f*H*P/((1+s)*(1+n*s)*(1+f*P/(1+s))^2)-dp #This is the partial derivative of the predator equation wrt the predator

##So what do I do with this? Input values for these 9 equations and find the eigenvalues

jmat=(matrix(c(aa,bb,cc,dd,ee,ff,gg,hh,ii), ncol=3, byrow=T))
eigen(jmat)$values
