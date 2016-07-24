calculate the LD decay fitted values
```
cal_LD<-function(distance,LD.data,n){
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
ld.df<-data.frame(distance,fpoints)
ld.df<-ld.df[order(ld.df$distance),]
plot(distance,LD.data,pch=20,cex=0.9,xlab="Distance (bp)",ylab=expression(r^2))
#axis(1,at=c(0,10000000,20000000,30000000),labels=c(0,10,20,30))
#axis(2,las=2);
lines(ld.df$distance,ld.df$fpoints,lty=1,lwd=1.2,col="red")
return(ld.df)
```
References:
Structure of linkage disequilibrium and phenotypic associations in the maize genome.PNAS,September 25, 2001.
