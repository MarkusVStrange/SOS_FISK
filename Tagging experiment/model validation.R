#models 
# M0: p = p0
# M1: p = p0 + PIT
# M2: p = p0 + PIT + year
# M3: p = p0 + PIT + year_month
# M4: p = p0 + PIT + year_month + release_location
# M5: p = p0 + PIT + Random(year_month)
# M6: p = p0 + PIT + Random(release_day)
# M7: p = p0 + PIT + year + Random(release_day)
# M8: p = p0 + PIT + year + Random(release_day) + release_location

# Run all models
# Compare models
G_M1.M0 <- 2*(M0$obj-M1$obj)
cat('Model 1 better than model 0; p-value =',round(pchisq(G_M1.M0,df=length(M1$par)-length(M0$par),lower.tail = FALSE),5),'\n') 

G_M2.M1 <- 2*(M1$obj-M2$obj)
cat('Model 2 better than model 1; p-value =',round(pchisq(G_M2.M1,df=length(M2$par)-length(M1$par),lower.tail = FALSE),5),'\n') 

G_M3.M2 <- 2*(M2$obj-M3$obj)
cat('Model 3 better than model 2; p-value =',round(pchisq(G_M3.M2,df=length(M3$par)-length(M2$par),lower.tail = FALSE),5),'\n') 

G_M4.M3 <- 2*(M3$obj-M4$obj)
cat('Model 4 not better than model 3; p-value =',round(pchisq(G_M4.M3,df=length(M4$par)-length(M3$par),lower.tail = FALSE),5),'\n') 

# Random effect models - NB likelihood ratio test not comparable as they are not sub-models
G_M5.M3 <- 2*(M5$obj-M3$obj)
cat('Model 3 better than model 5; p-value =',round(pchisq(G_M5.M3,df=length(M3$par)-length(M5$par),lower.tail = FALSE),5),'\n') 

G_M6.M5 <- 2*(M5$obj-M7$obj)
cat('Model 6 better than model 5; p-value =',round(pchisq(G_M6.M5,df=length(M6$par)-length(M5$par),lower.tail = FALSE),5),'\n') 

G_M3.M6 <- 2*(M6$obj-M3$obj)
cat('Model 3 better than model 5; p-value =',round(pchisq(G_M3.M6,df=length(M3$par)-length(M6$par),lower.tail = FALSE),5),'\n') 

# AICc
######
l0=M0[["objective"]]
k0=length(M0[["par"]])
AICc0 = 2*k0 + 2*l0 + 2*k0*(k0+1)/(length(dat$PIT)-k0-1)
l1=M1[["objective"]]
k1=length(M1[["par"]])
AICc1 = 2*k1 + 2*l1 + 2*k1*(k1+1)/(length(dat$PIT)-k1-1)
l2=M2[["objective"]]
k2=length(M2[["par"]])
AICc2 = 2*k2 + 2*l2 + 2*k2*(k2+1)/(length(dat$PIT)-k2-1)
l3=M3[["objective"]]
k3=length(M3[["par"]])
AICc3 = 2*k3 + 2*l3 + 2*k3*(k3+1)/(length(dat$PIT)-k3-1)
l4=M4[["objective"]]
k4=length(M4[["par"]])
AICc4 = 2*k4 + 2*l4 + 2*k4*(k4+1)/(length(dat$PIT)-k4-1)
l5=M5[["objective"]]
k5=length(M5[["par"]])
AICc5 = 2*k5 + 2*l5 + 2*k5*(k5+1)/(length(dat$PIT)-k5-1)
l6=M6[["objective"]]
k6=length(M6[["par"]])
AICc6 = 2*k6 + 2*l6 + 2*k6*(k6+1)/(length(dat$PIT)-k6-1)
l7=M7[["objective"]]
k7=length(M7[["par"]])
AICc7 = 2*k7 + 2*l7 + 2*k7*(k7+1)/(length(dat$PIT)-k7-1)
l8=M8[["objective"]]
k8=length(M8[["par"]])
AICc8= 2*k8 + 2*l8 + 2*k8*(k8+1)/(length(dat$PIT)-k8-1)

cat('Model 0 AICce =',round(AICc0,1),'\n') 
cat('Model 1 AICce =',round(AICc1,1),'\n') 
cat('Model 2 AICce =',round(AICc2,1),'\n') 
cat('Model 3 AICce =',round(AICc3,1),'\n') 
cat('Model 4 AICce =',round(AICc4,1),'\n') 
cat('Model 5 AICce =',round(AICc5,1),'\n') 
cat('Model 6 AICce =',round(AICc6,1),'\n') 
cat('Model 7 AICce =',round(AICc7,1),'\n')
cat('Model 8 AICce =',round(AICc8,1),'\n')

#####

# Plot model diagnostics
#####
par(mfrow=c(3,3))
qqnorm(res0$residual,main=paste("model 0: p = ",round(shapiro.test(res0$residual)$p.value,3)))
qqline(res0$residual)
qqnorm(res1$residual,main=paste("model 1: p = ",round(shapiro.test(res1$residual)$p.value,3)))
qqline(res1$residual)
qqnorm(res2$residual,main=paste("model 2: p = ",round(shapiro.test(res2$residual)$p.value,3)))
qqline(res2$residual)
qqnorm(res3$residual,main=paste("model 3: p = ",round(shapiro.test(res3$residual)$p.value,3)))
qqline(res3$residual)
qqnorm(res4$residual,main=paste("model 4: p = ",round(shapiro.test(res4$residual)$p.value,3)))
qqline(res4$residual)
qqnorm(res5$residual,main=paste("model 5: p = ",round(shapiro.test(res5$residual)$p.value,3)))
qqline(res5$residual)
qqnorm(res6$residual,main=paste("model 6: p = ",round(shapiro.test(res6$residual)$p.value,3)))
qqline(res6$residual)
qqnorm(res7$residual,main=paste("model 7: p = ",round(shapiro.test(res7$residual)$p.value,3)))
qqline(res7$residual)
qqnorm(res8$residual,main=paste("model 8: p = ",round(shapiro.test(res8$residual)$p.value,3)))
qqline(res8$residual)
#####