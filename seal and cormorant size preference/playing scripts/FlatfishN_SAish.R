
keyLogFpar0 <- -15.378	
keyLogFpar1 <- -11.702
keyLogFpar2 <- -10.329
keyLogFpar3 <- -10.441	
keyLogFpar4 <- -10.195
keyLogFpar5 <- -12.490
keyLogFpar6 <- -10.710
keyLogFpar7 <- -10.143
keyLogFpar8 <- -10.212

qQ1 <- c(keyLogFpar0,keyLogFpar1,mean(c(keyLogFpar2,keyLogFpar2,keyLogFpar2,
                                        keyLogFpar3,keyLogFpar4))) 
qQ4 <- c(keyLogFpar5,keyLogFpar6,mean(c(keyLogFpar7,keyLogFpar8,keyLogFpar8,
                                        keyLogFpar8,keyLogFpar8)))
plot(0:6,exp(c(qQ1,rep(qQ1[3],4))),pch=19,ylim=exp(c(min(c(qQ1,qQ4)),max(c(qQ1,qQ4)))))
points(0:6,exp(c(qQ4,rep(qQ4[3],4))),pch=19,col="red")
