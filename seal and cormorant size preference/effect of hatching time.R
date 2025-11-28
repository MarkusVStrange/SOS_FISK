

c1 <- as.data.frame(table(round(rnorm(exp(10),3,5))))
c1$Var1 <- as.numeric(as.character(c1$Var1))
c1$cohort <- 2000

c2 <- as.data.frame(table(round(rnorm(exp(9),20,5))))
c2$Var1 <- as.numeric(as.character(c2$Var1))
c2$cohort <- 1999

c3 <- as.data.frame(table(round(rnorm(exp(8),30,8))))
c3$Var1 <- as.numeric(as.character(c3$Var1))
c3$cohort <- 1998

c4 <- as.data.frame(table(round(rnorm(exp(7),47,10))))
c4$Var1 <- as.numeric(as.character(c4$Var1))
c4$cohort <- 1997

c5 <- as.data.frame(table(round(rnorm(exp(6),55,12))))
c5$Var1 <- as.numeric(as.character(c5$Var1))
c5$cohort <- 1996

c6 <- as.data.frame(table(round(rnorm(exp(5),58,13))))
c6$Var1 <- as.numeric(as.character(c6$Var1))
c6$cohort <- 1995

c7 <- as.data.frame(table(round(rnorm(exp(4),60,14))))
c7$Var1 <- as.numeric(as.character(c7$Var1))
c7$cohort <- 1994


library(ggplot2)
df <- rbind(c1,c2,c3,c4,c5,c6,c7)
df$cohort <- factor(df$cohort)
ggplot(df %>% filter(Var1>0), aes(x = Var1, y = Freq,fill=cohort)) +
  geom_bar(stat = "identity",width=1) +
  labs(x = "Length [cm]", y = "Frequency") +
  theme_minimal()

######################


c1 <- as.data.frame(table(round(rnorm(exp(10),12,5))))
c1$Var1 <- as.numeric(as.character(c1$Var1))
c1$cohort <- 2000

c2 <- as.data.frame(table(round(rnorm(exp(9),20,5))))
c2$Var1 <- as.numeric(as.character(c2$Var1))
c2$cohort <- 1999

c3 <- as.data.frame(table(round(rnorm(exp(8),30,8))))
c3$Var1 <- as.numeric(as.character(c3$Var1))
c3$cohort <- 1998

c4 <- as.data.frame(table(round(rnorm(exp(7),47,10))))
c4$Var1 <- as.numeric(as.character(c4$Var1))
c4$cohort <- 1997

c5 <- as.data.frame(table(round(rnorm(exp(6),55,12))))
c5$Var1 <- as.numeric(as.character(c5$Var1))
c5$cohort <- 1996

c6 <- as.data.frame(table(round(rnorm(exp(5),58,13))))
c6$Var1 <- as.numeric(as.character(c6$Var1))
c6$cohort <- 1995

c7 <- as.data.frame(table(round(rnorm(exp(4),60,14))))
c7$Var1 <- as.numeric(as.character(c7$Var1))
c7$cohort <- 1994


library(ggplot2)
df <- rbind(c1,c2,c3,c4,c5,c6,c7)
df$cohort <- factor(df$cohort)
ggplot(df %>% filter(Var1>0), aes(x = Var1, y = Freq,fill=cohort)) +
  geom_bar(stat = "identity",width=1) +
  labs(x = "Length [cm]", y = "Frequency") +
  theme_minimal()
