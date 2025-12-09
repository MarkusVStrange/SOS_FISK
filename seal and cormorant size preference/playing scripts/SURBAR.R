library(remotes)
#remotes::install_github("einarhjorleifsson/surbar")
library(surbar)
library(plyr)
library(reshape2)
library(ggplot2)
?surba_fit


##########################
# NS whiting
##########################
#####
Path <- "extdata"
sw.file <- "whi47d_sw.dat"
mat.file <- "whi47d_mat.dat"
survey.file <- "whi47d_ef_surba.dat"
zbarage.1 <- 2
zbarage.2 <- 4
ref.age <- 3
lambda <- 1
#wt.vec <- rep(1, 6)

# Read in stock weight data
tmp <- read.vpa.file(system.file(paste(Path,
                                       sw.file, sep = "/"), package = "surbar"))
sw <- data.frame(tmp$tab)
names(sw) <- c(tmp$a[1]:tmp$a[2])
rownames(sw) <- c(tmp$y[1]:tmp$y[2])
sw <- sw[, 2:6]
# Read in maturity data
tmp <- read.vpa.file(system.file(paste(Path,
                                       sw.file, sep = "/"), package = "surbar"))
mat <- data.frame(tmp$tab)
names(mat) <- c(tmp$a[1]:tmp$a[2])
rownames(mat) <- c(tmp$y[1]:tmp$y[2])
mat <- mat[, 2:6]
# Read in survey data
tmp <- read.survey.file(system.file(paste(Path,
                                          survey.file, sep = "/"), package = "surbar"))
idx <- tmp$idx
surveyNames <- unlist(lapply(idx, function(wk) {
  wk$name
}))
surveyTime <- unlist(lapply(idx, function(wk) {
  wk$rho
}))
u1 <- list(name = surveyNames[1], rho = surveyTime[1],
           tab = idx[[1]]$tab)
u2 <- list(name = surveyNames[2], rho = surveyTime[2],
           tab = idx[[2]]$tab)
dat <- list()
dat$ref.age <- ref.age
dat$zbarage.1 <- zbarage.1
dat$zbarage.2 <- zbarage.2
dat$lamda <- lambda
dat$sw <- sw
dat$mat <- mat
dat$u <- list(u1, u2)
Dat.nsWhiting <- dat
Surba.nsWhiting <- surba_setup(Dat.nsWhiting)


Fit.nsWhiting <- surba_fit(Surba.nsWhiting)
Retro.nsWhiting <- surba_retro(Surba.nsWhiting,
                               nYears = 10)



naa <- surba_rby(Fit.nsWhiting)
dim(Fit.nsWhiting$u)
head(Fit.nsWhiting$u)





Fit.nsWhiting$u[[1]]
Fit.nsWhiting$u[[2]]
# Results
FIT <- Fit.nsWhiting
RET <- Retro.nsWhiting


p <- ggplot(FIT$res, aes(year, value, colour = factor(age))) +
  facet_wrap(~survey, ncol = 1) + labs(x = "",
                                       y = "Log residuals", colour = "Age") +
  scale_color_brewer(palette = "Set3")
p + geom_line()



parF <- melt(FIT$boot$par$f)
names(parF) <- c("iter", "year", "value")
parF <- ddply(parF, c("year"), summarise,
              q05 = quantile(value, 0.05), q50 = quantile(value,0.5), q95 = quantile(value, 0.95),
              ave = mean(value))
parF$variable <- "Zmort"
ggplot(parF, aes(year)) + geom_ribbon(aes(ymin = q05,
                                          ymax = q95), fill = "grey") + geom_line(aes(y = q50)) +
  geom_point(aes(y = ave), colour = "red",
             size = 2) + facet_wrap(~variable,
                                    scale = "free_y") + labs(x = "", y = "")




rbyBootQ <- melt(FIT$boot$rby[, 1:5], id.vars = "year")
rbyBootQ <- ddply(rbyBootQ, c("year", "variable"),
                  summarise, q05 = quantile(value, 0.05),
                  q50 = quantile(value, 0.5), q95 = quantile(value,
                                                             0.95), ave = mean(value))
p <- ggplot(rbyBootQ, aes(year)) + geom_ribbon(aes(ymin = q05,
                                                   ymax = q95), fill = "grey80") + geom_line(aes(y = q50)) +
  geom_point(aes(y = ave), colour = "red",
             size = 2) + facet_wrap(~variable,
                                    scale = "free_y") + labs(x = "", y = "") +
  geom_point(data = melt(FIT$rby, id.var = "year"),
             aes(year, value), col = "green",
             shape = 3)
p
#####
Path <- paste(getwd(),"SURBAR data",sep="/")
sw.file <- "sw.dat"
mat.file <- "mo.dat"
survey.file <- "survey.dat"
zbarage.1 <- 1
zbarage.2 <- 3
ref.age <- 2
lambda <- 1
#wt.vec <- rep(1, 6)

# Read in stock weight data
tmp <- read.vpa.file(paste(Path,sw.file, sep = "/"))
sw <- data.frame(tmp$tab)
names(sw) <- c(tmp$a[1]:tmp$a[2])
rownames(sw) <- c(tmp$y[1]:tmp$y[2])
sw <- sw[15:40, 1:5]
# Read in maturity data
tmp <- read.vpa.file(paste(Path,mat.file, sep = "/"))
mat <- data.frame(tmp$tab)
names(mat) <- c(tmp$a[1]:tmp$a[2])
rownames(mat) <- c(tmp$y[1]:tmp$y[2])
mat <- mat[15:40, 1:5]
# Read in survey data
tmp <- read.survey.file(paste(Path,survey.file, sep = "/"))
idx <- tmp$idx
surveyNames <- unlist(lapply(idx, function(wk) {
  wk$name
}))
surveyTime <- unlist(lapply(idx, function(wk) {
  wk$rho
}))
u1 <- list(name = surveyNames[1], rho = surveyTime[1],
           tab = idx[[1]]$tab)
u2 <- list(name = surveyNames[2], rho = surveyTime[2],
           tab = idx[[2]]$tab)
dat <- list()
dat$ref.age <- ref.age
dat$zbarage.1 <- zbarage.1
dat$zbarage.2 <- zbarage.2
dat$lamda <- lambda
dat$sw <- sw
dat$mat <- mat
dat$u <- list(u1, u2)
Dat.nsWhiting <- dat
Surba.nsWhiting <- surba_setup(Dat.nsWhiting)


Fit.nsWhiting <- surba_fit(Surba.nsWhiting)
Retro.nsWhiting <- surba_retro(Surba.nsWhiting,
                               nYears = 10)

##########################
# WB cod
##########################
#####




naa <- surba_rby(Fit.nsWhiting)
dim(Fit.nsWhiting$u)
head(Fit.nsWhiting$u)





Fit.nsWhiting$u[[1]]
Fit.nsWhiting$u[[2]]
# Results
FIT <- Fit.nsWhiting
RET <- Retro.nsWhiting





s_names <- grep("^s", names(FIT$params0), value = TRUE)
s <- FIT$params0[s_names]



Z <- FIT$fit$Z







p <- ggplot(FIT$res, aes(year, value, colour = factor(age))) +
  facet_wrap(~survey, ncol = 1) + labs(x = "",
                                       y = "Log residuals", colour = "Age") +
  scale_color_brewer(palette = "Set3")
p + geom_line()



parF <- melt(FIT$boot$par$f)
names(parF) <- c("iter", "year", "value")
parF <- ddply(parF, c("year"), summarise,
              q05 = quantile(value, 0.05), q50 = quantile(value,0.5), q95 = quantile(value, 0.95),
              ave = mean(value))
parF$variable <- "Zmort"
ggplot(parF, aes(year)) + geom_ribbon(aes(ymin = q05,
                                          ymax = q95), fill = "grey") + geom_line(aes(y = q50)) +
  geom_point(aes(y = ave), colour = "red",
             size = 2) + facet_wrap(~variable,
                                    scale = "free_y") + labs(x = "", y = "")




rbyBootQ <- melt(FIT$boot$rby[, 1:5], id.vars = "year")
rbyBootQ <- ddply(rbyBootQ, c("year", "variable"),
                  summarise, q05 = quantile(value, 0.05),
                  q50 = quantile(value, 0.5), q95 = quantile(value,
                                                             0.95), ave = mean(value))
p <- ggplot(rbyBootQ, aes(year)) + facet_wrap(~variable,
                                    scale = "free_y") + labs(x = "", y = "") +
  geom_point(data = melt(FIT$rby, id.var = "year"),
             aes(year, value), col = "green",
             shape = 3)
p

#####

function (filename) 
{
  wk.n <- scan(filename, skip = 1, nlines = 1, quiet = TRUE) - 
    100
  wk.idx <- vector("list", length = wk.n)
  wk.start <- 3
  for (wk.k in 1:wk.n) {
    wk.idx[[wk.k]]$name <- paste(scan(filename, skip = wk.start - 
                                        1, nlines = 1, what = character(0), quiet = TRUE), 
                                 collapse = " ")
    wk.temp <- scan(filename, skip = wk.start, nlines = 1, 
                    quiet = TRUE)
    wk.idx[[wk.k]]$y1 <- wk.temp[1]
    wk.idx[[wk.k]]$y2 <- wk.temp[2]
    wk.idx[[wk.k]]$ny <- wk.temp[2] - wk.temp[1] + 1
    wk.temp <- scan(filename, skip = wk.start + 1, nlines = 1, 
                    quiet = TRUE)
    wk.idx[[wk.k]]$rho <- 0.5 * (wk.temp[4] + wk.temp[3])
    wk.temp <- scan(filename, skip = wk.start + 2, nlines = 1, 
                    quiet = TRUE)
    wk.idx[[wk.k]]$a1 <- wk.temp[1]
    wk.idx[[wk.k]]$a2 <- wk.temp[2]
    wk.idx[[wk.k]]$na <- wk.temp[2] - wk.temp[1] + 1
    wk.idx[[wk.k]]$tab <- read.table(filename, skip = wk.start + 
                                       3, nrows = wk.idx[[wk.k]]$ny)
    wk.temp <- wk.idx[[wk.k]]$tab[, 2:(wk.idx[[wk.k]]$na + 
                                         1)]
    wk.effort <- wk.idx[[wk.k]]$tab[, 1]
    wk.idx[[wk.k]]$tab <- data.frame(wk.temp/wk.effort)
    names(wk.idx[[wk.k]]$tab) <- wk.idx[[wk.k]]$a1:wk.idx[[wk.k]]$a2
    rownames(wk.idx[[wk.k]]$tab) <- wk.idx[[wk.k]]$y1:wk.idx[[wk.k]]$y2
    wk.start <- wk.start + 4 + wk.idx[[wk.k]]$ny
  }
  list(n = wk.n, idx = wk.idx)
}