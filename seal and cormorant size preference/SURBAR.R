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

dat <- Surba.nsWhiting
Fit <- dat$fit
y1 <- min(dat$years)
y2 <- max(dat$years)
ny <- y2 - y1 + 1
na <- length(dat$ages)
a1 <- min(dat$ages)
a2 <- max(dat$ages)
ref.age <- dat$ref.age
f <- dat$par$f
s <- dat$par$s
r <- dat$par$r
zbarage.1 <- dat$zbarage.1
zbarage.2 <- dat$zbarage.2
sw <- dat$sw
mat <- dat$mat
x.eigen <- eigen(Fit$hessian, only.values = TRUE)$values
if (length(x.eigen[x.eigen == 0]) > 0) {
  stop("At least one parameter cannot be estimated: \nCheck that there are data for each cohort.")
}
n.psim <- 1000
x.psim <- mvrnorm(n = n.psim, mu = Fit$par, vcov(Fit))
x.stock <- data.frame(array(NA, dim = c(ny, 5)))
names(x.stock) <- c("year", "rec", "ssb", "tsb", "meanz")
x.psim.stock <- vector("list", n.psim)
x.psim.stock <- lapply(x.psim.stock, function(wk) {
  wk <- x.stock
})
x.psim.s <- array(NA, dim = c(1000, na))
x.psim.s[, 1:(ref.age - 1)] <- x.psim[, 1:(ref.age - 1)]
x.psim.s[, ref.age] <- 1
x.psim.s[, (ref.age + 1):(na - 1)] <- x.psim[, ref.age:(na - 
                                                          2)]
x.psim.s[, na] <- x.psim.s[, na - 1]
x.psim.f <- array(NA, dim = c(1000, ny))
x.psim.f[, 1:(ny - 1)] <- x.psim[, (na - 1):(na + ny - 3)]
x.psim.f[, ny] <- apply(x.psim.f[, (ny - 3):(ny - 1)], 1, 
                        mean)
x.psim.r <- x.psim[, (na + ny - 2):length(Fit$par)]
x.s <- rep(NA, length = na)
x.f <- rep(NA, length = ny)
for (i in 1:n.psim) {
  x.s[1:(ref.age - 1)] <- x.psim[i, 1:(ref.age - 1)]
  x.s[ref.age] <- 1
  x.s[(ref.age + 1):(na - 1)] <- x.psim[i, ref.age:(na - 
                                                      2)]
  x.s[na] <- x.s[na - 1]
  x.f[1:(ny - 1)] <- x.psim[i, (na - 1):(na + ny - 3)]
  x.f[ny] <- mean(x.f[(ny - 3):(ny - 1)])
  x.r <- x.psim[i, (na + ny - 2):length(Fit$par)]
  zmort <- x.f %o% x.s
  lnn <- array(NA, dim = dim(zmort))
  lnn[1, ] <- rev(x.r[1:dim(lnn)[2]])
  lnn[2:dim(lnn)[1], 1] <- x.r[(dim(lnn)[2] + 1):length(x.r)]
  for (jj in 2:dim(zmort)[2]) {
    for (ii in 2:dim(zmort)[1]) {
      lnn[ii, jj] <- lnn[ii - 1, jj - 1] - zmort[ii - 
                                                   1, jj - 1]
    }
  }
  n <- exp(lnn)
  x.psim.stock[[i]]$year <- y1:y2
  x.psim.stock[[i]]$meanz <- apply(zmort[, zbarage.1:zbarage.2], 
                                   1, mean)
  x.psim.stock[[i]]$z <- zmort
  x.psim.stock[[i]]$rec <- exp(x.r[na:length(r)])
  x.psim.stock[[i]]$ssb <- apply(n * sw * mat, 1, sum)
  x.psim.stock[[i]]$tsb <- apply(n * sw, 1, sum)
}
rby <- do.call(rbind, x.psim.stock)[, 1:5]
rby$iter <- rep(1:n.psim, each = ny)
z <- do.call(rbind, x.psim.stock)[, 6]
z <- as.data.frame(z)
names(z) <- a1:a2
z$year <- rby$year
z$iter <- rep(1:n.psim, each = ny)
rownames(x.psim.s) <- 1:n.psim
colnames(x.psim.s) <- a1:a2
rownames(x.psim.f) <- 1:n.psim
colnames(x.psim.f) <- y1:y2
rownames(x.psim.r) <- 1:n.psim
colnames(x.psim.r) <- (y1 - (a2 - a1) - a1):(y2 - a1)
par <- list(s = x.psim.s, f = x.psim.f, r = x.psim.r)
ret <- list(rby = rby, z = z, par = par)
return(ret)



function (wk.p, dat = NULL) 
{
  x <- dat$u.std
  numk <- dat$numk
  na <- length(dat$ages)
  ny <- length(dat$years)
  ref.age <- dat$ref.age
  y1 <- min(dat$years)
  qval <- dat$qval
  rho <- dat$rho
  wt <- dat$wt
  params0 <- dat$params0
  lambda <- dat$lamda
  sw <- dat$sw
  mat <- dat$mat
  a1 <- min(dat$ages)
  a2 <- max(dat$ages)
  y2 <- max(dat$years)
  wk.nk <- numk
  wk.na <- na
  wk.ny <- ny
  wk.s <- rep(NA, length = wk.na)
  wk.s[1:(ref.age - 1)] <- wk.p[1:(ref.age - 1)]
  wk.s[ref.age] <- 1
  wk.s[(ref.age + 1):(wk.na - 1)] <- wk.p[ref.age:(wk.na - 
                                                     2)]
  wk.s[wk.na] <- wk.s[wk.na - 1]
  wk.f <- rep(NA, length = wk.ny)
  wk.f[1:(wk.ny - 1)] <- wk.p[(wk.na - 1):(wk.na + wk.ny - 
                                             3)]
  wk.f[wk.ny] <- mean(wk.f[(wk.ny - 3):(wk.ny - 1)])
  wk.r <- rep(NA, length = wk.na + wk.ny - 1)
  wk.r[1:(wk.na + wk.ny - 1)] <- wk.p[(wk.na + wk.ny - 2):length(wk.p)]
  wk.z <- wk.f %o% wk.s
  wk.n <- array(NA, dim = dim(wk.z))
  wk.n[1, ] <- rev(wk.r[1:wk.na])
  wk.n[2:wk.ny, 1] <- wk.r[(wk.na + 1):length(wk.r)]
  vecs <- array(NA, dim = c(wk.na * wk.ny, 4))
  vecs[, 1] <- matrix(data = wk.n, nrow = wk.na * wk.ny, ncol = 1)
  vecs[, 2] <- rep(1:wk.na, each = wk.ny)
  vecs[, 3] <- rep(y1:(y1 + wk.ny - 1), wk.na)
  vecs[, 4] <- vecs[, 3] - vecs[, 2]
  cz.list <- tapply(wk.z, vecs[, 4], cumsum)
  vecs.list <- lapply(levels(as.factor(vecs[, 4])), function(wk) {
    temp <- vecs[vecs[, 4] == wk, ]
    temp.rep <- dim(temp)[1]
    if (!is.null(temp.rep)) {
      temp[, 1] <- rep(temp[1], temp.rep)
    }
    temp
  })
  vecs.list <- lapply(vecs.list, function(wk) {
    temp.rep <- dim(wk)[1]
    if (!is.null(temp.rep)) {
      wk.zz <- unlist(cz.list[as.character(wk[1, 4])])
      wk <- cbind(wk, wk.zz)
      wk.a <- dim(wk)[1]
      wk[2:wk.a, 1] <- wk[1:(wk.a - 1), 1] - wk[1:(wk.a - 
                                                     1), 5]
    }
    else {
      wk.zz <- unlist(cz.list[as.character(wk[4])])
      wk <- c(wk, wk.zz)
    }
    wk
  })
  vecs.table <- do.call(rbind, vecs.list)
  vecs.table <- vecs.table[order(vecs.table[, 3], vecs.table[, 
                                                             2]), ]
  new.n <- matrix(vecs.table[, 1], nrow = wk.ny, ncol = wk.na, 
                  byrow = TRUE)
  wk.n <- exp(new.n)
  i.hat <- vector("list", length = wk.nk)
  i.dash.star <- vector("list", length = wk.nk)
  for (k in 1:wk.nk) {
    i.hat[[k]] <- wk.n * array(unlist(qval[[k]]), dim = dim(wk.n))
    i.dash.star[[k]] <- array(unlist(x[[k]] * exp(wk.z * 
                                                    rho[k])), dim = dim(wk.n))
  }
  res1 <- vector("list", length = wk.nk)
  out0 <- vector("list", length = wk.nk)
  for (k in 1:wk.nk) {
    res1[[k]] <- sqrt(wt[[k]]) * (log(i.dash.star[[k]]) - 
                                    log(i.hat[[k]]))
    out0[[k]] <- array(NA, dim = c(wk.na * wk.ny, 1))
    out0[[k]][1:(wk.na * wk.ny), 1] <- unlist(res1[[k]])
    out0[[k]] <- as.vector(na.exclude(out0[[k]]))
  }
  out1 <- as.vector(unlist(out0))
  f1 <- wk.f[1:(wk.ny - 2)]
  f2 <- wk.f[2:(wk.ny - 1)]
  res2 <- sqrt(lambda) * (f1 - f2)
  out2 <- as.vector(res2)
  as.numeric(c(out1, out2))
}