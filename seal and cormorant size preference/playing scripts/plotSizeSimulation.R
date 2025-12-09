
colors <- colorRampPalette(c("red", "blue"))(8)
faded_cols <- adjustcolor(colors, alpha.f = 0.2)  # 40% opacity
faded_cols


x <- (0:110)/10
c91 <- vbgrCod(x,1991)
c91sd <- vbgr.sdCod(x,1991,10000)

c92 <- vbgrCod(x,1992)
c92sd <- vbgr.sdCod(x,1992,10000)

c93 <- vbgrCod(x,1993)
c93sd <- vbgr.sdCod(x,1993,10000)

c94 <- vbgrCod(x,1994)
c94sd <- vbgr.sdCod(x,1994,10000)

c95 <- vbgrCod(x,1995)
c95sd <- vbgr.sdCod(x,1995,10000)

c96 <- vbgrCod(x,1996)
c96sd <- vbgr.sdCod(x,1996,10000)

c97 <- vbgrCod(x,1997)
c97sd <- vbgr.sdCod(x,1997,10000)

c98 <- vbgrCod(x,1998)
c98sd <- vbgr.sdCod(x,1998,10000)




plot(x+1990,c91,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[1],axes = FALSE,xlab="Year",ylab="Fish length [cm]")
axis(1, at = seq(1991,1999, by = 1))
axis(2)
box()
polygon(
  c(x+1990, rev(x+1990)),                          # x+1990 forward, then backward
  c(c91 - c91sd, rev(c91 + c91sd)),       # lower, then upper reversed
  col = faded_cols[1],                # semi-transparent red
  border = NA
)

lines(x+1990+1,c92,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[2])
polygon(
  c(x+1990+1, rev(x+1990+1)),                          # x+1990 forward, then backward
  c(c92 - c92sd, rev(c92 + c92sd)),       # lower, then upper reversed
  col = faded_cols[2],                # semi-transparent red
  border = NA
)

lines(x+1990+2,c93,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[3])
polygon(
  c(x+1990+2, rev(x+1990+2)),                          # x+1990 forward, then backward
  c(c93 - c93sd, rev(c93 + c93sd)),       # lower, then upper reversed
  col = faded_cols[3],                # semi-transparent red
  border = NA
)

lines(x+1990+3,c94,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[4])
polygon(
  c(x+1990+3, rev(x+1990+3)),                          # x+1990 forward, then backward
  c(c94 - c94sd, rev(c94 + c94sd)),       # lower, then upper reversed
  col = faded_cols[4],                # semi-transparent red
  border = NA
)

lines(x+1990+4,c95,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[5])
polygon(
  c(x+1990+4, rev(x+1990+4)),                          # x+1990 forward, then backward
  c(c95 - c95sd, rev(c95 + c95sd)),       # lower, then upper reversed
  col = faded_cols[5],                # semi-transparent red
  border = NA
)

lines(x+1990+5,c96,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[6])
polygon(
  c(x+1990+5, rev(x+1990+5)),                          # x+1990 forward, then backward
  c(c96 - c96sd, rev(c96 + c96sd)),       # lower, then upper reversed
  col = faded_cols[6],                # semi-transparent red
  border = NA
)

lines(x+1990+6,c97,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[7])
polygon(
  c(x+1990+6, rev(x+1990+6)),                          # x+1990 forward, then backward
  c(c97 - c97sd, rev(c97 + c97sd)),       # lower, then upper reversed
  col = faded_cols[7],                # semi-transparent red
  border = NA
)

lines(x+1990+7,c98,type='l',lwd=3,ylim=c(0,120),xlim=c(1991,1999),col=colors[8])
polygon(
  c(x+1990+7, rev(x+1990+7)),                          # x+1990 forward, then backward
  c(c98 - c98sd, rev(c98 + c98sd)),       # lower, then upper reversed
  col = faded_cols[8],                # semi-transparent red
  border = NA
)





















