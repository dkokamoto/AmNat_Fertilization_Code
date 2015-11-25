
### load necessary libraries 
library(deSolve);library(animation);library(grid);library(lattice)
library(plyr);library(reshape2);library(fields);library(grid);library(RColorBrewer)

####################################
### define compartmentalized model
####################################

dyn_fert <- function (time, state, pars, N1,N2, Ds.y, Ds.x, De, dx) {
  state_array <- array(state, dim = c(N1,N2,5))
  Sperm <- state_array[,,1]
  Eu <-   state_array[,,2]
  Ev  <- state_array[,,3]
  Ep  <- state_array[,,4]
  Em  <- state_array[,,5]
  
  with (as.list(pars), {
    ## sperm, egg, 
    dSperm  <- - beta* Sperm * (Eu + Ev + Ep + Em) - r * Sperm
    dEu  <- - gamma * beta * Sperm * Eu
    dEv   <- gamma * beta * Sperm * Eu - delta * Ev - gamma * beta * Sperm * Ev
    dEp   <- gamma * beta * Sperm * Ev
    dEm   <- delta * Ev
    
    zero1 <- rep(0, N1)
    zero2 <- rep(0, N2)
    
    ## 1. flux to the left ; zero fluxes near boundaries
    FluxSperm <- -Ds.x * rbind((Sperm[1:N1,] - rbind(zero2,Sperm[1:(N1-1),])), zero2)/dx
    FluxEu <- -De * rbind(zero2,(Eu[2:N1,] - Eu[1:(N1-1),]), zero2)/dx
    FluxEv <- -De * rbind(zero2,(Ev[2:N1,] - Ev[1:(N1-1),]), zero2)/dx
    FluxEp <- -De * rbind(zero2,(Ep[2:N1,] - Ep[1:(N1-1),]), zero2)/dx
    FluxEm <- -De * rbind(zero2,(Em[2:N1,] - Em[1:(N1-1),]), zero2)/dx
    
    ## Add flux gradient to rate of change
    dSperm    <- dSperm - (FluxSperm[2:(N1+1),] - FluxSperm[1:N1,])/dx
    dEu       <- dEu - (FluxEu[2:(N1+1),] - FluxEu[1:N1,])/dx
    dEv       <- dEv - (FluxEv[2:(N1+1),] - FluxEv[1:N1,])/dx
    dEp       <- dEp - (FluxEp[2:(N1+1),] - FluxEp[1:N1,])/dx
    dEm       <- dEm - (FluxEm[2:(N1+1),] - FluxEm[1:N1,])/dx
    
    ## 2. Fluxes in y-direction; zero fluxes near right boundary
    FluxSperm <- -Ds.y * cbind(zero1,(Sperm[,2:N2] - Sperm[,1:(N2-1)]), zero1)/dx
    FluxEu    <- -De * cbind(zero1,(Eu[,2:N2] - Eu[,1:(N2-1)]), zero1)/dx
    FluxEv    <- -De * cbind(zero1,(Ev[,2:N2] - Ev[,1:(N2-1)]), zero1)/dx
    FluxEp    <- -De * cbind(zero1,(Ep[,2:N2] - Ep[,1:(N2-1)]), zero1)/dx
    FluxEm    <- -De * cbind(zero1,(Em[,2:N2] - Em[,1:(N2-1)]), zero1)/dx
    
    ## Add flux gradient to rate of change
    dSperm    <- dSperm - (FluxSperm[,2:(N2+1)] - FluxSperm[,1:N2])/dx
    dEu       <- dEu - (FluxEu[,2:(N2+1)] - FluxEu[,1:N2])/dx
    dEv       <- dEv - (FluxEv[,2:(N2+1)] - FluxEv[,1:N2])/dx
    dEp       <- dEp - (FluxEp[,2:(N2+1)] - FluxEp[,1:N2])/dx
    dEm       <- dEm - (FluxEm[,2:(N2+1)] - FluxEm[,1:N2])/dx
    
    return(list(c(as.vector(dSperm), as.vector(dEu),as.vector(dEv),as.vector(dEp),as.vector(dEm))))
  })
}


####################################
## Apply the model 
####################################

### define model parameters
pars    <- c(beta  = 0.00259*0.45,    # /day, rate of ingestion
             gamma = 0.08,    # /day, growth rate of prey
             delta  = 1.84,   # /day, mortality rate of predator
             r = 0.00027)    # mmol/m3, carrying capacity

### define spatial parameters           
N1  <- 151
N2 <- 51
R  <-100
dx <- 1
Ds.y <- 0.5  # units/s, dispersion coefficient
Ds.x <- .5 # units/s, dispersion coefficient
De <- 0  
NN <- N1*N2   # total number of 
max_time = 7200
step_size= 30
n_times <- length(seq(from= 0, to = max_time, by = step_size))

## array of initial conditions
r=25+0.7071068
x0 <- round(N1/5)
x02 <- round(N1/5*3)
x03 <- round(N1/5*3)
y0 <- round(N2/2)

### define egg clouds
r2=0.6
a <- melt(outer(1:N1,1:N2))[,c(1,2)]
names(a) <- c("x","y") 
select  <- matrix(with(a,r^2 >= (x-x0)^2 + (y-y0)^2),ncol= N2)
select2 <- matrix(with(a,r^2 >= (x-x02)^2 + (y-y0)^2),ncol= N2)
select3 <- matrix(with(a,r2^2 >= (x-r-r2-x03)^2 + (y-y0)^2),ncol= N2)
image(select|select2|select3)
sperm_max <- 400000
egg_max <- 1

### define initial conditions 
yini_array <- array(0,dim=c(N1,N2,5))
yini_array[,,1][select3] <-sperm_max
yini_array[,,2][select] <- egg_max
yini_array[,,2][select2] <- egg_max
yini <- as.vector(yini_array)

## solve models (5000 state variables...  use Cash-Karp Runge-Kutta method
times   <- seq(0, max_time, by = step_size)
out <- ode.2D(y = yini, times = times, func = dyn_fert, parms = pars,
              dimens = c(N1, N2), names = c("Sperm", "Eu","Ev","Ep","Em"),
              N1 = N1,N2=N2, dx = dx, Ds.x = Ds.x,Ds.y=Ds.y,De=De, method = rkMethod("rk45ck"))

yini_array[,,2][select2] <- 0
yini_array[,,2][select] <- egg_max
yini2 <- as.vector(yini_array)
out2 <- ode.2D(y = yini2, times = times, func = dyn_fert, parms = pars,
              dimens = c(N1, N2), names = c("Sperm", "Eu","Ev","Ep","Em"),
              N1 = N1,N2=N2, dx = dx, Ds.x = Ds.x,Ds.y=Ds.y,De=De, method = rkMethod("rk45ck"))

diagnostics(out)
summary(out)

### get results
Spermout <- apply(subset(out, select = "Sperm", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))
Emout <- apply(subset(out, select = "Em", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))
Epout <- apply(subset(out, select = "Ep", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))
Euout <- apply(subset(out, select = "Eu", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))

Spermout2 <- apply(subset(out2, select = "Sperm", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))
Emout2 <- apply(subset(out2, select = "Em", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))
Epout2 <- apply(subset(out2, select = "Ep", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))
Euout2 <- apply(subset(out2, select = "Eu", arr = TRUE),c(1,2,3), function(x) ifelse(x<0,0,x))

rgb <- rgb.palette(100)
rgb2 <- tim.colors(10)
rgb.palette <- colorRampPalette(c("black",rgb2[2:10]), space = "rgb", bias = 1)
rgb3 <-rgb.palette(100)
sperm_breaks <- seq(from=-3, to= 2, by= 0.25)
sperm_labs <- c(expression(paste(10^-3)),
                expression(paste(10^-2)),
                expression(paste(10^-1)),
                expression(paste(10^0)),
                expression(paste(10^1)),
                expression(paste(10^2)),
                expression(paste(10^3)))

sperm_lab_breaks <- seq(from=-4, to= 2, by= 1)
iso_loc <- c(55,N2/2)

my.panel.levelplot2 <- function(...) {
  panel.levelplot(...)
  panel.arrows(x0=iso_loc[1],x1=iso_loc[1],y0=70,y1=iso_loc[2]+1, 
               col = "grey50",fill= "grey50",lwd= 2.5,
               length = 0.1,type= "closed")
}

my.panel.levelplot <- function(...) {
  panel.levelplot(...)
  panel.arrows(x0=iso_loc[1],x1=iso_loc[1],y0=70,y1=iso_loc[2]+1, 
               col = "black",fill= "black",lwd= 2.5,
               length = 0.1,type= "closed")
}

df1 <- data.frame(list(time = rep(times,2), 
                       ufert=c(Euout[iso_loc[1],iso_loc[2],],Euout2[iso_loc[1],iso_loc[2],]),
                       mfert=c(Emout[iso_loc[1],iso_loc[2],],Emout2[iso_loc[1],iso_loc[2],]), 
                       pfert= c(Epout[iso_loc[1],iso_loc[2],],Epout2[iso_loc[1],iso_loc[2],]),
                       src = rep(c("2 Egg Clouds","1 Egg Cloud "),each= n_times),
                       sperm= c(Spermout[iso_loc[1],iso_loc[2],],Spermout2[iso_loc[1],iso_loc[2],])))

options <- theme(
  panel.grid.minor = element_line(colour = NA),
  panel.grid.major = element_line(colour = NA),
  panel.margin= unit(.5,"lines"),
  legend.key = element_rect(colour=NA,fill=NA),
  panel.grid = element_line(colour = NA),
  panel.background=element_rect(fill= NA),
  plot.background=element_rect(fill= NA),
  legend.background=element_rect(fill= NA, colour= NA),
  legend.key.height= unit(1,"lines"),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.title.align = 0,
  axis.text.y = element_text(angle=90,hjust=0.5),
  legend.direction= "vertical",
  plot.title= element_text(size= 12),
  panel.border= element_rect(fill= NA),
  legend.key.width = unit(1,"lines"),
  plot.margin = unit(c(0.05,0.05,0.05,0.05),"inches"))

for (i in 1:n_times){
  png(filename=paste0("fifth_examplec",sprintf("%03d", i) ,".png"), width=700,height=400)
  a <- levelplot(log10(ifelse(Spermout2[ , ,i]<0.0011,0.0011,ifelse(Spermout2[ , ,i]>1000,1000,Spermout2[ , ,i]))),col.regions= rgb3,useRaster=T,interpolate= T,at= seq(from =-3, to =3, length.out= 100),colorkey=list(at= seq(from =-3, to = 3, length.out= 100),interpolate = T,labels= list(at=sperm_lab_breaks,labels=sperm_labs)),scales= list(draw= F,cex= 0),ylab= "",xlab= "",panel= my.panel.levelplot2,between = list(x = 1),
                 par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))
  Em <- apply(Emout2[ , ,i],c(1,2),function(x)min(max(x,0.0001),0.9999))
  Em[!select] <- NA 
  Ep <- apply(Epout[ , ,i],c(1,2),function(x)min(max(x,0.0001),0.9999))
  Ep[!select] <- NA 
  Eu <- apply(Euout2[ , ,i],c(1,2),function(x)min(max(x,0.0001),0.9999))
  Eu[!select] <- NA 
  b <- levelplot(Em,col.regions= rgb3, at= seq(from =0, to = 1, length.out= 100),useRaster=T,interpolate= T,scales= list(draw= F,cex= 0),colorkey= list(interpolate= T),ylab= "",xlab= "",xlab.top= "",panel= my.panel.levelplot2,between = list(x = 1),
                 par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0),layout.widths=list(left.padding=1, right.padding=2)))
  c <- levelplot(Ep,col.regions= rgb3, at= seq(from =0, to = 1, length.out= 100),useRaster=T,interpolate= T,scales= list(draw= F,cex= 0),colorkey= list(interpolate= T),ylab= "",xlab= "",xlab.top= "",panel= my.panel.levelplot2,between = list(x = 1),
                 par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0),layout.widths=list(left.padding=1, right.padding=2)))
  d <- levelplot(Eu,col.regions= rgb3, at= seq(from =0, to = 1, length.out= 100),useRaster=T,interpolate= T,scales= list(draw= F,cex= 0),colorkey= list(interpolate= T),ylab= "",xlab= "",xlab.top= "",panel= my.panel.levelplot2,between = list(x = 1),
                 par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0),layout.widths=list(left.padding=1, right.padding=2)))

  a2 <- levelplot(log10(ifelse(Spermout[ , ,i]<0.00011,0.00011,ifelse(Spermout[ , ,i]>1000,1000,Spermout[ , ,i]))),col.regions= rgb3,useRaster=T,interpolate= T,at= seq(from =-4, to =3, length.out= 100),colorkey=list(at= seq(from =-4, to = 3, length.out= 100),interpolate = T,labels= list(at=sperm_lab_breaks,labels=sperm_labs)),scales= list(draw= F,cex= 0),ylab= "Sperm",xlab= "",panel= my.panel.levelplot,between = list(x = 1),
 par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))
  Em <- apply(Emout[ , ,i],c(1,2),function(x)min(max(x,0.0001),0.9999))
  Em[!select&!select2] <- NA 
  Ep <- apply(Epout[ , ,i],c(1,2),function(x)min(max(x,0.0001),0.9999))
  Ep[!select&!select2] <- NA 
  Eu <- apply(Euout[ , ,i],c(1,2),function(x)min(max(x,0.0001),0.9999))
  Eu[!select&!select2] <- NA 
  b2 <- levelplot(Em,col.regions= rgb3, at= seq(from =0, to = 1, length.out= 100),useRaster= T,interpolate= T,scales= list(draw= F,cex= 0),xlab= "",ylab= "Normal Fert.",colorkey= list(interpolate= T),panel= my.panel.levelplot,between = list(x = 1),
                  par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0),layout.widths=list(left.padding=1, right.padding=2)))
  c2 <- levelplot(Ep,col.regions= rgb3, at= seq(from =0, to = 1, length.out= 100),useRaster=T,interpolate= T,scales= list(draw= F,cex= 0),colorkey= list(interpolate= T),xlab= "",xlab.top= "",ylab= "Abnormal Fert.",panel= my.panel.levelplot,between = list(x = 1),
                  par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0),layout.widths=list(left.padding=1, right.padding=2)))
  d2 <- levelplot(Eu,col.regions= rgb3, at= seq(from =0, to = 1, length.out= 100),useRaster=T,interpolate= T,scales= list(draw= F,cex= 0),colorkey= list(interpolate= T),xlab= "",xlab.top= "",ylab= "Unfert.",panel= my.panel.levelplot,between = list(x = 1),
                  par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0),layout.widths=list(left.padding=1, right.padding=2)))
  a3 <- c(a2,a,layout= c(2,1))
  d3 <- c(d2,d,layout= c(2,1))
  b3 <- c(b2,b,layout= c(2,1))
  c3 <- c(c2,c,layout= c(2,1))
 
  
  a <-ggplot(aes(time,sperm),data= subset(df1,time%in%times[1:i]))+
    geom_line(aes(colour= src),size= 1.5)+theme_bw()+
    xlab("Time (seconds)")+
    scale_y_continuous(limits= c(0,15), expand= c(0,0))+
    ylab(expression(paste("Sperm (#/",mu, "l)")))+
    scale_x_continuous(limits= c(0,max_time),expand= c(0,0))+
    options+scale_colour_manual(values= c("grey50","black"))+
    theme(axis.text.x= element_blank(),legend.position = "none",axis.title.x= element_blank(), plot.margin = unit(c(0.2,0.05,0.05,0.05),"inches"))
  
  b <-ggplot(aes(time,ufert),data= subset(df1,time%in%times[1:i]))+
    geom_line(aes(colour= src),size= 1.5 )+theme_bw()+
    xlab("time (s)")+
    scale_y_continuous(limits= c(0,1), expand= c(0.05,0),breaks= c(0,1))+
    ylab("% Unfert.")+
    scale_x_continuous(limits= c(0,max_time),expand= c(0,0))+
    options+scale_colour_manual(values= c("grey50","black"))+
    theme(legend.position = "none",axis.title.x= element_blank(),axis.text.x= element_blank())
  
  c <-ggplot(aes(time,mfert),data= subset(df1,time%in%times[1:i]))+
    geom_line(aes(colour= src),size= 1.5)+theme_bw()+
    xlab("time (s)")+
    scale_y_continuous(limits= c(0,1), expand= c(0.05,0),breaks= c(0,1))+
    ylab("% Normal Fert.")+
    scale_x_continuous(limits= c(0,max_time),expand= c(0,0))+
    options+scale_colour_manual(values= c("grey50","black"))+
    theme(legend.position = "none",axis.title.x= element_blank(),axis.text.x= element_blank())
  
  d <-ggplot(aes(time,pfert),data= subset(df1,time%in%times[1:i]))+
    geom_line(aes(colour= src),size= 1.5)+theme_bw()+
    xlab("time (s)")+
    scale_y_continuous(limits= c(0,1), expand= c(0.05,0),breaks= c(0,1))+
    ylab("% Abnormal Fert.")+
    scale_x_continuous(limits= c(0,max_time),expand= c(0,0))+
    options+scale_colour_manual(values= c("grey50","black"))+
    theme(legend.position= c(.4,.8))
 
 gg1 <- arrangeGrob(a,b,c,d,ncol= 1,nrow= 5,heights= c(1.05,0.95,0.95,1.1,0.05),just= "bottom")
 gg2 <- arrangeGrob(a3,d3,b3,c3,ncol= 1,nrow= 5,heights= c(1,1,1,1,0.25))
  gg.a <- grid.arrange(gg2,gg1,ncol= 2,widths=c(2,0.8))
 
  grid.rect(x= unit(0.32,"npc"),width=0.29, hjust=1,y= unit(.093,"npc"),vjust= 0, height= 0.88,
            gp=gpar(col="black",fill= "NA",lwd=5))
  grid.rect(x= unit(0.635,"npc"),width=0.295, hjust=1,y= unit(.093,"npc"),vjust= 0, height= 0.88,
                  gp=gpar(col="grey40",fill= "NA",lwd=5))
 grid.text("2 Egg Clouds", x=unit(0.18,"npc"),hjust=0.5,y= unit(.07,"npc"))
 grid.text("1 Egg Cloud", x=unit(0.48,"npc"),hjust=0.5,y= unit(.07,"npc"))
 grid.text("Solutions at Arrows", x=unit(0.87,"npc"),hjust=0.5,y= unit(.98,"npc"))
 grid.text("A", x=unit(0.3,"npc"),hjust=0.5,y= unit(.93,"npc"))
 grid.text("B", x=unit(0.3,"npc"),hjust=0.5,y= unit(.7,"npc"))
 grid.text("C", x=unit(0.3,"npc"),hjust=0.5,y= unit(.47,"npc"))
 grid.text("D", x=unit(0.3,"npc"),hjust=0.5,y= unit(.23,"npc"))
 grid.text("E", x=unit(0.61,"npc"),hjust=0.5,y= unit(.93,"npc"))
 grid.text("F", x=unit(0.61,"npc"),hjust=0.5,y= unit(.7,"npc"))
 grid.text("G", x=unit(0.61,"npc"),hjust=0.5,y= unit(.47,"npc"))
 grid.text("H", x=unit(0.61,"npc"),hjust=0.5,y= unit(.23,"npc"))
 grid.text("I", x=unit(0.98,"npc"),hjust=0.5,y= unit(.93,"npc"))
 grid.text("J", x=unit(0.98,"npc"),hjust=0.5,y= unit(.7,"npc"))
 grid.text("K", x=unit(0.98,"npc"),hjust=0.5,y= unit(.47,"npc"))
 grid.text("L", x=unit(0.98,"npc"),hjust=0.5,y= unit(.23,"npc"))
dev.off()
}

