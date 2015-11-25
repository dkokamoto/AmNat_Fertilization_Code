
library(deSolve);library(animation);library(grid);library(lattice)
library(plyr);library(reshape2);library(fields);library(grid);library(RColorBrewer)

dyn_fert <- function (t, y, parms) {
  state_vect <- y
  Sperm <- state_vect[1]
  Eu <-   state_vect[2]
  Ev  <- state_vect[3]
  Ep  <- state_vect[4]
  Em  <- state_vect[5]
  
  with (as.list(parms), {
    ## sperm, egg, 
    dSperm  <- - beta* Sperm * (Eu + Ev + Ep + Em) - r * Sperm
    dEu  <- - gamma * beta * Sperm * Eu
    dEv   <- gamma * beta * Sperm * Eu - delta * Ev - gamma * beta * Sperm * Ev
    dEp   <- gamma * beta * Sperm * Ev
    dEm   <- delta * Ev
    list(c(dSperm, dEu,dEv,dEp,dEm))
  })
}

loop_info <- data.frame(list(beta  = rep(c(0.00259*0.45,0.00259),each= 21),  
                           gamma = 0.08,    # /day, growth rate of prey
                           delta  = 1.84,   # /day, mortality rate of predator
                           r = 0.00027,
                           sperm=rep(10^seq(from= 1,to= 3,by= 0.1),2),
                           u=rep(c(1,1/64),each= 21),v=0,p=0,m=0))

loop_info$start <- do.call("rbind",lapply(exp(log(80)-0.25*log(loop_info$sperm)),function(x) max(x,1)))

row.names(loop_info) <-sprintf("%03d", 1:nrow(loop_info))

apply_fun <- function(list1){
  obj_fun <- function(x){
    maxtime <- x
    inits <- as.numeric(list1[c("sperm","u","v","p","m")])
    parms <- list1[c("beta","gamma","delta","r")]
    a <- ode(y= inits,parms= parms,func= dyn_fert,
             times= c(0.1,7200))
    b <- ode(y= inits,parms= parms,func= dyn_fert,
                      times= c(0.1,maxtime))
    abs(as.numeric((b[nrow(b),6]+b[nrow(b),5])/inits[2])-
          0.95*as.numeric((a[nrow(a),6]+a[nrow(a),5])/inits[2]))
  }
  start <- as.numeric(list1$start)
  optim(par= start,obj_fun,lower=0.101,upper=7199,method="L-BFGS-B")$par
}
apply_fun(loop_info[10,])
data_list <- split(loop_info,rownames(loop_info))
out_unfert <- lapply(data_list,apply_fun)
loop_info$t_sat <-as.numeric(out_unfert)
loop_info$Eggs <- factor(loop_info$u,labels= c("1/64","1"))
loop_info$Eggs <- factor(loop_info$Eggs,levels= c("1","1/64"))


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
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.title.align = 0,
  axis.text.y = element_text(hjust=0.5),
  legend.direction= "vertical",
  legend.position= c(0.8,0.8),
  plot.title= element_text(size= 12),
  panel.border= element_rect(fill= NA),
  legend.key.width = unit(1,"lines"))

ggplot(aes(sperm,t_sat/60),data= loop_info)+
  geom_line(aes(colour= Eggs))+
  theme_bw()+
  ylab("estimated minutes to reach 95% of max fertilization")+
  coord_cartesian(xlim= c(0,1000),ylim= c(0,1000/60))+
  scale_colour_manual(values= c(brewer.pal(5,"Blues")[5],"red"),name= expression(paste("Eggs ",mu,l^-1)))+
  options

+
  facet_grid(variable~init_sperm)


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


gg2 <- arrangeGrob(a3,d3,b3,c3,ncol= 1,nrow= 5,heights= c(1,1,1,1,0.25))