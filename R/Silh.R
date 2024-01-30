Silh <-function(p){
# %%%%Silhouette for probabilistic clustering methods
# %%%%%INPUT
# %p probability matrix
# %%%%%OUTPUT
# % Silhouettes
  n=nrow(p)
nc=ncol(p)

la=max.col(p)
ds=disS(p)

m=cbind(la,ds)
m2=m[sort.list(m[,2]), ]
m=m2[sort.list(m2[,1]), ]

ss=m[,2]

ll=table(la)
pl=cumsum(ll)
#pl=[0,pl(1:nc-1)];
ll=round(ll/2)
  pl=c(0,pl[1:(nc-1)])
tic=ll+pl
#tcks=ll+pl;
 # print(barplot(ss,space=0, main="Silhouette plot", horiz=TRUE,xlab='Silhouette Value', ylab='Cluster',xpd=F,axes = FALSE)
  ssdf=data.frame(ss=ss,x=1:length(ss))
 print( ggplot(ssdf, aes(x = 1:length(ss), y = ss)) + geom_col()+ 
    theme_bw()+ coord_flip()+ylab("Silhouette Value") +
    xlab("Cluster") + ggtitle('Silhouette plot')+  ggeasy::easy_center_title()+
    scale_x_continuous(
      breaks = as.numeric(tic),
      labels = (1:nc)
    ))
 
  
return(mean(ss))
}
