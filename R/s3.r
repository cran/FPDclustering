

plot.FPDclustering=function(x, ...){
  label=as.character(x$label)
  par(mfrow=c(1,1))
  Silh(x$probability)
  if(ncol(x$data)<10){
    par(mfrow=c(1,1))
 print(   ggparcoord(cbind(cluster= label,x$data),groupColumn = 1, columns=2:(ncol(x$data)+1))+ 
      theme_bw())
   print( pairs(x$data,col=label, ...))
  }else{
    ggparcoord(cbind(cluster=(label),x$data),groupColumn = 1, columns=2:11)+ 
      theme_bw()+ggtitle('First 10 varaibles')
  }
  if( exists('cont.loc',where=x)){
    if(length(x$cont.loc)>1){
    cont.loc=x$cont.loc
    if(ncol(x$data[,cont.loc])<10& ncol(x$data[,cont.loc])>1){
      print(     ggparcoord(cbind(cluster= label,x$data[,cont.loc]),groupColumn = 1, columns=2:(ncol(x$data[,cont.loc])+1))+
        theme_bw())

      print( pairs(x$data[,cont.loc],col=label, ...))

    }else if(ncol(x$data[,cont.loc])>=10){

      print(     ggparcoord(cbind(cluster=(label),x$data[,cont.loc]),groupColumn = 1, columns=2:11)+
        theme_bw()+ggtitle('First 10 varaibles'))
    }}
  
  
}}



summary.FPDclustering=function(object, ...){
  t=as.data.frame( table(object$label))
  names(t)[1]="Cluster"
  names(t)[2]="N. of elements"
  print(t)
  
}

