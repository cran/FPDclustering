


plot.FPDclustering=function(x,maxVar=30, ...){
  label=as.factor(x$label)
  par(mfrow=c(1,1))
  prob=apply(x$probability,1,max)
  Silh(x$probability)
  if( exists('cont.loc',where=x)){
    if(length(x$cont.loc)>1){
    cont.loc=x$cont.loc
    if(ncol(x$data[,cont.loc])<10& ncol(x$data[,cont.loc])>1){
      print(     ggparcoord(cbind(label= label,x$data[,cont.loc]),mapping=ggplot2::aes(color=as.factor(label)), columns=2:(ncol(x$data[,cont.loc])+1))+
        theme_bw()+scale_color_discrete("Clusters",labels=levels(label)))

      print( ggpairs(data.frame(x$data[,cont.loc]),mapping=aes(colour =label)))

    }else {if(ncol(x$data[,cont.loc])>=10){

      print(     ggparcoord(cbind(label=(label),x$data[,cont.loc]),mapping=ggplot2::aes(color=as.factor(label)), columns=2:11)+
        theme_bw()+scale_color_discrete("Clusters",labels=levels(label))+ggtitle('First 10 varaibles'))
    }}}} else{  if(ncol(x$data)<10){
      par(mfrow=c(1,1))
      print(   ggparcoord(cbind(label= label,x$data,prob=prob),mapping=ggplot2::aes(color=as.factor(label)),columns=2:(ncol(x$data)+1))+ 
                 theme_bw()+scale_color_discrete("Clusters",labels=levels(label)))
      print( ggpairs(data.frame(x$data),mapping=aes(colour =label)))
    }else{
      end=1
      repeat {
        begin=end+1
        end=min(begin+9,ncol(x$data))
        print( ggparcoord(cbind(label=(label),x$data),mapping=ggplot2::aes(color=as.factor(label)), columns=begin:end)+ 
                 theme_bw()+scale_color_discrete("Clusters",labels=levels(label)))
        
        if(end==ncol(x$data)|end>maxVar){break}}
    }}
  
  
}



summary.FPDclustering=function(object, ...){
  t=as.data.frame( table(object$label))
  names(t)[1]="Cluster"
  names(t)[2]="N. of elements"
  print(t)
  
}

