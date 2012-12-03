updateStruct<-function(
### Function incorporate data from SPkapp snapshot into the 3D model of the domain   
  struct=spine,##<< base model of the domain
  amend,##<< simulation results to be incorporated into the model
  show=TRUE,##<< logical value defining should be the new structure sent to RGL package for visualilation
  alpha=0.03##<< transparancy value for the replaced voxels
  ){
	dat<-struct
	dat$Alph<-alpha
	for(i in 1:dim(amend)[1]){
	 ind<-which(amend$X[i]==dat$X&amend$Y[i]==dat$Y&amend$Z[i]==dat$Z)
	 if(length(ind)>0){
	  dat[ind,c('Col','Alph')]<-amend[i,c('Col','Alph')]
	 }else{
	  dat<-rbind(dat,amend[i,])
	 }
	}
  if(show){
	showRGL(dat);
  }
  return(dat);
### data.frame with columns X,Y,Z,Col,Alph that contains coordinates (X,Y,Z) color (Col) and transparancy value (Alph) for voxels
}
showStruct<-function(
### Calculates and shows 3D structure of domain defined by logical function
struct=spine,##<< basic structure of the domain that have to be modified.
color='black',##<< color of spheres, matched by the function
alpha=1,##<< transparency of spheres, matched by the function
def=function(.x,...){.x$X>20&.x$Col=='green'},##<<logical function that defines coordinates of the voxels to be modified
show=FALSE,##<< logical flag, if show=TRUE new domain will be visualised by the 'rgl' library
...##<< additional parameters required by the logical function
){
  res<-struct;
  ind<-which(def(res,...))
  res[ind,c('Col','Alph')]<-list(Col=color,Alp=alpha)
  if(show){
	showRGL(res);
  }
  return(res);
  ### data.frame with columns X,Y,Z,Col,Alph that contains coordinates (X,Y,Z) color (Col) and transparancy value (Alph) for voxels
}
showRGL<-function(
###Visualisation function. Wrapper to RGL call.  
  res,##<<structure to be visualized
  col=res$Col,##<<vector of colors for voxels
  alpha=res$Alph,##<<vector of  transparancy values for voxels
  windowRect=c(0 ,0,900,700)##<<position and size of RGL window 
  ){
    if(!require(rgl)){
	  stop('Function is required package "rgl"');
    }
    open3d(windowRect=windowRect);
    spheres3d(res$X,res$Y,res$Z,radius=0.75,col=res$Col,alpha=res$Alph)
}

saveMovie<-function(
### wrapper for rgl function to create movie of spinning domain. 
  dir='snap'##<< directory to store snapshots
  ){
	system(paste('mkdir -p ',dir))
	for(i in seq(0, 3 * 360, by = 1)) {
	  rgl.viewpoint(theta = i, phi = 0)
	  rgl.snapshot(filename=paste(dir,'/snap.',i,'.png',sep=''),fmt='png')
	}
    }
    
prepareLegend<-function(
### Function create legend image with colored sheres and labels corresponding to the elements colorDef data.frame   
  colorDef,##<<Dictionary defining color of the voxel in dependence of its content.
  windowRect=c(0 ,0,900,700)##<<position and size of RGL window
  ){
    if(!require(rgl)){
	  stop('Function is required package "rgl"');
    }
    c<-colorDef[order(sapply(as.character(colorDef$kappa),nchar)),]
    L<-dim(c)[1]
    x<-rep(0,L)
    y<-rep(0,L)
    z<-10*seq(0,L-1,by=1)
    xt<-rep(10,L)
    yt<-y
    zt<-z
    open3d(windowRect=windowRect);
	spheres3d(x,y,z,col=c$color,alpha=1,radius=5)
	 text3d(xt,yt,zt,adj = c(0, 0.5),c$kappa)
	 spheres3d(350,0,5,radius=0.5)
}
readState<-function(
### main function of the package. Reads content of SPKappa snapshot file and defines color of identified voxels.
  file,##<<snapshot file to read
  anchor=c('GluR1_m:membrane','PSD95:cytosol'),##<<list of regular expressions to be identified in
  ### within content of the voxel. 
  colorDef=NA##<<predefined color set. If NA, new color definition will be generated from content of
  ### snapshot just read.
  ){
##note<< Complex definitions in anchor list processed sequentially in exclusive manner, so second complex define color
### of the voxel only in the case when voxel does not contains first complex and so on.
    if(!require(rkappa)){
	  stop('Function is required package "rkappa"');
    }
    if(!require(RColorBrewer)){
      stop('Function is required package "RColorBrewer"');
    }
    if(!require(rkappa)){
      stop('Function is required package "rkappa"');
    }
    ankoreg<-paste('^.*',anchor,'\\[([0-9]+)\\]\\[([0-9]+)\\]\\[([0-9]+)\\].*$',sep='')
    l<-readLines(file)
    l1<-gsub('%init: +[0-9]+ +','',gsub('\\[[0-9]+\\]','',gsub(':domainLink','',l)))
     bl1<-sapply(l1,makeBruttoStr)
     ubl1<-unique(bl1)
     colors<-colorRampPalette(brewer.pal(8,"Dark2"))(length(ubl1))
     if(any(is.na(colorDef))){
       colorDef<-data.frame(kappa=ubl1,color=colors,stringsAsFactors=FALSE)
     }
     l2<-l
     str<-data.frame(X=1,Y=1,Z=1,Col='green',Alph=1)[FALSE,]
     for(a in ankoreg){
       m<-regexpr(a,l2)
       ind<-which(m>=0)
       l3<-l2[ind]
       s<-substr(l3,m,attr(m,"match.length")+m)
       bl3<-sapply(gsub('%init: +[0-9]+ +','',gsub('\\[[0-9]+\\]','',gsub(':domainLink','',l3))),makeBruttoStr)
       cl3<-colorDef$color[sapply(bl3,pmatch,colorDef$kappa)]
       if(any(is.na(cl3))){
       	warning(bl3[is.na(cl3)])
       }
       coords<-t(data.frame(lapply(strsplit(sub(a,'\\1,\\2,\\3',l3),','),as.integer)))
       rownames(coords)<-NULL
       if(grepl('cytosol',a)){coords<-coords+1}
       str<-rbind(str,data.frame(X=coords[,1],Y=coords[,2],Z=coords[,3],Col=cl3,Alph=1,stringsAsFactors=FALSE))
       l2<-l2[m<0]
     }
     rownames(coords)<-NULL
     attr(str,'palette')<-colorDef
     return(str)
    ### Returns data.frame sutable for updateStruct amend argument. Attribute 'palette' of the
    ### result contains colorDef structure.
}
# makeBrutto<-function(kappa){
#  table(gsub('\\(.*$','',unlist(strsplit(kappa,'),',fixed=TRUE))))->brutto
#  return(brutto)
# }
# makeBruttoStr<-function(kappa){
#  brutto<-makeBrutto(kappa)
#  bstring<-paste(c(rbind(names(brutto),paste(brutto))),collapse='.')
#  return(bstring)
# }