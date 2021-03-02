#######################################################################################################################################
# NICEPLOT : function to build a graphical window with range, labels and title of X and Y axes as specified
# default values are for a single plot saved at resolution=150dpi, width=600pixels, height=600pixels
# (i.e. also works with res=300, width=1200, heigth=1200)
#
# "limX" and "limY": range of X and Y axes (min, max), see function "setrge" below to set limits given data
#
# "marg" : external margins width (in lines); bottom, left, top, right
#
# "tick" : length of tick marks, negative values for outside tick marks, positive for inner ones
#
# "nlab" : maximum number of ticks on each axis
# "labX" and "labY" : position of the labels on the axis
#						if "lab" argument remains empty, labels are set by default
#           			if "lab" argument is NA, no ticks and labels are added
# "nmlabX" and "nmlabY" : label names, default is nmlab=lab
# "lasX" and "lasY" : orientation of labels (1= horizontal, 3=vertical) 
# "lineX" and "lineY" : distance (in lines) between labels and axis
# "cexX" and "cexY" : size of labels characters
#
# "nmX" and "nmY": title of X and Y axes 
# lineXt" and "lineYt" : distance between an axis and its title
# "cexXt" and "cexYt" : size of axis title characters
#
#######################################################################################################################################

niceplot<-function(limX=c(0,1),limY=c(0,1), marg=c(4.5,4.5,2,2),
                      tick=-0.4,   nlab=7, labX=c(),labY=c(),    nmlabX=c(),nmlabY=c(),
                      lasX=1,lasY=1,    lineX=-0.1, lineY=0.1,   cexX=0.9, cexY=0.9,
                      nmX="X",nmY="Y",   lineXt=lineX+2, lineYt=lineY+2.5,   cexXt=1, cexYt=1   ) {
                      
par(mar=marg)   # margins
plot(limX,limY,type="n",axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=limX,ylim=limY) # window
rect(limX[1],limY[1],limX[2],limY[2])   # border
mtext(side=1,nmX,cex=cexXt,line=lineXt,font=2) # X title  
mtext(side=2,nmY,cex=cexYt,line=lineYt,font=2) # Y title 

# labels for X axis
if (is.na(sum(labX))==F) {
labx<-pretty(limX,n=nlab) ; labx<-labx[which(labx>=min(limX) & labx<=max(limX))] ; nmlabx<-labx  # default
if (length(labX)>0) { labx<-labX ; nmlabx<-nmlabX }                                                   # customized
axis(side=1, at=labx, labels=F, tcl=tick, pos=limY[1])  # ticks
mtext(side=1, nmlabx, at=labx, line=lineX, cex=cexX, las=lasX) # labels
                        } # end of if labels
# labels for Y axis
if (is.na(sum(labY))==F) {
laby<-pretty(limY,n=nlab) ; laby<-laby[which(laby>=min(limY) & laby<=max(limY))] ; nmlaby<-laby  # default
if (length(labY)>0) { laby<-labY ; nmlaby<-nmlabY }                                                   # customized
axis(side=2, at=laby, labels=F, tcl=tick, pos=limX[1]) # ticks 
mtext(side=2, nmlaby, at=laby, line=lineY, cex=cexY, las=lasY) # labels
 
                          } # end of if labels
 
} # end of function

##########################################################################################################################
##########################################################################################################################
# SETRGE : function to compute limits of an axis given data
# inputs: "x" a continuous variable (i.e. at least two values) and "p" a coefficent of extension (in %)
# output: limits of the axis = observed range extended by p (default 5%)
 
setrge<-function(x,p=5) {

# additional range 
plus<-(p/100)*(max(x,na.rm=T)-min(x,na.rm=T))

# lower limit
low<-min(x,na.rm=T)- plus

# upper limit
up<-max(x,na.rm=T)+ plus

return(c(low,up))
} # end of function

###########################################################################################################
###########################################################################################################
# SE : function to compute standard error of the mean for a continuous variable X
se<-function(x) { sd(x,na.rm=T)/sqrt(length(na.omit(x))) }

################################################################################################################################
################################################################################################################################
# MEANSEY : function to add points and bars representing mean and associated error (e.g. standard deviation or standard error)
#           for several modalities (X axis) and one variable (Y axis)
# inputs : 
#   x: positions on X axis at which points representing mean have to be plotted
#   meany and sey : vector of same length than X with mean and associated error values to be plotted 
#   pchp, colp, bgp, cexp: shape, border color, background color and size of points representing mean (single value or vector)
#   colb : color for error bars (single value or vector)
#   lgt : width of error bars ticks, in percentage of X range
meansey<-function(x,meany,sey,pchp=22,colp="black",bgp="black",cexp=1.5,colb="black",lgt=0.01) {
lgx<-lgt*(max(x,na.rm=T)-min(x,na.rm=T)) # tick width
segments(x,meany-sey,x,meany+sey,col=colb) # error bar
segments(x-lgx,meany-sey,x+lgx,meany-sey,col=colb) # top tick
segments(x-lgx,meany+sey,x+lgx,meany+sey,col=colb) # bottom tick
points(x,meany,pch=pchp,col=colp,bg=bgp,cex=cexp) # mean
} # end of meansey

################################################################################################################################
################################################################################################################################
# POINTBAR : generic function to add points and vertical bars for several modalities (x axis) and one variable (y axis)
#				e.g. median and 1st and 3rd quartile or mean and 95% confidence-interval      
# inputs : 
#   x: positions on X axis at which points have to be plotted
#   pointy : vector of same length than X with point position on Y axis (e.g. median)
#   boty : vector of same length than X with bottom limit of bar on Y axis, e.g. first quartile values 
#   topy : vector of same length than X with top limit of bar on Y axis, e.g. third quartlie values
#   pchp, colp, bgp, cexp: shape, border color, background color and size of points (single value or vector)
#   colb : color for error bars (single value or vector)
#   lgt : width of error bars ticks, in percentage of X range
pointbar<-function(x,pointy,boty,topy,pchp=22,colp="black",bgp="black",cexp=1.5,colb="black",lgt=0.01) {
lgx<-lgt*(max(x,na.rm=T)-min(x,na.rm=T)) # tick width
segments(x,boty,x,topy,col=colb) # error bar
segments(x-lgx,boty,x+lgx,boty,col=colb) # top tick
segments(x-lgx,topy,x+lgx,topy,col=colb) # bottom tick
points(x,pointy,pch=pchp,col=colp,bg=bgp,cex=cexp) # mean
} # end of pointbar

################################################################################################################################
################################################################################################################################
# MEANSEXY : function to add points and bars representing mean and associated error for two continuous variables for several objects
# inputs : 
#   meanxy : a matrix (k,2) with mean values for the k objects for the two variables (X=column 1, Y=column 2)
#   sexy : a 2-columns matrix (same rows than meanxy) with associated error values for the two variables (X=column 1, Y=column 2)
#   pchp, colp, bgp, cexp: shape, border color, background color and size of points representing means (single value or vector)
#   colb : color for error bars (single value or vector)
#   lgt : width of error bars ticks, in percentage of mean+-se ranges on the two axes

meansexy<-function(meanxy,sexy,pchp=22,colp="black",bgp="black",cexp=1.5,colb="black",lgt=0.01, lwd_seg=1) {
lgx<-lgt*(max(meanxy[,1]+sexy[,1],na.rm=T)-min(meanxy[,1]-sexy[,1],na.rm=T)) # x tick width
lgy<-lgt*(max(meanxy[,2]+sexy[,2],na.rm=T)-min(meanxy[,2]-sexy[,2],na.rm=T)) # y tick width

segments(meanxy[,1]-sexy[,1],meanxy[,2],meanxy[,1]+sexy[,1],meanxy[,2],col=colb, lwd=lwd_seg) # x error bar
segments(meanxy[,1],meanxy[,2]-sexy[,2],meanxy[,1],meanxy[,2]+sexy[,2],col=colb, lwd=lwd_seg) # y error bar
segments(meanxy[,1]-sexy[,1],meanxy[,2]-lgy,meanxy[,1]-sexy[,1],meanxy[,2]+lgy,col=colb, lwd=lwd_seg)
segments(meanxy[,1]+sexy[,1],meanxy[,2]-lgy,meanxy[,1]+sexy[,1],meanxy[,2]+lgy,col=colb, lwd=lwd_seg)
segments(meanxy[,1]-lgx,meanxy[,2]-sexy[,2],meanxy[,1]+lgx,meanxy[,2]-sexy[,2],col=colb, lwd=lwd_seg)
segments(meanxy[,1]-lgx,meanxy[,2]+sexy[,2],meanxy[,1]+lgx,meanxy[,2]+sexy[,2],col=colb, lwd=lwd_seg)
points(meanxy[,1],meanxy[,2],pch=pchp,col=colp,bg=bgp,cex=cexp) # mean
} # end of meansexy

###########################################################################################################
###########################################################################################################
# ADDSEGMENT : function to add a segment defined by an equation and the range on X and/or Y axes
# inputs : 
#  a, b : equation defining the segment (y=a+b*x)
# limX, limY : limit of the segment on both axes; limY is computed by default accoring to a,b and limX
# col_seg, lty_seg=1, lwd_seg : color, line type and width of the segment
addsegment<-function(a=0,b=1,limX,limY=sort(c(a+limX[1]*b,a+limX[2]*b)),col_seg="black",lty_seg=1,lwd_seg=1.5) {

x0<-limX[1] ; y0<-a+limX[1]*b
x1<-limX[2] ; y1<-a+limX[2]*b

if (y0<limY[1] ) {x0<-(limY[1]-a)/b  ; y0<-limY[1] }
if (y0>limY[2]) {x0<-(limY[2]-a)/b  ; y0<-limY[2] }

if (y1<limY[1] ) {x1<-(limY[1]-a)/b  ; y1<-limY[1] }
if (y1>limY[2]) {x1<-(limY[2]-a)/b  ; y1<-limY[2] }

segments(x0,y0,x1,y1,lty=lty_seg,lwd=lwd_seg,col=col_seg)

}  # end of function addsegment

###########################################################################################################
###########################################################################################################
# ROUND_TEXT: function to round a number (x) and get a character string with defined number of decimals (digits)

round_text<-function(x,digits=1) {

round_x<-round(x,digits)
length_round_x<-nchar(as.character(abs(round_x)))
length_int_round_x<-nchar(as.character(trunc(abs(round_x))))

# round x is an integer
if( length_round_x==length_int_round_x) {
	rep0<-""
	for (k in 1:digits)
		rep0<-paste(rep0,"0",sep="")
	res<-paste( trunc(round_x),".", rep0,sep="") } # end of if no decimals

# else
if( length_round_x>length_int_round_x) {
	nb_dec<-length_round_x - length_int_round_x -1
	
		if (nb_dec==digits) res<-as.character(round_x) # ok
	
		if (nb_dec<digits) {
								rep0<-""
								for (k in 1:(digits-nb_dec) )
									rep0<-paste(rep0,"0",sep="")
									res<-paste( round_x, rep0 ,sep="") } # end of if 0 missing
	}# end of if decimals
  
return(res)

} # end of function round_text

###########################################################################################################
###########################################################################################################
# PVALUE_STAR: function to convert a pvalue (p) in classical text notation: ns, *, **, ***
# by default if p>0.05, "ns" is returned but it can be changed with argument "ns" (e.g. to " ")
pvalue_star<-function(p,ns="ns") {
res<-ns
if(p<0.05) res<-"*"
if(p<0.01) res<-"**"
if(p<0.001) res<-"***"

return(res)	
	
} # end of function pvalue_star

###########################################################################################################
###########################################################################################################
# RESCOR : function to extract the results of a correlation test
# inputs : 
# - x, y : 2 continuous variables
# - meth : test to be performed ("spearman", "pearson")
# output :  a character vector with sample size ("n"), correlation coefficient ("r") and asociated pvalue (bilateral test) coded as ns,*,**,*** ("p")

rescor<-function(x,y,meth="spearman") {
cort<-cor.test(x,y, method=meth, exact=F )
n<-length(which(is.na(x+y)==F))
r<-round_text(cort$estimate,3)
p<-pvalue_star(cort$p.value)
res<-c(n,r,p) ; names(res)<-c("n","r","pvalue")
return(res)
} # end of rescor

###########################################################################################################
###########################################################################################################
# RESMANTEL : function to extract results o a Mantel correlation test
# inputs : 
# - x, y : 2 vectors or distance matrices (class dist) with pairwise distances
# - z1, z2: optional vector or dist object to consider in partial Mantel test
# output :  a character vector with sample size ("n"), correlation coefficient ("r") and asociated pvalue (bilateral test) coded as ns,*,**,*** ("p")
# NOTE: be sure that library "vegan" is detached to prevent from confusion with its function "mantel" ( use : detach("package:vegan") )

resmantel<-function(x,y,z1=c(),z2=c()) {
library(ecodist)

cor_mant<-mantel(y~x, nperm=10000 )
if(length(z1)!=0) cor_mant<-mantel(y~x+z1, nperm=10000 )
if(length(z1)!=0 & length(z2)!=0) cor_mant<-mantel(y~x+z1+z2, nperm=10000 )

n<-length(x)
r<-round_text(cor_mant["mantelr"],3)
p<-pvalue_star(cor_mant["pval3"])
res<-c(n,r,p) ; names(res)<-c("n","r","pvalue")
return(res)
} # end of resmantel

###########################################################################################################
###########################################################################################################
# EQUA_LM : function to write the equation of a linear regression 
# inputs : 
# - reslm: output from the linear regression done using the "lm" function
# - nmy= name of Y variable (e.g. Y axis title)
# - nmx= name of X varaible (e.g. X axis title)
# - prec= number of decimals of coefficients a and b in the character string
# output :  an expression with equation of the regression "Y=a+b*X, R²=r²"
#             which can be added to a plot using "title" or "text" function

equa_lm<-function(reslm,nmy="Y",nmx="X",prec=1) {

a<-round_text(reslm$coefficients[1],prec)
b<-round_text(abs(reslm$coefficients[2]),prec)
signb<-"+" ; if(reslm$coefficients[2]<0) signb<-"-"
r2<-round_text(summary(reslm)$r.squared,2)
pval<-pvalue_star(summary(reslm)$coefficients[2,4],ns="")

res<-substitute(y*"="*a*signb * b *" "* scriptstyle("X")*" "* x*", R²=" * r2 * pval,
		list(y=nmy,x=nmx,a=a,b=b,signb=signb,r2=r2,pval=pval) )
return(res)
} # end of equa_lm

###########################################################################################################
###########################################################################################################
# BUBBLE : function to draw circles at given positions with size (i.e. area) proportional to a third variable
# inputs : 
# - x, y : coordinates of the center of the circles
# - k : variable determining circle size
# - maxk, scle : the maximal value of the variable and a scale factor (used to scale size relative to plot) 
#             i.e. if k=max(k) then the circle has a radi of "scle"
#       by default the largest circle has a diameter of 20% of X axis range
# - col, fill : color for border and background of the circles

bubble<-function(x,y,k,maxk=max(k,na.rm=T),scle=0.1*(max(x,na.rm=T)-min(x,na.rm=T)),col="grey",fill="grey")  {
symbols(x,y,circles=scle*sqrt(k/maxk),inches=F,fg=col,bg=fill,add=T)
} # end of function bubble

###########################################################################################################
###########################################################################################################
# PROP_BIQUALI : function to plot distribution of samples among modalities of two categorical variables
# inputs : 
# - conting: contingency table with number of samples in each combination of modalities of the two variables
# - nmmodX and nmmodY: names of respective modalities for the two variables ;  row.names and colnames of "conting" will be used  by default
# - titX and titY: names of the two variables
# proportion is scaled so that 100% filled a unit square

prop_biquali<-function(conting, nmmodX=row.names(conting), nmmodY=colnames(conting), titX="X", titY="Y") {

niceplot(nmX=titX,nmY=titY,limX=c(0.5,length(nmmodX)+0.5),limY=c(0.5,length(nmmodY)+0.5),labX=1:length(nmmodX),labY=1:length(nmmodY),nmlabX=nmmodX,nmlabY=nmmodY) 

segments(0.5,(2:length(nmmodY))-0.5,length(nmmodX)+0.5,(2:length(nmmodY))-0.5,col="grey70")
segments((2:length(nmmodX))-0.5,0.5,(2:length(nmmodX))-0.5,length(nmmodY)+0.5,col="grey70")

conting_sc<-sqrt(conting/sum(conting)) # scaling
for (i in 1:length(nmmodX) )
for (j in 1:length(nmmodY) )
{
rect(i-conting_sc[i,j]/2,j-conting_sc[i,j]/2,i+conting_sc[i,j]/2,j+conting_sc[i,j]/2,col="grey80")
text(i,j,round(conting_sc[i,j]^2*100),font=2)   } # end of i,j
} # end of function biquali
 
###########################################################################################################
###########################################################################################################

dimfig <- function() {
p <- par("usr")
f <- par("plt")
dx <- (p[2] - p[1])/(f[2] - f[1])
dy <- (p[4] - p[3])/(f[4] - f[3])
rect(xleft = p[1] - f[1] * dx, ybottom = p[3] - f[3] * dy, xright = p[2] +
(1 - f[2]) * dx, ytop = p[4] + (1 - f[4]) * dy, col = rgb(0,
0, 1, 0.8), xpd = NA)
arrows(x0 = p[1] - f[1] * dx, y0 = p[3] + (p[4] - p[3])/2, x1 = p[2] +
(1 - f[2]) * dx, y1 = p[3] + (p[4] - p[3])/2, xpd = NA,
code = 3)
text(x = p[1] + (p[2] - p[1])/2, y = p[3] + (p[4] - p[3])/2,
labels = paste(round(par("fin")[1], 1), "pouces"), pos = 3,
xpd = NA)
arrows(x0 = p[1] + (p[2] - p[1])/3, y0 = p[3] - f[3] * dy, x1 = p[1] +
(p[2] - p[1])/3, y1 = p[4] + (1 - f[4]) * dy, xpd = NA,
code = 3)
text(x = p[1] + (p[2] - p[1])/3, y = p[3] + (p[4] - p[3])/2.1,
labels = paste(round(par("fin")[2], 1), "pouces"), pos = 2,
srt = 90, xpd = NA)
}


showmeplot <- function() {
p <- par("usr")
midx <- (p[1] + p[2])/2
midy <- (p[3] + p[4])/2
rect(xleft = p[1], ybottom = p[3], xright = p[2], ytop = p[4],
col = rgb(1, 0, 0, 0.8))
arrows(x0 = p[1], y0 = midy, x1 = p[2], y1 = midy, code = 3)
text(x = midx, y = midy, labels = paste(round(par("pin")[1],
2), "pouces"), xpd = NA, pos = 3)
arrows(x0 = midx, y0 = p[3], x1 = midx, y1 = p[4], code = 3)
text(x = midx, y = p[3] + (p[4] - p[3])/2, labels = paste(round(par("pin")[2],
2), "pouces"), xpd = NA, pos = 2, srt = 90)
}

#	-------------------------------------------------------------------------------------
	# Fonction qui calcule le R²

erdeu=function(model,y){
			#model: c'est le modèle dont tu veux le r²
			#y: ce sont les variables à prédire (un vecteur, ça doit marcher)
			#exemple:  erdeu(nls_all,data[,2])
			Res=sum((summary(model)$residuals)^2)
			r=1-(Res/sum((y-mean(y))^2))
	
			a=list(r,Res)
			names(a)=c("r2","RSS")
			return(a)
		}

#	-------------------------------------------------------------------------------------
	# Fonction qui enlève les niveaux des lites

flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}                     
 
 
 
 roundUpNice <- function(x, nice=c(0.1,1,1.5,2,2.5,4,4.5,5,5.5,6,7,8,9,10)) {
   if(length(x) != 1) stop("x must be of length 1")
   10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
                  
 error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
   arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
 }
 
 
 pointbar<-function(x,pointy,boty,topy,pchp=22,colp="black",bgp="black",
                    cexp=1.5,colb="black",lgt=0.01, linwd = 1) {
   lgx<-lgt*(max(x,na.rm=T)-min(x,na.rm=T)) # tick width
   segments(x,boty,x,topy,col=colb, lwd = linwd) # error bar
   segments(x-lgx,boty,x+lgx,boty,col=colb, lwd = linwd) # top tick
   segments(x-lgx,topy,x+lgx,topy,col=colb, lwd = linwd) # bottom tick
   points(x,pointy,pch=pchp,col=colp,bg=bgp,cex=cexp) # mean
 } # end of pointbar
 
 