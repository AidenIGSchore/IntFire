#Corresponding Author: Aiden I. G. Schore, contact aschore2(at)illinois.edu
##############################################################################
#Citation: Schore, A. I. G. et al., Plant functional types improve satellite-
#derived fire severity assessments in interior Alaska. Environmental Research:
#Ecology. In prep.
##############################################################################
#This code takes ground-based fire severity measurements and compares them to
#remotely sensed images. After determining the strongest relationship for each
#plant functional type, it generates equations to use in Google Earth Engine
#code linked at the appropriate spot. The code then continues using the outputs
#from Google Earth Engine to upscale the PlanetScope products and then runs a
#validation of the performance of each.
##############################################################################

require(outliers)

# Initial/PlanetScope Run ------------------------------------------------------
#This section uses the PlanetScope data provided at https://arcticdata.io/catalog/view/doi:10.18739/A25H7BW3S to determine what burn index
#best corresponds to each plant functional type (PFT).

#If you want to run this code, make sure the spreadsheets are located within your working directory.

#Read in the ground data points and fire index values at those points
ConDat=read.csv("GroundPoints_Planet.csv") #When using your own data, replace the filename here.
#Other code may need to be updated if the input data is not in the same format as the provided data.

#Remove ground point located in a topographic shadow
which(ConDat$NAME=="NC023") #51
ConDat=ConDat[c(1:50,52:135),]

#Eliminate the points with CBI values of 0 and separate them into PFTs. Then remove outliers.
ConDat2=ConDat[which(ConDat$CBI>0),]
Conifer=ConDat2[which(ConDat2$ValVegClass=="Coniferous"|ConDat2$ShrubBreak=="Coniferous"),]
#The ShrubBreak column splits the removed shrub classification into Woody coniferous and woody deciduous classes,
#which is idiosyncratic to this data. For other data you'll likely want to remove the | term
plot(Conifer$CBI) #There may be an outlier
grubbs.test(Conifer$CBI) #Lowest value is an outlier
which.min(Conifer$CBI) #29
Conifer2=Conifer[c(1:28,30:51),]
plot(Conifer2$CBI) #There may still be an outlier
grubbs.test(Conifer2$CBI) #Lowest value is an outlier
which.min(Conifer2$CBI) #11
Conifer2=Conifer2[c(1:10,12:50),]
plot(Conifer2$CBI) #Test for more outliers
grubbs.test(Conifer2$CBI) #No more outliers
Deciduous=ConDat2[which(ConDat2$ValVegClass=="Deciduous"|ConDat2$ShrubBreak=="Deciduous"),]
plot(Deciduous$CBI) #There may be an outlier
grubbs.test(Deciduous$CBI) #Lowest value is an outlier
which.min(Deciduous$CBI) #23
Deciduous2=Deciduous[c(1:22,24:27),]
plot(Deciduous2$CBI) #More outliers look unlikely
grubbs.test(Deciduous2$CBI) #No more outliers
Graminoid=ConDat2[which(ConDat2$ValVegClass=="Graminoid"),]
plot(Graminoid$CBI) #Outlier looks unlikely
grubbs.test(Graminoid$CBI) #no outliers

#Because CBI caps out asymptotically at 3, we'll be running some logistic regressions. Those columns need names.
lognames=rep(NA,18)
for(i in 1:18){
    lognames[i]=paste("log", names(ConDat2)[i+6], sep="_")
  }

#Set up our regression data
compDat=array(data=NA, dim=c(36,3,3), dimnames = list(c(lognames,names(ConDat2)[7:24]), c("Coniferous", "Deciduous", "Graminoid"), c("R2", "p_var", "n")))
for(i in 1:18){
  j=i+6 #j values are based on the setup of the input tables
  
  #Run logarithmic regressions for the Coniferous class
  if(length(which(is.na(Conifer2[,j])==FALSE)>2)){
    k=which(is.na(Conifer2[,j])==FALSE)
    #Log can only be run on positive values, so we'll modify any negative values by subtracting twice the minimum value to shift the range. A minimum value
    #of exactly 0 is unlikely enough that I haven't written code for the scenario
    if(min(Conifer2[k,j])<0){
      x=summary(lm(Conifer2[k,6] ~ log(Conifer2[k,j]-2*min(Conifer2[k,j]))))
    } else{
      x=summary(lm(Conifer2[k,6] ~ log(Conifer2[k,j])))
    }
    compDat[i,1,1]=x$r.squared
    compDat[i,1,2]=x$coefficients[8]
    compDat[i,1,3]=length(which(is.na(Conifer2[,j])==FALSE))
    
    #Now the linear regressions
    y=summary(lm(Conifer2[k,6]~Conifer2[k,j]))
    l=i+18
    compDat[l,1,1]=y$r.squared
    compDat[l,1,2]=y$coefficients[8]
    compDat[l,1,3]=length(which(is.na(Conifer2[,j])==FALSE))
  }
  
  #Repeat for Deciduous class
  if(length(which(is.na(Deciduous2[,j])==FALSE)>2)){
    k=which(is.na(Deciduous2[,j])==FALSE)
    if(min(Deciduous2[k,j])<0){
      x=summary(lm(Deciduous2[k,6] ~ log(Deciduous2[k,j]-2*min(Deciduous2[k,j]))))
    } else{
      x=summary(lm(Deciduous2[k,6] ~ log(Deciduous2[k,j])))
    }
    compDat[i,2,1]=x$r.squared
    compDat[i,2,2]=x$coefficients[8]
    compDat[i,2,3]=length(which(is.na(Deciduous2[,j])==FALSE))
    y=summary(lm(Deciduous2[k,6]~Deciduous2[k,j]))
    l=i+18
    compDat[l,2,1]=y$r.squared
    compDat[l,2,2]=y$coefficients[8]
    compDat[l,2,3]=length(which(is.na(Deciduous2[,j])==FALSE))
  }
  
  #Now Graminoid
  if(length(which(is.na(Graminoid[,j])==FALSE)>2)){
    k=which(is.na(Graminoid[,j])==FALSE)
    if(min(Graminoid[k,j])<0){
      x=summary(lm(Graminoid[k,6] ~ log(Graminoid[k,j]-2*min(Graminoid[k,j]))))
    } else{
      x=summary(lm(Graminoid[k,6] ~ log(Graminoid[k,j])))
    }
    compDat[i,3,1]=x$r.squared
    compDat[i,3,2]=x$coefficients[8]
    compDat[i,3,3]=length(which(is.na(Graminoid[,j])==FALSE))
    y=summary(lm(Graminoid[k,6]~Graminoid[k,j]))
    l=i+18
    compDat[l,3,1]=y$r.squared
    compDat[l,3,2]=y$coefficients[8]
    compDat[l,3,3]=length(which(is.na(Graminoid[,j])==FALSE))
  }
}

View(compDat[,,1])
View(compDat[,,2])
#View(compDat[,,3]) #If you want to see n

#I prefer to manually look at the table to see additional information, but if you want to just get the top performers, uncomment the code below
# for(i in 1:3){
#  j=which.max(compDat[,i,1])
#  print(paste0("Best Performing Index for ", names(compDat[1,,1])[i], ": ", rownames(compDat[,,1])[j]))
# }

#Based on the top performing regressions, we get the necessary coefficients
ConiferEq_PS=lm(CBI~dGNDVI, data=Conifer2)
DeciduousEq_PS=lm(CBI~NDVI, data=Deciduous2) #The top performer had an outlier in the x domain, so we're using the next best fire index (R-square diff. < 0.06)
GraminoidEq_PS=lm(CBI~SR, data=Graminoid)

summary(ConiferEq_PS)
summary(DeciduousEq_PS)
summary(GraminoidEq_PS)

#Now take those coefficients and move to the Google Earth Engine code to generate the first set of unmixing values
#GEE Code: https://code.earthengine.google.com/1a53e974627becd00a53aeee700a20d5
#All necessary assets for the Earth Engine code are public

# Harmonized Landsat/Sentinel Run ----------------------------------------------------------
#Now we test the performance using Harmonized Landsat/Sentinel (HSL) imagery
#The spreadsheet for this section is available at https://arcticdata.io/catalog/view/doi:10.18739/A25H7BW3S.

#Load datasets generated by Google Earth Engine
ConDat3=read.csv("GroundPoints_HLS.csv") #When using your own data, replace the filename here.
#Other code may need to be updated if the input data is not in the same format as the provided data.

#Some modifications to the spreadsheet
dNBR2=(ConDat3$dNBR)^2 #Create a squared dNBR column for our polynomial
ConDat3=cbind(ConDat3, dNBR2)
ConDat3$RdNBR=ConDat3$RdNBR/1000 #The values in the spreadsheets are scaled by a factor of 1000, that's easy enough to fix here
ConDat3$Class30m[which(ConDat3$Class30m=="")]<-NA

#Remove a data point in a topographic shadow
which(ConDat3$NAME=="NC023") #51
ConDat3=ConDat3[c(1:50,52:135),]

#Remove ponts with 0 CBI
ConDat4=ConDat3[which(ConDat3$CBI>0),]

#Separate out PFTs to calculate equations for discrete 30m HSL imagery
Conifer4=ConDat4[which(ConDat4$Class30m=="Coniferous"),]
plot(Conifer4$CBI) #There may be an outlier
grubbs.test(Conifer4$CBI) #Lowest value is an outlier
which.min(Conifer4$CBI) #25
Conifer4=Conifer4[c(1:24,26:48),]
plot(Conifer4$CBI) #More outliers look unlikely
grubbs.test(Conifer4$CBI) #No more outliers
Deciduous4=ConDat4[which(ConDat4$Class30m=="Deciduous"),]
plot(Deciduous4$CBI) #There may be an outlier
grubbs.test(Deciduous4$CBI) #No statistically significant outliers
Graminoid4=ConDat4[which(ConDat4$Class30m=="Graminoid"),]
plot(Graminoid4$CBI) #Outliers look unlikely
grubbs.test(Graminoid4$CBI) #No outliers

#Update the column names for the logarithmic regressions for the variables we're testing on this run
lognames2=rep(NA,20)
for(i in 1:20){
  lognames2[i]=paste("log", names(ConDat4)[i+5], sep="_")
}

compDat2=array(data=NA, dim=c(41,3,3), dimnames = list(c(lognames2, names(ConDat4)[6:25], "dNBR Polynomial"), c("Coniferous", "Deciduous", "Graminoid"), c("R2", "p_var", "n")))
for(i in 1:20){
  j=i+5 #j values are based on the setup of the input tables
  
  #Run regressions for the Coniferous class
  if(length(which(is.na(Conifer4[,j])==FALSE)>2)){
    k=which(is.na(Conifer4[,j])==FALSE)
    #Log can only be run on positive values, so we'll modify any negative values by subtracting twice the minimum value to shift the range. Values of
    #exactly 0 are unlikely enough that I haven't written code for them
    if(min(Conifer4[k,j])<0){
      x=summary(lm(Conifer4[k,5] ~ log(Conifer4[k,j]-2*min(Conifer4[k,j]))))
    } else{
      x=summary(lm(Conifer4[k,5] ~ log(Conifer4[k,j])))
    }
    compDat2[i,1,1]=x$r.squared
    compDat2[i,1,2]=x$coefficients[8]
    compDat2[i,1,3]=length(which(is.na(Conifer4[,j])==FALSE))
    
    #Now the linear regressions
    y=summary(lm(Conifer4[k,5]~Conifer4[k,j]))
    l=i+20
    compDat2[l,1,1]=y$r.squared
    compDat2[l,1,2]=y$coefficients[8]
    compDat2[l,1,3]=length(which(is.na(Conifer4[,j])==FALSE))
  }
  
  #Repeat for Deciduous class
  if(length(which(is.na(Deciduous4[,j])==FALSE)>2)){
    k=which(is.na(Deciduous4[,j])==FALSE)
    if(min(Deciduous4[k,j])<0){
      x=summary(lm(Deciduous4[k,5] ~ log(Deciduous4[k,j]-2*min(Deciduous4[k,j]))))
    } else{
      x=summary(lm(Deciduous4[k,5] ~ log(Deciduous4[k,j])))
    }
    compDat2[i,2,1]=x$r.squared
    compDat2[i,2,2]=x$coefficients[8]
    compDat2[i,2,3]=length(which(is.na(Deciduous4[,j])==FALSE))
    y=summary(lm(Deciduous4[k,5]~Deciduous4[k,j]))
    l=i+20
    compDat2[l,2,1]=y$r.squared
    compDat2[l,2,2]=y$coefficients[8]
    compDat2[l,2,3]=length(which(is.na(Deciduous4[,j])==FALSE))
  }
  
  #Now Graminoid
  if(length(which(is.na(Graminoid4[,j])==FALSE)>2)){
    k=which(is.na(Graminoid4[,j])==FALSE)
    if(min(Graminoid4[k,j])<0){
      x=summary(lm(Graminoid4[k,5] ~ log(Graminoid4[k,j]-2*min(Graminoid4[k,j]))))
    } else{
      x=summary(lm(Graminoid4[k,5] ~ log(Graminoid4[k,j])))
    }
    compDat2[i,3,1]=x$r.squared
    compDat2[i,3,2]=x$coefficients[8]
    compDat2[i,3,3]=length(which(is.na(Graminoid4[,j])==FALSE))
    y=summary(lm(Graminoid4[k,5]~Graminoid4[k,j]))
    l=i+20
    compDat2[l,3,1]=y$r.squared
    compDat2[l,3,2]=y$coefficients[8]
    compDat2[l,3,3]=length(which(is.na(Graminoid4[,j])==FALSE))
  }
}

#Tack on the polynomial terms that didn't fit in the loops
#First the one based on the given dNBR values
x=summary(lm(CBI ~ dNBR+dNBR2, data=Conifer4))
compDat2[41,1,1]=x$r.squared
compDat2[41,1,2]=x$coefficients[12]
compDat2[41,1,3]=length(which(is.na(Conifer4[,24])==FALSE))

x=summary(lm(CBI ~ dNBR+dNBR2, data=Deciduous4))
compDat2[41,2,1]=x$r.squared
compDat2[41,2,2]=x$coefficients[12]
compDat2[41,2,3]=length(which(is.na(Deciduous4[,24])==FALSE))

x=summary(lm(CBI ~ dNBR+dNBR2, data=Graminoid4))
compDat2[41,3,1]=x$r.squared
compDat2[41,3,2]=x$coefficients[12]
compDat2[41,3,3]=length(which(is.na(Graminoid4[,24])==FALSE))

View(compDat2[,,1])
View(compDat2[,,2])
#View(compDat2[,,3]) #If you want to see n

#I prefer to manually look at the table to see additional information, but if you want to just get the top performers, uncomment the code below
# for(i in 1:3){
#  j=which.max(compDat2[,i,1])
#  print(paste0("Best Performing Index for ", names(compDat2[1,,1])[i], ": ", rownames(compDat2[,,1])[j]))
# }

#Based on the top performing regressions, we get the necessary coefficients
ConiferEq_HLS=lm(CBI~GNDVI, data=Conifer4)
DeciduousEq_HLS=lm(CBI~log(dNBR), data=Deciduous4) #The top performer had an outlier in the x domain, so we're using the next best fire index (R-square diff. < 0.06)
GraminoidEq_HLS=lm(CBI~SR, data=Graminoid4)

summary(ConiferEq_HLS)
summary(DeciduousEq_HLS)
summary(GraminoidEq_HLS)

#Now take those coefficients and move to the Google Earth Engine code to generate the second set of unmixing values
#GEE Code: https://code.earthengine.google.com/88c730cdfb38c808b4ef791fd36d8ed7
#All necessary assets for the Earth Engine code are public
#Come back after generating the unmixed IntFire values


# Unmixing Evaluation -----------------------------------------------------
#Come back here after generating the unmixed IntFire in Google Earth Engine.
#For this code, the generation is already complete and can be found at https://arcticdata.io/catalog/view/doi:10.18739/A25H7BW3S.
#Here we evaluate the performance of IntFire calculated at 3m and 30m compared to
#the spectrally unmixed HLS imagery based on each set of relationships.

#Read in the unmixed data
unmixed=read.csv("GroundPoints_Unmixed.csv")

#Reset the previous tables so we can properly bind the unmixed values
ConDat3=read.csv("GroundPoints_HLS.csv") #When using your own data, replace the filename here.
#Other code may need to be updated if the input data is not in the same format as the provided data.
dNBR2=(ConDat3$dNBR)^2 #Create a squared dNBR column for our polynomial
ConDat3=cbind(ConDat3, dNBR2)
ConDat3$RdNBR=ConDat3$RdNBR/1000 #The values in the spreadsheets are scaled by a factor of 1000, that's easy enough to fix here
ConDat3$Class30m[which(ConDat3$Class30m=="")]<-NA
ConDat3=cbind(ConDat3, unmixed[,2:3])

#Remove a data point in a topographic shadow
which(ConDat3$NAME=="NC023") #51
ConDat3=ConDat3[c(1:50,52:135),]

#Remove ponts with 0 CBI
ConDat4=ConDat3[which(ConDat3$CBI>0),]

#Apply the equations we generated earlier to the discrete HSL pixels
IntFire_HLS_Coarse=rep(NA,96)
for(i in 1:nrow(ConDat4)){
  if(is.na(ConDat4$Class30m[i])==FALSE){
  if(ConDat4$Class30m[i]=="Coniferous"){
    IntFire_HLS_Coarse[i]=1.54527*ConDat4$dGNDVI[i]+2.45477
  } else if(ConDat4$Class30m[i]=="Deciduous"){
    IntFire_HLS_Coarse[i]=-1.7445*ConDat4$NDVI[i]+2.8281
  } else if(ConDat4$Class30m[i]=="Graminoid"){
    IntFire_HLS_Coarse[i]=-0.1986*ConDat4$SR[i]+2.58694
  }
  }}

IntFire_HLS_Discrete=rep(NA,96)
for(i in 1:nrow(ConDat4)){
  if(is.na(ConDat4$Class30m[i])==FALSE){
    if(ConDat4$Class30m[i]=="Coniferous"){
      IntFire_HLS_Discrete[i]=-4.4527*ConDat4$GNDVI[i]+4.7336
    } else if(ConDat4$Class30m[i]=="Deciduous"){
      IntFire_HLS_Discrete[i]=0.6418*log(ConDat4$dNBR[i])+2.7639
    } else if(ConDat4$Class30m[i]=="Graminoid"){
      IntFire_HLS_Discrete[i]=-0.22705*ConDat4$SR[i]+2.87169
    }
  }}

ConDat4=cbind(ConDat4, IntFire_HLS_Coarse, IntFire_HLS_Discrete)
#Remove the outlying points we excluded from the initial analysis
which(ConDat4$CBI==0.088) #68 - conifer outlier 1
which(ConDat4$CBI==0.895) #13 - conifer outlier 2
which(ConDat4$CBI==0.0888) #67 - deciduous outlier
ConDat4=ConDat4[c(1:12,14:66,69:96),]

#Generate a table of what indices best correspond to the ground-based CBI data
compDat3=array(data=NA, dim=c(45,3), dimnames = list(c(names(ConDat4)[6:25], names(ConDat4)[28:31], lognames2, "dNBR Polynomial"), c("R2", "p_var", "n")))
#Linear regressions
for(i in 1:24){
  if(i<21){
    j=i+5 #j values for burn severity indices
  } else{
    j=i+7 #j values for IntFire variants
  }
  
  if(length(which(is.na(ConDat4[,j])==FALSE)>2)){
    x=summary(lm(ConDat4[,5] ~ ConDat4[,j]))
    compDat3[i,1]=x$r.squared
    compDat3[i,2]=x$coefficients[8]
    compDat3[i,3]=length(which(is.na(ConDat4[,j])==FALSE))
  }
}

#Logarithmic regressions
for(i in 1:20){
  j=i+5
  if(length(which(is.na(ConDat4[,j])==FALSE)>2)){
    k=which(is.na(ConDat4[,j])==FALSE)
    if(min(ConDat4[k,j])<0){
      x=summary(lm(ConDat4[k,5] ~ log(ConDat4[k,j]-2*min(ConDat4[k,j]))))
    } else{
      x=summary(lm(ConDat4[k,5] ~ log(ConDat4[k,j])))
    }
    l=i+24
    compDat3[l,1]=x$r.squared
    compDat3[l,2]=x$coefficients[8]
    compDat3[l,3]=length(which(is.na(ConDat4[,j])==FALSE))
  }
}

#Tack on the polynomial term and the discrete pixel calculations that didn't fit in the loops
#First the one based on the given dNBR values
x=summary(lm(CBI ~ dNBR+dNBR2, data=ConDat4))
compDat3[45,1]=x$r.squared
compDat3[45,2]=x$coefficients[12]
compDat3[45,3]=length(which(is.na(ConDat4[,24])==FALSE))

View(compDat3) #Discrete IntFire R2=0.469, GNDVI R2=0.414, 13.2% improvement

# PlanetScope Evaluation ------------------------------------------------------------
#Here we evaluate the relationships generated in the first section of code to see the
#overall performance of IntFire.

IntFire_fine=rep(NA,96)
for(i in 1:nrow(ConDat2)){
  if(ConDat2$ValVegClass[i]=="Coniferous"){
    IntFire_fine[i]=1.54527*ConDat2$dGNDVI[i]+2.45477
  } else if(ConDat2$ValVegClass[i]=="Deciduous"){
    IntFire_fine[i]=-1.7445*ConDat2$NDVI[i]+2.8281
  } else if(ConDat2$ValVegClass[i]=="Shrub"){
    if(ConDat2$ShrubBreak[i]=="Coniferous"){
      IntFire_fine[i]=1.54527*ConDat2$dGNDVI[i]+2.45477
    } else if(ConDat2$ShrubBreak[i]=="Deciduous"){
      IntFire_fine[i]=-1.7445*ConDat2$NDVI[i]+2.8281
    }
  } else{
    IntFire_fine[i]=-0.1986*ConDat2$SR[i]+2.58694
  }
}

IntFire_PS_Coarse=rep(NA,96)
for(i in 1:nrow(ConDat2)){
  if(is.na(ConDat2$Class30m[i])==FALSE){
    if(ConDat2$Class30m[i]=="Coniferous"){
      IntFire_PS_Coarse[i]=1.54527*ConDat2$dGNDVI[i]+2.45477
    } else if(ConDat2$Class30m[i]=="Deciduous"){
      IntFire_PS_Coarse[i]=-1.7445*ConDat2$NDVI[i]+2.8281
    } else if(ConDat2$Class30m[i]=="Graminoid"){
      IntFire_PS_Coarse[i]=-0.1986*ConDat2$SR[i]+2.58694
    }
  }}

ConDat2=cbind(ConDat2, IntFire_fine, IntFire_PS_Coarse)
which(ConDat2$CBI==0.088) #68 - conifer outlier
which(ConDat2$CBI==0.895) #13 - conifer outlier 2
which(ConDat2$CBI==0.0888) #67 - deciduous outlier
ConDat5=ConDat2[c(1:12,14:66,69:96),]

#Run regressions of each fire index (including CBIm) and ground-based CBI
compDat4=array(data=NA, dim=c(38,3), dimnames = list(c(names(ConDat5)[7:24], lognames, "IntFire","Coarse_IntFire"), c("R2", "p_var", "n")))
#Linear Regressions
for(i in 1:18){
  j=i+6
  
  if(length(which(is.na(ConDat5[,j])==FALSE)>2)){
    x=summary(lm(ConDat5[,6] ~ ConDat5[,j]))
    compDat4[i,1]=x$r.squared
    compDat4[i,2]=x$coefficients[8]
    compDat4[i,3]=length(which(is.na(ConDat5[,j])==FALSE))
  }
}
#Logarithmic Regressions
for(i in 1:18){
  j=i+6
  if(length(which(is.na(ConDat5[,j])==FALSE)>2)){
    k=which(is.na(ConDat5[,j])==FALSE)
    if(min(ConDat5[k,j])<0){
      x=summary(lm(ConDat5[k,6] ~ log(ConDat5[k,j]-2*min(ConDat5[k,j]))))
    } else{
      x=summary(lm(ConDat5[k,6] ~ log(ConDat5[k,j])))
    }
    l=i+18
    compDat4[l,1]=x$r.squared
    compDat4[l,2]=x$coefficients[8]
    compDat4[l,3]=length(which(is.na(ConDat5[,j])==FALSE))
  }
}

#Add on IntFire and Coarse IntFire
x=summary(lm(CBI ~ IntFire_fine, data=ConDat5))
compDat4[37,1]=x$r.squared
compDat4[37,2]=x$coefficients[8]
compDat4[37,3]=length(which(is.na(ConDat5[,27])==FALSE))

x=summary(lm(CBI ~ IntFire_PS_Coarse, data=ConDat5))
compDat4[38,1]=x$r.squared
compDat4[38,2]=x$coefficients[8]
compDat4[38,3]=length(which(is.na(ConDat5[,28])==FALSE))


View(compDat4) #CBIm R2=0.591, SR R2=0.323, 82.6% improvement
