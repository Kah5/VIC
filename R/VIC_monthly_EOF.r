##This is the working code for working with VIC output data from the model runst
#that Chun-Mei set up over the Prairie Peninsula region

#the Flux files used in this analysis reside in thing2.crc.nd.edu
# under /tmp/kheilman_SM/RUN/RESULTS

R #Open R in the crc
library(sp)

##Pacific Oscillation index data for later use
#this data is from http://www.cpc.ncep.noaa.gov/data/teledoc/pna.shtml
PNA <- read.table("PNA.txt", header=TRUE, sep="\t")
NAO <- read.table("NAO.txt", header=TRUE, sep="\t")
#########################################3
#Loading Monthly data from chun-mei's VIC run
################################
# list the files to process that are in the current directory with the pattern beginng with "FLUX_"
filesToProcess <- dir(pattern = "Monthly_Flux_.*")

file_list<-filesToProcess
##extract column of interest from all flux files into columns using cbind
dataset<-do.call("cbind", lapply(file_list,FUN=function(files){read.table(files)$V3}))
colnames(dataset)<-substr(file_list, 14, 30) #this takes the latlon from each file makes it the header for each row

latlong<-cbind(substr(file_list,14,21),substr(file_list,23,30)) #create latlong vector from file names


######################
##EOF Analysis
#####################
#express the space-time data as a 2D matrix with "M" rows (space) by "N" columns(time)
#transpose our data so that time is represented in the columns and space in rows

A <- t(as.matrix(dataset))  #convert to matrix and transpose

#A <- matrix(c(2,4,-6,8,1,2,-3,4,4,1,3,2), nrow=3, byrow=TRUE) #testing with dummy matrix
##preprocess the data by subtracting the long-term mean (mean of each row) from each row
##This creates an "anomoly" matrix
ltmeans <- apply(A, 1, mean)   #calculate the longterm means for each row and puts it in a list
X <- A - ltmeans   #subtract longterm means from the matrix A

## Calculate the eigenvectors (E) and Eigen Values EV of the covariance matrix C
esv <- svd(X) #svd analysis of a matrix returns a lise with d, u, and v objections
diag.s <- diag(esv$d) #this is equivalent to S in matlab
V <- esv$v #this should be columns of the EOF in time dimension
U <-esv$u

## calculate the z matrix
z <- diag.s %*% t(V) # z is the matrix multiplcation of S and V'
#rows of the z matrix should be the Principle components associated with each EOF

## z can be confirmed another way:
uta <- t(U) %*% X

## calculate the orginal data based on the first 2 EOFS
A <- U[,1:2] %*% z[1:2,] # this gives the first 2 rows in the first column as non-zero, but rounding, they are zero

#the Principle components here now represent time  series of length N
plot(z[1,], type='l')  #time series explaining the most variation in the data
plot(z[2,], type='l')  #time series for 2nd EOF
plot(z[3,], type='l')

#create separate objects for eofs 1-3
eof1 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,1]))
colnames(eof1) <- c("lat", "lon","eof1")
eof2 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,2]))
colnames(eof2) <- c("lat", "lon","eof2")
eof3 <- data.frame(cbind(as.numeric(latlong[,1]), as.numeric(latlong[,2]), U[,3]))
colnames(eof3) <- c("lat", "lon","eof3")

write.table(U,"precipU.txt",sep='\t')
write.table(Z, "precipPCs.txt", sep='\t')

##make the data spatial points using coordinates() from the sp package
coordinates(eof1) <-~lon+lat
coordinates(eof2) <-~lon+lat
coordinates(eof3) <-~lon+lat

#this visualizes the EOF's in space
precip.eof1 <- spplot(eof1, main ="EOF1 Precipitation")
precip.eof2 <- spplot(eof2, main = "EOF2 Precipitation")
precip.eof3 <- spplot(eof3, main = "EOF3 Precipitation")


##the next steps:
#1. do this for precip? and soil moisture data
#2. correlate the PC's associated with EOFs to the zonal indices



########################################################
#EOF analysis for evaporation
###########################################################
dataset<-do.call("cbind", lapply(file_list,FUN=function(files){read.table(files)$V4}))
colnames(dataset)<-substr(file_list, 14, 30) #this takes the latlon from each file makes it the header for each row

latlong<-cbind(substr(file_list,14,21),substr(file_list,23,30)) #create latlong vector from file names


######################
##EOF Analysis
#####################
#express the space-time data as a 2D matrix with "M" rows (space) by "N" columns(time)
#transpose our data so that time is represented in the columns and space in rows

A <- t(as.matrix(dataset))  #convert to matrix and transpose

#A <- matrix(c(2,4,-6,8,1,2,-3,4,4,1,3,2), nrow=3, byrow=TRUE) #testing with dummy matrix
##preprocess the data by subtracting the long-term mean (mean of each row) from each row
##This creates an "anomoly" matrix
ltmeans <- apply(A, 1, mean)   #calculate the longterm means for each row and puts it in a list
X <- A - ltmeans   #subtract longterm means from the matrix A

## Calculate the eigenvectors (E) and Eigen Values EV of the covariance matrix C
esv <- svd(X) #svd analysis of a matrix returns a lise with d, u, and v objections
diag.s <- diag(esv$d) #this is equivalent to S in matlab
V <- esv$v #this should be columns of the EOF in time dimension
U <-esv$u

## calculate the z matrix
z <- diag.s %*% t(V) # z is the matrix multiplcation of S and V'
#rows of the z matrix should be the Principle components associated with each EOF

## z can be confirmed another way:
uta <- t(U) %*% X

## calculate the orginal data based on the first 2 EOFS
A <- U[,1:2] %*% z[1:2,] # this gives the first 2 rows in the first column as non-zero, but rounding, they are zero

#the Principle components here now represent time  series of length N
plot(z[1,], type='l')  #time series explaining the most variation in the data
plot(z[2,], type='l')  #time series for 2nd EOF
plot(z[3,], type='l')

#create separate objects for eofs 1-3
eof1 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,1]))
colnames(eof1) <- c("lat", "lon","eof1")
eof2 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,2]))
colnames(eof2) <- c("lat", "lon","eof2")
eof3 <- data.frame(cbind(as.numeric(latlong[,1]), as.numeric(latlong[,2]), U[,3]))
colnames(eof3) <- c("lat", "lon","eof3")

write.table(U,"evapU.txt",sep='\t')
write.table(Z, "evapPCs.txt", sep='\t')

##make the data spatial points using coordinates() from the sp package
coordinates(eof1) <-~lon+lat
coordinates(eof2) <-~lon+lat
coordinates(eof3) <-~lon+lat

#this visualizes the EOF's in space
evap.eof1 <- spplot(eof1, main ="EOF1 Evapotranspiration")
evap.eof2 <- spplot(eof2, main = "EOF2 Evapotranspiration")
evap.eof3 <- spplot(eof3, main = "EOF3 Evapotranspiration")


#######################################################

dataset<-do.call("cbind", lapply(file_list,FUN=function(files){read.table(files)$V8}))
colnames(dataset)<-substr(file_list, 14, 30) #this takes the latlon from each file makes it the header for each row

latlong<-cbind(substr(file_list,14,21),substr(file_list,23,30)) #create latlong vector from file names

##Now that data is formatted properly, can do analyses: EOF 

#timeseriesplot for 1 gridcell
#soil1.time<-ts(dataset[,1], freq = 12, start=c(1950,1))
#plot.ts(soil1.time, main="Time Series of Soil Moisture layer 1")

######################
##EOF Analysis
#####################
#express the space-time data as a 2D matrix with "M" rows (space) by "N" columns(time)

#transpose our data so that time is represented in the columns and space in rows

A <- t(as.matrix(dataset))  #convert to matrix and transpose

#A <- matrix(c(2,4,-6,8,1,2,-3,4,4,1,3,2), nrow=3, byrow=TRUE) #testing with dummy matrix
##preprocess the data by subtracting the long-term mean (mean of each row) from each row
##This creates an "anomoly" matrix
ltmeans <- apply(A, 1, mean)   #calculate the longterm means for each row and puts it in a list
X <- A - ltmeans   #subtract longterm means from the matrix A

#dont need this if we use svd()
##Calculate the covariance matrix, C for the anomoly matrix X
#C_1 <- X %*% t(X) #by hand
#C <- cov(X) #using function cov

## Calculate the eigenvectors (E) and Eigen Values EV of the covariance matrix C
esv <- svd(X) #svd analysis of a matrix returns a lise with d, u, and v objections
diag.s <- diag(esv$d) #this is equivalent to S in matlab
V <- esv$v #this should be columns of the EOF in time dimension
U <-esv$u


## calculate the z matrix
z <- diag.s %*% t(V) # z is the matrix multiplcation of S and V'
#rows of the z matrix should be the Principle components associated with each EOF

## z can be confirmed another way:
uta <- t(U) %*% X

## calculate the orginal data based on the first 2 EOFS
A <- U[,1:2] %*% z[1:2,] # this gives the first 2 rows in the first column as non-zero, but rounding, they are zero


#the Principle components here now represent time  series of length N
plot(z[1,],type='l')  #time series explaining the most variation in the data
plot(z[2,],type='l')  #time series for 2nd EOF

#create separate objects for eofs 1-3
eof1 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,1]))
colnames(eof1) <- c("lat", "lon","eof1")
eof2 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,2]))
colnames(eof2) <- c("lat", "lon","eof2")
eof3 <- data.frame(cbind(as.numeric(latlong[,1]), as.numeric(latlong[,2]), U[,3]))
colnames(eof3) <- c("lat", "lon","eof3")
write.table(U, "soilm1U.txt", sep="\t") # write a table with EOFS
write.table(z, "soilm1Z.txt", sep="\t")
##make the data spatial points using coordinates() from the sp package
coordinates(eof1) <-~lon+lat
coordinates(eof2) <-~lon+lat
coordinates(eof3) <-~lon+lat

#this visualizes the EOF's in space
soilm1.eof1 <- spplot(eof1, main ="EOF1 Soil moisture layer1")
soilm1.eof2 <- spplot(eof2, main = "EOF2 Soil moisture layer1")
soilm1.eof3 <- spplot(eof3, main = "EOF3 Soil moisture layer1")

##Hypothesis: some soil moisture dipole effects will be correlated to zonal index (likely EOF2, with a dipole like effect)
#to test this, I will use the PC1, PC2, and PC3 associated with the EOF1, EOF2, and EOF3
#PC1 is z[,1]
plot(z[1,],type='l') 
PNO.1950<-PNO[1:768,] #takes only 1950-2013 data from PNO index

#find correlation between PNO and Z[2,]
cor(z[1,], PNO.1950$INDEX)
cor(z[2,],PNO.1950$INDEX)
cor(z[3,], PNO.1950$INDEX)

#create linear model fits for this
PC2.lm<-lm(z[2,]~PNO.1950$INDEX)
##################################################3
#soil moisture layer 2
####################################################
dataset<-do.call("cbind", lapply(file_list,FUN=function(files){read.table(files)$V9}))
colnames(dataset)<-substr(file_list, 14, 30) #this takes the latlon from each file makes it the header for each row

latlong<-cbind(substr(file_list,14,21),substr(file_list,23,30)) #create latlong vector from file names

##Now that data is formatted properly, can do analyses: EOF 

#timeseriesplot for 1 gridcell
#soil1.time<-ts(dataset[,1], freq = 12, start=c(1950,1))
#plot.ts(soil1.time, main="Time Series of Soil Moisture layer 1")

######################
##EOF Analysis
#####################
#express the space-time data as a 2D matrix with "M" rows (space) by "N" columns(time)

#transpose our data so that time is represented in the columns and space in rows

A <- t(as.matrix(dataset))  #convert to matrix and transpose

#A <- matrix(c(2,4,-6,8,1,2,-3,4,4,1,3,2), nrow=3, byrow=TRUE) #testing with dummy matrix
##preprocess the data by subtracting the long-term mean (mean of each row) from each row
##This creates an "anomoly" matrix
ltmeans <- apply(A, 1, mean)   #calculate the longterm means for each row and puts it in a list
X <- A - ltmeans   #subtract longterm means from the matrix A

#dont need this if we use svd()
##Calculate the covariance matrix, C for the anomoly matrix X
#C_1 <- X %*% t(X) #by hand
#C <- cov(X) #using function cov

## Calculate the eigenvectors (E) and Eigen Values EV of the covariance matrix C
esv <- svd(X) #svd analysis of a matrix returns a lise with d, u, and v objections
diag.s <- diag(esv$d) #this is equivalent to S in matlab
V <- esv$v #this should be columns of the EOF in time dimension
U <-esv$u

## calculate the z matrix
z <- diag.s %*% t(V) # z is the matrix multiplcation of S and V'
#rows of the z matrix should be the Principle components associated with each EOF

## z can be confirmed another way:
uta <- t(U) %*% X

## calculate the orginal data based on the first 2 EOFS
A <- U[,1:2] %*% z[1:2,] # this gives the first 2 rows in the first column as non-zero, but rounding, they are zero

#the Principle components here now represent time  series of length N
plot(z[1,],type='l')  #time series explaining the most variation in the data
plot(z[2,],type='l')  #time series for 2nd EOF

write.table(U, "soilm2U.txt", sep="\t")
write.table(z, "soilm2Z.txt", sep="\t")
#create separate objects for eofs 1-3
eof1 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,1]))
colnames(eof1) <- c("lat", "lon","eof1")
eof2 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,2]))
colnames(eof2) <- c("lat", "lon","eof2")
eof3 <- data.frame(cbind(as.numeric(latlong[,1]), as.numeric(latlong[,2]), U[,3]))
colnames(eof3) <- c("lat", "lon","eof3")

##make the data spatial points using coordinates() from the sp package
coordinates(eof1) <-~lon+lat
coordinates(eof2) <-~lon+lat
coordinates(eof3) <-~lon+lat

#this visualizes the EOF's in space
soilm2.eof1 <- spplot(eof1, main ="EOF1 Soil moisture layer2")
soilm2.eof2 <- spplot(eof2, main = "EOF2 Soil moisture layer2")
soilm2.eof3 <- spplot(eof3, main = "EOF3 Soil moisture layer2")



####################################
#soil moisture layer 3
####################################
dataset<-do.call("cbind", lapply(file_list,FUN=function(files){read.table(files)$V10}))
colnames(dataset)<-substr(file_list, 14, 30) #this takes the latlon from each file makes it the header for each row

latlong<-cbind(substr(file_list,14,21),substr(file_list,23,30)) #create latlong vector from file names

##Now that data is formatted properly, can do analyses: EOF 

#timeseriesplot for 1 gridcell
#soil1.time<-ts(dataset[,1], freq = 12, start=c(1950,1))
#plot.ts(soil1.time, main="Time Series of Soil Moisture layer 1")

######################
##EOF Analysis
#####################
#express the space-time data as a 2D matrix with "M" rows (space) by "N" columns(time)

#transpose our data so that time is represented in the columns and space in rows

A <- t(as.matrix(dataset))  #convert to matrix and transpose

#A <- matrix(c(2,4,-6,8,1,2,-3,4,4,1,3,2), nrow=3, byrow=TRUE) #testing with dummy matrix
##preprocess the data by subtracting the long-term mean (mean of each row) from each row
##This creates an "anomoly" matrix
ltmeans <- apply(A, 1, mean)   #calculate the longterm means for each row and puts it in a list
X <- A - ltmeans   #subtract longterm means from the matrix A

#dont need this if we use svd()
##Calculate the covariance matrix, C for the anomoly matrix X
#C_1 <- X %*% t(X) #by hand
#C <- cov(X) #using function cov

## Calculate the eigenvectors (E) and Eigen Values EV of the covariance matrix C
esv <- svd(X) #svd analysis of a matrix returns a lise with d, u, and v objections
diag.s <- diag(esv$d) #this is equivalent to S in matlab
V <- esv$v #this should be columns of the EOF in time dimension
U <-esv$u

## calculate the z matrix
z <- diag.s %*% t(V) # z is the matrix multiplcation of S and V'
#rows of the z matrix should be the Principle components associated with each EOF

## z can be confirmed another way:
uta <- t(U) %*% X

## calculate the orginal data based on the first 2 EOFS
A <- U[,1:2] %*% z[1:2,] # this gives the first 2 rows in the first column as non-zero, but rounding, they are zero

#the Principle components here now represent time  series of length N
plot(z[1,],type='l')  #time series explaining the most variation in the data
plot(z[2,],type='l')  #time series for 2nd EOF

write.table(U, "soilm3U.txt", sep="\t")
write.table(z, "soilm3Z.txt", sep="\t")
#create separate objects for eofs 1-3
eof1 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,1]))
colnames(eof1) <- c("lat", "lon","eof1")
eof2 <- data.frame(cbind(as.numeric(latlong[,1]),as.numeric(latlong[,2]),U[,2]))
colnames(eof2) <- c("lat", "lon","eof2")
eof3 <- data.frame(cbind(as.numeric(latlong[,1]), as.numeric(latlong[,2]), U[,3]))
colnames(eof3) <- c("lat", "lon","eof3")

##make the data spatial points using coordinates() from the sp package
coordinates(eof1) <-~lon+lat
coordinates(eof2) <-~lon+lat
coordinates(eof3) <-~lon+lat

#this visualizes the EOF's in space
soilm3.eof1 <- spplot(eof1, main ="EOF1 Soil moisture layer3")
soilm3.eof2 <- spplot(eof2, main = "EOF2 Soil moisture layer3")
soilm3.eof3 <- spplot(eof3, main = "EOF3 Soil moisture layer3")



###save the soil moisture EOF's to a pdf

pdf("EOFmaps.pdf")
precip.eof1
precip.eof2
precip.eof3
soilm1.eof1
soilm1.eof2
soilm1.eof3
soilm2.eof1
soilm2.eof2
soilm2.eof3
soilm3.eof1
soilm3.eof2
soilm3.eof3
dev.off()
