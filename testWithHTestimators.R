
library( MBHdesign)
library( terra)

data( volcano)

n <- 500

volcano <- rast( t( volcano))
Total.actual <- sum( values( volcano))

IP <- rast( volcano, nlyrs=2, names=c("even","gradient"), vals=NA)
values( IP$even) <- n / ncell( IP)
values( IP$gradient) <- rowFromCell( IP, 1:ncell(IP)) + colFromCell( IP, 1:ncell( IP))
values( IP$gradient) <- n * values( IP$gradient) / sum( values( IP$gradient))
#plot( IP)

### using quasiSamp.raster

#sample volcano with even IP
set.seed( 727)
samp <- quasiSamp.raster( n=n, inclusion.probs=IP$even)  #correct inclusion probs
y <- extract( volcano, as.matrix( samp[,1:2]))
ip <- extract( IP$even, as.matrix( samp[,1:2]))
t.hat.even <- sum( y / ip)

#sample volcano with uneven IP
set.seed( 727)
samp <- quasiSamp.raster( n=n, inclusion.probs=IP$gradient)  #correct inclusion probs
y <- extract( volcano, as.matrix( samp[,1:2]))
ip <- extract( IP$gradient, as.matrix( samp[,1:2]))
t.hat.gradient <- sum( y / ip)

print( c(Total.actual, t.hat.even, t.hat.gradient))

#### using quasiSamp

#sample volcano with even IP
set.seed( 727)
tmp <- t( as.matrix( ext( IP)))
my.area <- matrix( c( tmp[,1], rev( tmp[,1]), rep( tmp[,2], each=2)), ncol=2)
samp <- quasiSamp( n=n, study.area=my.area, potential.sites=crds( IP), inclusion.probs=values( IP$even))  #correct inclusion probs
y <- extract( volcano, as.matrix( samp[,1:2]))
ip <- extract( IP$even, as.matrix( samp[,1:2]))
t.hat.even <- sum( y / ip)

#sample volcano with uneven IP
set.seed( 727)
samp <- quasiSamp( n=n, study.area=my.area, potential.sites=crds( IP), inclusion.probs=values( IP$gradient))  #correct inclusion probs
y <- extract( volcano, as.matrix( samp[,1:2]))
ip <- extract( IP$even, as.matrix( samp[,1:2]))
t.hat.gradient <- sum( y / ip)

print( c(Total.actual, t.hat.even, t.hat.gradient))

#### using quasiSamp.cluster

#sample volcano with even IP
set.seed(727)
samp <- quasiSamp.cluster( nCluster=10, clusterSize=5, clusterRadius=4, inclusion.probs=IP$even)
#plot( volcano)
#points( samp[,1:2], pch=20, col='orange')
y <- extract( volcano, as.matrix( samp[,1:2]))
ip <- extract( IP$even, as.matrix( samp[,1:2]))
t.hat.even <- sum( y / ip)

#sample volcano with uneven IP
set.seed( 727)
samp <- quasiSamp.cluster( nCluster=10, clusterSize=5, clusterRadius=4, inclusion.probs=IP$gradient)
#plot( volcano)
#points( samp[,1:2], pch=20, col='orange')
y <- extract( volcano, as.matrix( samp[,1:2]))
ip <- extract( IP$even, as.matrix( samp[,1:2]))
t.hat.gradient <- sum( y / ip)

print( c(Total.actual, t.hat.even, t.hat.gradient))


