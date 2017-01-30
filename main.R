library(XML)
library(tidyverse)

# Einlesen der Daten
patientsNS<-getNodeSet(xmlParse(file="./data/data.xml"),"//patients/patient")
data.df<-data.frame(matrix(ncol=3,nrow=0))
colnames(data.df)<-c("patID","t","inr")
for(index in 1:length(patientsNS)){
  patient<-patientsNS[[index]]
  
  t<-as.integer(sapply(getNodeSet(patient,"./medications/medication/timeDiff"),xmlValue))
  inr<-as.numeric(sapply(getNodeSet(patient,"./medications/medication/inr"),xmlValue))
  patID<-rep(as.integer(xmlAttrs(patient)[["patientID"]]),length.out=length(t))
  
  data.df<-rbind(data.df,data.frame(patID,t,inr))
  
}

# Loglike-Funktionen zur Schätzung der Hyperparameter

# Berechnen der loglikelihood für einen einzelnen Datensatz
loglike<-function(param,x, y) {
  
  k <- calcSigma(x, x, param[2]) + diag(x = param[1], length(x), length(x))
  
  ll <- -.5*t(y - param[3]) %*% solve(k) %*% (y - param[3]) - .5*log(det(k)) - .5*length(x)*log(2*pi)
  return(ll)
}

# Berechnen der kumulativen loglikelihood
loglike.sum<-function(param,data.df){
  
  ll<-0
  patientIDs<-unique(data.df$patID)
  for (patID in patientIDs){
    t<-data.df$t[data.df$patID==patID]
    inr<-data.df$inr[data.df$patID==patID]
    ll<-ll-1*loglike(param,t,inr)
  }
  
  return(ll)
}

# End Helper-Funktionen

# function calcHyperParams
calcHyperParams<-function(data.df,lower){
  
  hyperParams<-optim(par=c(2,4,2), loglike.sum, data=data.df, method = "L-BFGS-B", lower = lower, upper = c(Inf, Inf, Inf))
  
  hyperParams<-list(sigmasq=hyperParams$par[1],
                    l=hyperParams$par[2],
                    mu=hyperParams$par[3],
                    loglike.value=hyperParams$value,
                    message=as.character(hyperParams$message)
  )
  
  
  return(hyperParams)
}

# Funktionen zur Berechnug des Gauss-Prozesses

# Berechnen der Kovarianzmatrix
# parameter: 
# x1 Datenvektor 1
# x2 Datenvektor 2
# l  Längenparameter 
#
# Rückgabewert:
# Kovarianzmatrix

calcSigma<-function(x1,x2,l=1){
  
  sigma<-matrix(rep(0,length(x1)*length(x2)),nrow=length(x1))
  
  for(row in 1:nrow(sigma)){
    for(col in 1:ncol(sigma)){
      sigma[row,col]<-exp(-(x1[row] - x2[col])^2/(2*l^2))
    }
  }
  
  return(sigma)
} # calcSigma


# Berechnung des Gauss-Prozesses
#
# Parameter:
# x.test    Testwerte
# x.train   x-Trainingswerte = Zeitpunkte der Messungen
# y.train   y-Trainingswerte = gemessene inr-Werte
# l         Längenparamter des Gaussprozesses
# sigmasq   
# mu        Mittelwert
gpReg<-function(x.test,x.train, y.train,l=1, sigmasq=0, mu=0){
  
  
  k.train.train.inv<-solve(calcSigma(x1=x.train,x2=x.train,l=l))+diag(x=sigmasq,nrow=length(x.train),ncol=length(x.train))
  
  k.train.test<-calcSigma(x1=x.train,x2=x.test,l=l)
  
  k.test.train<-calcSigma(x1=x.test,x2=x.train,l=l)
  
  k.test.test<-calcSigma(x1=x.test,x2=x.test,l=l)
  
  vcov <- k.test.test - k.test.train %*% k.train.train.inv %*% k.train.test
  diag(vcov)[diag(vcov) < 0] <- 0 # Elimination von Werten Nahe 0
  mean <- mu  + k.test.train %*% k.train.train.inv %*% (y.train - mu)
  
  return(list(predictions = data.frame(x = x.test, mean = mean, lo95 = mean - 1.96*sqrt(diag(vcov)), hi95 = mean + 1.96*sqrt(diag(vcov))), vcov = vcov))
  
}


# Render GP-Plot
renderGPPlot<-function(gp){
  
  lo95=gp$predictions$lo95
  hi95=gp$predictions$hi95
  
  gp.plot<-ggplot(data=gp$predictions,aes(x, mean))+
    geom_line()
    #geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = .25)+ 
    #geom_point(aes(x = t, y = inr), data = data.frame(t,inr))
  
  return(gp.plot)
  
  
}

# Festlegen für welchen Patienten der Prozess berechnet werden soll:
patID<-28
# Schritt 1: Schätzung der Hyperparameter
# Entweder kumulativ aus dem Datensatz aller Patienten:
# hyperParams<-calcHyperParams(data.df,lower=c(0.1,90,0))
#
# oder aus einem einzelnen Patientendatensatz:
  hyperParams<-calcHyperParams(data.df[data.df$patID==patID,],lower=c(0.2,90,2))
  print(paste("Hyperparameter:","sigmasq:",hyperParams$sigmasq,"l:",hyperParams$l,"mu",hyperParams$mu))
# Schritt 2: Festlegung der Daten zum Patienten
  patient<-data.df[data.df$patID==patID,]
  x.train<-patient$t
  y.train<-patient$inr
  x.test<-seq(from=min(x.train)+1,to=max(x.train)-1,length.out=50)
  
#Schritt 3: Berechnung des Gauss Prozesses
  
gp<-gpReg(x.test=x.test,
          x.train=x.train,
          y.train=y.train,
          l=hyperParams$l,
          sigmasq=hyperParams$sigmasq,
          mu=hyperParams$mu)

print(gp$predictions$mean)