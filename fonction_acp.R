fonctionacp <- function(X, nb_composantes){
  
  nbre_lig=dim(X)[1]
  nbre_col=dim(X)[2]
  
  #Centrer et reduire la matrice
  Xcr <- scale(X, center = T, scale = T)*sqrt((nbre_lig)/(nbre_lig-1))

  
  #Calcul de la matrice de variance-covariance
  matcovariance <- cov(Xcr)
  
  #Calcul de la matrice de correlation
  matcorrelation <- cor(Xcr)
  
  #Matrice d'inertie, on prend Q = I_p et D = (1/n)*I_n
  Si <- t(Xcr) %*% Xcr * (1/nbre_lig)

  #Calcul de l'inertie
  Ig <- sum(diag(Si))
  
  #Pour controler la valeur de l'inertie
  WD <- Xcr %*% t(Xcr) * (1/nbre_lig)
  sumdiag <- eigen(WD)
  sum(sumdiag$values)
  
  #Calcul des valeurs propres de la matrice d'inertie 
  diagonaliser=eigen(Si)
  VP=round(diagonaliser$values, 3)
  #Stockage des valeurs propres non nulles dans un autre vecteur
  VP2=VP[-which(VP==0)]
  
  #Calcul des pourcentages cumules
  Pourcentage=round(VP2*100/Ig, 3)
  PourcentageCum=cumsum(VP2*100/Ig)
  
  Inertie=data.frame(AXE=1:length(VP2), VP2, Pourcentage , PourcentageCum)
  
  #Composantes Principales, chaque colonne correspond aux coord selon un axe
  
  Fmat = matrix(Xcr %*% diagonaliser$vectors[,1:length(VP2)], nrow = nbre_lig)
  Gmat = t( sqrt(diagonaliser$values[1:length(VP2)]) * t(diagonaliser$vectors[,1:length(VP2)]) )
  
  #Individus
  CTRi= round(((Fmat*Fmat*100)/nbre_lig)%*%diag(1/diagonaliser$values[1:length(VP2)]) ,3)
  QLTi= t(t(Fmat*Fmat) %*% diag(1/rowSums(Fmat*Fmat)) )
  
  
  #Variables
  CTRv =  round(((Gmat*Gmat*100))%*%diag(1/diagonaliser$values[1:length(VP2)]) ,3)
  QLTv= t(t(Gmat*Gmat) %*% diag(1/rowSums(Gmat*Gmat)) )
  
  return(list(
    VaP = VP2,
    iner=Inertie,
    Mcovar = matcovariance,
    Mcorre = matcorrelation,
    Miner = Si,
    coorindiv = Fmat[,1:nb_composantes],
    coorvar = Gmat[,1:nb_composantes] ,
    indcontri = CTRi[,1:nb_composantes],
    indqlt = QLTi[,1:nb_composantes],
    varcontri = CTRv[,1:nb_composantes],
    varqlt = QLTv[,1:nb_composantes]
    ))
}


#name                   Description                          
#1  "$VaP"             "Valeurs propres"                        
#2  "$iner"            "Data frame contenant inerties et interties cumulees"          
#3  "$Mcovar"          "Matrice de variance-covariance"           
#4  "$Mcorre"          "Matrice des correlations"
#5  "$Miner"           "Matrice d'inertie"             
#6  "$coorindiv"       "Coordonnees des individus"     
#7  "$coorvar"         "Coordonnees des variables"        
#8  "$indcontri"       "Contributions des individus"         
#9  "$indqlt"          "Qualites des individus"           
#10 "$varcontri"       "Qontributions des variables"   
#11 "$varqlt"          "Qualites des variables"  
 
