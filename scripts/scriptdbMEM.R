
# load data ---------------------------------------------------------------
# Communidade
spp<-read.csv("data/densi.csv", header=T, sep=";", dec = ",", row.names=1)  

# Ambiente: topografia
env.1<-read.csv("data/topo.csv", header=T, sep=";", dec = ",", row.names=1)

# Ambiente: quimica do solo no periodo de seca ou chuva
# env.2<-read.csv("data/solo_sc.csv", header= T, sep= ";", dec= ",", row.names=1) 
env.2<-read.csv("data/solo_ch.csv", header= T, sep= ";", dec= ",", row.names=1)

# Espaço (latitude, longitude)
spa<-read.csv("data/space.csv", header=T, sep=";", dec = ",", row.names=1) 


# load packages -----------------------------------------------------------

library(vegan)     
library(stats) # contem a função poly
# ?poly
library(packfor) 
library(PCNM)

# Standardization  --------------------------------------------------------
# padronizacao dos dados da comunidade segundo a soma dos valores registrados na unidade amostral

spp.t<-decostand(spp, "total")   
spp.bray<-vegdist(spp.t, "horn")  
spp.pcoa<-cmdscale(spp.bray, k=nrow(spp)-1, eig=TRUE, add=TRUE)    

# transformando variavel aspecto
trf.asp<-cbind(sin(env.1[,2]), cos(env.1[,2]))   
colnames(trf.asp)<-c("asp1", "asp2")    
env <- cbind(env.2, env.1, trf.asp)   
str(env)


# Compute Orthogonal Polynomials and name rows and columns ----------------
# nomear para não ter problemas depois porque as colunas ficam com "nomes" 1, 2, 1, 2, ...

envP<-cbind(poly(env$elev,2),    
            poly(env$MO,2),      
            poly(env$pH,2),      
            poly(env$K,2),        
            poly(env$Ca,2),      
            poly(env$Mg,2),
            poly(env$SB,2),
            poly(env$CTC,2),   
            poly(env$SLT,2),
            poly(env$ARG,2),
            poly(env$N,2),
            poly(env$sol,2),
            poly(env$P,2),
            poly(env$H_Al,2),
            poly(env$V,2),
            poly(env$asp1,2),
            poly(env$asp2,2),
            poly(env$AF,2),
            poly(env$slop,2))

colnames(envP)<-c("elev", "elev^2", "MO", "MO^2", "pH", "pH^2", "K", "K^2", "Ca", 
                  "Ca^2", "Mg", "Mg^2", "SB", "SB^2", "CTC", "CTC^2", "SLT", "SLT^2",   
                  "ARG", "ARG^2", "N", "N^2", "sol", "sol^2", "P", "P^2", "H_Al", "H_Al^2", 
                  "V", "V^2", "asp1", "asp1^2", "asp2", "asp2^2", "AF", "AF^2", "slop", "slop^2")  

# nomeacao das linhas do objeto "envP" de acordo com as linhas do objeto "env"
rownames(envP)<-rownames(env)   


# Factor Selection --------------------------------------------------------
# disociación de factores entre cuadráticos y lineares
# para faciltiar rodar la función fordwar.sel

envP.l<-subset(envP, select = - c (`elev^2`, `MO^2`, `pH^2`, `K^2`, `Ca^2`, `Mg^2`, 
                                   `SB^2`, `CTC^2`, `SLT^2`, `ARG^2`, `N^2`, `sol^2`, `P^2`, 
                                   `H_Al^2`, `V^2`, `asp1^2`, `asp2^2`, `AF^2`, `slop^2`))     
                                                                                                                            
envP.q<-subset(envP, select = - c (elev, MO, pH, K, Ca, Mg, SB, CTC, SLT, ARG, 
                                   N, sol, P, H_Al, V, asp1, asp2, AF, slop))   
                                                                                                                                
                                                                                                                                                                                                                                                    # 
# Criação de objeto com correlaçao entre fatores ambientais para teste de correlção

cor_envP.l<-cor(envP.l, envP.l)      

ifelse(abs(cor_envP.l)>0.7, "OUT", "-")  #selecionando correlaçoes maiores do 0.7

# Ajuste y seleccion de variables de los datos lineares

calib.adj.l<-RsquareAdj(rda(spp.pcoa$points ~  elev + H_Al + K + Mg + SLT + N +                                             
                                               sol + P + asp1 + asp2 + AF + slop, 
                            data= as.data.frame(envP.l)))$adj.r.squared    

envP.l2<-subset(envP.l, select= - c (MO, pH, CTC, SB, ARG, Ca, V))    

forward.sel(spp.pcoa$points, envP.l2,  adjR2thresh=calib.adj.l , alpha = 0.1)   

RsquareAdj(rda(spp.pcoa$points ~ elev + P + H_Al + K + N + sol + SLT, 
               data= as.data.frame(envP.l)))$adj.r.squared                                                                   

cor_envP.q<-cor(envP.q, envP.q)
ifelse(abs(cor_envP.q)>0.7, "OUT", "-")

# Ajuste y seleccion de variables de los datos cuadráticos
calib.adj.q<-RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `SLT^2` + 
                              `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `V^2` + `asp2^2` + `AF^2` +
                              `slop^2`, data= as.data.frame(envP.q)))$adj.r.squared                               

envP.q2<-subset(envP.q, select= - c (`SB^2`, `Ca^2`, `Mg^2`, `asp1^2`))   

forward.sel(spp.pcoa$points, envP.q2,  adjR2thresh= 0.2, alpha = 0.1)     

# aqui calculado o R2 final segundo o conjunto de fatores
RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `N^2` + `AF^2` + `K^2` , 
               data= as.data.frame(envP.q)))$adj.r.squared                                                          

# Creating the model: putting variables together again after fordw --------

envP_FIM<-cbind(envP.l2[,1], envP.l2[,8:9], envP.l2[,2], envP.l2[,6:7], 
                envP.l2[,4], envP.q2[,1], envP.q2[,9], envP.q2[,15], envP.q2[,4])                                                           

colnames(envP_FIM)<-c("Elev", "P", "H_Al", "K", "N", "Sol", "Silt", "Elev^2", "N^2", "AF^2", "K^2") 

# conferindo se por acaso existe ainda alguma correlaçao entre fatores
# Nada esta correlacionado. Beleza!

cor_envP_FIM<-cor(envP_FIM, envP_FIM)     
ifelse(abs(cor_envP_FIM)>0.7, "OUT", "-") 


# Testing model's significance --------------------------------------------

# centralizar dados espaciais
spa_cent<-scale(spa, scale=F, center=TRUE)                              
spa.rda<-capscale(spp.t ~ ., distance="horn", data = as.data.frame(spa_cent))  
anova(spa.rda)    

# Se nao for significativo, pode parar por aqui 
# e n?o precisa incluir os fatores espaciais lineares
# (tamb?m conhecidos como latitude e longitude...). 
# Como aqui ? significativo, damos continuidade a rotina.

# Second Factor's Selection -----------------------------------------------
## fatores espaciais lineares (X e Y = longitude e latitude) 

spa_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = as.data.frame(spa_cent)))$adj.r.squared   
forward.sel(spp.pcoa$points, spa_cent, adjR2thresh = spa_R2a, alpha = 0.05)  
spa.dist<-dist(spa_cent)                                   

# Distance matrix creation and dbMEM --------------------------------------

spa.PCNM<-PCNM(spa.dist, dbMEM=T, moran=T, all=TRUE)    

## selecionar os fatores espaciais refinados (dbMeMs) positivos
select<-which(spa.PCNM$Moran_I$Positive) 

## Fatores espaciais refinados submetida ao procedimento de sele?ao de variaveis
PCNMs<-as.data.frame(spa.PCNM$vectors)[,select] 

# testar significancia do modelo
PCNM.rda<-capscale(spp.t ~ ., distance="horn", data = PCNMs)  
anova.cca(PCNM.rda)    

# seleçao de fatores e criação do dataframe com as variaveis do espaço refinado
PCNM_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = PCNMs))$adj.r.squared     
forward.sel(spp.pcoa$points, PCNMs, adjR2thresh=PCNM_R2a, alpha=0.05)   
PCNM_FIM<-as.data.frame(cbind(PCNMs[,1]))   
colnames(PCNM_FIM)<-c("PCNM_1")   
rownames(PCNM_FIM)<-rownames(spa_cent)  


# Variance Partition ------------------------------------------------------

varpart(spp.pcoa$points, envP_FIM, spa_cent, PCNM_FIM)

# juntando as variaveis
fts.exp<-as.data.frame(cbind(envP_FIM, spa_cent, PCNM_FIM))   
# testando significancia do modelo completo
anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                           + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct")                 
summary(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`    
                         + X + Y + PCNM_1, distance= "horn", fts.exp))                                   

# Significancia dos termos que compoem o modelo
anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`     
                          + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="terms")

# Significancia dos eixos que compoem o modelo
anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`     
                   + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="axis")


# Disentangling factors ---------------------------------------------------
# Contribuição pura do ambiente
#STEP 1a
pure_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                        + Condition (X + Y + PCNM_1), distance= "horn", fts.exp)                              
anova.cca(pure_E, step= 10000, perm= 10000, model="reduced")                                        
# R2 ajustado do modelo referente a contribui?ao pura do ambiente
RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                + Condition (X + Y + PCNM_1), fts.exp))                                       


#STEP 1b
# Pure linear space
pure_S.l<-capscale(spp.t ~ X + Y + 
                     Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+
                                  PCNM_1), distance= "horn", fts.exp)                           
anova.cca(pure_S.l, step= 10000, perm= 10000, model= "reduced")                                               
RsquareAdj(rda(spp.pcoa$points ~ X + Y + 
                 Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` 
                            + `N^2` + `AF^2` + `K^2` + PCNM_1), fts.exp))                                      


#STEP 1c
# Pure refined space
pure_S.p<-capscale(spp.t ~  PCNM_1 + 
                     Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` +
                                  `N^2` + `AF^2` + `K^2`+ X + Y), distance= "horn", fts.exp)                        
anova.cca(pure_S.p, step= 100000, perm= 100000, model="reduced")                                              
RsquareAdj(rda(spp.pcoa$points ~  PCNM_1 + 
                 Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` +
                              `N^2` + `AF^2` + `K^2`+ X + Y), fts.exp))                                     


#STEP 2a 
# Linear space + Environment
S.l_j_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                         + X + Y + Condition (PCNM_1), distance= "horn", fts.exp)         

anova.cca(S.l_j_E, step=1000, perm= 1000, model="reduced")                                
RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                + X + Y + Condition (PCNM_1), fts.exp))                   


#STEP 2b
# Linear space + refined space
S.l_j_S.p<-capscale(spp.t ~   X + Y + PCNM_1 + 
                      Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` +`N^2` + 
                                   `AF^2` + `K^2`), distance= "horn", fts.exp)        
anova.cca(S.l_j_S.p, step= 1000, perm= 1000, model= "reduced")               
RsquareAdj(rda(spp.pcoa$points ~   X + Y + PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + 
                                                              `Elev^2` + `N^2` + `AF^2` + `K^2`), fts.exp))                    # linhas de codigo para definir o R2 ajustado do modelo referente a contribui?ao do espa?o linear junto ao espa?o refinado

#STEP 2c
# Refined space + Environment

E_j_S.p<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + 
                          PCNM_1 + Condition (X + Y), distance= "horn", fts.exp)                               
anova.cca(E_j_S.p, step= 1000, perm= 1000, model="reduced")                                
RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ 
                                 PCNM_1 + Condition (X + Y), fts.exp))                     


# dbRDA  ------------------------------------------------------------------

palm.rda<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` +
                     `K^2` + X + Y + PCNM_1, distance= "horn", fts.exp)   
plot(palm.rda, scaling=1, main="Triplot RDA spe.hor ~ allfcts - scaling 1 - wa scores")     

palm.sp.sc<-scores(palm.rda, choices= 1:2, scaling=1, display="sp")
arrows(0, 0, palm.sp.sc[,1], palm.sp.sc[,2], length= 0, lty=1, col="red")        

plot(palm.rda, scaling=2, main="Triplot RDA spe.hor ~ allfcts - scaling 2 - wa scores")
palm.sp.sc2<-scores(palm.rda, choices= 1:2, scaling=2, display="sp")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")    

plot(palm.rda, scaling=2, main="Triplot RDA spe.hor ~ allfcts - scaling 2 - wa scores")
palm.sp.sc2<-scores(palm.rda, choices= 1:2, scaling=2, display="sp")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")   


plot(palm.rda, scaling=2, display=c("sp","lc","cn"), 
     main="Triplot RDA spe.hor ~ allfcts - scaling 2 - sp+lc+cn ")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")  

plot(palm.rda, scaling=2, display=c("sp","lc"), main="Biplot RDA spe.hor ~ allfcts - scaling 2 - sp+lc")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red") 

plot(palm.rda, scaling=2, display=c("sp","cn"), main="Biplot RDA spe.hor ~ allfcts - scaling 2 - sp+cn")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")  

plot(palm.rda, scaling=2, display=c("lc","cn"), main="Biplot RDA spe.hor ~ allfcts - scaling 2 - lc+cn")



