# Indica a pasta onde estão os arquivos das analises neste script 
setwd("C:/Users/Sara/Desktop/R/R/Doctorado/First")

# chamada dos dados referentes as spp da comunidade: abundancia x plot
spp<-read.csv("densi.csv", header=T, sep=";", dec = ",", row.names=1)   

# iniciação/chamada do pacote vegan
library(vegan)  

# padronização dos dados 
?decostand
spp.t<-decostand(spp, "total")  

# calculo da matriz de distancias
spp.distmx<-vegdist(spp.t, "horn")  

# Criaçao de matriz resposta da comunidade segundo a ordenaçao PCoA. 
?cmdscale 
spp.pcoa<-cmdscale(spp.distmx, k=nrow(spp)-1, eig=TRUE, add=TRUE)   

# topografia - que serve tanto para chuva quanto para a seca)
env.topog<-read.csv("topo.csv", header=T, sep=";", dec = ",", row.names=1)   

# transformaçao do fator aspecto,  
trf.asp<-cbind(sin(env.topog[,2]), cos(env.topog[,2]))   

# nomeação  das colunas do obejto "trf.asp"
colnames(trf.asp)<-c("asp1", "asp2")    

# dados ambientais (solo)
#env.soil<-read.csv("solo_sc.csv", header= T, sep= ";", dec= ",", row.names=1) 
env.soil<-read.csv("solo_ch.csv", header= T, sep= ";", dec= ",", row.names=1) 

# criaçao de objeto contendo todas os fatores ambientais usando a funçao "cbind"
env<-cbind(env.soil, env.topog, trf.asp)   

# iniciação/chamada do pacote stats, que contem a funçao "poly"
library(stats)   

# criaçao de objeto/matriz com todos os fatores ambientais 
#e seus respectivos valores quadráticos 
?poly
envP<-cbind(poly(env$elev,2),    
            poly(env$MO,2),      
            poly(env$pH,2),     
            poly(env$K,2),       
            poly(env$Ca,2),      
            poly(env$Mg,2),
            poly(env$SB,2),
            poly(env$CTC,2),
            #poly(env$AG,2),
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

# nomeação das colunas do objeto criado com auxílio da funçao "poly". Se não nomear
# terá problemas depois porque as colunas ficam com "nomes" 1, 2, 1, 2, ...

colnames(envP)<-c("elev", "elev^2", "MO", "MO^2", "pH", "pH^2", "K", "K^2", "Ca", "Ca^2",     
                  "Mg", "Mg^2", "SB", "SB^2", "CTC", "CTC^2", "SLT", "SLT^2",   
                  "ARG", "ARG^2", "N", "N^2", "sol", "sol^2", "P", "P^2", "H_Al", "H_Al^2", "V", "V^2"
                  , "asp1", "asp1^2", "asp2", "asp2^2", "AF", "AF^2", "slop", "slop^2")  

# nomeaçao das linhas do objeto "envP" de acordo com as linhas do objeto "env"
rownames(envP)<-rownames(env)   

# Dissociaçao das informaçoes/fatores ambientais lineares daqueles quadráticos. 

envP.l<-subset(envP, select = - c (`elev^2`, `MO^2`, `pH^2`, `K^2`, `Ca^2`, `Mg^2`, `SB^2`, `CTC^2`, `SLT^2`,       
                                   `ARG^2`, `N^2`, `sol^2`, `P^2`, `H_Al^2`, `V^2`, `asp1^2`, `asp2^2`, `AF^2`, `slop^2`))  

# Dissociaçao  das informaçoes/fatores ambientais quadraticos daqueles lineares. 

envP.q<-subset(envP, select = - c (elev, MO, pH, K, Ca, Mg, SB, CTC, SLT, ARG, N, sol, P, H_Al, V, asp1, asp2, AF, slop))        # 

# Criação de objeto com correlaçao entre fatores ambientais
cor_envP.l<-cor(envP.l, envP.l)   

# solicitaçao de matriz na tela de operações do R,
#respondendo quais correlaçoes sao maiores do 0.7
ifelse(abs(cor_envP.l)>0.7, "OUT", "-")  

# COMBINAÇÔES 

# OPÇAO 1: MANTER elev	+pH+	K+					ARG+	N	+sol+	P	+		asp1	+asp2	+AF+	slop
# OPÇAO 2: MANTER elev		pH	K					SLT	ARG		sol	P			asp1	asp2	AF	slop
# OPÇAO 3: MANTER elev				Ca					ARG	N	sol	P			asp1	asp2	AF	slop
# OPÇAO 4: MANTER elev					Mg				ARG	N	sol	P			asp1	asp2	AF	slop
# OPÇAO 5: MANTER elev						SB			ARG	N	sol	P			asp1	asp2	AF	slop  
# OPÇAO 6: MANTER MO	pH	K						ARG		sol	P			asp1	asp2	AF	slop   
# OPÇAO 7: MANTER MO			Ca					ARG	N	sol	P			asp1	asp2	AF	slop
# OPÇAO 8: MANTER  MO				Mg				ARG	N	sol	P			asp1	asp2	AF	slop
# OPÇAO 9: MANTER MO					SB			ARG	N	sol	P			asp1	asp2	AF	slop
# OPÇAO 10: MANTER  CTC	SLT	ARG		sol	P		V	asp1	asp2	AF	slop 
# OPÇAO 11: MANTER ARG	N	sol	P	H_Al	V	asp1	asp2	AF	slop  
# OPÇAO 12: MANTER pH	K						ARG	N	sol	P	H_Al		asp1	asp2	AF	slop
# OPÇAO 13: MANTER  pH	K				CTC	SLT	ARG		sol	P			asp1	asp2	AF	slop 
# OPÇAO 14: MANTER  SLT	ARG	N	sol	P	H_Al	V	asp1	asp2	AF	slop    
# OPÇAO 15: MANTER    elev	MO	pH	K	Ca	Mg	SB				N			H_Al	


# Criaçao de objeto para execuçao da seleçao de fatores.

calib.adj.l<-RsquareAdj(rda(spp.pcoa$points ~ elev +	pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	asp2 +	AF + slop, 
                            data= as.data.frame(envP.l)))$adj.r.squared  
calib.adj.l

envP.l2<-subset(envP.l, select= - c (MO, Ca,	Mg,	SB,	CTC, N, H_Al,	V))    

library(packfor)    
forward.sel(spp.pcoa$points, envP.l2,  adjR2thresh=calib.adj.l, alpha = 0.05)   
RsquareAdj(rda(spp.pcoa$points ~ SLT+P+elev, data= as.data.frame(envP.l)))$adj.r.squared   

# Conjuntos de fatores ambientais lineares selecionados 
# segundo cada opçao de combinaçao de fatores montada com base nas correlações.
# Escolher o maior R2. 

# OPÇAO 1: FICAM elev + P + K + sol R2 0.1696979        16,96
# OPÇAO 2: FICAM elev + P + SLT R2 0.1988946            19,88
# OPÇAO 3: FICAM elev + P + sol + Ca R2 0.1688526       16.88
# OPÇAO 4: FICAM elev + P + sol R2  0.1547667           15.47
# OPÇAO 5: FICAM elev + P + +SB + sol R2 0.1698414      16.98
# OPÇAO 6: FICAM  MO + P + K + sol+slop R2 0.159928     15.99
# OPÇAO 7: FICAM MO + P + Ca + sol + slop R2 0.1610222  16.10
# OPÇAO 8: FICAM MO + P + sol + slop R2 0.1461174       14.61
# OPÇAO 9: FICAM MO + P + sol R2 0.1304828              13.04
# OPÇAO 10: FICAM SLT+P+V+slop r2 0.1923694             19.23
# OPÇAO 11: FICAM N+P+V+slop R2 0.1368387               13.68
# OPÇAO 12: FICAM N+P+sol+slop R2 0.1326203             13.26
# OPÇAO 13: FICAM SLT+P+pH+slop R2 0.1871225            18.71
 

# Criação de objeto com correlaçao entre fatores ambientais quadradas
# solicitaçao de matriz na tela de operações do R, respondendo quais correlaçoes 
#sao maiores do 0.7

cor_envP.q<-cor(envP.q, envP.q)
ifelse(abs(cor_envP.q)>0.7, "OUT", "-")

#OPÇÃO 1 `elev^2` +	`MO^2`+	`K^2`	+	`Mg^2`+	`SB^2`+	`SLT^2`+	`ARG^2`	+`N^2`+	`sol^2`+	`P^2`	+	`asp2^2`+	`AF^2`+	`slop^2`
#OPÇÃO 2 `elev^2` +	`MO^2`+	`K^2`	+	`Mg^2` + `SB^2`	+	`SLT^2` + `ARG^2` +	`N^2`	+ `sol^2` + `P^2`	+	`asp1^2` + `AF^2`+	`slop^2`
#OPÇÃO 3 `elev^2` +	`MO^2` +	`pH^2` +	`K^2` +	`Ca^2` +	`Mg^2` + `SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` + `asp2^2` +	`AF^2` + `slop^2`
#OPÇÃO 4 `elev^2` +	`MO^2` +	`pH^2` +	`K^2` +	`Ca^2` +	`Mg^2` + `SLT^2` +	`ARG^2` +	`N^2`+	`sol^2` +	`P^2` +	`asp1^2` + `AF^2` +	`slop^2`
#OPÇÃO 5 `elev^2` +	`K^2` +	`Mg^2` + `SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` +	`H_Al^2` +	`V^2` +	`asp2^2` +	`AF^2` +	`slop^2`
#OPÇÃO 6 `elev^2`+	`K^2` +		`Mg^2` + `SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` +	`H_Al^2` +	`V^2` +	`asp1^2` + `AF^2` +	`slop^2`
#OPÇÃO 7 `elev^2` +	`MO^2` +	`K^2` +	`Mg^2` +	`SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` +	`V^2` +	`asp1^2` +	`AF^2` +	`slop^2`
#OPÇÃO 8 `elev^2` +	`MO^2` +	`K^2` +	`Mg^2` + `SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` +	`V^2` +	`asp2^2` +	`AF^2` +	`slop^2`
#OPÇÃO 9 `elev^2` +	`pH^2` +	`K^2` +	`Mg^2` +	`SB^2` +	`CTC^2` +	`SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` +	`asp2^2` +	`AF^2` +	`slop^2`
#OPÇÃO 10 `elev^2`+		`K^2` +		`Mg^2`+	`SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` +	`H_Al^2` +	`V^2`	+ `asp1^2` +	`AF^2`+	`slop^2`
#OPÇÃO 11 `elev^2` +	`K^2` +	`Mg^2` + `SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` +	`H_Al^2` +	`V^2` + `asp2^2` +	`AF^2` +	`slop^2`


# Seleçao de variaveis dos fatores ambientais quadráticos.

calib.adj.q<-RsquareAdj(rda(spp.pcoa$points ~ `elev^2` +	`MO^2` +	`pH^2` +	`K^2` +	`Ca^2` +	`Mg^2` + `SLT^2` +	`ARG^2` 
                            +	`N^2` +	`sol^2` +	`P^2` + `asp2^2` +	`AF^2` + `slop^2`, data= as.data.frame(envP.q)))$adj.r.squared    
calib.adj.q

#elev^2	MO^2	pH^2	K^2	Ca^2	Mg^2	SB^2	CTC^2	SLT^2	ARG^2	N^2	sol^2	P^2	H_Al^2	V^2	asp1^2	asp2^2	AF^2	slop^2

envP.q2<-subset(envP.q, select= - c (`SB^2`, `CTC^2`, `H_Al^2` ,	`V^2`,	`asp1^2`))   

# alpha menos restritivo, de 0.10. este foi menos restritivo do adj.q para poder tener algo
forward.sel(spp.pcoa$points, envP.q2,  adjR2thresh= 0.2, alpha = 0.1)     

# R2 final segundo o conjunto de fatores ambientais quadráticos selecionado
RsquareAdj(rda(spp.pcoa$points ~ `elev^2` +	`AF^2` +`ARG^2`+`K^2`+ `MO^2`+`P^2`+`Ca^2`,	data= as.data.frame(envP.q)))$adj.r.squared    

# Fatores ambientais quadraticos que foi selecionado em cada uma das opçoes

# OPÇAO 1:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+`P^2`+`MO^2`+`SB^2` R2= 0.1954707         19,54
# OPÇAO 2:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+`P^2`+`MO^2`+`SB^2` R2= 0.1954707         19,54
# OPÇAO 3:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+`P^2`+`MO^2`+`Ca^2` R2= 0.1960484         19,60
# OPÇAO 4:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+`P^2`+`MO^2`+`Ca^2` R2= 0.1960484         19,60
# OPÇAO 5:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+`P^2`+`SLT^2`       R2= 0.1790101         17,90
# OPÇAO 6:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+`P^2`+`SLT^2`       R2= 0.1790101         17,90
# OPÇAO 7:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+ `MO^2`+`P^2`+`V^2` R2= 0.1952378         19,52
# OPÇAO 8:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+ `MO^2`+`P^2`+`V^2` R2= 0.1952378         19,52
# OPÇAO 9:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+ `SLT^2`+`P^2`      R2 0.1790101          17,90
# OPÇAO 10:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+ `SLT^2`+`P^2`      R2 0.1790101          17,90
# OPÇAO 11:  `elev^2` +	`AF^2` +`ARG^2`+`K^2`+ `SLT^2`+`P^2`      R2 0.1790101          17,90



# ficou elev + P + SLT + `elev^2` +	`AF^2` +`ARG^2`+`K^2`+`P^2`+`MO^2`+`Ca^2


# matriz de fatores ambientais (quadráticos e lineares) final 

envP_FIM<-cbind(envP.l2[,1], envP.l2[,4], envP.l2[,7], envP.q2[,1:2],envP.q2[,8], envP.q2[,13], envP.q2[,4:5],  envP.q2[,11])   
# Nomeando as colunas da matriz final de fatores ambientais construida na linha de comando acima
colnames(envP_FIM)<-c("elev", "SLT", "P", "Elev^2", "MO^2", "ARG^2", "AF^2", "K^2", "Ca^2", "P^2")   
envP_FIM
# conferindo alguma correlaçao entre os fatores ambientais lineares e os quadráticos, 
# Se Nada está correlacionado. Beleza!

cor_envP_FIM<-cor(envP_FIM, envP_FIM) 
ifelse(abs(cor_envP_FIM)>0.7, "OUT", "-")

# Dados espaciais brutos (latitude e longitude, em graus decimais, que é como está na planilha base)
spa<-read.csv("space.csv", header=T, sep=";", dec = ",", row.names=1)   

# Dados espaciais/geograficos centrados         
spa_cent<-scale(spa, scale=F, center=TRUE)       

# db_rda (com funçao capscale) para testar se o modelo é significativo
spa.rda<-capscale(spp.t ~ ., distance="horn", data = as.data.frame(spa_cent))  

# aqui eu testo o modelo de db-RDA 
anova(spa.rda)    

# O procedimento é o mesmo que o feito com os fatores ambientais. 
# Threshold de referencia para a seleçao de fatores. 

spa_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = as.data.frame(spa_cent)))$adj.r.squared 

# Seleçao de variaveis para o conjunto de fatores espaciais lineares (X e Y = longitude e latitude)
forward.sel(spp.pcoa$points, spa_cent, adjR2thresh = spa_R2a, alpha = 0.05)   

# Modelar os fatores espaciais refinados (distance based moran eigenvector maps) e a respectiva rotina para
# selecionar aqueles que ficarão nas análises finais. Essa linha de comando aqui cria uma matriz de distancia
# entre as unidades amostrais com base nos fatores espaciais lineares.
spa.dist<-dist(spa_cent)                                    

# aqui chamamos/abrimos o pacote que permite rodar a pcnm (que é a analise que determina os fatores espaciais lineares)
library(PCNM)   

# depois de aberto o pacote, criamos um objeto com o resultado da pcnm rodada nessa linha de comando 
# usando a matriz de distancia criada acima.
spa.PCNM<-PCNM(spa.dist, dbMEM=T, moran=T, all=TRUE)   

# seleciona os fatores espaciais refinados (dbMeMs) positivos
# cria um objeto/matriz com fatores espaciais refinados que sera 
# submetida a seguir ao procedimento de seleçao de variaveis

select<-which(spa.PCNM$Moran_I$Positive)   
PCNMs<-as.data.frame(spa.PCNM$vectors)[,select] 

# como feito para os fatores espaciais lineares, primeiro rodamos uma dbRDA para checar se o modelo é significativo 
PCNM.rda<-capscale(spp.t ~ ., distance="horn", data = PCNMs) 

# aqui checamos se o modelo é significativo. Como é, damos continuidade a rotina.
anova.cca(PCNM.rda)    

# aqui é a criaçao do threshold usado na seleçao de fatores espaciais refinados. mesma lógica que descrita mais acima
PCNM_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = PCNMs))$adj.r.squared  

# seleçao dos fatores espaciais refinados que serão mantidos nas análises
forward.sel(spp.pcoa$points, PCNMs, adjR2thresh=PCNM_R2a, alpha=0.05)   

# criaçao do objeto/matriz (que nesse caso é uma coluna só mesmo..) com o fator espacial refinado 
PCNM_FIM<-as.data.frame(cbind(PCNMs[,1]))   

# nomeaçao da coluna da matriz criada na linha de codigo acima
colnames(PCNM_FIM)<-c("PCNM_1")   

# nomeaçao das linhas da matriz criada duas linhas de comando acima  
rownames(PCNM_FIM)<-rownames(spa_cent) 


# PARTIÇÃO DE VARIANCIA final 
varpart(spp.pcoa$points, envP_FIM, spa_cent, PCNM_FIM)   

# aqui crio um objeto que é uma matriz juntando todos os fatores preditores. 
fts.exp<-as.data.frame(cbind(envP_FIM, spa_cent, PCNM_FIM))   


# Significancia do modelo completo
anova.cca(capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`
                   + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct")                 

# Resultados gerais da db_RDA
# Eu uso para olhar o escores dos preditores nos eixos da db-RDA
summary(capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`        
                 + X + Y + PCNM_1, distance= "horn", fts.exp))                                  

# Linhas de código para pedir a significancia dos termos que compoem o modelo
anova.cca(capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`     
                   + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="terms")

# linhas de codigo para pedir a significancia dos eixos que compoem o modelo
anova.cca(capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`     
                   + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="axis")



# Modelo referente a contribuiçao pura do ambiente
#STEP 1a
pure_E<-capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`    
                 + Condition (X + Y + PCNM_1), distance= "horn", fts.exp)      

# teste do modelo que define a contribuiçao pura do ambiente
anova.cca(pure_E, step= 10000, perm= 10000, model="reduced")    

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do ambiente
RsquareAdj(rda(spp.pcoa$points ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2` 
               + Condition (X + Y + PCNM_1), fts.exp))                                       


# Modelo referente a contribuiçao pura do espaço linear
#STEP 1b
pure_S.l<-capscale(spp.t ~ X + Y + Condition (elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`
                                              + PCNM_1), distance= "horn", fts.exp)

# teste do modelo que define a contribuiçao pura do espaço linear
anova.cca(pure_S.l, step= 10000, perm= 10000, model= "reduced")  

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do espaço linear
RsquareAdj(rda(spp.pcoa$points ~ X + Y + Condition (elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`
                                                    + PCNM_1), fts.exp))                                      

# Modelo referente a contribuiçao pura do espaço refinado
#STEP 1c
pure_S.p<-capscale(spp.t ~  PCNM_1 + Condition (elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`
                                                + X + Y), distance= "horn", fts.exp)                          

# teste do modelo que define a contribuiçao pura do espaço refinado
anova.cca(pure_S.p, step= 100000, perm= 100000, model="reduced") 

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do espaço refinado
RsquareAdj(rda(spp.pcoa$points ~  PCNM_1 + Condition (elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`
                                                      + X + Y), fts.exp))                                     


# Modelo referente a contribuiçao do ambiente junto do espaço linear
#STEP 2a
S.l_j_E<-capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`
                  + X + Y + Condition (PCNM_1), distance= "horn", fts.exp) 

# teste do modelo que define a contribuiçao do ambiente junto do espaço linear
anova.cca(S.l_j_E, step=1000, perm= 1000, model="reduced")      

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do ambiente junto do espaço linear
RsquareAdj(rda(spp.pcoa$points ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`
               + X + Y + Condition (PCNM_1), fts.exp))                   


# Modelo referente a contribuiçao do espaço linear junto ao espaço refinado
#STEP 2b
S.l_j_S.p<-capscale(spp.t ~   X + Y + PCNM_1 + Condition (elev + SLT + P + `Elev^2` + `MO^2` 
                                                          + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`), distance= "horn", fts.exp)        

# teste do modelo que define a contribuiçao do espaço linear junto ao espaço refinado
anova.cca(S.l_j_S.p, step= 1000, perm= 1000, model= "reduced")                               

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do espaço linear junto ao espaço refinado
RsquareAdj(rda(spp.pcoa$points ~   X + Y + PCNM_1 + Condition (elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + 
                                                                 `K^2` + `Ca^2` + `P^2`), fts.exp))                    


# Modelo referente a contribuiçao do ambiente junto ao espaço refinado
#STEP 2c
E_j_S.p<-capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2` +
                    PCNM_1 + Condition (X + Y), distance= "horn", fts.exp)  

# teste do modelo que define a contribuiçao do ambiente junto ao espaço refinado
anova.cca(E_j_S.p, step= 1000, perm= 1000, model="reduced")                                

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do ambiente junto ao espaço refinado
RsquareAdj(rda(spp.pcoa$points ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2`+
                 PCNM_1 + Condition (X + Y), fts.exp))                    

# linha de codigo para rodar dbRDA que servirá de base para as figuras
palm.rda<-capscale(spp.t ~ elev + SLT + P + `Elev^2` + `MO^2` + `ARG^2` + `AF^2` + `K^2` + `Ca^2` + `P^2` + X + Y + PCNM_1, distance= "horn", fts.exp)   


# basicamente, daqui para baixo o que eu fiz foi elaborar diferentes plots 
# com os dados da db-RDA gerada na linha de comando acima. Da uma olhada no livro do
# Legendre & Legendre 2012 (NUmerical Ecology pg 638 ou próximo) e no Borcard et al. 2011 (Numerical Ecology
# with R) para entender as diferenças entre os tipos de scaling e os tipos de escores 
# mostrados no plot... ok?
#PLOTS
plot(palm.rda, scaling=1, main="Wet Triplot RDA spe.hor ~ allfcts - scaling 1 - wa scores")    
palm.sp.sc<-scores(palm.rda, choices= 1:2, scaling=1, display="sp")
arrows(0, 0, palm.sp.sc[,1], palm.sp.sc[,2], length= 0, lty=1, col="red")        


plot(palm.rda, scaling=2, main="Wet Triplot RDA spe.hor ~ allfcts - scaling 2 - wa scores")
palm.sp.sc2<-scores(palm.rda, choices= 1:2, scaling=2, display="sp")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")    

plot(palm.rda, scaling=2, display=c("sp","lc","cn"), main="Wet Triplot RDA spe.hor ~ allfcts - scaling 2 - sp+lc+cn ")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")  

plot(palm.rda, scaling=2, display=c("sp","lc"), main="Wet Biplot RDA spe.hor ~ allfcts - scaling 2 - sp+lc")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red") 

plot(palm.rda, scaling=2, display=c("sp","cn"), main="Wet Biplot RDA spe.hor ~ allfcts - scaling 2 - sp+cn")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")  

plot(palm.rda, scaling=2, display=c("lc","cn"), main="Wet Biplot RDA spe.hor ~ allfcts - scaling 2 - lc+cn")



#RDA DE TUDO pega somente todos os dados das opoções lineares e quadraticas escolhidas, sem a segunda seleção.

#elev +	pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	asp2 +	AF + slop, `elev^2` +	`MO^2` +	`pH^2` +	`K^2` +	`Ca^2` +	`Mg^2` + `SLT^2` +	`ARG^2` 
#+	`N^2` +	`sol^2` +	`P^2` + `asp2^2` +	`AF^2` + `slop^2`


RsquareAdj(rda(spp.pcoa$points ~  elev + pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	asp2 +	AF + slop, data= as.data.frame(envP.l)))$adj.r.squared    

#Criaçao de objeto contendo valor de referencia como threshold para execuçao da seleçao de fatores

plot(rda(spp.pcoa$points ~  elev + pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	
           asp2 +	AF + slop+ `elev^2` +	`MO^2` +	`pH^2` +	`K^2` +	`Ca^2` +	
           `Mg^2` + `SLT^2` +	`ARG^2` +	`N^2` +	`sol^2` +	`P^2` + `asp2^2` +	`AF^2` + 
         	  `slop^2`, data= as.data.frame(envP)),main="Wet RDA")

summary(rda(spp.pcoa$points ~  elev + pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	asp2 +	AF + slop, data= as.data.frame(envP)))




#RDA DE TUDO pega somente todos os dados das opoções lineares e quadraticas escolhidas, sem a segunda seleção.

#elev +	pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	asp2 +	AF + slop, `elev^2` +	`MO^2` +	`pH^2` +	`K^2` +	`Ca^2` +	`Mg^2` + `SLT^2` +	`ARG^2` 
#+	`N^2` +	`sol^2` +	`P^2` + `asp2^2` +	`AF^2` + `slop^2`


RsquareAdj(rda(spp.pcoa$points ~  elev + pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	asp2 +	AF + slop, data= as.data.frame(envP.l)))$adj.r.squared    

#Criaçao de objeto contendo valor de referencia como threshold para execuçao da seleçao de fatores

plot(rda(spp.pcoa$points ~  elev + pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	
           asp2 +	AF + slop, data= as.data.frame(envP)),main="Wet RDA")

summary(rda(spp.pcoa$points ~  elev + pH +	K	+	SLT	+ ARG	+	sol	+ P	+	asp1 +	asp2 +	AF + slop, data= as.data.frame(envP)))

coef(envP)
