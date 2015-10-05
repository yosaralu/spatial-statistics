
# Indica a pasta onde estão os arquivos das analises neste script 
setwd("C:/Users/Sara/Desktop/R/R/Doctorado/First")

# chamada dos dados referentes as spp da comunidade: abundancia x plot
spp<-read.csv("densi.csv", header=T, sep=";", dec = ",", row.names=1)   

# iniciação/chamada do pacote vegan
library(vegan)  

# padronização dos dados da comunidade segundo a soma dos valores registrados na unidade amostral
# decostand is a function that provides some popular (and effective) standardization methods for community ecologists
# total: divide by margin total (default MARGIN = 1).
?decostand
spp.t<-decostand(spp, "total")  

# calculo da matriz de distancias
#(usei a medida de morisita-horn/também rodei com bray-curtis mas obtive uma % de explicaçao um pouco menor)
# a mesma medida de distancia deve ser usada ao avaliar os dados da seca e da chuva ok?? 

spp.distmx<-vegdist(spp.t, "horn")  

# Criaçao de matriz resposta da comunidade segundo a ordenaçao PCoA. 
# Ao usarmos a matriz da ordenaçao para representar a comunidade
# estamos rodando uma distance-based RDA e não a RDA mais clássica ok? 
# O principal beneficio dessa abordagem, no caso 
# desses nossos dados, é que podemos usar variáveis (spp) com baixa ocorrência. 
# A RDA clássica depende do registro de pelo menos cinco pontos para as variáveis... 
# A db-RDA também permite usar qualquer medida de distancia para 
# modelar os dados, enquanto a RDA classica é "fechada" nesse aspecto 
# (salvo algumas transformações, ver Legendre & Gallagher 2001 ou Gallagher & Legendre 2001).
# Classical multidimensional scaling of a data matrix. Also known as principal coordinates analysis (Gower, 1966).
?cmdscale 
spp.pcoa<-cmdscale(spp.distmx, k=nrow(spp)-1, eig=TRUE, add=TRUE)   

# chamada/emtrada dos dados ambientais (parte referente a topografia
# - que serve tanto para chuva quanto para a seca)
env.topog<-read.csv("topo.csv", header=T, sep=";", dec = ",", row.names=1)   

# transformaçao do fator aspecto, dentro dos dados de topografia, pelo seu seno e coseno (1 fator viram dois mesmo). 
# Nos trabalhos da Punchi-Manage et al. (2013, 2014) eles citam esta transformaçao para fatores 
# que sejam de dados circulares, como no caso dos 360 graus usados para registrar o fator aspecto...

trf.asp<-cbind(sin(env.topog[,2]), cos(env.topog[,2]))   

# nomeação  das colunas do obejto "trf.asp"
colnames(trf.asp)<-c("asp1", "asp2")    

# chamada/entrada dos dados ambientais (parte referente a quimica do solo no periodo de seca)
env.soil<-read.csv("solo_sc.csv", header= T, sep= ";", dec= ",", row.names=1) 

# criaçao de objeto contendo todas os fatores ambientais usando a funçao "cbind"
env<-cbind(env.soil, env.topog, trf.asp)   

# iniciação/chamada do pacote stats, que contem a funçao "poly"
library(stats)   

# criaçao de objeto/matriz com todos os fatores ambientais e seus respectivos valores quadráticos (ao quadrado)
# segundo a funçao "poly". Essa funçao já transforma e padroniza os fatores ambientais, então não é preciso qualquer
# outra trasnformação desses dados. A funçao poly calcula os valores quadráticos e garante que eles sejam ortogonais 
# entre si. Se você calcular direto o valor quadrático, os fatores terão uma correlação de 1...
# Por isso essa função é tão legal, ela padroniza os dados e garante que não sejam lineares.
?poly
envP<-cbind(poly(env$elev,2),    
            poly(env$MO,2),      
            poly(env$pH,2),     
            poly(env$K,2),       
            poly(env$Ca,2),      
            poly(env$Mg,2),
            poly(env$SB,2),
            poly(env$CTC,2),
            poly(env$AG,2),
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
                  "Mg", "Mg^2", "SB", "SB^2", "CTC", "CTC^2", "AG", "AG^2", "SLT", "SLT^2",   
                  "ARG", "ARG^2", "N", "N^2", "sol", "sol^2", "P", "P^2", "H_Al", "H_Al^2", "V", "V^2"
                  , "asp1", "asp1^2", "asp2", "asp2^2", "AF", "AF^2", "slop", "slop^2")  

# nomeaçao das linhas do objeto "envP" de acordo com as linhas do objeto "env"
rownames(envP)<-rownames(env)   

# Dissociaçao das informaçoes/fatores ambientais lineares daqueles quadráticos. 
# Aqui usei a funçao subset para selecionar uma parte das colunas da matriz de dados ambientais "envP", 
# criando o obejto "envP.l". 

envP.l<-subset(envP, select = - c (`elev^2`, `MO^2`, `pH^2`, `K^2`, `Ca^2`, `Mg^2`, `SB^2`, `CTC^2`, `AG^2`, `SLT^2`,       
                                   `ARG^2`, `N^2`, `sol^2`, `P^2`, `H_Al^2`, `V^2`, `asp1^2`, `asp2^2`, `AF^2`, `slop^2`))  

# Dissociaçao  das informaçoes/fatores ambientais quadraticos daqueles lineares. 
# A razão pela qual separei essas informaçoes foi para executar a seleçao de fatores ambientais com base em correlacao. 
# ou seja, para diminuir a colinearidade. Fiz isso porque sao muitos fatores
# e ficava melhor para visualizar as correlações e para rodar a função "forward.sel",
# que seleciona os fatores a serem usados na análise final (depois de lidar com a multicolinearidade...). 
# Repara que se eu não separasse seria uma matriz com 40 fatores.      

envP.q<-subset(envP, select = - c (elev, MO, pH, K, Ca, Mg, SB, CTC, AG, SLT, ARG, N, sol, P, H_Al, V, asp1, asp2, AF, slop))        # 

# Criação de objeto com correlaçao entre fatores ambientais
cor_envP.l<-cor(envP.l, envP.l)   

# solicitaçao de matriz na tela de operações do R, respondendo quais correlaçoes sao maiores do 0.7
ifelse(abs(cor_envP.l)>0.7, "OUT", "-")  

# Com base na tabela de correlações maiores ou iguais a 0.7, bolei diferentes combinações de fatores ambientais, 
# segundo a eliminaçao de diferentes fatores... 
# Essas combinaçoes foram então usadas na seleçao de variáveis ("forward.sel") usando a matriz de spp criada
# com a ordenaçao PCoA como resposta. Essas foram as opçoes que bolei. Existem outras, mas acho que cobri as principais para 
# esse conjunto de dados. Você precisará atualizar isso para os dados de chuva ok?

# OPÇAO 1: MANTER MO + SB + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OPÇAO 2: MANTER elev + SB + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop
# OPÇAO 3: MANTER MO + SB + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop   
# OPÇAO 4: MANTER MO + pH + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OPÇAO 5: MANTER MO + pH + K + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OPÇAO 6: MANTER MO + pH + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop     
# OPÇAO 7: MANTER MO + Ca + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 8: MANTER MO + Ca + K + CTC + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 9: MANTER MO + Ca + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 10: MANTER MO + H_Al + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 11: MANTER MO + H_Al + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 12: MANTER elev + SB + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop   
# OPÇAO 13: MANTER elev + pH + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 14: MANTER elev + pH + K + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OPÇAO 15: MANTER elev + pH + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop     
# OPÇAO 16: MANTER elev + Ca + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 17: MANTER elev + Ca + K + CTC + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 18: MANTER elev + Ca + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 19: MANTER elev + H_Al + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OPÇAO 20: MANTER elev + H_Al + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop

#######################------------------------------
#RDA DE TUDO pega somente todos os dados das opoções lineares e quadraticas escolhidas, sem a segunda seleção.
#DEIXA PARA O FINAL DA MELHOR OPCAO DA LINEAR COM A MELHOR DOS QUADRATICOS

#RsquareAdj(rda(spp.pcoa$points ~  elev + H_Al + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop, data= as.data.frame(envP.l)))$adj.r.squared    

#Criaçao de objeto contendo valor de referencia como threshold para execuçao da seleçao de fatores

#plot(rda(spp.pcoa$points ~  elev + H_Al + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop +
#                    `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + 
#                   `V^2` + `asp2^2` + `AF^2` + `slop^2`, data= as.data.frame(envP)))

#summary(rda(spp.pcoa$points ~  elev + H_Al + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop +
#       `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + 
#      `V^2` + `asp2^2` + `AF^2` + `slop^2`, data= as.data.frame(envP)))
#######################----------------------------

# Criaçao de objeto contendo valor de referencia como threshold para execuçao da seleçao de fatores.
# A seleçao de fatores com a funçao "forward.sel" deve ser usada com dois critérios de "stop". Um é o valor 
# de alpha (0.05 geralmente, mas eu usei um menos restritivo, de 0.10) e o outro é o R2 da rda 
# usando todos os fatores ambientais (depois da eliminaçao) segundo está aqui nessa linha (dupla) de código.

calib.adj.l<-RsquareAdj(rda(spp.pcoa$points ~  elev + H_Al + K + Mg + SLT + AG + N +                                            
                              sol + P + asp1 + asp2 + AF + slop, data= as.data.frame(envP.l)))$adj.r.squared  
calib.adj.l

# Aqui estou criando um objeto com o conjunto de fatores ambientais que não entraram na opçao 
# Novamente usando a funçao "subset" para criar um subconjunto de fatores
# daqueles existentes na matriz de fatores ambientais lineares. Isso precisa ser feito porque 
# a função "forward.sel" não roda como a regressao, onde especificamos os fatores que queremos
# usar. Ela pede que forneça um objeto contendo todos os fatores... e por isso ele precisa ser criado.

envP.l2<-subset(envP.l, select= - c (MO, pH, CTC, SB, AG, Ca, V))    

# chamada/abertura do pacote "packfor"
library(packfor)    

# execuçao da seleçao de variaveis para o conjunto de dados ambientais lineares
# segundo a opção de combinaçao de fatores escolhida. 
forward.sel(spp.pcoa$points, envP.l2,  adjR2thresh=calib.adj.l, alpha = 0.05)   

# determinação do R2 segundo o conjunto de fatores ambientais lineares selecionados na seleção de variáveis feita acima. 
RsquareAdj(rda(spp.pcoa$points ~ elev + P + H_Al + K + N + sol + SLT, data= as.data.frame(envP.l)))$adj.r.squared   

# abaixo estão os conjuntos de fatores ambientais lineares selecionados 
# segundo cada opçao de combinaçao de fatores montada com base nas correlações. 
# Cada conjunto de fatores selecionados é seguido de seu respectivo R2 ao ser usado na dbRDA
# Esse R2 serve de referencia para o conjunto final de fatores selecionado e que sera usado no modelo final. 
# Eu escolhi o conjunto de fatores que propiciou o maior R2. 

# OPÇAO 1: FICAM MO + sol + slop + SB + P + ARG + N    R2 = 15,79
# OPÇAO 2: FICAM elev + P + SB + sol + N + SLT    R2 = 16,32
# OPÇAO 3: FICAM MO + sol + slop + SB + P + ARG + N    R2 = 15,79
# OPÇAO 4: FICAM MO + sol + slop + SLT + P + pH + K + N   R2 = 19.01    
# OPÇAO 5: FICAM MO + sol + slop + SLT + P + pH + K + N   R2 = 19.01
# OPÇAO 6: FICAM MO + sol + slop + SLT + P + pH + K + N   R2 = 19.01 
# OPÇAO 7: FICAM MO + sol + slop + Ca + P + K + SLT + N   R2 = 19.62 
# OPÇAO 8: FICAM MO + sol + slop + Ca + P + K + SLT + N   R2 = 19.62 
# OPÇAO 9: FICAM MO + sol + slop + Ca + P + K + ARG + SLT  R2 = 19.51
# OPÇAO 10: FICAM MO + sol + slop + SLT + P + K + H_Al + N  R2 = 19.55
# OPÇAO 11: FICAM MO + sol + slop + SLT + P + K + H_Al + N  R2 = 19.55
# OPÇAO 12: FICAM elev + P + SB + sol + N + SLT  R2 = 16,32
# OPÇAO 13: FICAM elev + P + sol + pH + K + N + SLT  R2 = 17.91
# OPÇAO 14: FICAM elev + P + sol + pH + K + N + SLT  R2 = 17.91
# OPÇAO 15: FICAM elev + P + sol + pH + K + N + SLT  R2 = 17.91
# OPÇAO 16: FICAM elev + P + Ca + K + sol + N + SLT  R2 = 19.66
# OPÇAO 17: FICAM elev + P + Ca + K + sol + N + SLT  R2 = 19.66
# OPÇAO 18: FICAM elev + P + Ca + K + sol + N + SLT  R2 = 19.66
# OPÇAO 19: FICAM elev + P + H_Al + K + N + sol + SLT  R2 = 19.93
# OPÇAO 20: FICAM elev + P + H_Al + K + N + sol + SLT  R2 = 19.93


# Criação de objeto com correlaçao entre fatores ambientais quadradas
# solicitaçao de matriz na tela de operações do R, respondendo quais correlaçoes sao maiores do 0.7

cor_envP.q<-cor(envP.q, envP.q)
ifelse(abs(cor_envP.q)>0.7, "OUT", "-")
# OPÇAO 1: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`   
# OPÇAO 2: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`              
# OPÇAO 3: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`         
# OPÇAO 4: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`
# OPÇAO 5: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OPÇAO 6: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `V^2` + `asp1^2` + `AF^2` + `slop^2`
# OPÇAO 7: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OPÇAO 8: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OPÇAO 9: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OPÇAO 10: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `V^2` + `asp2^2` + `AF^2` + `slop^2`


# Aqui eu aplico de novo o procedimento para a seleçao de variaveis dos fatores ambientais quadráticos.

calib.adj.q<-RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + 
                              `V^2` + `asp2^2` + `AF^2` + `slop^2`, data= as.data.frame(envP.q)))$adj.r.squared    
calib.adj.q

#criaçao da matriz segundo o conjunto de fatores ambientais que serão submetidos
# a seleçao de fatores ambientais a seguir (segundo uma das opçoes de combinação de fatores 
# ambientais quadráticos com base nas correlaçoes). 
envP.q2<-subset(envP.q, select= - c (`SB^2`, `Ca^2`, `Mg^2`, `asp1^2`))   

# aqui está a seleção de fatores ambientais quadraticos com o threshold criado na linha 
# de código acima e um alpha menos restritivo, de 0.10. este foi menos restritivo do adj.q para poder tener algo
forward.sel(spp.pcoa$points, envP.q2,  adjR2thresh= 0.2, alpha = 0.1)     

# aqui é calculado o R2 final segundo o conjunto de fatores ambientais quadráticos selecionado
RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `N^2` + `AF^2` + `ARG^2` + `K^2` + `V^2`, data= as.data.frame(envP.q)))$adj.r.squared    

# abaixo segue o conjunto de fatores ambientais quadraticos que foi selecionado em cada uma das opçoes
# de combinaçao dos fatores segundo as correlações e os respectivos R2 das db-RDA usado esses fatores.
# Para os fatores ambientais quadráticos, o conjunto final de fatores selecionado foi sempre o mesmo. 

# OPÇAO 1: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 2: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 3: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 4: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 5: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 6: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 7: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 8: FICAM `elev^2` + `N^2` + `AF^2`   R2= 14.03
# OPÇAO 9: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OPÇAO 10: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03

# aqui está a construçao, com auxílio da função "cbind" da matriz de fatores ambientais (quadráticos
# e lineares) final a ser usada nos modelos finais e na partiçao de variancia.

envP_FIM<-cbind(envP.l2[,1], envP.l2[,8:9], envP.l2[,2], envP.l2[,6:7], envP.l2[,4], envP.q2[,1], envP.q2[,9], envP.q2[,15], envP.q2[,4], envP.q2[,8],envP.q2[,13])   

# Aqui estou nomeando as colunas da matriz final de fatores ambientais construida na linha de comando acima
colnames(envP_FIM)<-c("Elev", "P", "H_Al", "K", "N", "Sol", "Silt", "Elev^2", "N^2", "AF^2", "K^2","ARG^2", "V^2")   

# conferindo se por acaso existe ainda alguma correlaçao entre os fatores ambientais lineares e os quadráticos, 
# visto que estes foram selecionados separadamente. Nada está correlacionado. Beleza!

cor_envP_FIM<-cor(envP_FIM, envP_FIM) 
ifelse(abs(cor_envP_FIM)>0.7, "OUT", "-")

# Aqui entramos com os dados espaciais brutos (latitude e longitude, em graus decimais, que é como está na planilha base)
spa<-read.csv("space.csv", header=T, sep=";", dec = ",", row.names=1)   

# Antes de começar a trabalhar os dados espaciais/geograficos deve-se centra-los, como é feito nesta linha de comando         
spa_cent<-scale(spa, scale=F, center=TRUE)       

# Aqui eu crio um objeto que é uma db_rda (com funçao capscale) para testar se o modelo é significativo
spa.rda<-capscale(spp.t ~ ., distance="horn", data = as.data.frame(spa_cent))  

# aqui eu testo o modelo de db-RDA criado na linha de comando acima. 
# Se não for significativo, pode parar por aqui e não precisa incluir os fatores espaciais lineares
# (também conhecidos como latitude e longitude...). Como aqui é significativo, damos continuidade a rotina. 

anova(spa.rda)    

# O procedimento é o mesmo que o feito com os fatores ambientais. Sendo assim, nessa linha 
# de comando é criado um valor threshold de referencia para a seleçao de fatores. 

spa_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = as.data.frame(spa_cent)))$adj.r.squared 

# E aqui é executada a seleçao de variaveis para o conjunto de fatores espaciais lineares (X e Y = longitude e latitude)
# Sim, é para rodar para o conjunto de apeans 2 fatores sim... pelo menos e isso que recomenda o livro 
# do Borcard et al. (2011 - numerical ecology with R). Acabaram sendo selecionados os dois fatores (X e Y) mesmo.
forward.sel(spp.pcoa$points, spa_cent, adjR2thresh = spa_R2a, alpha = 0.05)   

# Depois de lidar com os fatores espaciais lineares, damos inicio nessa linha de codigo aqui a
# modelar os fatores espaciais refinados (distance based moran eigenvector maps) e a respectiva rotina para
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


# PARTIÇÃO DE VARIANCIA final usando a matriz de spp segundo a ordenaçao PCoA,
# e cada um dos conjuntos de dados preditores criados ao longo das linhas de código acima. 
varpart(spp.pcoa$points, envP_FIM, spa_cent, PCNM_FIM)   

# aqui crio um objeto que é uma matriz juntando todos os fatores preditores. 
# Faço isso para rodar as linhas de código a seguir, onde é testada significancia 
# do modelo inteiro, dos termos, dos eixos e das partições 
# (step 1a, step 1b, step 1c, ..., step 2c - da uma olhada na tabela do arquivo word que passei na primeira
# leva de resultados que se chama partição_RDA que tem os "steps"). 
fts.exp<-as.data.frame(cbind(envP_FIM, spa_cent, PCNM_FIM))   


# essas duas linhas de código (poderia ser uma so, pus em duas para 
# facilitar a visualizaçao) testam a significancia do modelo completo
anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ `ARG^2` + `V^2`  
                   + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct")                 

# Essas duas linhas de código são para pedir os resultados gerais da db_RDA
# Eu uso para olhar o escores dos preditores nos eixos da db-RDA
summary(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ `ARG^2` + `V^2`        
                 + X + Y + PCNM_1, distance= "horn", fts.exp))                                  

# Linhas de código para pedir a significancia dos termos que compoem o modelo
anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ `ARG^2` + `V^2`     
                   + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="terms")

# linhas de codigo para pedir a significancia dos eixos que compoem o modelo
anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ `ARG^2` + `V^2`     
                   + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="axis")




# Modelo referente a contribuiçao pura do ambiente
#STEP 1a
pure_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ + `ARG^2` + `V^2`   
                 + Condition (X + Y + PCNM_1), distance= "horn", fts.exp)      

# teste do modelo que define a contribuiçao pura do ambiente
anova.cca(pure_E, step= 10000, perm= 10000, model="reduced")    

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do ambiente
RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + `ARG^2` + `V^2`
               + Condition (X + Y + PCNM_1), fts.exp))                                       


# Modelo referente a contribuiçao pura do espaço linear
#STEP 1b
pure_S.l<-capscale(spp.t ~ X + Y + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  + `ARG^2` + `V^2`
                                              + PCNM_1), distance= "horn", fts.exp)

# teste do modelo que define a contribuiçao pura do espaço linear
anova.cca(pure_S.l, step= 10000, perm= 10000, model= "reduced")  

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do espaço linear
RsquareAdj(rda(spp.pcoa$points ~ X + Y + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  +`ARG^2` + `V^2`
                                                    + PCNM_1), fts.exp))                                      

# Modelo referente a contribuiçao pura do espaço refinado
#STEP 1c
pure_S.p<-capscale(spp.t ~  PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  + `ARG^2` + `V^2`
                                                + X + Y), distance= "horn", fts.exp)                          

# teste do modelo que define a contribuiçao pura do espaço refinado
anova.cca(pure_S.p, step= 100000, perm= 100000, model="reduced") 

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do espaço refinado
RsquareAdj(rda(spp.pcoa$points ~  PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  + `ARG^2` + `V^2`
                                                      + X + Y), fts.exp))                                     


# Modelo referente a contribuiçao do ambiente junto do espaço linear
#STEP 2a
S.l_j_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  + `ARG^2` + `V^2`
                  + X + Y + Condition (PCNM_1), distance= "horn", fts.exp) 

# teste do modelo que define a contribuiçao do ambiente junto do espaço linear
anova.cca(S.l_j_E, step=1000, perm= 1000, model="reduced")      

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do ambiente junto do espaço linear
RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  + + `ARG^2` + `V^2`
               + X + Y + Condition (PCNM_1), fts.exp))                   


# Modelo referente a contribuiçao do espaço linear junto ao espaço refinado
#STEP 2b
S.l_j_S.p<-capscale(spp.t ~   X + Y + PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + 
                                                            `Elev^2` + `N^2` + `AF^2` + `K^2` + `ARG^2` + `V^2`), distance= "horn", fts.exp)        

# teste do modelo que define a contribuiçao do espaço linear junto ao espaço refinado
anova.cca(S.l_j_S.p, step= 1000, perm= 1000, model= "reduced")                               

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do espaço linear junto ao espaço refinado
RsquareAdj(rda(spp.pcoa$points ~   X + Y + PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + 
                                                                 `Elev^2` + `N^2` + `AF^2` + `K^2`+ `ARG^2` + `V^2`), fts.exp))                    


# Modelo referente a contribuiçao do ambiente junto ao espaço refinado
#STEP 2c
E_j_S.p<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + `ARG^2` + `V^2`+
                    PCNM_1 + Condition (X + Y), distance= "horn", fts.exp)  

# teste do modelo que define a contribuiçao do ambiente junto ao espaço refinado
anova.cca(E_j_S.p, step= 1000, perm= 1000, model="reduced")                                

# linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do ambiente junto ao espaço refinado
RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+  `ARG^2` + `V^2`+
                 PCNM_1 + Condition (X + Y), fts.exp))                    




# linha de codigo para rodar dbRDA que servirá de base para as figuras
palm.rda<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ `ARG^2` + `V^2` + X + Y + PCNM_1, distance= "horn", fts.exp)   


# basicamente, daqui para baixo o que eu fiz foi elaborar diferentes plots 
# com os dados da db-RDA gerada na linha de comando acima. Da uma olhada no livro do
# Legendre & Legendre 2012 (NUmerical Ecology pg 638 ou próximo) e no Borcard et al. 2011 (Numerical Ecology
# with R) para entender as diferenças entre os tipos de scaling e os tipos de escores 
# mostrados no plot... ok?
#PLOTS
plot(palm.rda, scaling=1, main="Triplot RDA spe.hor ~ allfcts - scaling 1 - wa scores")    
palm.sp.sc<-scores(palm.rda, choices= 1:2, scaling=1, display="sp")
arrows(0, 0, palm.sp.sc[,1], palm.sp.sc[,2], length= 0, lty=1, col="red")        


plot(palm.rda, scaling=2, main="Triplot RDA spe.hor ~ allfcts - scaling 2 - wa scores")
palm.sp.sc2<-scores(palm.rda, choices= 1:2, scaling=2, display="sp")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")    

plot(palm.rda, scaling=2, main="Triplot RDA spe.hor ~ allfcts - scaling 2 - wa scores")
palm.sp.sc2<-scores(palm.rda, choices= 1:2, scaling=2, display="sp")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")   


plot(palm.rda, scaling=2, display=c("sp","lc","cn"), main="Triplot RDA spe.hor ~ allfcts - scaling 2 - sp+lc+cn ")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")  

plot(palm.rda, scaling=2, display=c("sp","lc"), main="Biplot RDA spe.hor ~ allfcts - scaling 2 - sp+lc")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red") 

plot(palm.rda, scaling=2, display=c("sp","cn"), main="Biplot RDA spe.hor ~ allfcts - scaling 2 - sp+cn")
arrows(0, 0, palm.sp.sc2[,1], palm.sp.sc2[,2], length= 0, lty=1, col="red")  

plot(palm.rda, scaling=2, display=c("lc","cn"), main="Biplot RDA spe.hor ~ allfcts - scaling 2 - lc+cn")



