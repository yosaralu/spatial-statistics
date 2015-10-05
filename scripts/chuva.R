
setwd("C:/Users/Sara/Desktop/R/R/Doctorado/First")

spp<-read.csv("densi.csv", header=T, sep=";", dec = ",", row.names=1)   # chamada dos dados referentes as spp da comunidade

library(vegan)     # iniciação/chamada do pacote vegan

spp.t<-decostand(spp, "total")    # padronização dos dados da comunidade segundo a soma dos valores registrados na unidade amostral

spp.bray<-vegdist(spp.t, "horn")  
spp.pcoa<-cmdscale(spp.bray, k=nrow(spp)-1, eig=TRUE, add=TRUE)    


env.1<-read.csv("topo.csv", header=T, sep=";", dec = ",", row.names=1)   
trf.asp<-cbind(sin(env.1[,2]), cos(env.1[,2]))   

colnames(trf.asp)<-c("asp1", "asp2")    # nomeação  das colunas do obejto "trf.asp"

#env.2<-read.csv("solo_sc.csv", header= T, sep= ";", dec= ",", row.names=1) # chamada/entrada dos dados ambientais (parte referente a quimica do solo no periodo de seca)
env.2<-read.csv("solo_ch.csv", header= T, sep= ";", dec= ",", row.names=1)



env<-cbind(env.2, env.1, trf.asp)   # criaçao de objeto contendo todas os fatores ambientais usando a funçao "cbind"

library(stats)    # iniciação/chamada do pacote stats, que contem a funçao "poly"

envP<-cbind(poly(env$elev,2),    # criaçao de objeto/matriz com todos os fatores ambientais e seus respectivos valores quadráticos (ao quadrado)
            poly(env$MO,2),      # segundo a funçao "poly". Essa funçao já transforma e padroniza os fatores ambientais, então não é preciso qualquer
            poly(env$pH,2),      # outra trasnformação desses dados. A funçao poly calcula os valores quadráticos e garante que eles sejam ortogonais entre
            poly(env$K,2),       # si. Se você calcular direto o valor quadrático, os fatores terão uma correlação de 1...Por isso essa função é tão legal, ela 
            poly(env$Ca,2),      # padroniza os dados e garante que não sejam lineares.
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

colnames(envP)<-c("elev", "elev^2", "MO", "MO^2", "pH", "pH^2", "K", "K^2", "Ca", "Ca^2",     # nomeação das colunas do objeto criado com auxílio da funçao "poly". Se não nomear
                  "Mg", "Mg^2", "SB", "SB^2", "CTC", "CTC^2", "SLT", "SLT^2",   # terá problemas depois porque as colunas ficam com "nomes" 1, 2, 1, 2, ...
                  "ARG", "ARG^2", "N", "N^2", "sol", "sol^2", "P", "P^2", "H_Al", "H_Al^2", "V", "V^2"
                 , "asp1", "asp1^2", "asp2", "asp2^2", "AF", "AF^2", "slop", "slop^2")  

rownames(envP)<-rownames(env)   # nomeaçao das linhas do objeto "envP" de acordo com as linhas do objeto "env"

envP.l<-subset(envP, select = - c (`elev^2`, `MO^2`, `pH^2`, `K^2`, `Ca^2`, `Mg^2`, `SB^2`, `CTC^2`, `SLT^2`,       # Dissociaçao das informaçoes/fatores ambientais lineares daqueles
                                   `ARG^2`, `N^2`, `sol^2`, `P^2`, `H_Al^2`, `V^2`, `asp1^2`, `asp2^2`, `AF^2`, `slop^2`))  # quadráticos. Aqui usei a funçao subset para selecionar uma parte     
                                                                                                                            # das colunas da matriz de dados ambientais "envP", criando o obejto
                                                                                                                            # "envP.l". 

envP.q<-subset(envP, select = - c (elev, MO, pH, K, Ca, Mg, SB, CTC, SLT, ARG, N, sol, P, H_Al, V, asp1, asp2, AF, slop))   # Dissociaçao  das informaçoes/fatores ambientais quadraticos daqueles
                                                                                                                                # lineares. A razão pela qual separei essas informaçoes foi para 
                                                                                                                                # executar a seleçao de fatores ambienatias com base em correlacao. 
                                                                                                                                # ou seja, para diminuir a colinearidade. Fiz isso porque sao muitos fatores
                                                                                                                                # e ficava melhor para visualizar as correlações e para rodar a função "forward.sel",
                                                                                                                                # que seleciona os fatores a serem usados na análise final (depois de lidar
                                                                                                                                # com a multicolinearidade...). Repara que se eu não separasse seria uma matriz
                                                                                                                                # com 40 fatores.                                                                                                                    # 


cor_envP.l<-cor(envP.l, envP.l)      # Criação de objeto com correlaçao entre fatores ambientais

ifelse(abs(cor_envP.l)>0.7, "OUT", "-")  # solicitaçao de matriz na tela de operações do R, respondendo quais correlaçoes sao maiores do 0.7

# Com base na tabela de correlações maiores ou iguais a 0.7, bolei diferentes combinações de fatores ambientais, segundo a eliminaçao de
# diferentes fatores... Essas combinaçoes foram então usadas na seleçao de variáveis ("forward.sel") usando a matriz de spp criada
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


calib.adj.l<-RsquareAdj(rda(spp.pcoa$points ~  elev + H_Al + K + Mg + SLT + N +                                             # Criaçao de objeto contendo valor de referencia
                                               sol + P + asp1 + asp2 + AF + slop, data= as.data.frame(envP.l)))$adj.r.squared    # como threshold para execuçao da seleçao de fatores.
                                                                                                                                 # A seleçao de fatores com a funçao "forward.sel" deve
                                                                                                                                 # ser usada com dois critérios de "stop". Um é o valor 
                                                                                                                                 # de alpha (0.05 geralmente, mas eu usei um menos restritivo, de 0.10) 
                                                                                                                                 # e o outro é o R2 da rda 
                                                                                                                                 # usando todos os fatores ambientais (depois da eliminaçao)
                                                                                                                                 # segundo está aqui nessa linha (dupla) de código.

envP.l2<-subset(envP.l, select= - c (MO, pH, CTC, SB, ARG, Ca, V))    # Aqui estou criando um objeto com o conjunto de fatores ambientais que serão submetidos
                                                                     # a seleçao de fatores ambientais a seguir (segundo uma das opçoes de combinação de fatores
                                                                     # com base nas correlaçoes). Novamente usando a funçao "subset" para criar um subconjunto de fatores
                                                                     # daqueles existentes na matriz de fatores ambientais lineares. Isso precisa ser feito porque 
                                                                     # a função "forward.sel" não roda como a regressao, onde especificamos os fatores que queremos
                                                                     # usar. Ela pede que forneça um objeto contendo todos os fatores... e por isso ele precisa ser criado.

library(packfor)    # chamada/abertura do pacote "packfor"

forward.sel(spp.pcoa$points, envP.l2,  adjR2thresh=calib.adj.l , alpha = 0.1)   # execuçao da seleçao de variaveis para o conjunto de dados ambientais lineares
                                                                                # segundo a opção de combinaçao de fatores escolhida. 

RsquareAdj(rda(spp.pcoa$points ~ elev + P + H_Al + K + N + sol + SLT, data= as.data.frame(envP.l)))$adj.r.squared   # determinação do R2 segundo o conjunto de fatores ambientais lineares
                                                                                                                    # selecionados na seleção. 

# abaixo estão os conjuntos de fatores ambientais lineares selecionados segundo cada opçao de combinaçao de fatores montada 
# com base nas correlações. Cada conjunto de fatores selecionados é seguido de seu respectivo R2 ao ser usado na dbRDA
# Esse R2 serve de referencia para o conjunto final de fatores selecionado e que sera usado no modelo final. Eu escolhi
# o conjunto de fatores que propiciou o maior R2. 

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


calib.adj.q<-RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + 
                                              `V^2` + `asp2^2` + `AF^2` + `slop^2`, data= as.data.frame(envP.q)))$adj.r.squared    # Aqui eu aplico de novo o procedimento para a seleçao de variaveis ambientais
                                                                                                                                   # porém agora é os fatores ambientais quadráticos... Sendo assim, essa linha de comando
                                                                                                                                   # cria o criterio de corte (threshold) do R2 com a db-RDA usando todos os fatores
                                                                                                                                   # ambientais quadraticos apos a eliminaçao com base nas correlaçoes.                            

envP.q2<-subset(envP.q, select= - c (`SB^2`, `Ca^2`, `Mg^2`, `asp1^2`))   # aqui é a criaçao da matriz segundo o conjunto de fatores ambientais que serão submetidos
                                                                          # a seleçao de fatores ambientais a seguir (segundo uma das opçoes de combinação de fatores
                                                                          # ambientais quadráticos com base nas correlaçoes). 

forward.sel(spp.pcoa$points, envP.q2,  adjR2thresh= 0.2, alpha = 0.1)     # aqui está a seleção de fatores ambientais quadraticos com o threshold criado na linha 
                                                                          # de código acima e um alpha menos restritivo, de 0.10. 

RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `N^2` + `AF^2` + `K^2` , data= as.data.frame(envP.q)))$adj.r.squared    # aqui é calculado o R2 final segundo o conjunto de fatores
                                                                                                                    # ambientais quadráticos selecionado

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


envP_FIM<-cbind(envP.l2[,1], envP.l2[,8:9], envP.l2[,2], envP.l2[,6:7], envP.l2[,4], envP.q2[,1], envP.q2[,9], envP.q2[,15], envP.q2[,4])   # aqui está a construçao, com auxílio da função
                                                                                                                                            # "cbind" da matriz de fatores ambientais (quadráticos
                                                                                                                                            # e lineares) final a ser usada nos modelos finais e
                                                                                                                                            # e na partiçao de variancia.
                                                                                                                                      
colnames(envP_FIM)<-c("Elev", "P", "H_Al", "K", "N", "Sol", "Silt", "Elev^2", "N^2", "AF^2", "K^2")   # Aqui estou nomeando as colunas da matriz final de fatores ambientais construida na linha de comando acima

cor_envP_FIM<-cor(envP_FIM, envP_FIM)     # conferindo se por acaso existe ainda alguma correlaçao entre os fatores ambientais lineares e os quadráticos, visto que estes foram 
ifelse(abs(cor_envP_FIM)>0.7, "OUT", "-") # selecionados separadamente. Nada está correlacionado. Beleza!



spa<-read.csv("space.csv", header=T, sep=";", dec = ",", row.names=1)   # Aqui entramos com os dados espaciais brutos (latitude e longitude, em graus decimais, que é como está na planilha base)

spa_cent<-scale(spa, scale=F, center=TRUE)                              # Antes de começar a trabalhar os dados espaciais/geograficos deve-se centra-los, como é feito nesta linha de comando         

spa.rda<-capscale(spp.t ~ ., distance="horn", data = as.data.frame(spa_cent))  # Aqui eu crio um objeto que é uma db_rda (com funçao capscale) para testar se o modelo é significativo

anova(spa.rda)    # aqui eu testo o modelo de db-RDA criado na linha de comando acima. Se não for significativo, pode parar por aqui e não precisa incluir os fatores espaciais lineares
                  # (também conhecidos como latitude e longitude...). Como aqui é significativo, damos continuidade a rotina. 

spa_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = as.data.frame(spa_cent)))$adj.r.squared   # O procedimento é o mesmo que o feito com os fatores ambientais. Sendo assim, nessa linha 
                                                                                              # de comando é criado um valor threshold de referencia para a seleçao de fatores. 

forward.sel(spp.pcoa$points, spa_cent, adjR2thresh = spa_R2a, alpha = 0.05)   # E aqui é executada a seleçao de variaveis para o conjunto de fatores espaciais lineares (X e Y = longitude e latitude)
                                                                              # Sim, é para rodar para o conjunto de apeans 2 fatores sim... pelo menos e isso que recomenda o livro 
                                                                              # do Borcard et al. (2011 - numerical ecology with R). Acabaram sendo selecionados os dois fatores (X e Y) mesmo.


spa.dist<-dist(spa_cent)                                    # Depois de lidar com os fatores espaciais lineares, damos inicio nessa linha de codigo aqui a
                                                            # modelar os fatores espaciais refinados (distance based moran eigenvector maps) e a respectiva rotina para
                                                            # selecionar aqueles que ficarão nas análises finais. Essa linha de comando aqui cria uma matriz de distancia
                                                            # entre as unidade amostrais com base nos fatores espaciais lineares.

library(PCNM)   # aqui chamamos/abrimos o pacote que permite rodar a pcnm (que é a analise que determina os fatores espaciais lineares)

spa.PCNM<-PCNM(spa.dist, dbMEM=T, moran=T, all=TRUE)    # depois de aberto o pacote, criamos um objeto com o resultado da pcnm rodada nessa linha de comando usando a matriz
                                                        # de distancia criada acima.

select<-which(spa.PCNM$Moran_I$Positive)   # seleciona os fatores espaciais refinados (dbMeMs) positivos
PCNMs<-as.data.frame(spa.PCNM$vectors)[,select] # cria um objeto/matriz com fatores espaciais refinados que sera submetida a seguir ao procedimento de seleçao de variaveis


PCNM.rda<-capscale(spp.t ~ ., distance="horn", data = PCNMs)   # como feito para os fatores espaciais lineares, primeiro rodamos uma dbRDA para checar se o modelo é significativo 

anova.cca(PCNM.rda)    # aqui checamos se o modelo é significativo. Como é, damos continuidade a rotina.

PCNM_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = PCNMs))$adj.r.squared      # aqui é a criaçao do threshold usado na seleçao de fatores espaciais refinados. mesma lógica que descrita mais acima

forward.sel(spp.pcoa$points, PCNMs, adjR2thresh=PCNM_R2a, alpha=0.05)   # seleçao dos fatores espaciais refinados que serão mantidos nas análises

PCNM_FIM<-as.data.frame(cbind(PCNMs[,1]))   # criaçao do objeto/matriz (que nesse caso é uma coluna só mesmo..) com o fator espacial refinado 

colnames(PCNM_FIM)<-c("PCNM_1")   # nomeaçao da coluna da matriz criada na linha de codigo acima

rownames(PCNM_FIM)<-rownames(spa_cent) # nomeaçao das linhas da matriz criada duas linhas de comando acima  



varpart(spp.pcoa$points, envP_FIM, spa_cent, PCNM_FIM)   # aqui é a linha de comando com a partiçao de variancia final usando a matriz de spp segundo a ordenaçao PCoA,
                                                         # e cada um dos conjuntos de dados preditores criados ao longo das linhas de código acima. 


fts.exp<-as.data.frame(cbind(envP_FIM, spa_cent, PCNM_FIM))   # aqui crio um objeto que é uma matriz juntando todos os fatores preditores. Faço isso para rodar as linhas de 
                                                              # código a seguir, onde é testada significancia do modelo inteiro, dos termos, dos eixos e das 
                                                              # partições (step 1a, step 1b, step 1c, ..., step 2c - da uma olhada na tabela do arquivo word que passei na primeira
                                                              # leva de resultados que se chama partição_RDA que tem os "steps"). 


anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                           + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct")                 # essas duas linhas de código (poderia ser uma so, pus em duas para 
                                                                                                         # facilitar a visualizaçao) testam a significancia do modelo completo

summary(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`       # Essas duas linhas de código são para pedir os resultados gerais da db_RDA
                         + X + Y + PCNM_1, distance= "horn", fts.exp))                                   # Eu uso para olhar o escores dos preditores nos eixos da db-RDA

anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`     # Linhas de código para pedir a significancia dos termos que compoem o modelo
                          + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="terms")

anova.cca(capscale(spp.t ~ Silt +  `Alt^2` + AG + `SB^2` + SOL + V + `Arg^2` + `K^2`                     # linhas de codigo para pedir a significancia dos eixos que compoem o modelo
                          + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="axis")





#STEP 1a
pure_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                        + Condition (X + Y + PCNM_1), distance= "horn", fts.exp)                              # linhas de código para rodar o modelo referente a contribuiçao pura do ambiente

anova.cca(pure_E, step= 10000, perm= 10000, model="reduced")                                                  # teste do modelo que define a contribuiçao pura do ambiente

RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                + Condition (X + Y + PCNM_1), fts.exp))                                       # linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do ambiente


#STEP 1b
pure_S.l<-capscale(spp.t ~ X + Y + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                              + PCNM_1), distance= "horn", fts.exp)                           # linhas de código para rodar o modelo referente a contribuiçao pura do espaço linear

anova.cca(pure_S.l, step= 10000, perm= 10000, model= "reduced")                                               # teste do modelo que define a contribuiçao pura do espaço linear

RsquareAdj(rda(spp.pcoa$points ~ X + Y + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                                    + PCNM_1), fts.exp))                                      # linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao pura do espaço linear



#STEP 1c
pure_S.p<-capscale(spp.t ~  PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                                + X + Y), distance= "horn", fts.exp)                          # linhas de código para rodar o modelo referente a contribuiçao pura do espaço refinado

anova.cca(pure_S.p, step= 100000, perm= 100000, model="reduced")                                              # teste do modelo que define a contribuiçao pura do espaço refinado

RsquareAdj(rda(spp.pcoa$points ~  PCNM_1 + Condition (Silt +  `Alt^2` + ARG + `SB^2` + SOL + V + `Arg^2` + `K^2` + X + Y), fts.exp))                                     


#STEP 2a
S.l_j_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                         + X + Y + Condition (PCNM_1), distance= "horn", fts.exp)         # linhas de código para rodar o modelo referente a contribuiçao do ambiente junto do espaço linear

anova.cca(S.l_j_E, step=1000, perm= 1000, model="reduced")                                # teste do modelo que define a contribuiçao do ambiente junto do espaço linear

RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                + X + Y + Condition (PCNM_1), fts.exp))                   # linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do ambiente junto do espaço linear



#STEP 2b
S.l_j_S.p<-capscale(spp.t ~   X + Y + PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + 
                                                         `Elev^2` + `N^2` + `AF^2` + `K^2`), distance= "horn", fts.exp)        # linhas de código para rodar o modelo referente a contribuiçao do espaço linear junto ao espaço refinado

anova.cca(S.l_j_S.p, step= 1000, perm= 1000, model= "reduced")                                # teste do modelo que define a contribuiçao do espaço linear junto ao espaço refinado

RsquareAdj(rda(spp.pcoa$points ~   X + Y + PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + 
                                                              `Elev^2` + `N^2` + `AF^2` + `K^2`), fts.exp))                    # linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do espaço linear junto ao espaço refinado



#STEP 2c
E_j_S.p<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + 
                          PCNM_1 + Condition (X + Y), distance= "horn", fts.exp)                                # linhas de código para rodar o modelo referente a contribuiçao do ambiente junto ao espaço refinado

anova.cca(E_j_S.p, step= 1000, perm= 1000, model="reduced")                                # teste do modelo que define a contribuiçao do ambiente junto ao espaço refinado

RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ 
                                 PCNM_1 + Condition (X + Y), fts.exp))                     # linhas de codigo para definir o R2 ajustado do modelo referente a contribuiçao do ambiente junto ao espaço refinado






palm.rda<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + X + Y + PCNM_1, distance= "horn", fts.exp)   # linha de codigo para rodar dbRDA que servirá de base para as figuras

plot(palm.rda, scaling=1, main="Triplot RDA spe.hor ~ allfcts - scaling 1 - wa scores")     # basicamente, daqui para baixo o que eu fiz foi elaborar diferentes plots 
                                                                                            # com os dados da db-RDA gerada na linha de comando acima. Da uma olhada no livro do
                                                                                            # Legendre & Legendre 2012 (NUmerical Ecology pg 638 ou próximo) e no Borcard et al. 2011 (Numerical Ecology
                                                                                            # with R) para entender as diferenças entre os tipos de scaling e os tipos de escores 
                                                                                            # mostrados no plot... ok?

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



