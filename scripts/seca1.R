
setwd("C:/Users/Sara/Desktop/R/R/Doctorado/First")

#setwd("C:\\Users\\user\\Documents\\PROJETOS_MANUSCRITOS\\Palmeiras_RJ\\Analise_atual")  # Indica a pasta onde est�o os arquivos das analises neste script 


spp<-read.csv("densi.csv", header=T, sep=";", dec = ",", row.names=1)   # chamada dos dados referentes as spp da comunidade

library(vegan)     # inicia��o/chamada do pacote vegan

spp.t<-decostand(spp, "total")    # padroniza��o dos dados da comunidade segundo a soma dos valores registrados na unidade amostral

spp.bray<-vegdist(spp.t, "horn")  # calculo da matriz de distancias (usei a medida de morisita-horn/tamb�m rodei com bray-curtis mas obtive uma % de explica�ao um pouco menor)
                                  # a mesma medida de distancia deve ser usada ao avaliar os dados da seca e da chuva ok?? N�o sei se sabe, mas para ver como uma dada fun��o
                                  # funciona � so voc� colocar uma interroga��o antes da fun�ao, selecionar a interroga�ao junto com a fun�ao e apertar crtl+R ok?

spp.pcoa<-cmdscale(spp.bray, k=nrow(spp)-1, eig=TRUE, add=TRUE)   # Cria�ao de matriz resposta da comunidade segundo a ordena�ao PCoA. Para ver o que significa cada argumento, 
                                                                  # da uma olhada no help da fun�ao cmdscale. AO usarmos a matriz da ordena�ao para representar a comunidade
                                                                  # estamos rodadndo uma distance-based RDA e n�o a RDA mais cl�ssica ok? O principal beneficio dessa abordagem, no caso 
                                                                  # desses nossos dados, � que podemos usar vari�veis (spp) com baixa ocorr�ncia. A RDA cl�ssica depende do registro                                                                   # 
                                                                  # de pelo menos cinco pontos para as vari�veis... A db-RDA tamb�m permite usar qualquer medida de distancia para 
                                                                  # modelar os dados, enquanto a RDA classica � "fechada" nesse aspecto (salvo algumas transforma��es 
                                                                  # (ver Legendre & Gallagher 2001 ou Gallagher & Legendre 2001 - n�o tenho certeza da ordem dos autores...).


env.1<-read.csv("topo.csv", header=T, sep=";", dec = ",", row.names=1)   # chamada/emtrada dos dados ambientais (parte referente a topografia - que serve tanto para chuva quanto para a seca)
trf.asp<-cbind(sin(env.1[,2]), cos(env.1[,2]))   # transforma�ao do fator aspecto, dentro dos dados de topografia, pelo seu seno e coseno (1 fator viram dois mesmo). Nos trabalhos 
                                                 # da Punchi-Manage et al. (2013, 2014) eles citam esta transforma�ao para fatores que sejam de dados circulares, como no caso dos 360 graus
                                                 # usados para registrar o fator aspecto...

colnames(trf.asp)<-c("asp1", "asp2")    # nomea��o  das colunas do obejto "trf.asp"

env.2<-read.csv("solo_sc.csv", header= T, sep= ";", dec= ",", row.names=1) # chamada/entrada dos dados ambientais (parte referente a quimica do solo no periodo de seca)

env<-cbind(env.2, env.1, trf.asp)   # cria�ao de objeto contendo todas os fatores ambientais usando a fun�ao "cbind"

library(stats)    # inicia��o/chamada do pacote stats, que contem a fun�ao "poly"

envP<-cbind(poly(env$elev,2),    # cria�ao de objeto/matriz com todos os fatores ambientais e seus respectivos valores quadr�ticos (ao quadrado)
            poly(env$MO,2),      # segundo a fun�ao "poly". Essa fun�ao j� transforma e padroniza os fatores ambientais, ent�o n�o � preciso qualquer
            poly(env$pH,2),      # outra trasnforma��o desses dados. A fun�ao poly calcula os valores quadr�ticos e garante que eles sejam ortogonais entre
            poly(env$K,2),       # si. Se voc� calcular direto o valor quadr�tico, os fatores ter�o uma correla��o de 1...Por isso essa fun��o � t�o legal, ela 
            poly(env$Ca,2),      # padroniza os dados e garante que n�o sejam lineares.
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

colnames(envP)<-c("elev", "elev^2", "MO", "MO^2", "pH", "pH^2", "K", "K^2", "Ca", "Ca^2",     # nomea��o das colunas do objeto criado com aux�lio da fun�ao "poly". Se n�o nomear
                  "Mg", "Mg^2", "SB", "SB^2", "CTC", "CTC^2", "AG", "AG^2", "SLT", "SLT^2",   # ter� problemas depois porque as colunas ficam com "nomes" 1, 2, 1, 2, ...
                  "ARG", "ARG^2", "N", "N^2", "sol", "sol^2", "P", "P^2", "H_Al", "H_Al^2", "V", "V^2"
                 , "asp1", "asp1^2", "asp2", "asp2^2", "AF", "AF^2", "slop", "slop^2")  

rownames(envP)<-rownames(env)   # nomea�ao das linhas do objeto "envP" de acordo com as linhas do objeto "env"

envP.l<-subset(envP, select = - c (`elev^2`, `MO^2`, `pH^2`, `K^2`, `Ca^2`, `Mg^2`, `SB^2`, `CTC^2`, `AG^2`, `SLT^2`,       # Dissocia�ao das informa�oes/fatores ambientais lineares daqueles
                                   `ARG^2`, `N^2`, `sol^2`, `P^2`, `H_Al^2`, `V^2`, `asp1^2`, `asp2^2`, `AF^2`, `slop^2`))  # quadr�ticos. Aqui usei a fun�ao subset para selecionar uma parte     
                                                                                                                            # das colunas da matriz de dados ambientais "envP", criando o obejto
                                                                                                                            # "envP.l". 

envP.q<-subset(envP, select = - c (elev, MO, pH, K, Ca, Mg, SB, CTC, AG, SLT, ARG, N, sol, P, H_Al, V, asp1, asp2, AF, slop))   # Dissocia�ao  das informa�oes/fatores ambientais quadraticos daqueles
                                                                                                                                # lineares. A raz�o pela qual separei essas informa�oes foi para 
                                                                                                                                # executar a sele�ao de fatores ambienatias com base em correlacao. 
                                                                                                                                # ou seja, para diminuir a colinearidade. Fiz isso porque sao muitos fatores
                                                                                                                                # e ficava melhor para visualizar as correla��es e para rodar a fun��o "forward.sel",
                                                                                                                                # que seleciona os fatores a serem usados na an�lise final (depois de lidar
                                                                                                                                # com a multicolinearidade...). Repara que se eu n�o separasse seria uma matriz
                                                                                                                                # com 40 fatores.                                                                                                                    # 


cor_envP.l<-cor(envP.l, envP.l)      # Cria��o de objeto com correla�ao entre fatores ambientais

ifelse(abs(cor_envP.l)>0.7, "OUT", "-")  # solicita�ao de matriz na tela de opera��es do R, respondendo quais correla�oes sao maiores do 0.7

# Com base na tabela de correla��es maiores ou iguais a 0.7, bolei diferentes combina��es de fatores ambientais, segundo a elimina�ao de
# diferentes fatores... Essas combina�oes foram ent�o usadas na sele�ao de vari�veis ("forward.sel") usando a matriz de spp criada
# com a ordena�ao PCoA como resposta. Essas foram as op�oes que bolei. Existem outras, mas acho que cobri as principais para 
# esse conjunto de dados. Voc� precisar� atualizar isso para os dados de chuva ok?

# OP�AO 1: MANTER MO + SB + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OP�AO 2: MANTER elev + SB + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop
# OP�AO 3: MANTER MO + SB + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop   
# OP�AO 4: MANTER MO + pH + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OP�AO 5: MANTER MO + pH + K + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OP�AO 6: MANTER MO + pH + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop     
# OP�AO 7: MANTER MO + Ca + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 8: MANTER MO + Ca + K + CTC + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 9: MANTER MO + Ca + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 10: MANTER MO + H_Al + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 11: MANTER MO + H_Al + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 12: MANTER elev + SB + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop   
# OP�AO 13: MANTER elev + pH + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OP�AO 14: MANTER elev + pH + K + Mg + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop     
# OP�AO 15: MANTER elev + pH + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop     
# OP�AO 16: MANTER elev + Ca + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 17: MANTER elev + Ca + K + CTC + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 18: MANTER elev + Ca + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 19: MANTER elev + H_Al + K + CTC + SLT + ARG + N + sol + P + asp1 + asp2 + AF + slop  
# OP�AO 20: MANTER elev + H_Al + K + Mg + SLT + AG + N + sol + P + asp1 + asp2 + AF + slop


calib.adj.l<-RsquareAdj(rda(spp.pcoa$points ~  elev + H_Al + K + Mg + SLT + AG + N +                                             # Cria�ao de objeto contendo valor de referencia
                                               sol + P + asp1 + asp2 + AF + slop, data= as.data.frame(envP.l)))$adj.r.squared    # como threshold para execu�ao da sele�ao de fatores.
                                                                                                                                 # A sele�ao de fatores com a fun�ao "forward.sel" deve
                                                                                                                                 # ser usada com dois crit�rios de "stop". Um � o valor 
                                                                                                                                 # de alpha (0.05 geralmente, mas eu usei um menos restritivo, de 0.10) 
                                                                                                                                 # e o outro � o R2 da rda 
                                                                                                                                 # usando todos os fatores ambientais (depois da elimina�ao)
                                                                                                                                 # segundo est� aqui nessa linha (dupla) de c�digo.

envP.l2<-subset(envP.l, select= - c (MO, pH, CTC, SB, AG, Ca, V))    # Aqui estou criando um objeto com o conjunto de fatores ambientais que ser�o submetidos
                                                                     # a sele�ao de fatores ambientais a seguir (segundo uma das op�oes de combina��o de fatores
                                                                     # com base nas correla�oes). Novamente usando a fun�ao "subset" para criar um subconjunto de fatores
                                                                     # daqueles existentes na matriz de fatores ambientais lineares. Isso precisa ser feito porque 
                                                                     # a fun��o "forward.sel" n�o roda como a regressao, onde especificamos os fatores que queremos
                                                                     # usar. Ela pede que forne�a um objeto contendo todos os fatores... e por isso ele precisa ser criado.

library(packfor)    # chamada/abertura do pacote "packfor"

forward.sel(spp.pcoa$points, envP.l2,  adjR2thresh=calib.adj.l , alpha = 0.1)   # execu�ao da sele�ao de variaveis para o conjunto de dados ambientais lineares
                                                                                # segundo a op��o de combina�ao de fatores escolhida. 

RsquareAdj(rda(spp.pcoa$points ~ elev + P + H_Al + K + N + sol + SLT, data= as.data.frame(envP.l)))$adj.r.squared   # determina��o do R2 segundo o conjunto de fatores ambientais lineares
                                                                                                                    # selecionados na sele��o. 

# abaixo est�o os conjuntos de fatores ambientais lineares selecionados segundo cada op�ao de combina�ao de fatores montada 
# com base nas correla��es. Cada conjunto de fatores selecionados � seguido de seu respectivo R2 ao ser usado na dbRDA
# Esse R2 serve de referencia para o conjunto final de fatores selecionado e que sera usado no modelo final. Eu escolhi
# o conjunto de fatores que propiciou o maior R2. 

# OP�AO 1: FICAM MO + sol + slop + SB + P + ARG + N    R2 = 15,79
# OP�AO 2: FICAM elev + P + SB + sol + N + SLT    R2 = 16,32
# OP�AO 3: FICAM MO + sol + slop + SB + P + ARG + N    R2 = 15,79
# OP�AO 4: FICAM MO + sol + slop + SLT + P + pH + K + N   R2 = 19.01    
# OP�AO 5: FICAM MO + sol + slop + SLT + P + pH + K + N   R2 = 19.01
# OP�AO 6: FICAM MO + sol + slop + SLT + P + pH + K + N   R2 = 19.01 
# OP�AO 7: FICAM MO + sol + slop + Ca + P + K + SLT + N   R2 = 19.62 
# OP�AO 8: FICAM MO + sol + slop + Ca + P + K + SLT + N   R2 = 19.62 
# OP�AO 9: FICAM MO + sol + slop + Ca + P + K + ARG + SLT  R2 = 19.51
# OP�AO 10: FICAM MO + sol + slop + SLT + P + K + H_Al + N  R2 = 19.55
# OP�AO 11: FICAM MO + sol + slop + SLT + P + K + H_Al + N  R2 = 19.55
# OP�AO 12: FICAM elev + P + SB + sol + N + SLT  R2 = 16,32
# OP�AO 13: FICAM elev + P + sol + pH + K + N + SLT  R2 = 17.91
# OP�AO 14: FICAM elev + P + sol + pH + K + N + SLT  R2 = 17.91
# OP�AO 15: FICAM elev + P + sol + pH + K + N + SLT  R2 = 17.91
# OP�AO 16: FICAM elev + P + Ca + K + sol + N + SLT  R2 = 19.66
# OP�AO 17: FICAM elev + P + Ca + K + sol + N + SLT  R2 = 19.66
# OP�AO 18: FICAM elev + P + Ca + K + sol + N + SLT  R2 = 19.66
# OP�AO 19: FICAM elev + P + H_Al + K + N + sol + SLT  R2 = 19.93
# OP�AO 20: FICAM elev + P + H_Al + K + N + sol + SLT  R2 = 19.93




cor_envP.q<-cor(envP.q, envP.q)
ifelse(abs(cor_envP.q)>0.7, "OUT", "-")
# OP�AO 1: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`   
# OP�AO 2: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`              
# OP�AO 3: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`         
# OP�AO 4: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp1^2` + `AF^2` + `slop^2`
# OP�AO 5: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OP�AO 6: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `V^2` + `asp1^2` + `AF^2` + `slop^2`
# OP�AO 7: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OP�AO 8: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `SB^2` + `Mg^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OP�AO 9: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `Ca^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `asp2^2` + `AF^2` + `slop^2`
# OP�AO 10: MANTER `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + `V^2` + `asp2^2` + `AF^2` + `slop^2`


calib.adj.q<-RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `MO^2` + `pH^2` + `K^2` + `CTC^2` + `AG^2` + `SLT^2` + `ARG^2` + `N^2` + `sol^2` + `P^2` + `H_Al^2` + 
                                              `V^2` + `asp2^2` + `AF^2` + `slop^2`, data= as.data.frame(envP.q)))$adj.r.squared    # Aqui eu aplico de novo o procedimento para a sele�ao de variaveis ambientais
                                                                                                                                   # por�m agora � os fatores ambientais quadr�ticos... Sendo assim, essa linha de comando
                                                                                                                                   # cria o criterio de corte (threshold) do R2 com a db-RDA usando todos os fatores
                                                                                                                                   # ambientais quadraticos apos a elimina�ao com base nas correla�oes.                            

envP.q2<-subset(envP.q, select= - c (`SB^2`, `Ca^2`, `Mg^2`, `asp1^2`))   # aqui � a cria�ao da matriz segundo o conjunto de fatores ambientais que ser�o submetidos
                                                                          # a sele�ao de fatores ambientais a seguir (segundo uma das op�oes de combina��o de fatores
                                                                          # ambientais quadr�ticos com base nas correla�oes). 

forward.sel(spp.pcoa$points, envP.q2,  adjR2thresh= 0.2, alpha = 0.1)     # aqui est� a sele��o de fatores ambientais quadraticos com o threshold criado na linha 
                                                                          # de c�digo acima e um alpha menos restritivo, de 0.10. 

RsquareAdj(rda(spp.pcoa$points ~ `elev^2` + `N^2` + `AF^2` + `K^2` , data= as.data.frame(envP.q)))$adj.r.squared    # aqui � calculado o R2 final segundo o conjunto de fatores
                                                                                                                    # ambientais quadr�ticos selecionado

# abaixo segue o conjunto de fatores ambientais quadraticos que foi selecionado em cada uma das op�oes
# de combina�ao dos fatores segundo as correla��es e os respectivos R2 das db-RDA usado esses fatores.
# Para os fatores ambientais quadr�ticos, o conjunto final de fatores selecionado foi sempre o mesmo. 

# OP�AO 1: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 2: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 3: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 4: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 5: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 6: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 7: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 8: FICAM `elev^2` + `N^2` + `AF^2`   R2= 14.03
# OP�AO 9: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03
# OP�AO 10: FICAM `elev^2` + `N^2` + `AF^2` + `K^2`   R2= 14.03


envP_FIM<-cbind(envP.l2[,1], envP.l2[,8:9], envP.l2[,2], envP.l2[,6:7], envP.l2[,4], envP.q2[,1], envP.q2[,9], envP.q2[,15], envP.q2[,4])   # aqui est� a constru�ao, com aux�lio da fun��o
                                                                                                                                            # "cbind" da matriz de fatores ambientais (quadr�ticos
                                                                                                                                            # e lineares) final a ser usada nos modelos finais e
                                                                                                                                            # e na parti�ao de variancia.
                                                                                                                                      
colnames(envP_FIM)<-c("Elev", "P", "H_Al", "K", "N", "Sol", "Silt", "Elev^2", "N^2", "AF^2", "K^2")   # Aqui estou nomeando as colunas da matriz final de fatores ambientais construida na linha de comando acima

cor_envP_FIM<-cor(envP_FIM, envP_FIM)     # conferindo se por acaso existe ainda alguma correla�ao entre os fatores ambientais lineares e os quadr�ticos, visto que estes foram 
ifelse(abs(cor_envP_FIM)>0.7, "OUT", "-") # selecionados separadamente. Nada est� correlacionado. Beleza!



spa<-read.csv("space.csv", header=T, sep=";", dec = ",", row.names=1)   # Aqui entramos com os dados espaciais brutos (latitude e longitude, em graus decimais, que � como est� na planilha base)

spa_cent<-scale(spa, scale=F, center=TRUE)                              # Antes de come�ar a trabalhar os dados espaciais/geograficos deve-se centra-los, como � feito nesta linha de comando         

spa.rda<-capscale(spp.t ~ ., distance="horn", data = as.data.frame(spa_cent))  # Aqui eu crio um objeto que � uma db_rda (com fun�ao capscale) para testar se o modelo � significativo

anova(spa.rda)    # aqui eu testo o modelo de db-RDA criado na linha de comando acima. Se n�o for significativo, pode parar por aqui e n�o precisa incluir os fatores espaciais lineares
                  # (tamb�m conhecidos como latitude e longitude...). Como aqui � significativo, damos continuidade a rotina. 

spa_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = as.data.frame(spa_cent)))$adj.r.squared   # O procedimento � o mesmo que o feito com os fatores ambientais. Sendo assim, nessa linha 
                                                                                              # de comando � criado um valor threshold de referencia para a sele�ao de fatores. 

forward.sel(spp.pcoa$points, spa_cent, adjR2thresh = spa_R2a, alpha = 0.05)   # E aqui � executada a sele�ao de variaveis para o conjunto de fatores espaciais lineares (X e Y = longitude e latitude)
                                                                              # Sim, � para rodar para o conjunto de apeans 2 fatores sim... pelo menos e isso que recomenda o livro 
                                                                              # do Borcard et al. (2011 - numerical ecology with R). Acabaram sendo selecionados os dois fatores (X e Y) mesmo.


spa.dist<-dist(spa_cent)                                    # Depois de lidar com os fatores espaciais lineares, damos inicio nessa linha de codigo aqui a
                                                            # modelar os fatores espaciais refinados (distance based moran eigenvector maps) e a respectiva rotina para
                                                            # selecionar aqueles que ficar�o nas an�lises finais. Essa linha de comando aqui cria uma matriz de distancia
                                                            # entre as unidade amostrais com base nos fatores espaciais lineares.


library(PCNM)   # aqui chamamos/abrimos o pacote que permite rodar a pcnm (que � a analise que determina os fatores espaciais lineares)

spa.PCNM<-PCNM(spa.dist, dbMEM=T, moran=T, all=TRUE)    # depois de aberto o pacote, criamos um objeto com o resultado da pcnm rodada nessa linha de comando usando a matriz
                                                        # de distancia criada acima.

select<-which(spa.PCNM$Moran_I$Positive)   # seleciona os fatores espaciais refinados (dbMeMs) positivos
PCNMs<-as.data.frame(spa.PCNM$vectors)[,select] # cria um objeto/matriz com fatores espaciais refinados que sera submetida a seguir ao procedimento de sele�ao de variaveis


PCNM.rda<-capscale(spp.t ~ ., distance="horn", data = PCNMs)   # como feito para os fatores espaciais lineares, primeiro rodamos uma dbRDA para checar se o modelo � significativo 

anova.cca(PCNM.rda)    # aqui checamos se o modelo � significativo. Como �, damos continuidade a rotina.

PCNM_R2a<-RsquareAdj(rda(spp.pcoa$points ~ ., data = PCNMs))$adj.r.squared      # aqui � a cria�ao do threshold usado na sele�ao de fatores espaciais refinados. mesma l�gica que descrita mais acima

forward.sel(spp.pcoa$points, PCNMs, adjR2thresh=PCNM_R2a, alpha=0.05)   # sele�ao dos fatores espaciais refinados que ser�o mantidos nas an�lises

PCNM_FIM<-as.data.frame(cbind(PCNMs[,1]))   # cria�ao do objeto/matriz (que nesse caso � uma coluna s� mesmo..) com o fator espacial refinado 

colnames(PCNM_FIM)<-c("PCNM_1")   # nomea�ao da coluna da matriz criada na linha de codigo acima

rownames(PCNM_FIM)<-rownames(spa_cent) # nomea�ao das linhas da matriz criada duas linhas de comando acima  



varpart(spp.pcoa$points, envP_FIM, spa_cent, PCNM_FIM)   # aqui � a linha de comando com a parti�ao de variancia final usando a matriz de spp segundo a ordena�ao PCoA,
                                                         # e cada um dos conjuntos de dados preditores criados ao longo das linhas de c�digo acima. 


fts.exp<-as.data.frame(cbind(envP_FIM, spa_cent, PCNM_FIM))   # aqui crio um objeto que � uma matriz juntando todos os fatores preditores. Fa�o isso para rodar as linhas de 
                                                              # c�digo a seguir, onde � testada significancia do modelo inteiro, dos termos, dos eixos e das 
                                                              # parti��es (step 1a, step 1b, step 1c, ..., step 2c - da uma olhada na tabela do arquivo word que passei na primeira
                                                              # leva de resultados que se chama parti��o_RDA que tem os "steps"). 


anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                           + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct")                 # essas duas linhas de c�digo (poderia ser uma so, pus em duas para 
                                                                                                         # facilitar a visualiza�ao) testam a significancia do modelo completo

summary(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`       # Essas duas linhas de c�digo s�o para pedir os resultados gerais da db_RDA
                         + X + Y + PCNM_1, distance= "horn", fts.exp))                                   # Eu uso para olhar o escores dos preditores nos eixos da db-RDA

anova.cca(capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`     # Linhas de c�digo para pedir a significancia dos termos que compoem o modelo
                          + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct" , by="terms")

anova.cca(capscale(spp.t ~ Silt +  `Alt^2` + AG + `SB^2` + SOL + V + `Arg^2` + `K^2`                     # linhas de codigo para pedir a significancia dos eixos que compoem o modelo
                          + X + Y + PCNM_1, distance= "horn", fts.exp), model="direct", by="axis")





#STEP 1a
pure_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                        + Condition (X + Y + PCNM_1), distance= "horn", fts.exp)                              # linhas de c�digo para rodar o modelo referente a contribui�ao pura do ambiente

anova.cca(pure_E, step= 10000, perm= 10000, model="reduced")                                                  # teste do modelo que define a contribui�ao pura do ambiente

RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                + Condition (X + Y + PCNM_1), fts.exp))                                       # linhas de codigo para definir o R2 ajustado do modelo referente a contribui�ao pura do ambiente


#STEP 1b
pure_S.l<-capscale(spp.t ~ X + Y + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                              + PCNM_1), distance= "horn", fts.exp)                           # linhas de c�digo para rodar o modelo referente a contribui�ao pura do espa�o linear

anova.cca(pure_S.l, step= 10000, perm= 10000, model= "reduced")                                               # teste do modelo que define a contribui�ao pura do espa�o linear

RsquareAdj(rda(spp.pcoa$points ~ X + Y + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                                    + PCNM_1), fts.exp))                                      # linhas de codigo para definir o R2 ajustado do modelo referente a contribui�ao pura do espa�o linear



#STEP 1c
pure_S.p<-capscale(spp.t ~  PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + X + Y), distance= "horn", fts.exp)                          # linhas de c�digo para rodar o modelo referente a contribui�ao pura do espa�o refinado

anova.cca(pure_S.p, step= 100000, perm= 100000, model="reduced")                                              # teste do modelo que define a contribui�ao pura do espa�o refinado

RsquareAdj(rda(spp.pcoa$points ~  PCNM_1 + Condition (Silt +  `Alt^2` + AG + `SB^2` + SOL + V + `Arg^2` + `K^2`  
                                                      + X + Y), fts.exp))                                     # linhas de codigo para definir o R2 ajustado do modelo referente a contribui�ao pura do espa�o refinado
 


#STEP 2a
S.l_j_E<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                         + X + Y + Condition (PCNM_1), distance= "horn", fts.exp)         # linhas de c�digo para rodar o modelo referente a contribui�ao do ambiente junto do espa�o linear

anova.cca(S.l_j_E, step=1000, perm= 1000, model="reduced")                                # teste do modelo que define a contribui�ao do ambiente junto do espa�o linear

RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`  
                                + X + Y + Condition (PCNM_1), fts.exp))                   # linhas de codigo para definir o R2 ajustado do modelo referente a contribui�ao do ambiente junto do espa�o linear



#STEP 2b
S.l_j_S.p<-capscale(spp.t ~   X + Y + PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + 
                                                         `Elev^2` + `N^2` + `AF^2` + `K^2`), distance= "horn", fts.exp)        # linhas de c�digo para rodar o modelo referente a contribui�ao do espa�o linear junto ao espa�o refinado

anova.cca(S.l_j_S.p, step= 1000, perm= 1000, model= "reduced")                                # teste do modelo que define a contribui�ao do espa�o linear junto ao espa�o refinado

RsquareAdj(rda(spp.pcoa$points ~   X + Y + PCNM_1 + Condition (Elev +  P + H_Al + K + N + Sol + Silt + 
                                                              `Elev^2` + `N^2` + `AF^2` + `K^2`), fts.exp))                    # linhas de codigo para definir o R2 ajustado do modelo referente a contribui�ao do espa�o linear junto ao espa�o refinado



#STEP 2c
E_j_S.p<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + 
                          PCNM_1 + Condition (X + Y), distance= "horn", fts.exp)                                # linhas de c�digo para rodar o modelo referente a contribui�ao do ambiente junto ao espa�o refinado

anova.cca(E_j_S.p, step= 1000, perm= 1000, model="reduced")                                # teste do modelo que define a contribui�ao do ambiente junto ao espa�o refinado

RsquareAdj(rda(spp.pcoa$points ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2`+ 
                                 PCNM_1 + Condition (X + Y), fts.exp))                     # linhas de codigo para definir o R2 ajustado do modelo referente a contribui�ao do ambiente junto ao espa�o refinado






palm.rda<-capscale(spp.t ~ Elev +  P + H_Al + K + N + Sol + Silt + `Elev^2` + `N^2` + `AF^2` + `K^2` + X + Y + PCNM_1, distance= "horn", fts.exp)   # linha de codigo para rodar dbRDA que servir� de base para as figuras

plot(palm.rda, scaling=1, main="Triplot RDA spe.hor ~ allfcts - scaling 1 - wa scores")     # basicamente, daqui para baixo o que eu fiz foi elaborar diferentes plots 
                                                                                            # com os dados da db-RDA gerada na linha de comando acima. Da uma olhada no livro do
                                                                                            # Legendre & Legendre 2012 (NUmerical Ecology pg 638 ou pr�ximo) e no Borcard et al. 2011 (Numerical Ecology
                                                                                            # with R) para entender as diferen�as entre os tipos de scaling e os tipos de escores 
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



