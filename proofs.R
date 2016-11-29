#Feito para ~1 e ~rx LDM (ver fun??o survLDM3)-> atualizar para PLDM-KMW-IPCW
#Colocar class para o resultado -> utilizar em plot.surv
#library(doParallel)

f <- survCOND(survCS(time1, event1, Stime, event)~age, x = 365, y = 730,
              z.value=48, data = colonCS, conf = TRUE)

fit <- survCOND(survCS(time1, event1, Stime, event)~age, x = 365, y = 365*1:8, z.value=48, data = colonCS)

fit <- survCOND(survCS(time1, event1, Stime, event)~nodes, x = 365, y = 730, z.value=1, data = colonCS)

survCOND(survCS(time1, event1, Stime, event)~1, x = 365, y = 730, data = colonCS, method = "LDM", conf=TRUE)





fit <- survCOND(survCS(time1, event1, Stime, event)~rx, x = 365, data = colonCS, method = "LDM", conf = TRUE)




fit1 <- survCOND(survCS(time1, event1, Stime, event)~1, x = 365, data = colonCS, method = "KMW", conf = TRUE)


f <- survCOND(survCS(time1, event1, Stime, event)~factor(sex), x = 365, data = colonCS, method = "LDM")

fit <- survCOND(survCS(time1, event1, Stime, event)~factor(sex), x = 365, data = colonCS, method = "KMW", conf=TRUE)



f <-survCOND(survCS(time1, event1, Stime, event)~rx, x = 365, y = 730, data = colonCS, method = "LDM")

fit <- survCOND(survCS(time1, event1, Stime, event)~rx, x = 365, data = colonCS, method = "LDM", conf=T)

fit2 <- survCOND(survCS(time1, event1, Stime, event)~rx, x = 365, y = seq(365,2999,222), data = colonCS, method = "LDM", conf=T)

matplot(x=fit$y, y=fit$estimates, type="s")

fit <- survCOND(survCS(time1, event1, Stime, event)~1, x = c(365), data = colonCS, method = "LDM", conf=T)



plot.surv(fit, type="s", conftype="s")


#FALTA INTRODUZIR ESTA FUN??O
survCOND2(survCS(time1, event1, Stime, event)~1, x = 365, y = 730, data = colonCS, method = "LDM", presmooth = TRUE)

survCOND2(survCS(time1, event1, Stime, event)~rx, x = 365, y = 730, data = colonCS, method = "LDM", presmooth = TRUE)

#TESTAR PARA OS DADOS DE BLADDER

#TESTAR COM OS COMANDO DO PAPER

#FUN??O PLOT

#Manual do package

#
