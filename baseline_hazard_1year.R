library(survival)
basehaz_noneprom <- vector()
for(i in 1:40){
  basehaz_noneprom[i] <- basehaz(non_ePROM_md[["analyses"]][[i]])[179,1]
}
mean(basehaz_noneprom)

basehaz_eprom1 <- vector()
for(i in 1:40){
  basehaz_eprom1[i] <- basehaz(eprom_md1[["analyses"]][[i]])[179,1]
}
mean(basehaz_eprom1)
basehaz_eprom2 <- vector()
for(i in 1:40){
  basehaz_eprom2[i] <- basehaz(eprom_md2[["analyses"]][[i]])[179,1]
}
mean(basehaz_eprom2)

basehaz_eprom3 <- vector()
for(i in 1:40){
  basehaz_eprom3[i] <- basehaz(eprom_md3[["analyses"]][[i]])[179,1]
}
mean(basehaz_eprom3)
