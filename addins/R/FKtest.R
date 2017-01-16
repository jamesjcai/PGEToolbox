setwd("C:/biodata/GEO/evQTL_Mapping/model3popgeno2/r_fun");
DATA<-read.table("input.txt");
gen<-DATA[,1];
y<-DATA[,2];


f<-fligner.test(y,gen);
write.table(f$p.value, file="output.txt", sep = "\t", col.names = FALSE, row.names = FALSE, qmethod = "double",append=FALSE)


