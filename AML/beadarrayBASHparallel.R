library("synapseClient")
library("illuminaHumanv4.db")
library(beadarray)
library(parallel)

synapseLogin()
tmpAML = synGet("syn5521819") 


mclapply(as.list(list.files("RadichAMLmRNAData/", full.names = T)), function(x){system(paste("gunzip", x))})
files = as.list(list.files("RadichAMLmRNAData/", full.names = T))

write.csv("",file = "temp.txt", row.names= F)

loadNbash = function(x)
{ 
  sectionName = gsub("RadichAMLmRNAData//|\\.txt$","",x); 
  bld         = readIllumina(dir="RadichAMLmRNAData", sectionNames = sectionName); 
  bashd       = BASH(bld, array = 1); 
  
  write.csv(sectionName,file = "temp.txt", append = T)
  
  return(bashd)
}

bashList = mclapply(files,loadNbash ,mc.cores = 35)

