out_root <-  "/home/joel/sim/qmap/CORDEX/output/EUR-44/final_results/"
freq='day'
outpath='/home/joel/sim/qmap/'
domain='EUR-44'
indir=paste0(outpath,"CORDEX/output/", domain,"/")
vars=c('tas', 'pr', 'ps', 'hurs', 'rsds','rlds', 'uas', 'vas')

#===============================================================================
# finally construct TS files outside of loop
#===============================================================================
if(freq=='day'){strsplitPar <- 23}
if(freq=='3hr'){strsplitPar <- 31}
files = list.files(path=indir, pattern="_ll.nc$", recursive=T, full.name=T)

base = substring(files, 1,nchar(files)-strsplitPar)
base2 = unique(base)

#===============================================================================
# checksums timeperiods
#===============================================================================
correct_length <- 12

NtimePeriods <- c()

for (i in 1:length(base2)){
  NtimePeriods[i] <- length(which(base==base2[i]))
  }

mylist <- strsplit(base2, "/")
filenames <- sapply(mylist, "[", 18)
df <- (data.frame(filenames, NtimePeriods))

#===============================================================================
# checksums variables
#===============================================================================
mylist2 <- strsplit(filenames, paste0("_",domain))
correct_nvars=length(vars)
basenames <- (sapply(mylist2, "[", c(2)))
basenames_unique <- unique(sapply(mylist2, "[", c(2)))

Nvars <- c()

for (i in 1:length(basenames_unique)){
  Nvars[i] <- length(which(basenames==basenames_unique[i]))
  }

df_vars <- (data.frame(basenames_unique, correct_nvars -Nvars))

#===============================================================================
# merge datasets by time
#===============================================================================

for (mybase in base2){
#mybase2 = unlist(strsplit(mybase,'/'))[[length(unlist(strsplit(mybase,'/')))]]
#base3 = list.files(path=indir, pattern=mybase2, recursive=T, full.name=T)
system(paste0("cdo -b F64 -f nc2 mergetime ",mybase,"* ", mybase,"_TS.nc"))
print(paste0("concatenated file: ", mybase,"_TS.nc"))
}

# analyse
final_ts = list.files(path=indir, pattern="_TS.nc$", recursive=T, full.name=T)
mylist <- strsplit(final_ts, '/')
myvars <- unique(sapply(mylist, "[", 18))
myinst <-  unique(sapply(mylist, "[", 11))
mymodel <-  unique(sapply(mylist, "[", 12))
myexper <-  unique(sapply(mylist, "[", 13))





#===============================================================================
# merge all TS files by variable
#===============================================================================

final_ts = list.files(path=indir, pattern="_TS.nc$", recursive=T, full.name=T)
mylist <- strsplit(final_ts, "/")
filenames <- sapply(mylist, "[", 18)
mylist2 <- strsplit(filenames, paste0("_",domain))
basenames <- (sapply(mylist2, "[", c(2)))
basenames_unique <- unique(sapply(mylist2, "[", c(2)))

for (i in 1:length(basenames_unique)){
  var_index <- which(basenames==basenames_unique[i])
  system(paste(c("cdo merge " ,final_ts[var_index], paste0(out_root,basenames_unique[i],"_",length(var_index),"VARS.nc") ), collapse=' '))
}
