#obtaining data from motus towers using R API

library(motus)
library(lubridate)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)


Sys.setenv(TZ = "UTC") #sets the working environment to be in UTC to match motus

#tower at Western UWO - SSC (#11645)
#receiver CTT-7292E3C5E490

# =================================================
# downloading the data
#==================================================
proj.num <- "CTT-7292E3C5E490"
sql_motus <- tagme(projRecv = proj.num, 
                   new = TRUE,
                   skipNodes = TRUE, 
                   skipDeprecated = TRUE)

sql_motus <- tagme(proj.num, 
                   skipNodes = TRUE, 
                   skipDeprecated = TRUE)
metadata(sql_motus) #needed to get the metadata

## testing Sauble Beach Tower

sql_motus.sauble <- tagme(projRecv = "CTT-V3023D10185C",
                          new = TRUE) 
metadata(sql_motus.sauble)

# ================================================
# getting detection data
# ================================================

library(DBI)
library(RSQLite)

dbListTables(sql_motus) #to view tables
dbListFields(sql_motus, "alltags") #to view fields of tables

alltags.tbl <- tbl(sql_motus, "alltags") #convert to table
deps.tbl <- tbl(sql_motus, "tagDeps")
sp.tbl <- tbl(sql_motus, "species")
tagProp.tbl <- tbl(sql_motus, "tagProps")
recvDep.tbl <- tbl(sql_motus, "recvDeps")

alltags.df <- alltags.tbl %>% 
  collect() %>% 
  as.data.frame

write.csv (alltags.df, "UWO-SSCAllTags.csv")

deps.df <- deps.tbl %>% 
  collect() %>% 
  as.data.frame

write.csv (deps.df, "UWO-SSCDeps.csv")

recvDeps.df <- recvDep.tbl %>% 
  collect() %>% 
  as.data.frame

#====================================================
# adding additional bird metadata
#====================================================

#adding sex and age where possible
alltags.deps <- alltags.df %>%
  left_join(
    deps.df %>% 
      select(tagID, speciesID, markerNumber, sex, age),
    by = c(
      "motusTagID" = "tagID",        
      "speciesID" = "speciesID",    
      "markerNumber" = "markerNumber"  
    )
  )

##adding in receiver data
alltags.recvs <- alltags.deps %>% 
#Drop the empty columns from the main data so there's no naming conflict
  select(-recvSiteName, -recvDeployLat, -recvDeployLon) %>% 
  #now join
  left_join(
    recvDeps.df %>% 
      select(serno, name, stationLat, stationLon) %>% 
      rename(
        recvSiteName  = name,
        recvDeployLat = stationLat,
        recvDeployLon = stationLon
      ), 
    by = c("recv" = "serno")
  )
  
