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
dbListFields(sql_motus, "recvs") #to view fields of tables

#creating tbls
alltags.tbl <- tbl(sql_motus, "alltags") #convert to table
tags.tbl <- tbl(sql_motus, "tags")
deps.tbl <- tbl(sql_motus, "tagDeps")
sp.tbl <- tbl(sql_motus, "species")
tagProp.tbl <- tbl(sql_motus, "tagProps")
recvDep.tbl <- tbl(sql_motus, "recvDeps")
antDep.tbl <-tbl(sql_motus, "antDeps")
runs.tbl <- tbl (sql_motus, "allruns")
ambigs.tbl <- tbl(sql_motus, "allambigs")

#converting tbls into data frames
alltags.df <- alltags.tbl %>% 
  collect() %>% 
  as.data.frame

alltags.df2 <- alltags.df %>%
mutate(time = as_datetime(ts))

write.csv (alltags.df, "UWO-SSCAllTags.csv")

ambigs.df <- ambigs.tbl %>% 
  collect () %>% 
  as.data.frame

tags.df <- tags.tbl %>% 
  collect() %>% 
  as.data.frame

deps.df <- deps.tbl %>% 
  collect() %>% 
  as.data.frame

write.csv (deps.df, "UWO-SSCDeps.csv")

recvDeps.df <- recvDep.tbl %>% 
  collect() %>% 
  as.data.frame

antDep.df <- antDep.tbl %>% 
  select(deployID, port, antennaType, bearing, heightMeters) %>%
  collect() %>%
  as.data.frame()

species.df <- sp.tbl %>% 
  collect () %>% 
  as.data.frame

runs.df <- runs.tbl %>% 
  collect() %>% 
  as.data.frame

# checking if there are multiple deployments for any dags
alltags.df %>%
  select(motusTagID, tagDeployID) %>%
  filter(!(is.na(tagDeployID))) %>% # remove NA tagDeployIDs
  distinct() %>%
  group_by(motusTagID) %>%
  mutate(n = n()) %>%
  filter(n > 1)

#there are multiple tag IDs so creating motudTagDepID to have unique IDs

alltags.df <- alltags.df %>%
  mutate(motusTagDepID = paste(motusTagID, tagDeployID, sep = "."))
#same but for metadata
deps.df <- deps.df %>%
  mutate(motusTagDepID = paste(tagID, deployID, sep = "."))

#fix time in receiver deployments

recvDeps.df <- recvDeps.df %>% 
  mutate(timeStart = as_datetime(tsStart),
         timeEnd = as_datetime(tsEnd))

#adding in receiver and antenna metadata

recvDeps.df <- recvDeps.df %>% 
  select(deployID, receiverType, deviceID, name, latitude, longitude, 
         isMobile, timeStart, timeEnd, projectID, elevation) 

stationDeps.df <- left_join(recvDeps.df, antDep.df, by = "deployID")

# adding additional bird metadata
head(alltags.df) #fill columns in this one
alltagsspecies.df <- species.df %>% 
  select(speciesID = id, 
         speciesEN = english, 
         speciesSci = scientific, 
         speciesGroup = group) %>%  #to get species name
  right_join(alltags.df)

alltags.deps <- deps.df %>%
  select(motusTagDepID, sex, age ) %>% 
  right_join(alltagsspecies.df) %>% 
  mutate(time = as_datetime(ts))

runs.df <- runs.df %>% 
  mutate(RunBegin = as_datetime(tsBegin),
         RunEnd = as_datetime(tsEnd))

write.csv(runs.df, "UWO-SSCruns.csv")

alltags.runs <- runs.df %>% 
  select(runID, RunBegin, RunEnd) %>% 
  right_join(alltags.deps) %>% 
  distinct()

write.csv(alltags.deps, "UWO-SSCDetectionsRaw.csv")

