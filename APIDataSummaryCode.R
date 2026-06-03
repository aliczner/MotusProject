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

## testing Sauble Beach Tower

sql_motus.sauble <- tagme(projRecv = "CTT-V3023D10185C",
                          new = TRUE) #started at 3:22 pm


# ================================================
# getting detection data
# ================================================

library(DBI)
library(RSQLite)

dbListTables(sql_motus) #to view tables
dbListFields(sql_motus, "alltags") #to view fields of tables
dbListFields(sql_motus, "species")

alltags.tbl <- tbl(sql_motus, "alltags") #convert to table
species.table <- tbl(sql_motus, "species")

alltags.df <- alltags.tbl %>% 
  collect() %>% 
  as.data.frame

species.df <- species.table %>% 
  collect() %>% 
  as.data.frame



write.csv(alltags.df, "alltagstest.csv")
