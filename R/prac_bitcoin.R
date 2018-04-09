rm(list=ls())

# load --------------------------------------------------------------------
library(bigrquery)

project <- "prac-dataviz"

sql <- "SELECT *
FROM `bigquery-public-data.bitcoin_blockchain.transactions`"

todo_copies <- query_exec(sql, project = project, useLegacySql = FALSE)
