library(tidyverse)

file <- "experiments/sachs/sachs_data_pooled.csv"
newfile <- "experiments/sachs/sachs_data.csv"

.col_types <- cols('AKT inh'=col_integer(), 'G0076'=col_integer(),
                   'LY294002'=col_integer(), 'PMA/beta2CAMP + noAlphaCD3/28'=col_integer(), 
                   'Psitectorigenin'=col_integer(), 'U0126'=col_integer(),
                   'experiment'=col_integer())

sachs_data_pooled <- read_csv(file, col_types = .col_types)
combined_intervention <- sachs_data_pooled[['PMA/beta2CAMP + noAlphaCD3/28']]

PMA <- sapply(combined_intervention, function (x) (x == 1)*1)
sachs_data_pooled[['PMA']] <- PMA

b2CAMP <- sapply(combined_intervention, function (x) (x == 2)*1)
sachs_data_pooled[['b2CAMP']] <- b2CAMP

drop <- c('PMA/beta2CAMP + noAlphaCD3/28')
sachs_data_pooled <- sachs_data_pooled[, !(names(sachs_data_pooled) %in% drop)]

write.csv(sachs_data_pooled,file=newfile,row.names=FALSE)
