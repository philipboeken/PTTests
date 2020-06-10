library(usethis)

col_types <- readr::cols('experiment' = readr::col_integer(), 'CD3/28' = readr::col_integer(),
                         'AKT inh' = readr::col_integer(), 'G0076' = readr::col_integer(),
                         'Psitectorigenin' = readr::col_integer(), 'U0126' = readr::col_integer(), 
                         'LY294002' = readr::col_integer(), 'PMA + noCD3/28' = readr::col_integer(), 
                         'b2CAMP + noCD3/28' = readr::col_integer())

sachs_data <- readr::read_csv("data-raw/sachs_data.csv", col_types = col_types)

usethis::use_data(sachs_data, compress = "xz", overwrite = T)
