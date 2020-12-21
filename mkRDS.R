
setwd("D:/2020-2/Zarathu/Lithium")

library(readxl)
library(data.table)

## Sheet name
name.sheet <- excel_sheets("201221 anonymized data.xlsx")

## Warning when reading Sheet 1
list.data <- parallel::mclapply(name.sheet, function(i){
  if (i == "clinical data"){
    data.table(read_excel("201221 anonymized data.xlsx", sheet = i, col_types = "text"))
  } else{
    data.table(read_excel("201221 anonymized data.xlsx", sheet = i))
  }
})

## Remove column with all NA
list.data[[4]] <- list.data[[4]][, 1:11]

## list name: same to excel sheet
names(list.data) <- name.sheet

## Save to RDS
saveRDS(list.data, "lithium.RDS")
