
#setwd("~/ShinyApps/jihyunbaek/lithium")
setwd("/home/heejooko/ShinyApps/lithium")

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

## ICD data
ICD.data <- data.table(read_excel("ICD data.xlsx",skip=2))
saveRDS(ICD.data,"ICD_data.RDS")

# W21
W210226 <- data.table(read_excel("W210226-1(환자번호 지움) (1).xlsx", skip = 2))
W210216 <- data.table(read_excel("W210216-4(환자번호 지움) (1).xlsx", skip = 2))

saveRDS(W210226, "W210226.RDS");saveRDS(W210216, "W210216.RDS")