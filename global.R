library(data.table);library(magrittr);library(DT);library(jstable);library(dplyr);library(stats)

# setwd("~/ShinyApps/jihyunbaek/lithium")
setwd("/home/heejooko/ShinyApps/lithium")
lithium <- readRDS("lithium.RDS")


# CKD-EPI----------------------------------------
CKDEPI<-function(scr,age,sex){
  if(sex=="F"){ k<-0.7; alpha<-(-0.329); const<-1.018 }
  else{ k<-0.9; alpha<-(-0.411); const<-1 }
  GFR <- 141 * min(scr/k,1)^alpha * max(scr/k,1)^(-1.209) * 0.993^age * const
  return(GFR)
}

df <- lithium$MEDI[, c("NO","처방일","처방명","일수","횟수")] 
names(df) <- c("NO","date","drug","day","times")
df[, drug := ifelse(drug == "Lithium carbonate 300mg", "Lithium", "Valproate")]
df <- unique(df)
df <- df[order(-day), .(maxday = max(day, na.rm = T), maxnotqd = day[which(times != 1)[1]]), by=c("NO","date","drug")]
df <- df[, .(maxday, qd = ifelse(is.na(maxnotqd), maxday , maxday - maxnotqd)), by=c("NO","date","drug")]
df <- df[, .(totDay = sum(maxday, na.rm = T), qd = sum(qd, na.rm = T)/sum(maxday, na.rm = T)),by = c("NO", "drug")]

df.long <- dcast(df, NO ~ drug, value.var = c("totDay", "qd"))

#left_join 위해서 NO의 class를 맞춰주기----------------------------------------
lithium$`clinical data`$NO <- lithium$`clinical data`$NO %>% as.numeric() %>% as.character()
df.long$NO <- df.long$NO %>% as.character()
lithium$`clinical data` <- merge(lithium$`clinical data`, df.long, by = "NO")
## Dx group
lithium$`clinical data`[, group_bipolar_schizoaffective_other := factor(ifelse(grepl("Bipolar|bipolar", lithium$`clinical data`$주상병명), "Bipolar disorder",
                                                                               ifelse(grepl("Schizoaffective|schizoaffective", lithium$`clinical data`$주상병명), "Schizoaffective disorder", "vOthers")))]


# Data inclusion----------------------------------------
a <- lithium$`clinical data`[,
                             .(NO,성별,생년월일,totDay_Lithium,totDay_Valproate,qd_Lithium,qd_Valproate, 
                               HTN = factor(as.integer(!is.na(`고혈압 여부`))), DM = factor(as.integer(!is.na(`당뇨 여부`))), group_bipolar_schizoaffective_other)]

N_profile<-cbind("전체",NA,nrow(a),NA,NA)

## NO. list: TDM both
NO.tdmboth <- lithium$`renal function & TDM`[, NO := as.character(NO)][`세부검사명` %in% c("Lithium", "Valproic Acid")][, c("NO", "세부검사명")][, unique(`세부검사명`), by = "NO"][, names(which(table(NO) == 2))]

a <- a[xor(is.na(totDay_Lithium),is.na(totDay_Valproate)),,][!(NO %in% c("2250", NO.tdmboth))]

N_profile<-rbind(N_profile,cbind("Li+Valp combination",as.integer(N_profile[nrow(N_profile),3])-nrow(a),nrow(a),NA,NA))

a[, drug := factor(ifelse(is.na(totDay_Lithium), 0, 1))]

ICD_data <- readRDS("ICD_data.RDS")
setnames(ICD_data,c("개인정보동의여부","정렬순서"),c("Privacy Consent","NO"))
ICD_data$NO<-ICD_data$NO %>% as.character()
ICD_data<-ICD_data[`Privacy Consent`=="Y",,]

a<-merge(a,ICD_data[,.(NO),],by="NO")
N_profile<-rbind(N_profile,cbind("개인정보사용미동의",as.integer(N_profile[nrow(N_profile),3])-nrow(a),nrow(a),a[drug==0,.N,],a[drug==1,.N,]))

a<-a[(totDay_Lithium>180 | totDay_Valproate>180),,]
N_profile<-rbind(N_profile,cbind("총처방일수 180일 초과",as.integer(N_profile[nrow(N_profile),3])-nrow(a),nrow(a),a[drug==0,.N,],a[drug==1,.N,]))


## Date age----------------------------------------
df <- lithium$MEDI[, NO := as.character(NO)][,.SD,]
setnames(df,c("처방일","함량단위투여량","일수"),c("date","dose","day"))

data.main <- a %>% 
  merge(df[,.(firstPrescriptionDay=min(date, na.rm = T)), by = "NO"], by = "NO",all.x = T) %>% 
  merge(df[,.(lastPrescriptionDay=max(date, na.rm = T)), by = "NO"], by = "NO",all.x = T) %>% 
  merge(df[, .(avgDose_1day = sum(dose * day)/sum(day)), by = "NO"], by = "NO",all.x = T)

data.main[, Age := floor(as.numeric(as.Date(firstPrescriptionDay) - as.Date(`생년월일`))/365.25)]


## NEW HTN/DM
HTN_DM <- apply(ICD_data, 1, function(x){
  HTN <- substr(x, 1, 3) %in% c(paste0("I", 10:16), 401:405)
  HTN_which <- which(HTN == T)[1]
  HTN_yn <- as.integer(!is.na(HTN_which))
  HTN_date <- ifelse(HTN_yn == 0, NA, x[HTN_which - 1])
  
  DM <- (substr(x, 1, 3) %in% paste0("E", c("08", "09", 10:13))) | (substr(x, 1, 4) %in% paste0("250.", 0:9)) 
  DM_which <- which(DM == T)[1]
  DM_yn <- as.integer(!is.na(DM_which))
  DM_date <- ifelse(DM_yn == 0, NA, x[DM_which - 1])
  
  return(c(HTN2 = HTN_yn, HTN_date = HTN_date, DM2 = DM_yn, DM_date = DM_date))
}) %>% t %>% data.table

HTN_DM$NO <- ICD_data$NO

HTN_DM_info <- merge(data.main[, c("NO", "firstPrescriptionDay")], HTN_DM, by = "NO")
HTN_DM_info[, `:=`(HTN = factor(as.integer(HTN2 == 1 & as.Date(HTN_date) <= as.Date(firstPrescriptionDay))),
                   DM = factor(as.integer(DM2 == 1 & as.Date(DM_date) <= as.Date(firstPrescriptionDay))))]
## merge
data.main <- merge(data.main[, -c("HTN", "DM")], HTN_DM_info[, c("NO", "HTN", "DM")], by = "NO")


## F code --------------------------------------

W210226 <- readRDS("W210226.RDS")
setnames(W210226,c("개인정보동의여부","정렬순서","진단코드"),c("Privacy Consent","NO","dcode"))
W210226$NO<-W210226$NO %>% as.character()
W210226<-W210226[`Privacy Consent`=="Y",,]
W210226<-merge(W210226,data.main[,.(NO),],by="NO")

W210226<-W210226[dcode %like% "F2" | dcode %like% "F30" | dcode %like% "F31" | dcode %like% "F32" | dcode %like% "F33",.SD,]


data.main <- merge(data.main,W210226[,.(NO),],by="NO")
N_profile<-rbind(N_profile,cbind("주진단코드 F20~F33",as.integer(N_profile[nrow(N_profile),3])-data.main[,.N,],data.main[,.N,],data.main[drug==0,.N,],data.main[drug==1,.N,]))

W210226[, Fcode := factor(ifelse(dcode %like% "F2", "F20-29", ifelse(dcode %like% "F30|F31", "F30-F31", ifelse(dcode %like% "F32|F33", "F32-F33", "Others"))))]
W210226[, Fcode2 := substr(dcode, 1, 3)]

data.main <- data.main[Age>=18,,]
N_profile<-rbind(N_profile,cbind("첫처방일기준 만 18세 이상",as.integer(N_profile[nrow(N_profile),3])-data.main[,.N,],data.main[,.N,],data.main[drug==0,.N,],data.main[drug==1,.N,]))


#W210226[,schizo:=factor(as.integer((dcode %like% "F2"))),]
#W210226[,mood:=factor(as.integer((dcode %like% "F3"))),]
#W210226[,bipolar:=factor(as.integer(((dcode %like% "F30")|(dcode %like% "F31")))),]
#W210226[,depressive:=factor(as.integer(((dcode %like% "F32")|(dcode %like% "F33")))),]



data.main<-merge(data.main,W210226[,c("NO","Fcode", "Fcode2"),],by="NO")
# LithiumToxicity----------------------------------------

df <- lithium$`renal function & TDM`[, NO := as.character(NO)][] %>% 
  merge(data.main[, .(NO, firstPrescriptionDay, lastPrescriptionDay)], by = "NO", all.x = T)
setnames(df,c("세부검사명","결과","시행일시"),c("test","result","testdate"))

data.main <- data.main %>% 
  merge(df[test=="Lithium" & as.numeric(result) > 1.0 & (testdate - firstPrescriptionDay >= 0) & (lastPrescriptionDay - testdate  >= 0),  .(LithiumToxicity1.0 = .N), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Lithium" & as.numeric(result) > 0.8 & (testdate - firstPrescriptionDay >= 0) & (lastPrescriptionDay - testdate  >= 0), .(LithiumToxicity0.8 = .N), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Lithium" & as.numeric(result) > 1.2 & (testdate - firstPrescriptionDay >= 0) & (lastPrescriptionDay - testdate  >= 0), .(LithiumToxicity1.2 = .N), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Lithium" & (testdate - firstPrescriptionDay >= 0), .(avgTDM_Lithium = mean(as.numeric(result), na.rm = T)), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Valproic Acid" & (testdate - firstPrescriptionDay >= 0), .(avgTDM_Valproate = mean(as.numeric(result), na.rm = T)), by="NO"], by="NO", all.x = T)

# for (v in c("LithiumToxicity1.0", "LithiumToxicity0.8", "LithiumToxicity1.2")){
#   data.main[[v]] <- factor(ifelse(is.na(data.main[[v]]), 0, data.main[[v]]))
# }

for (v in c("LithiumToxicity1.0", "LithiumToxicity0.8", "LithiumToxicity1.2")){
  data.main[[v]] <- factor(ifelse(is.na(data.main[[v]]), 0, 1))
}

df<-lithium$`renal function & TDM`
df$NO <- as.character(df$NO)
df <- merge(df, lithium$`clinical data`[,.(NO,`성별`,`생년월일`),], by="NO", mult=all)
df$`결과`<- as.numeric(df$`결과`) 
df$시행일시 <- as.Date(df$시행일시); df$생년월일 <- as.Date(df$생년월일)
setnames(df,c("세부검사명","시행일시","생년월일","결과","성별"),c("test", "testDate", "birthDate", "result", "sex"))
df[,age:=as.numeric(testDate-birthDate)/365.25]
df[,eGFR:=ifelse(test=="Creatinine",CKDEPI(result,age,sex),NA),by=seq_len(nrow(df))]


## data for figure 1----------------------------------------
data.f1 <- df[!is.na(eGFR),.(NO,testDate,eGFR)]
setnames(data.f1, "testDate", "date")


## Main data----------------------------------------

data.main <- merge(data.main, df[eGFR < 60, .(eGFRbelow60Date = min(testDate)), by = "NO"], all.x = TRUE) %>% 
  merge(df[test == "Creatinine", .(testNum = .N), by="NO"], by="NO", all.x=TRUE) %>% 
  .[!is.na(testNum) & (is.na(eGFRbelow60Date) | as.Date(firstPrescriptionDay) < as.Date(eGFRbelow60Date))]

data.main <- merge(data.main, df[eGFR < 30, .(eGFRbelow30Date = min(testDate)), by = "NO"], all.x = TRUE)
data.main <- merge(data.main, df[eGFR < 15, .(eGFRbelow15Date = min(testDate)), by = "NO"], all.x = TRUE)


# baseline에서 6개월 / 30일 이상 데이터

data.f1<-merge(data.f1,data.f1[,.(baselineDate=min(date)),by="NO"],by="NO")
data.f1[,dateDiff:=(date-baselineDate),]
data.f1[dateDiff>=180,firstPostBaselineDate:=min(date),by="NO"]
data.f1[,dateDiff2:=(date-firstPostBaselineDate)]
data.f1[dateDiff2>=30,secondPostBaselineDate:=min(date),by="NO"]

data.main<-merge(data.main,data.f1[,.(PostBaselineDate1=max(firstPostBaselineDate,na.rm=T),PostBaselineDate2=max(secondPostBaselineDate,na.rm=T)),by="NO"],by="NO")
data.main<-data.main[!is.na(as.POSIXlt(PostBaselineDate1)) & !is.na(as.POSIXlt(PostBaselineDate2))]
N_profile<-rbind(N_profile,cbind("baseline + 2개 이상의 post-baseline eGFR (post-baseline 두 eGFR 사이 최소 30일 이상 기간)",as.integer(N_profile[nrow(N_profile),3])-data.main[,.N,],data.main[,.N,],data.main[drug==0,.N,],data.main[drug==1,.N,]))


# 맨 뒤<60 && 30일 차이나는 그 다음 맨 뒤<60 --> eGFRbelow60

data.f1[,lastTestDate:=max(date,na.rm=T),by="NO"]
data.f1[,dateDiff_fromlast:=(lastTestDate-date),]

data.f1[dateDiff_fromlast>=30,secondLastTestDate:=max(date),by="NO"]

# data.f1[date==secondLastTestDate,.N,by="NO"][(N)!=1,]  같은 날짜의 lastTestDate이나 secondLastTestDate에 측정값 여러개인 것이 있어서 mean 취함
# data.f1[NO=="7934"]

library(zoo)

df.lastEGFRs<-merge(data.f1[date==lastTestDate,.(lasteGFR=mean(eGFR,na.rm=T), lastTestDate),by="NO"],data.f1[date==secondLastTestDate,.(secondLasteGFR=mean(eGFR,na.rm=T), secondLastTestDate),by="NO"],by="NO",all=T)
df.lastEGFRs<-df.lastEGFRs[!duplicated(df.lastEGFRs$NO)]

data.main<-merge(data.main,df.lastEGFRs,by="NO")

data.main[,eGFRbelow60:=factor(ifelse(lasteGFR<60 & secondLasteGFR<60,1,0)),]
data.main[,eGFRbelow60Date:=ifelse(eGFRbelow60==1,as.Date(lastTestDate),NA)]
data.main[, duration60 := ifelse(eGFRbelow60==0,as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay), as.Date(eGFRbelow60Date) - as.Date(firstPrescriptionDay))]
data.main[, year_FU := duration60/365.25]

data.main[,eGFRbelow30:=factor(ifelse(lasteGFR<30 & secondLasteGFR<30,1,0)),]
data.main[,eGFRbelow30Date:=ifelse(eGFRbelow30==1,as.Date(lastTestDate),NA)]
data.main[, duration30 := ifelse(eGFRbelow30==0,as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay), as.Date(eGFRbelow30Date) - as.Date(firstPrescriptionDay))]
data.main[, year_FU_eGFRbelow30 := duration30/365.25]

data.main[,eGFRbelow15:=factor(ifelse(lasteGFR<15 & secondLasteGFR<15,1,0)),]
data.main[,eGFRbelow15Date:=ifelse(eGFRbelow15==1,as.Date(lastTestDate),NA)]
data.main[, duration15 := ifelse(eGFRbelow15==0,as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay), as.Date(eGFRbelow15Date) - as.Date(firstPrescriptionDay))]
data.main[, year_FU_eGFRbelow15 := duration15/365.25]

data.f1<-data.f1[,-c("baselineDate","dateDiff","firstPostBaselineDate","dateDiff2","secondPostBaselineDate","lastTestDate","dateDiff_fromlast","secondLastTestDate")]

# data.main[, duration60 := ifelse(is.na(eGFRbelow60Date),as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay), as.Date(eGFRbelow60Date) - as.Date(firstPrescriptionDay))]
# data.main[, eGFRbelow60 := factor(as.integer(!is.na(eGFRbelow60Date)))]
# data.main[, year_FU := duration60/365.25]
# 
# data.main[, duration30 := ifelse(is.na(eGFRbelow30Date),as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay), as.Date(eGFRbelow30Date) - as.Date(firstPrescriptionDay))]
# data.main[, eGFRbelow30 := factor(as.integer(!is.na(eGFRbelow30Date)))]
# data.main[, year_FU_eGFRbelow30 := duration30/365.25]
# 
# data.main[, duration15 := ifelse(is.na(eGFRbelow15Date),as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay), as.Date(eGFRbelow15Date) - as.Date(firstPrescriptionDay))]
# data.main[, eGFRbelow15 := factor(as.integer(!is.na(eGFRbelow15Date)))]
# data.main[, year_FU_eGFRbelow15 := duration15/365.25]

data.main[, `:=`(totYear_Lithium = totDay_Lithium/365.25, totYear_Valproate = totDay_Valproate/365.25)]
# duration Full
data.main[, year_FU_full := as.numeric(as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay))/365.25]

setnames(data.main, "성별", "Sex")
data.main[, Sex := factor(Sex)]

data.main <- data.main[, .SD, .SDcols = -c("생년월일", "firstPrescriptionDay", "lastPrescriptionDay", "duration60", "duration30", "duration15", "totDay_Valproate", "totDay_Lithium","testNum")]

## Figure 1 data----------------------------------------

# 처방 정보

df <- lithium$MEDI[, c("NO","처방일","처방명","일수","횟수")]
names(df) <- c("NO","date","drug","day","times")
df[, drug := factor(ifelse(drug == "Lithium carbonate 300mg", 1, 0))]
df <- unique(df)[, `:=`(NO = as.character(NO), date = as.Date(date))][]
df <- df[, .(maxday = max(day, na.rm = T)), by=c("NO","date","drug")]

# 일단 합친 후 cumsum
data.f1 <- rbindlist(list(data.f1, df),use.names = TRUE, fill=TRUE)[order(NO,date)][, maxday:=ifelse(is.na(maxday),0,maxday)][]
data.f1[, cumulativePrescriptionDay := cumsum(maxday),by=.(NO)]

data.f1 <- data.f1[!is.na(eGFR), !c("maxday","drug")]
data.f1 <- merge(data.f1, data.main[,.(NO,drug),], by="NO")

## eGFR linear regression ---------------------------------------

data.f1[,cumulativePrescriptionYear:=cumulativePrescriptionDay/365.25,]
data.main<-merge(data.main,data.f1[,.(eGFRchange=coef(lm(eGFR~cumulativePrescriptionYear))[2]),by=NO])


## Base eGFR, GFR change----------------------------------------

data.GFRchange <- data.f1[cumulativePrescriptionDay<365.25,.(year0GFR=mean(eGFR,na.rm=T)),by="NO"] %>% 
  merge(.,data.f1[, .(base_eGFR = ifelse(any(as.integer(cumulativePrescriptionDay) == 0), eGFR[last(which(cumulativePrescriptionDay == 0))], eGFR[1])), by ="NO"], by= "NO", all = T) %>% ## base eGFR
  merge(.,data.f1[(365.25*3)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*4),.(year3GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*5)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*6),.(year5GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*7)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*8),.(year7GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*10)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*11),.(year10GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*12)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*13),.(year12GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*15)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*16),.(year15GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*20)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*21),.(year20GFR=mean(eGFR,na.rm=T)),by="NO"],all=T)

data.main<-merge(data.main,data.GFRchange,all=T)

data.main<-data.main[base_eGFR>=60,,]
N_profile<-rbind(N_profile,cbind("baseline eGFR<60",as.integer(N_profile[nrow(N_profile),3])-data.main[,.N,],data.main[,.N,],data.main[drug==0,.N,],data.main[drug==1,.N,]))

## ----------------------------------------

ICD_data<-merge(ICD_data,data.main[,.(NO,drug),],by="NO")

ICD_data<-ICD_data[,alldiagnosis:=Reduce(paste,.SD),.SDcols=grep("진단코드",colnames(ICD_data))][,c("NO","drug","alldiagnosis"),]
ICD_data<-ICD_data[!(alldiagnosis %like% "N0|N1" & !(alldiagnosis %like% "N09")),.SD,]
N_profile<-rbind(N_profile,cbind("ICD N00-N08 or N10-N19",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,],ICD_data[drug==0,.N,],ICD_data[drug==1,.N,]))

ICD_data<-ICD_data[!(alldiagnosis %like% "T86.1"),.SD,]
N_profile<-rbind(N_profile,cbind("ICD T86.1",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,],ICD_data[drug==0,.N,],ICD_data[drug==1,.N,]))

ICD_data<-ICD_data[!(alldiagnosis %like% "Z94.0"),.SD,]
N_profile<-rbind(N_profile,cbind("ICD Z94.0",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,],ICD_data[drug==0,.N,],ICD_data[drug==1,.N,]))

data.main <- merge(data.main,ICD_data[,.(NO),],by="NO")

## 복용년수별 n수 ----------------------------------------

Year_N<-data.frame(Year=0:26,
                   Lithium_N=sapply(0:26,function(x) data.main[totYear_Lithium>x,.N,]),
                   Valproate_N=sapply(0:26,function(x) data.main[totYear_Valproate>x,.N,]))


## 해당 연차에 eGFR<60 된 n수 ----------------------------------------

data.f1<-merge(data.f1,data.main[,.(NO),],all.y=TRUE)

data.f1<-merge(data.f1,data.main[,.(NO,base_eGFR),],by="NO",all.x=TRUE)
data.f1<-data.f1[!(cumulativePrescriptionDay==0 & eGFR!=base_eGFR),,]
data.f1<-data.f1[,-c("base_eGFR"),]

#data.f1<-data.f1[eGFR>=30,,]

colnames(N_profile)<-c("조건","제외","N","Valproate","Lithium")

eGFRbelow60ratio<-
  lapply(0:26,function(x){
    NthYear<-unique(data.f1[(365.25*x)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*(x+1)),.(NthYeareGFR=mean(eGFR,na.rm=T),drug),by="NO"])
    nth<-merge(NthYear[NthYeareGFR<60,.(below60=.N),by=drug],NthYear[,.N,by=drug],by="drug",all=TRUE)
    if(NthYear[drug==1,.N,]==0){  nth<-rbind(nth,data.table(drug=1,below60="NA",N="NA"))  }
    if(NthYear[drug==0,.N,]==0){  nth<-rbind(nth,data.table(drug=0,below60="NA",N="NA"))  }
    nth[,yn:=paste(below60,N,sep = "/"),]  
    nth<-data.table::transpose(nth[,4,])
    return(nth)}) %>% Reduce(rbind,.)

eGFRbelow60ratio<-rbind(
  merge(data.main[base_eGFR<60,.(aa=.N),by=drug],data.main[,.(bb=.N),by=drug],by="drug",all=TRUE) %>% 
    .[,.(cc=paste(aa,bb,sep="/")),by="drug"] %>% .[,.(cc),] %>% transpose,
  eGFRbelow60ratio)

eGFRbelow60ratio<-as.data.frame(eGFRbelow60ratio)
colnames(eGFRbelow60ratio)<-c("Valproate","Lithium")
rownames(eGFRbelow60ratio)<-c("baseline",unlist(lapply(0:26,function(x){paste0("Year ",x)})))

## eGFR<60 최초발생일 연차별 n수

findCumYear<-function(ID,eGFRbelow60Date){
  return(data.f1[NO==ID & date==eGFRbelow60Date,cumulativePrescriptionYear])
}

dt<-unique(data.main[eGFRbelow60==1,.(yy=findCumYear(NO,eGFRbelow60Date),drug),by=c("NO","eGFRbelow60Date")])

eGFRbelow60Years<-
  lapply(0:20,function(x){
    nth<-dt[x<=yy & yy<x+1,.N,by=drug]
    if(dt[x<=yy & yy<x+1 & drug==0,.N,]==0){  nth<-rbind(nth,data.table(drug=0, N=0))  }
    if(dt[x<=yy & yy<x+1 & drug==1,.N,]==0){  nth<-rbind(nth,data.table(drug=1, N=0))  }
    nth<-data.table::transpose(nth[order(drug)])[2]
  }) %>% Reduce(rbind,.)



eGFRbelow60Years <- as.data.frame(eGFRbelow60Years)

colnames(eGFRbelow60Years)<-c("Valproate","Lithium")
rownames(eGFRbelow60Years)<-unlist(lapply(0:20,function(x){paste0("Year ",x)}))

eGFRbelow60Years$Valproate <- as.integer(eGFRbelow60Years$Valproate)
eGFRbelow60Years$Lithium <- as.integer(eGFRbelow60Years$Lithium)

eGFRbelow60Years<-rbind(eGFRbelow60Years,
                        data.frame(row.names="Sum",Valproate=sum(eGFRbelow60Years$Valproate),Lithium=sum(eGFRbelow60Years$Lithium)))

## ----------------------------------------

data.main <- data.main[, -c("eGFRbelow60Date","eGFRbelow30Date","eGFRbelow15Date","lastTestDate","secondLastTestDate", "PostBaselineDate1", "PostBaselineDate2")]  ## NO 제외
data.main$eGFR_last_minus_base <- data.main$lasteGFR - data.main$base_eGFR

outlierNO<-data.f1[eGFR<60,NO] %>% unique
outlier.data.f1<-data.f1[NO %in% outlierNO,,]


label.main <- jstable::mk.lev(data.main)

label.main[variable == "eGFRbelow60", `:=`(var_label = "eGFR < 60", val_label = c("No", "Yes"))]
# label.main[variable == "eGFRbelow30", `:=`(var_label = "eGFR < 30", val_label = c("No", "Yes"))]
# label.main[variable == "eGFRbelow15", `:=`(var_label = "eGFR < 15", val_label = c("No", "Yes"))] ## eGFRbelow15 Yes 없어서 에러나서 지움
label.main[variable == "drug", `:=`(var_label = "Drug", val_label = c("Valproate", "Lithium"))]
label.main[variable == "DM", `:=`(var_label = "DM", val_label = c("No", "Yes"))]
label.main[variable == "HTN", `:=`(var_label = "HTN", val_label = c("No", "Yes"))]
label.main[variable == "LithiumToxicity1.0", `:=`(var_label = "Lithium > 1.0 발생 여부", val_label = c("No", "Yes"))]
label.main[variable == "LithiumToxicity1.2", `:=`(var_label = "Lithium > 1.2 발생 여부", val_label = c("No", "Yes"))]
label.main[variable == "LithiumToxicity0.8", `:=`(var_label = "Lithium > 0.8 발생 여부", val_label = c("No", "Yes"))]
label.main[variable == "avgDose_1day", `:=`(var_label = "Average 1day dose")]
label.main[variable == "totYear_Lithium", `:=`(var_label = "Cumulative Lithium year")]
label.main[variable == "totYear_Valproate", `:=`(var_label = "Cumulative Valproate year")]
label.main[variable == "qd_Lithium", `:=`(var_label = "Lithium QD proportion")]
label.main[variable == "qd_Valproate", `:=`(var_label = "Valproate QD proportion")]

label.main[variable == "year0GFR", `:=`(var_label = "복용 1년 이내 GFR")]
label.main[variable == "year3GFR", `:=`(var_label = "복용 3년차 GFR")]
label.main[variable == "year5GFR", `:=`(var_label = "복용 5년차 GFR")]
label.main[variable == "year7GFR", `:=`(var_label = "복용 7년차 GFR")]
label.main[variable == "year10GFR", `:=`(var_label = "복용 10년차 GFR")]
label.main[variable == "year12GFR", `:=`(var_label = "복용 12년차 GFR")]
label.main[variable == "year15GFR", `:=`(var_label = "복용 15년차 GFR")]
label.main[variable == "year20GFR", `:=`(var_label = "복용 20년차 GFR")]
label.main[variable == "Fcode", `:=`(var_label = "Primary psychiatric diagnoses, ICD-10")]

## variable order : 미리 만들어놓은 KM, cox 모듈용

varlist_kmcox <- list(variable = c("eGFRbelow60", "year_FU", "drug", setdiff(names(data.main), c("eGFRbelow60", "year_FU", "drug", "NO" ))))