# Readme and TODO comments ------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.10.
#' 
#' Updated from 2018.06.27.
#' 
#' Create and clean datasets for the analysis of cancer risk and survival
#' as a function of co-occurring birth defects in children diagnosed with 
#' Down syndrome.
#' 
#' I've verified that the existing number of total defects and binary 
#' organ systems defects variables are suitable for this analysis.
#' 
#' Several previously outstanding questions were answered through meetings
#' with Philip and/or Karen.
#' 1 We will include other leukemias in the analysis.  We will structure 
#'   the leukemia variables such that we study "any leukemia" and then
#'   pull out ALL and AML cases for more specific analyses.
#' 2 For the (planned, eventual) analysis of specific defects, we will 
#'   include defects that occur in 5% or more of the DS cases.  This is 
#'   equivalent to a defect occuring in more than 462 kids with DS.
#'   Because they are inherently of interest, we will also include:
#'   Tetralogy of Fallot, ASD, VSD, AVSD, and PDA regardless of whether they
#'   satisfy this criterion. 
#' 3 We will include all DS cases as a single unit.  We will not consider
#'   children with mosaicism or translocation karyotypes separately.
#' 4 For the few NC children with invalid BPA codes (758.008 or 758.098),
#'   we will do nothing different.  These children phenotypically resemble
#'   the other DS cases in terms of comorbid birth defects and the finer
#'   distinctions based on karyotype are not of interest to us in this 
#'   analysis.
#' 5 We will include Michigan children contingent on locating BD codes in 
#'   MI cancer cases.  As of today, these children are not included in the 
#'   data.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' TODO: Consider whether it is appropriate that the reference group of 
#'       non-DS kids only includes ALL and AML as opposed to any leukemias.
#'       The DS kids include those diagnosed with any leukemia.
#' TODO: Consider that the birth defects variables need revision.  Rather than comparing
#'       DS cases with the index anomaly to DS cases with no co-occurring
#'       anomalies, we should be comparing them to DS cases without the 
#'       index anomaly, yes?
#' TODO: Should the NA values for minordefect.total be replaced with zeroes?
#' TODO: Clean MI birth defects codes and recalculate number of defects in this state.
#'       Include only defects captured by ICD9 codes that map to something in the BPA 
#'       system.  This should exclude the additional conditions MI collects data on.



# Generate the core analytic dataset --------------------------------------

require(dplyr)

load('Z:/Jeremy/GOBACK/Datasets/goback.v20180711.rdata')

#' Extract DS kids and non-DS kids with ALL or AML.  Exclude AR children.
ds.leuk <- filter(goback, down.syndrome == 1 | all == 1 | aml == 1)
ds.leuk <- filter(ds.leuk, state != 'AR')

#' Exclude DS cases with non-leukemia cancers listed as the first primary.
exclude.cancer1 <- filter(ds.leuk, down.syndrome == 1 & cancer == 1 & !(cancer1 %in% c('all','aml','leu.other')))

ds.leuk <- filter(ds.leuk, !(studyid %in% exclude.cancer1$studyid))

#' 289 of the remaining children are marked NA for Down syndrome.
#' These children have a cancer diagnosis, do not have Down syndrome, and do have some other defect(s).
#' For purposes of this analysis, these children should be coded 0 with respect to Down syndrome.
ds.leuk$down.syndrome <- ifelse(is.na(ds.leuk$down.syndrome), 0, ds.leuk$down.syndrome)

#' Compute a three-level factor variable grouping children into the categories described in point 1 in the Readme.
ds.leuk$ds.leuk <- factor(ifelse(ds.leuk$down.syndrome == 0 & (ds.leuk$aml == 1 | ds.leuk$all == 1), 0,
                                 ifelse(ds.leuk$down.syndrome == 1 & ds.leuk$leu.any == 0, 1,
                                        ifelse(ds.leuk$down.syndrome == 1 & ds.leuk$leu.any == 1, 2, NA))),
                          levels = c(0:2),
                          labels = c('Non-DS ALL/AML', 'DS Non-Leukemia', 'DS Leukemia'))

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' The organ system level birth defects dummy variables are valid for this 
#' analysis.
#' 
#' The number of birth defects variable is as well, although we should 
#' substract 1 to account for the fact that all kids have a DS diagnosis.
#' 
#' We will need to compute a variable for any co-occurring birth defect.
#' This basically could be 1 if the value (defect.total - 1) > 0.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' Compute number of comorbid defects for DS kids.
#' Results in negative number of birth defects for non-DS ALL/AML cases with no birth defects.
ds.leuk$number.of.comorbid.defects <- ds.leuk$defect.total - 1
ds.leuk$number.of.comorbid.defects <- ifelse(ds.leuk$number.of.comorbid.defects < 0, 0, ds.leuk$number.of.comorbid.defects)

#' Compute dummy variable for any comorbid defect.
ds.leuk$any.comorbid.defect <- ifelse(ds.leuk$number.of.comorbid.defects >= 1 & ds.leuk$ds.leuk != 'Non-DS ALL/AML', 1,
                                      ifelse(ds.leuk$ds.leuk == 'Non-DS ALL/AML', NA, 0))

#' Replace with NA's with 0's for unaffected children across all the organ system dummy variables.
#' We will just be comparing children with or without a defect in that organ system regardless of their other defects.
orgsys.vars <- c(22,30,35,38,61,65,68,76,83,94)

for (i in orgsys.vars){
  print(names(ds.leuk[i]))
  print(table(ds.leuk[,i]))
  ds.leuk[,i] <- ifelse(is.na(ds.leuk[,i]), 0, ds.leuk[,i])
  print(table(ds.leuk[,i]))
}

#' There's an obvious overabundance of craniosynostosis diagnoses in TX.
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata')

#' This is caused, I think, by considering kids with 'other specified skull and face bone anomalies' to have craniosynostosis.
#' You can see that here.
match <- grep('756.0', names(ds.leuk.bd.codes.txnc.transpose))

tmp <- select(ds.leuk.bd.codes.txnc.transpose, studyid, match)
for (i in 2:41){
  tmp[,i] <- ifelse(is.na(tmp[,i]), 0, tmp[,i])
}
tmp$has.code <- ifelse(rowSums(tmp[2:41]) >= 1, 1, 0)

#' These are the known and presumed codes for craniosynostosis in TX and NC.
#' There are technically no codes for craniosynostosis in MI: they all map to 756.0, congenital anomalies of the skull and face bones.
cranio.codes <- c('756.000','756.005','756.006','756.010','756.020','756.030', # These are the valid BPA codes for craniosynostosis per the TX DHHS handbook.
                  '756.008','756.011','756.012','756.013','756.014','756.018', '756.032','756.034' # These are the presumed synonymous codes sporadically recored in NC.
)

ds.leuk.bd.codes.txnc.transpose <- select(ds.leuk.bd.codes.txnc.transpose, studyid, cranio.codes)

for (i in 2:15){
  ds.leuk.bd.codes.txnc.transpose[,i] <- ifelse(is.na(ds.leuk.bd.codes.txnc.transpose[,i]), 0, ds.leuk.bd.codes.txnc.transpose[,i])
}

ds.leuk.bd.codes.txnc.transpose$has.code <- ifelse(rowSums(ds.leuk.bd.codes.txnc.transpose[2:15]) >= 1, 1, 0)
ds.leuk.bd.codes.txnc.transpose <- filter(ds.leuk.bd.codes.txnc.transpose, has.code == 1)

cranio.ids <- c(ds.leuk.bd.codes.txnc.transpose$studyid)

ds.leuk$craniosynostosis <- ifelse(ds.leuk$studyid %in% cranio.ids, 1, 0)

save(ds.leuk, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180718.rdata')

rm(list = ls()); gc()



# Generate variables for any major defect in each organ sys ---------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' There are major differences in the rates at which minor defects are 
#' reported by state, with states other than Texas appearing not to have 
#' captured them to any appreciable extent.
#' 
#' With that in mind, the organ system level dummy variables need to be 
#' revised.  We will compute variables for whether there is any MAJOR 
#' defect in that organ system, with major defects defined according to the 
#' NBDPN document in the data dictionaries directory.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(gmodels)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180718.rdata')
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata')
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.mi.transpose.v20180712.rdata')

#' Replace NAs with 0's for all structural birth defects among children without the index birth defect.
for (i in 22:94){ds.leuk[,i] <- ifelse(is.na(ds.leuk[,i]), 0, ds.leuk[,i])}

ds.leuk$cns.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(23,24,26,28)]) >= 1, 1, 0)
ds.leuk$eye.any.major.anomaly <- ifelse(rowSums(ds.leuk[31:32]) >= 1, 1, 0)
ds.leuk$ear.any.major.anomaly <- ifelse(rowSums(ds.leuk[36]) >= 1, 1, 0)
ds.leuk$cardiovascular.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(49,55,43,48,40,52,46,47,51,56,41,44,42,57,54)]) >= 1, 1, 0)
ds.leuk$orofacial.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(62,65)]) >= 1, 1, 0)
ds.leuk$gi.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(73,69,70,74)]) >= 1, 1, 0)

#' Must compute variable for congenital posterior urethral valves and cloacal exstrophy to proceed.
cpuv.ids.txnc <- filter(ds.leuk.bd.codes.txnc.transpose, `753.600` == 1 | `753.604` == 1 | `753.608` == 1 | `753.609` == 1)
cpuv.ids.mi <- filter(ds.leuk.bd.codes.mi.transpose, `753.6` == 1)
cpuv.ids <- c(cpuv.ids.txnc$studyid, cpuv.ids.mi$studyid)

cloaca.ids.txnc <- filter(ds.leuk.bd.codes.txnc.transpose, `751.550` == 1 | `756.790` == 1)
cloaca.ids.mi <- filter(ds.leuk.bd.codes.mi.transpose, `751.5` == 1)
cloaca.ids <- c(cloaca.ids.mi$studyid, cloaca.ids.txnc$studyid)

ds.leuk$congenital.posterior.urethral.vales <- ifelse(ds.leuk$studyid %in% cpuv.ids, 1, 0)
ds.leuk$cloacal.exstrophy <- ifelse(ds.leuk$studyid %in% cloaca.ids, 1, 0)

ds.leuk$genitourinary.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(78,80,77,166,167)]) >= 1, 1, 0)
ds.leuk$musculoskeletal.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(84,87,88,90,91,92)] >= 1), 1, 0)

ds.leuk <- select(ds.leuk, -cloacal.exstrophy, -congenital.posterior.urethral.vales, -any.comorbid.defect)

ds.leuk$any.major.comorbid.anomaly <- ifelse(rowSums(ds.leuk[159:166]) >= 1, 1, 0)

save(ds.leuk, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180718.2.rdata')



# Generate BD  (TX and NC) and cancer (all states) datasets ---------------

require(dplyr)

#' Need a clean environment for the references to ls() below to work.
rm(list = ls())

#' As of now, only TX and NC children in the ds.leuk dataset.
load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180710.rdata')
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")

data.frames.names <- paste0('ds.leuk.',ls()[1:3])

data.frames <- list(#bd.codes.mi, 
                    #bd.codes.mi.transpose, 
                    bd.codes.txnc, 
                    bd.codes.txnc.transpose, 
                    cancer.codes)

names(data.frames) <- data.frames.names

#' Today's teachable moment: if you want to supply the name of the object to be saved as a character vector, 
#' you must use the 'list =' argument. 
for (i in 1:length(data.frames)){
  
  tmp <- c(ds.leuk$studyid)
  
  print(names(data.frames[i]))
  print(dim(data.frames[[i]]))
  data.frames[[i]] <- filter(data.frames[[i]], data.frames[[i]]$studyid %in% tmp)
  print(dim(data.frames[[i]]))
  
  tmp <- data.frames[[i]]
  
  assign(data.frames.names[i],tmp)
  
  path <- paste0('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/',data.frames.names[i],'.v20180711.rdata')
  
  save(list = data.frames.names[i], file = path)
  
}

rm(list = ls()); gc()



# Generate birth defects codes data frames: MI ----------------------------

require(readstata13); require(dplyr)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.transpose.v20180712.rdata')

ds.leuk <- filter(ds.leuk, state == 'MI')
ds.leuk <- c(ds.leuk$studyid)

#' Results in 3302 of 3307 rows returned.
ds.leuk.bd.codes.mi <- filter(bd.codes.mi, studyid %in% ds.leuk)
ds.leuk.bd.codes.mi.transpose <- filter(bd.codes.mi.transpose, studyid %in% ds.leuk)

save(ds.leuk.bd.codes.mi, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.mi.v20180712.rdata')
save(ds.leuk.bd.codes.mi.transpose, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.mi.transpose.v20180712.rdata')

