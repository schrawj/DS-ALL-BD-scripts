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

#' TODO: Find BD codes for MI cancer cases.
#' TODO: Once found, update birth defects codes data frames if necessary
#'       with the new set of MI children to be included.
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

#' Replace 99 with NA for missing paternal age values.
ds.leuk$f.age <- ifelse(ds.leuk$f.age == 99, NA, ds.leuk$f.age)

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

save(ds.leuk, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.rdata')

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

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.rdata')

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

save(ds.leuk, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')

for (i in 159:167){
  print(names(ds.leuk[i]))
  print(CrossTable(ds.leuk[,i], ds.leuk$state, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE))
}

for (i in c(84,87,88,90,91,92)){
  print(names(ds.leuk[i]))
  print(CrossTable(ds.leuk[,i], ds.leuk$state, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE))
}

craniosynostosis.ids <- filter(ds.leuk, down.syndrome == 1 & craniosynostosis == 1)
craniosynostosis.ids <- c(craniosynostosis.ids$studyid)

craniosynostosis.ids <- filter(ds.leuk.bd.codes.txnc, studyid %in% craniosynostosis.ids)
craniosynostosis.ids <- filter(ds.leuk.bd.codes.txnc.transpose, studyid %in% craniosynostosis.ids)

cranio.cases <- filter(ds.leuk.bd.codes.txnc.transpose, rowSums(ds.leuk.bd.codes.txnc.transpose[1873:1891], na.rm = TRUE)  >= 1)

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

# Investigate the spectrum of co-occuring birth defects -------------------

#' One question I've been pondering today is whether we should use the computed birth defects variables, or return to the codes themselves
#' for a more comprehensive picture.  One way to make that decision might be to look at the defects that occur in a non-trivial fraction of 
#' individuals with DS (let's say 5% or more) and see how well captured these are by the defects variables we've already computed.

require(dplyr); require(xlsx)

load("Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180627.rdata")
load("Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.rdata")

ds.ids <- select(filter(ds.leuk, down.syndrome == 1), studyid)
ds.ids <- c(ds.ids$studyid)

#' All elements of ds.ids map successfully to a row in the bd.codes data frame.
ds.leuk.bd.codes.txnc <- filter(ds.leuk.bd.codes.txnc, studyid %in% ds.ids)

#' Generate a data frame with one column for defect codes, one column for the count of the number of times they're observed.
running.total <- c()

for (i in 1:nrow(ds.leuk.bd.codes.txnc)){
  
  new.elements <- as.character(ds.leuk.bd.codes.txnc[i,2:67])
  new.elements <- subset(new.elements, !is.na(new.elements))
  
  running.total <- c(running.total, new.elements)
  
}

tab <- table(running.total)

tab.out <- data.frame(bpa.code = names(tab), defect.count = tab)
tab.out <- select(rename(tab.out, freq = defect.count.Freq), bpa.code, freq)
tab.out <- arrange(tab.out, desc(freq))

#' We can put names to these by reading in the dx_map file Peter Langlois created.
#' Do so, keep only rows with at least 1 occurrence in a child with DS.
dx.map <- read.xlsx(file = 'Z:/Jeremy/GOBACK/Data dictionaries and data definitions/dxmap_current.xlsx',
                    sheetName = 'Current', colIndex = c(1,2,13,14), colClasses = c('character','character','character','character'))
dx.map <- rename(dx.map, bpa.code = BPA_number, bpa.name = BPANAME, icd9.code = ICD9_code, icd9.name = ICD9_name)

ds.leuk.defects <- left_join(dx.map, tab.out, by = 'bpa.code')
ds.leuk.defects <- arrange(filter(ds.leuk.defects, !is.na(ds.leuk.defects$freq)), desc(freq))

#' It's obvious from reviewing the results for common defets that the computed defects variables will not
#' capture the spectrum of co-occurring BDs in children with DS well.  Unless they are of special interest, 
#' we will want to use a new set of defects for this analysis.  Discuss with Philip.
save(ds.leuk.defects, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.comorbid.defect.counts.v20180628.rdata')

write.xlsx(ds.leuk.defects, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.comorbid.defect.counts.xlsx',
           row.names = FALSE)



# Address discrepancy in DS cases by computed variable and codes ----------

#' NOTE: The cause of this issue is now known.  
#' The numbers are different because a small
#' number of NC DS cases have codes ending in .xx8 instead of .xx0.
require(dplyr); require(xlsx)

load("Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180627.rdata")
load("Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.rdata")
load("Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.rdata")

ds.ids.computed <- filter(ds.leuk, state %in% c('TX','NC') & down.syndrome == 1)
ds.ids.computed <- c(ds.ids.computed$studyid)

ds.ids.bycode <- filter(ds.leuk.bd.codes.txnc.transpose, studyid %in% ds.ids.computed)
ds.ids.bycode <- filter(ds.ids.bycode,
                        `758.000` == 1 | `758.010` == 1 | `758.020` == 1 | 
                        `758.030` == 1 | `758.040` == 1 | `758.090` == 1 |
                        `758.008` == 1 | `758.098` == 1) # I can't find evidence that these two are valid codes, but if they are/were they fall in the range for DS.
                                                         # They're also sufficient to account for the difference in the computed vs. coded number of DS cases.

#' Searching these children in the bd.codes data frames reveals that they do in fact have co-occurring defects that are common in the DS cases with valid codes.
#' I'll mention this issue to Philip but my inclination is that these really are DS cases, and should probably be treated as their nearest valid code.
questionable.ids <- filter(ds.leuk.bd.codes.txnc.transpose,
                           `758.008` == 1 | `758.098` == 1)
questionable.ids <- c(questionable.ids$studyid)

tmp <- filter(ds.leuk.bd.codes.txnc, studyid %in% questionable.ids)
tmp[1:100, 1:8]

write.xlsx(tmp, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/defects.in.kids.with.unusal.DS.codes.xlsx',row.names = FALSE, showNA = FALSE)

#' In TX and NC, how many kids have each code?
codes <- c('758.000', '758.010', '758.020', '758.030', '758.040', '758.090', '758.008', '758.098')

for (i in codes){
  print(i)
  tmp <- filter(ds.leuk.bd.codes.txnc.transpose, ds.leuk.bd.codes.txnc.transpose[,i] == 1)
  tmp <- c(tmp$studyid)
  tmp <- filter(ds.leuk, studyid %in% tmp)
#  print(nrow(tmp))
  print(table(tmp$state, tmp$leu.any))
}


load("Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.mi.transpose.rdata")

mi.codes <- c("758", "758.000", "758.008", "758.010", "758.020", "758.030", "758.040", "758.090", "758.098")

for (i in mi.codes){
  print(i)
  tmp <- filter(ds.leuk.bd.codes.mi.transpose, ds.leuk.bd.codes.mi.transpose[,i] == 1)
  tmp <- c(tmp$studyid)
  tmp <- filter(ds.leuk, studyid %in% tmp)
  #  print(nrow(tmp))
  print(table(tmp$all, useNA = 'ifany'))
}

# Investigate other leukemias in DS kids ----------------------------------

#' What are the specific cancer diagnoses in DS kids with other leukemias?
#' Would any of these be admissible to roll into the DS-ALL/AML group?
require(dplyr); require(gmodels); require(xlsx)

load("Z:/Jeremy/GOBACK/Datasets/goback.v20180611.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata")

leu.other <- filter(goback, down.syndrome == 1 & leu.other == 1)

leu.other.ids <- c(leu.other$studyid)
leu.other.ids <- select(filter(cancer.codes, studyid %in% leu.other.ids),
                        studyid, morph31, site_code1)

#' Generate some aggregate info about the number of times each ICD-O-3 morphology code appears.
tab <- table(leu.other.ids$morph31)

leu.other.freq <- rename(data.frame(morph.codes = names(tab), freq = tab),
                         freq = freq.Freq)
leu.other.freq <- leu.other.freq[,c(1,3)]

write.xlsx(leu.other.freq, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/leu.other.morph.codes.in.ds.kids.xlsx',
           row.names = FALSE)

#' Karen says the ages at diagnosis can be informative for determining whether some of these other leukemias are 
#' actually classifiable as ALL or AML.
#' Generate an individual-level report including the description of the morphology code (which I manually added to the file above, outside of R)
#' and the age at diagnosis for that child.

morph.descriptions <- read.xlsx(file = 'Z:/Jeremy/DS-ALL BD project/R outputs/leu.other.morph.codes.in.ds.kids.xlsx', sheetIndex = 1,
                                stringsAsFactors = FALSE)
morph.descriptions <- rename(morph.descriptions, morph31 = morph.codes)
morph.descriptions$morph31 <- as.numeric(morph.descriptions$morph31)

leu.other.ids <- left_join(select(leu.other.ids, studyid, morph31), 
                           select(morph.descriptions, morph31, description), 
                           by = 'morph31')
leu.other.ids <- left_join(leu.other.ids,
                           select(
                                  filter(goback, down.syndrome == 1 & leu.other == 1), studyid, person.yrs),
                           by = 'studyid')

write.xlsx(leu.other.ids, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/ds.other.leukemia.cases.xlsx', row.names = FALSE)

rm(list = ls()); gc()

# Verify accuracy of existing organ system level variables ------------

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180710.rdata')
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata')

#' Some of the initial analyses will focus on identifying whether there is a difference in the average number of birth defects 
#' among DS cases accoding to cancer status.
#' Because the defects that co-occur in this population don't overlap entirely with the ones we have analyzed before, there will
#' be some computation of new dummy variables.  It also strikes me as being worth verifying that the existing organ-system level variables
#' account for all codes in the range of that organ system, and not just the major defects we investigated previously.
cns.patterns <- c('740.','741.','742.')

cns.columns <- data.frame(column.names = names(ds.leuk.bd.codes.txnc.transpose),
                          column.index = 1:length(names(ds.leuk.bd.codes.txnc.transpose)))
cns.columns$column.names <- as.character(cns.columns$column.names)
cns.columns <- subset(cns.columns, substr(cns.columns$column.names, 1, 4) %in% cns.patterns)

tmp <- ds.leuk.bd.codes.txnc.transpose[, c(1, cns.columns$column.index)]
tmp$any.cns.defect <- rowSums(tmp[2:length(names(tmp))], na.rm = TRUE)
table(tmp$any.cns.defect)

#' This agrees exactly with the existing conganomalies.cns variable.
#' Let's try it for one more organ system to make sure.
efn.patterns <- '744.' 

efn.columns <- data.frame(column.names = names(ds.leuk.bd.codes.txnc.transpose),
                          column.index = 1:length(names(ds.leuk.bd.codes.txnc.transpose)))
efn.columns$column.names <- as.character(efn.columns$column.names)
efn.columns <- subset(efn.columns, substr(efn.columns$column.names, 1, 4) %in% efn.patterns)

tmp <- ds.leuk.bd.codes.txnc.transpose[, c(1, efn.columns$column.index)]
tmp$any.efn.defect <- rowSums(tmp[2:length(names(tmp))], na.rm = TRUE)
table(tmp$any.efn.defect >= 1)

#' Also matches.  I'm satisfied.
rm(list = ls()); gc()



# Verify accuracy of existing number of birth defects variable ------------

require(dplyr)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180710.rdata')
load("Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata")
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.v20180711.rdata')

table(ds.leuk$majordefect.total, useNA = 'ifany')

#' NAs may need replaced with zero.
table(ds.leuk$minordefect.total, useNA = 'ifany')

table(ds.leuk$defect.total, useNA = 'ifany')

#' Total is indeed sum of major and minor.
tmp <- select(ds.leuk, studyid, defect.total, majordefect.total, minordefect.total, down.syndrome)
tmp <- filter(tmp, defect.total > 0)

ds.leuk.bd.codes.txnc.transpose$defect.total <- rowSums(ds.leuk.bd.codes.txnc.transpose[2:length(names(ds.leuk.bd.codes.txnc.transpose))], na.rm = TRUE)

table(ds.leuk.bd.codes.txnc.transpose$defect.total, useNA = 'ifany')
table(tmp$defect.total, useNA = 'ifany')

olddefecttotal <- select(tmp, studyid, defect.total, down.syndrome)
newdefecttotal <- select(ds.leuk.bd.codes.txnc.transpose, studyid, defect.total)
defecttotal <- filter(left_join(olddefecttotal, newdefecttotal, by = 'studyid'),
                      defect.total.x != defect.total.y)
defecttotal$delta <- defecttotal$defect.total.x - defecttotal$defect.total.y
head(defecttotal, 100)

#' Old variable always counts more defects than the new.
table(defecttotal$delta)

ids <- c(defecttotal$studyid)

tmp <- filter(ds.leuk.bd.codes.txnc, studyid %in% ids)

#' This is due to the fact that the new transposed birth defects data frames do not account for the fact that the child may have the same code appear more than
#' once when that code possibly describes more than one defect.
head(defecttotal, 10)
head(tmp, 10)

#' Compare old against new computed from the non-transposed data frame.
tmp <- ds.leuk.bd.codes.txnc

for (i in 2:length(names(tmp))){
  tmp[, i] <- ifelse(!is.na(tmp[,i]), 1, 0)
}

tmp$defect.total <- rowSums(tmp[2:length(names(tmp))], na.rm = TRUE)
newdefecttotal <- select(tmp, studyid, defect.total)
defecttotal <- filter(left_join(olddefecttotal, newdefecttotal, by = 'studyid'))
defecttotal$delta <- defecttotal$defect.total.x - defecttotal$defect.total.y

#' Again, an exact match.
table(defecttotal$delta)

rm(list = ls()); gc()

# Verify accuracy of BD variables in MI kids with cancer ------------------

require(readstata13); require(dplyr)

load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.transpose.v20180614.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180606.rdata")
load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.rdata')

#' A file supplied by Tiffany with the actual ICD9 birth defects codes in MI children with cancer.
mi.bd.cancercases <- read.dta13('Z:/Jeremy/DS-ALL BD project/Datasets/Raw datasets/mi_cancer_long_jeremy_071118.dta', convert.underscore = TRUE)
mi.bd.cancercases <- subset(mi.bd.cancercases, !duplicated(mi.bd.cancercases$studyid))

#' Extract BD codes and some demographic variables.
mi.bd.cancercases <- mi.bd.cancercases[ , c(2,6,10,109:132)]
mi.bd.cancercases$row.sum <- rowSums(mi.bd.cancercases[4:27], na.rm = TRUE)
mi.bd.cancercases <- select(filter(mi.bd.cancercases, row.sum > 0), -row.sum)
mi.bd.cancercases$studyid <- paste0('mi',mi.bd.cancercases$studyid)

tmp <- filter(ds.leuk, state == 'MI', down.syndrome == 1, cancer == 1)
tmp <- tmp[, c(1,7,6,19,22,30,35,38,61,65,68,76,83,94)]
tmp <- arrange(tmp, studyid)

check <- left_join(tmp, mi.bd.cancercases, by = 'studyid')
check <- subset(check, !duplicated(check$studyid))

#' Many codes are recorded that do not map to defects assessed in Texas.
#' These include the cancer itself in some instances.
#' Any analysis of number of defects or presence or absence of specific defects will need to address this.
#' Probably we will have to wipe records of the codes that don't map to anything in TX and NC.
#' However, the current organ system level variables appear to be valid.
check[1:60, c(1, 5, 17:30)]
check[1:60, c(1, 6, 17:30)]
check[1:60, c(1, 11, 17:30)]

#' Direct conversion of these variables to character will not work well.
check$ICD9COD1 <- as.character(check$ICD9COD1)






# Which defects are captured equally well across states? ------------------

require(xlsx); require(dplyr); require(gmodels); require(stringr)

#' Read in co-occurring defects and filter at ~5% frequency. 
ds.defects <- read.xlsx('Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.comorbid.defect.counts.xlsx', sheetIndex = 1, stringsAsFactors = FALSE,
                        colIndex = c(1,2,5))
ds.defects <- filter(ds.defects, freq >= 594)

bpa.defects <- read.xlsx('Z:/Jeremy/GOBACK/Data dictionaries and data definitions/dxmap_current.xlsx', sheetName = 'Current', header = TRUE, 
                         stringsAsFactors = FALSE, colIndex = 1:2)

load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata')
load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')

common.defects <- ds.defects$bpa.code

all.defects <- colnames(ds.leuk.bd.codes.txnc.transpose[2:2362])

#' Only include DS children from TX or NC in this analysis.
ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML' & state %in% c('TX','NC'))

defect.freq <- data.frame(bpa.code = as.character(),
                          tx.count = as.numeric(),
                          tx.percent = as.numeric(),
                          nc.count = as.numeric(),
                          nc.percent = as.numeric(),
                          p.chisq = as.numeric(),
                          p.fisher = as.numeric(),
                          similar.nc.codes = as.character())

for (i in 1:length(common.defects)){
  
  defect <- common.defects[i]
  
  #'   tmp <- subset(ds.leuk.bd.codes.txnc.transpose, ds.leuk.bd.codes.txnc.transpose[, defect] == 1)
  #'   tmp <- c(tmp$studyid)
  #'   tmp <- filter(ds.leuk, studyid %in% tmp)
  #'   tmp <- c(tmp$studyid)
  
  #'   tmp.ds.leuk <- ds.leuk
  
  #'   tmp.ds.leuk[, defect] <- ifelse(tmp.ds.leuk$studyid %in% tmp, 1, 0)
  
  #'   tab <- CrossTable(tmp.ds.leuk[, defect], tmp.ds.leuk$state, chisq = TRUE, fisher = TRUE)
  
  pattern <- substr(defect, 1, 6)
  
  possible.nc.synonyms <- which(substr(all.defects, 1, 6) == pattern) + 1
  possible.nc.synonyms <- names(ds.leuk.bd.codes.txnc.transpose[possible.nc.synonyms])
  possible.nc.synonyms <- subset(possible.nc.synonyms, !(possible.nc.synonyms %in% bpa.defects$BPA_number))
  
  #' TODO: Figure out how to concatenate all NC codes into a single column.
  #'   new.defect <- data.frame(bpa.code = defect,
  #'                            tx.freq = tab$t[1,2],
  #'                            nc.freq = tab$t[1,1],
  #'                            similar.nc.codes = c(possible.nc.synonyms),
  #'                            p.chisq = tab$chisq$p.value,
  #'                            p.fisher = tab$fisher.gt$p.value)  
  
  #' Another thought I had: generate a vector of the code plus all likely nc synonyms.
  #' Select those columns from the transposed BD data frame.
  #' Filter rows with rowSums >= 1 (i.e., tagged with at least one of the codes).
  #' Extract studyid as a character vector.
  #' Filter ds.leuk for these ids, then creata table for state and record frequencies for this.
  #' Should capture NC codes better.
  
  code.plus.synonyms <- c(defect, possible.nc.synonyms)
  
  tmp.transpose <- ds.leuk.bd.codes.txnc.transpose[, c('studyid',code.plus.synonyms)]
  for (i in 2:length(names(tmp.transpose))){
    tmp.transpose[,i ] <- ifelse(is.na(tmp.transpose[,i]), 0, tmp.transpose[,i])
  }
  tmp.transpose$has.code <- ifelse(rowSums(tmp.transpose[2:length(names(tmp.transpose))]) >= 1, 1, 0)
  #tmp.transpose <- filter(tmp.transpose, has.code == 1)
  
  tmp.transpose$state <- substr(tmp.transpose$studyid, 1, 2)
  
  #tmp.ds.leuk <- filter(ds.leuk, studyid %in% tmp.transpose$studyid)
  
  tab <- CrossTable(tmp.transpose$state, tmp.transpose$has.code, chisq = TRUE, fisher = TRUE)
  
  new.defect <- data.frame(bpa.code = defect,
                           tx.count = tab$t[2,2],
                           tx.percent = (tab$t[2,2]/7766)*100,
                           nc.count = tab$t[1,2],
                           nc.percent = (tab$t[1,2]/1479)*100, 
                           p.chisq = tab$chisq$p.value,
                           p.fisher = tab$fisher.gt$p.value,
                           similar.nc.codes = paste(code.plus.synonyms, collapse = ", "))    
  
  defect.freq <- rbind(defect.freq, new.defect)
  

  

}

#' TODO: Bind in BPA names.
#' TODO: Compute percentage of affected children in each state.