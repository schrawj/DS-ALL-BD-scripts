# Readme and TODO comments ------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.18.
#' 
#' Pasted this code in from the 'generating analytic datasets' script.
#' 
#' It deals with investigating the validity and utility of variables
#' rather than creating the analytic datasets per se.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



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





# Are common comorbid defects captured equally well across states? --------

require(xlsx); require(dplyr); require(gmodels); require(stringr)

#' Read in co-occurring defects. 
ds.defects <- read.xlsx('Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.comorbid.defect.counts.xlsx', sheetIndex = 1, stringsAsFactors = FALSE,
                        colIndex = c(1,2,5))

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
                          codes.used.to.dx = as.character())

for (i in 1:length(common.defects)){
  
  defect <- common.defects[i]
  
  pattern <- substr(defect, 1, 6)
  
  #' Returns integers mapping to indices of the columns in the transposed BD data frame that match the first five digits of the index BPA code.
  #' Then fetches the names of those columns and keeps the ones that are not valid codes.
  #' The idea being that if the first five digits match but the code is not listed in the TX DHHS codes, it is 
  #' one of these odd NC codes that is probably a synonym for the index BPA code.
  possible.nc.synonyms <- which(substr(all.defects, 1, 6) == pattern) + 1
  possible.nc.synonyms <- names(ds.leuk.bd.codes.txnc.transpose[possible.nc.synonyms])
  possible.nc.synonyms <- subset(possible.nc.synonyms, !(possible.nc.synonyms %in% bpa.defects$BPA_number))
  
  #' A list of the code plus all likely NC synonyms.
  code.plus.synonyms <- c(defect, possible.nc.synonyms)
  
  #' Get those columns.  Set NA values to 0.  Compute a dummy variable to flag rows where at least one code diagnosed.
  tmp.transpose <- ds.leuk.bd.codes.txnc.transpose[, c('studyid',code.plus.synonyms)]
  for (i in 2:length(names(tmp.transpose))){
    tmp.transpose[,i ] <- ifelse(is.na(tmp.transpose[,i]), 0, tmp.transpose[,i])
  }
  tmp.transpose$has.code <- ifelse(rowSums(tmp.transpose[2:length(names(tmp.transpose))]) >= 1, 1, 0)
  
  #' Reconstitute a variable for state; compute n and percent of DS cases in each state with the index code or a synonym.
  tmp.transpose$state <- substr(tmp.transpose$studyid, 1, 2)
  
  tab <- CrossTable(tmp.transpose$state, tmp.transpose$has.code, chisq = TRUE, fisher = TRUE)
  
  #' Save output, including a list of the additional codes considered to indicate the diagnosis.
  new.defect <- data.frame(bpa.code = defect,
                           tx.count = tab$t[2,2],
                           tx.percent = (tab$t[2,2]/7766)*100,
                           nc.count = tab$t[1,2],
                           nc.percent = (tab$t[1,2]/1479)*100, 
                           p.chisq = tab$chisq$p.value,
                           p.fisher = tab$fisher.ts$p.value,
                           codes.used.to.dx = paste(code.plus.synonyms, collapse = ", "))    
  
  defect.freq <- rbind(defect.freq, new.defect)
  
}

defect.freq <- left_join(defect.freq, bpa.defects, by = c('bpa.code' = "BPA_number"))
defect.freq <- defect.freq[, c(1,9,2:8)]

write.xlsx(defect.freq, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.comorbid.defect.counts.by.state.xlsx', row.names = FALSE)

rm(list = ls()); gc()



# Are NDBPN major defects captured equally well across states? ------------

require(dplyr); require(gmodels)

load('Z:/Jeremy/GOBACK/Datasets/goback.chrom.v20180711.rdata')

down <- filter(goback.chrom, down.syndrome == 1)

rm(goback.chrom); gc()

#' A vector of the specific major structural birth defects we've computed.
defect.names <- names(down[22:94])
defect.names <- subset(defect.names, !grepl('conganomalies.', defect.names))

defect.positions <- which(names(down) %in% defect.names) 

for (i in defect.positions){
  down[,i] <- ifelse(is.na(down[,i]), 0, down[,i])
}

#' Initialize an empty data frame to hold output.
defect.freq <- data.frame(name = as.character(),
                          tx.count = as.numeric(),
                          tx.percent = as.numeric(),
                          nc.count = as.numeric(),
                          nc.percent = as.numeric(), 
                          mi.count = as.numeric(),
                          mi.percent = as.numeric(),
                          ar.count = as.numeric(),
                          ar.percent = as.numeric(),
                          p.chisq = as.numeric())   

for (i in 1:length(defect.names)){
  
  defect <- defect.names[i]
  
  tab <- CrossTable(down$state, down[, defect], chisq = TRUE)
  
  #' Save output, including a list of the additional codes considered to indicate the diagnosis.
  new.defect <- data.frame(name = defect,
                           tx.count = tab$t[4,2],
                           tx.percent = (tab$t[4,2]/7775)*100,
                           nc.count = tab$t[3,2],
                           nc.percent = (tab$t[3,2]/1482)*100, 
                           mi.count = tab$t[2,2],
                           mi.percent = (tab$t[2,2]/3141)*100,
                           ar.count = tab$t[1,2],
                           ar.percent = (tab$t[1,2]/723)*100,
                           p.chisq = tab$chisq$p.value)
  
  defect.freq <- rbind(defect.freq, new.defect)
  
}

defect.freq <- arrange(defect.freq, desc(tx.count))

write.xlsx(defect.freq, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.major.defect.counts.by.state.xlsx', row.names = FALSE, sheetName = 'all_states')

rm(list = ls()); gc()

#' Repeat, sans Arkansas.
#' Not all names defects were provided to us by that state.
load('Z:/Jeremy/GOBACK/Datasets/goback.chrom.v20180711.rdata')

down <- filter(goback.chrom, down.syndrome == 1 & state != 'AR')

rm(goback.chrom); gc()

#' A vector of the specific major structural birth defects we've computed.
defect.names <- names(down[22:94])
defect.names <- subset(defect.names, !grepl('conganomalies.', defect.names))

defect.positions <- which(names(down) %in% defect.names) 

for (i in defect.positions){
  down[,i] <- ifelse(is.na(down[,i]), 0, down[,i])
}

defect.freq <- data.frame(name = as.character(),
                          tx.count = as.numeric(),
                          tx.percent = as.numeric(),
                          nc.count = as.numeric(),
                          nc.percent = as.numeric(), 
                          mi.count = as.numeric(),
                          mi.percent = as.numeric(),
                          p.chisq = as.numeric())   

for (i in 1:length(defect.names)){
  
  defect <- defect.names[i]
  
  tab <- CrossTable(down$state, down[, defect], chisq = TRUE)
  

  new.defect <- data.frame(name = defect,
                           tx.count = tab$t[3,2],
                           tx.percent = (tab$t[3,2]/7775)*100,
                           nc.count = tab$t[2,2],
                           nc.percent = (tab$t[2,2]/1482)*100, 
                           mi.count = tab$t[1,2],
                           mi.percent = (tab$t[1,2]/3141)*100,
                           p.chisq = tab$chisq$p.value)
  
  defect.freq <- rbind(defect.freq, new.defect)
  
}

defect.freq <- arrange(defect.freq, desc(tx.count))

write.xlsx(defect.freq, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.major.defect.counts.by.state.xlsx', row.names = FALSE, append = TRUE, sheetName = 'wo_AR')

rm(list = ls()); gc()



# Are organ system level variables captured equally well across st --------

require(dplyr); require(gmodels)

load('Z:/Jeremy/GOBACK/Datasets/goback.chrom.v20180711.rdata')

down <- filter(goback.chrom, down.syndrome == 1 & state != 'AR')

orgsys.vars <- grep('conganomalies.', names(down), value = TRUE)

exclude <- grep('.other', orgsys.vars, value = TRUE)

orgsys.vars <- c(orgsys.vars[!(orgsys.vars %in% exclude)], 'oral.clefts')

for (i in orgsys.vars){
  down[,i] <- ifelse(is.na(down[,i]), 0, down[,i])
}

#' Initialize an empty data frame to hold output.
defect.freq <- data.frame(name = as.character(),
                          tx.count = as.numeric(),
                          tx.percent = as.numeric(),
                          nc.count = as.numeric(),
                          nc.percent = as.numeric(), 
                          mi.count = as.numeric(),
                          mi.percent = as.numeric(),
                          p.chisq = as.numeric())   

for (i in 1:length(orgsys.vars)){
  
  defect <- orgsys.vars[i]
  
  tab <- CrossTable(down$state, down[, defect], chisq = TRUE)
  
  #' Save output, including a list of the additional codes considered to indicate the diagnosis.
  new.defect <- data.frame(name = defect,
                           tx.count = tab$t[3,2],
                           tx.percent = (tab$t[3,2]/7775)*100,
                           nc.count = tab$t[2,2],
                           nc.percent = (tab$t[2,2]/1482)*100, 
                           mi.count = tab$t[1,2],
                           mi.percent = (tab$t[1,2]/3141)*100,
                           p.chisq = tab$chisq$p.value)
  
  defect.freq <- rbind(defect.freq, new.defect)
  
}

defect.freq <- arrange(defect.freq, desc(tx.count))

write.xlsx(defect.freq, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.orgsys.defect.counts.by.state.xlsx', row.names = FALSE, sheetName = 'wo_AR')



# Summary statistics pertaining to Karen's suggestions --------------------

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180723.rdata')

#' At the same time Karen returned thoughts on the other luekemia diagnoses,
#' she suggested stratifying analyses of AML by those who are 4 years of age or less 
#' vs those who are more than 4 years of age.  Relationship to TMD is presumed to differ
#' among these strata.
#' Similarly, she suggested looking at AML cases <= 4 years with TMD vs without.
tmp <- filter(ds.leuk, aml == 1)
table(tmp$ds.leuk, floor(tmp$person.yrs))
table(tmp$conganomalies.digestivesystem.other, tmp$transient.myeloproliferative.disorder, useNA = 'ifany')
