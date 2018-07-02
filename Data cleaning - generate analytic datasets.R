#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.06.27.
#' 
#' Create and clean datasets for the analysis of cancer risk and survival
#' as a function of co-occurring birth defects in children diagnosed with 
#' Down syndrome.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' TODO: Do we want to include other leukemias?  
#' TODO: Which defects do we want to include in the analysis?  I say those that occur in some non-trivial portion of the kids.
#' TODO: Get Philip's thoughts on what to do about mosaic cases - these are identified by the 758.040 code in the BPA system.
#'       Mosaic and complete cases cannot be distinguished from ICD9 codes.  Should we include them?  
#'       For reference, there are 135 such cases in TX and NC.  This is almost exactly the expected 1%.  
#'       3 of them eventually developed leukemias.
#'       A further 4% should have translocation trisomy - where a full or partial copy of 21 is attached to another chromosome.
#'       They still have 46 chromosomes, but have the symptoms of Down syndrome.  
#'       I don't know whether their expressivitiy is comparable.
#' TODO: What do we do about children with 758.008 and 758.098 DS codes?  These are not valid codes as far as I can tell, but these
#'       kids do have defects that are characteristic of the other DS cases.  I say, treat them as their nearest valid neighbor:
#'       758.008 becomes 758.000.
#'       758.098 becomes 758.090.

# Generate the core analytic dataset --------------------------------------

require(dplyr); require(gmodels); require(tictoc)

load('Z:/Jeremy/GOBACK/Datasets/goback.v20180611.rdata')

#' Very few solid tumors were diagnosed.
#' There are 48 other leukemia cases among children with Down syndrome.
tmp <- filter(goback, down.syndrome == 1 & cancer == 1 & state != 'AR')
table(tmp$cancer1)

#' Generate a list of DS children with cancers other than ALL or AML.
#' We will exclude these from the analysis unless Philip says otherwise.
exclusions <- filter(goback, down.syndrome == 1 & cancer == 1 & all == 0 & aml == 0)
exclusions <- c(exclusions$studyid)

ds.leuk <- filter(filter(goback, down.syndrome == 1 | all == 1 | aml == 1),
                  state != 'AR')
ds.leuk <- subset(ds.leuk, !(ds.leuk$studyid %in% exclusions))

#' 289 of these children are marked NA for Down syndrome.
#' These children have a cancer diagnosis, do not have Down syndrome, and do have some other defect(s).
#' For purposes of this analysis, these children should be coded 0 with respect to Down syndrome.
ds.leuk$down.syndrome <- ifelse(is.na(ds.leuk$down.syndrome), 0, ds.leuk$down.syndrome)

#' Compute a three-level factor variable grouping children into the categories described above.
ds.leuk$ds.leuk <- factor(ifelse(ds.leuk$down.syndrome == 0 & (ds.leuk$aml == 1 | ds.leuk$all == 1), 0,
                                 ifelse(ds.leuk$down.syndrome == 1 & ds.leuk$all == 0 & ds.leuk$aml == 0, 1,
                                        ifelse(ds.leuk$down.syndrome == 1 & (ds.leuk$all == 1 | ds.leuk$aml == 1), 2, NA))),
                          levels = c(0:2),
                          labels = c('Non-DS ALL/AML', 'DS Non-ALL/AML', 'DS ALL/AML'))

#' Check validity of new computed variable against the old ones.
CrossTable(ds.leuk$down.syndrome, ds.leuk$ds.leuk)
CrossTable(ds.leuk$leu.any, ds.leuk$ds.leuk)
CrossTable(ds.leuk$cancer, ds.leuk$ds.leuk)

#' Group membership counts by state.
CrossTable(ds.leuk$state, ds.leuk$ds.leuk)

save(ds.leuk, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180627.rdata')

rm(list = ls()); gc()



# Generate birth defects and cancer codes datasets ------------------------

require(dplyr)

#' Need a clean environment for the references to ls() below to work.
rm(list = ls())

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180627.rdata')
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.transpose.v20180614.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180606.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")

data.frames.names <- paste0('ds.leuk.',ls()[1:5])

data.frames <- list(bd.codes.mi, bd.codes.mi.transpose, bd.codes.txnc, bd.codes.txnc.transpose, cancer.codes)
names(data.frames) <- data.frames.names

#' Today's teachable moment: if you want to supply the name of the object to be saved as a character vector, 
#' you must use the 'list =' argument. 
for (i in 1:length(data.frames)){
  
  print(names(data.frames[i]))
  print(dim(data.frames[[i]]))
  data.frames[[i]] <- filter(data.frames[[i]], data.frames[[i]]$studyid %in% tmp)
  print(dim(data.frames[[i]]))
  
  tmp <- data.frames[[i]]
  
  assign(data.frames.names[i],tmp)
  
  path <- paste0('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/',data.frames.names[i],'.rdata')
  
  save(list = data.frames.names[i], file = path)
  
}

rm(list = ls()); gc()

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

#' What are the specific diagnoses in DS kids with other leukemias?
#' Would any of these be admissible to roll into the DS-ALL/AML group?
require(dplyr); require(gmodels); require(xlsx)

load("Z:/Jeremy/GOBACK/Datasets/goback.v20180611.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata")

leu.other <- filter(goback, down.syndrome == 1 & leu.other == 1)

leu.other.ids <- c(leu.other$studyid)
leu.other.ids <- filter(cancer.codes, studyid %in% leu.other.ids)
print(arrange(leu.other.ids, morph31))

tab <- table(leu.other.ids$morph31)

leu.other.freq <- rename(data.frame(morph.codes = names(tab), freq = tab),
                         freq = freq.Freq)
leu.other.freq <- leu.other.freq[,c(1,3)]

write.xlsx(leu.other.freq, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/leu.other.morph.codes.in.ds.kids.xlsx',
           row.names = FALSE)

rm(list = ls()); gc()
