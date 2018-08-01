# Readme and TODO comments ------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.08.01.
#' 
#' Updated from 2018.07.10.
#' 
#' Create and clean datasets for the analysis of cancer risk and survival
#' as a function of co-occurring birth defects in children diagnosed with 
#' Down syndrome.
#' 
#' Continued discussion and exploration of the data has led to further
#' revisions to our analytic strategy.
#' 
#' 1. We will restric the analysis to DS ALL and non-DS ALL cases.  This 
#'    was decided due to the thorny issue of the hepatomegaly association
#'    among AML children, likely secondary to TMD that was not reported 
#'    to/by the cancer registry.  
#' 2. We will compare ALL risk according to the presence/absence of NBDPN
#'    major birth defects.  We will compare burden of defects between 
#'    states on a relative basis, by computing the state-specific median
#'    number of birth defects among children with DS an no cancer, then 
#'    dichotomizing children from that state based on whether they are 
#'    above or below this figure.  This is done to skirt the problem of 
#'    differential reporting of co-occurring defects across states.  The 
#'    logic is that comparisons of the absolute number of defects will
#'    be confounded by the differences in the propensity for these defects
#'    to show up in the BD registries, but that comparisons children who have
#'    "more" vs. "fewer" defects based on the typical number of defects 
#'    in their state of birth will be less affected. In contrast to minor 
#'    defects, we feel that the NBDPN-defined major birth defects are 
#'    comparatively well captured by all registries involved in this project. 
#' 
#' This approach has the side effect of solving a few other issues: 
#' 
#' 1. Whether we consider TMD a cancer diagnosis or an exposure (previously 
#' an unresolved TODO).
#' 2. What sorts of leukemias should be in the non-DS cancer group for any 
#' analysis involving these kids (now there will be none). 
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' TODO: Should the NA values for minordefect.total be replaced with zeroes?
#' TODO: Clean MI birth defects codes and recalculate number of defects in this state.
#'       Include only defects captured by ICD9 codes that map to something in the BPA 
#'       system.  This should exclude the additional conditions MI collects data on.
#' TODO: Fix craniosynostosis in Texas.



# Generate the core analytic dataset --------------------------------------

require(dplyr)

load('Z:/Jeremy/GOBACK/Datasets/goback.v20180711.rdata')

#' Extract DS kids.  Exclude AR children.
ds.leuk <- filter(goback, down.syndrome == 1 & state != 'AR')

#' Exclude DS kids with any cancer except ALL as the first primary.
exclude.cancer1 <- filter(ds.leuk, cancer1 != 'all')
ds.leuk <- filter(ds.leuk, !(studyid %in% exclude.cancer1$studyid))

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Compute the relevant birth defects variables.
#' 
#' 1. Replace NAs with 0s for major defects.
#' 
#' 2. Compute state-specific median number of birth defects among DS kids
#'    without cancer.  Within each state, dichotomize children based on 
#'    whether they are above or below this number.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

for (i in 22:94){
  ds.leuk[,i] <- ifelse(is.na(ds.leuk[,i]), 0, ds.leuk[,i])
}

#' Compute dichotomous variables for number of total, major, and minor defects.
#' Note that because the median number of minor defects in MI and NC is 0, the variable is essentially
#' "any minor defect" in those states.
print(subset(aggregate(defect.total ~ state + all, data = ds.leuk, median), all == 0))
ds.leuk$defect.total.high.or.low <-  ifelse(ds.leuk$state == "TX" & ds.leuk$defect.total > 10, 1,
                                     ifelse(ds.leuk$state != 'TX' & ds.leuk$defect.total > 4, 1, 0))

print(subset(aggregate(minordefect.total ~ state + all, data = ds.leuk, median), all == 0))
ds.leuk$defect.minor.high.or.low <-  ifelse(ds.leuk$state == "TX" & ds.leuk$minordefect.total > 5, 1,
                                     ifelse(ds.leuk$state != 'TX' & ds.leuk$minordefect.total > 0, 1, 0))

print(subset(aggregate(majordefect.total ~ state + all, data = ds.leuk, median), all == 0))
ds.leuk$defect.major.high.or.low <-  ifelse(ds.leuk$state == "TX" & ds.leuk$majordefect.total > 4, 1,
                                            ifelse(ds.leuk$state != 'TX' & ds.leuk$majordefect.total > 3, 1, 0))

save(ds.leuk, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180801.1.rdata')

rm(list = ls()); gc()



# Compute organ system level dummy variables ------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' Compute the relevant birth defects variables.
#' 
#' 3. Compute a set of organ system level dummy variables that take value
#'    1 only if the child has one or more NBDPN-defined major birth defects
#'    in that organ system and 0 otherwise.
#' 
#' There are a few prerequisites to computing these, such as repairing 
#' craniosynostosis and computing defects for CPUV and cloacal exstrophy.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180801.1.rdata')
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata')
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.mi.transpose.v20180712.rdata')

#' Must compute variables for congenital posterior urethral valves and cloacal exstrophy to proceed.
cpuv.ids.txnc <- filter(ds.leuk.bd.codes.txnc.transpose, `753.600` == 1 | `753.604` == 1 | `753.608` == 1 | `753.609` == 1)
cpuv.ids.mi <- filter(ds.leuk.bd.codes.mi.transpose, `753.6` == 1)
cpuv.ids <- c(cpuv.ids.txnc$studyid, cpuv.ids.mi$studyid)

cloaca.ids.txnc <- filter(ds.leuk.bd.codes.txnc.transpose, `751.550` == 1 | `756.790` == 1)
cloaca.ids.mi <- filter(ds.leuk.bd.codes.mi.transpose, `751.5` == 1)
cloaca.ids <- c(cloaca.ids.mi$studyid, cloaca.ids.txnc$studyid)

ds.leuk$congenital.posterior.urethral.vales <- ifelse(ds.leuk$studyid %in% cpuv.ids, 1, 0)
ds.leuk$cloacal.exstrophy <- ifelse(ds.leuk$studyid %in% cloaca.ids, 1, 0)

#' There's an obvious overabundance of craniosynostosis diagnoses in TX.
#' This is evident from reviewing the statewise counts of comorbid defects.
#' These are the correct codes for craniosynostosis, including the additional NC versions.
#' Fixing this problem is a prerequisite to computing the new organ system variable for musculoskeletal defects.
cranio.codes <- c('756.000','756.005','756.006','756.010','756.020','756.030','756.011', '756.012', '756.013', 
                  '756.014', '756.018', '756.032', '756.034', '756.008')

ds.leuk.bd.codes.txnc.transpose <- select(ds.leuk.bd.codes.txnc.transpose, studyid, cranio.codes)

for (i in 2:15){
  ds.leuk.bd.codes.txnc.transpose[,i] <- ifelse(is.na(ds.leuk.bd.codes.txnc.transpose[,i]), 0, ds.leuk.bd.codes.txnc.transpose[,i])
}

ds.leuk.bd.codes.txnc.transpose$has.code <- ifelse(rowSums(ds.leuk.bd.codes.txnc.transpose[2:15]) >= 1, 1, 0)
ds.leuk.bd.codes.txnc.transpose <- filter(ds.leuk.bd.codes.txnc.transpose, has.code == 1)

ds.leuk$craniosynostosis <- ifelse(ds.leuk$studyid %in% ds.leuk.bd.codes.txnc.transpose$studyid, 1, 0)

#' Compute organ system level variables for any major defect.
#' These defects were chosen because they are the NBDPN reported birth defects.
#' See the 'NBDPN major birth defects from population-based birht defects surveillance programs in the united states' 
#' file in this project's data dictionaries directory.
ds.leuk$cns.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(23,24,26,28)]) >= 1, 1, 0)
ds.leuk$eye.any.major.anomaly <- ifelse(rowSums(ds.leuk[31:32]) >= 1, 1, 0)
ds.leuk$ear.any.major.anomaly <- ifelse(rowSums(ds.leuk[36]) >= 1, 1, 0)
ds.leuk$cardiovascular.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(49,55,43,48,40,52,46,47,51,56,41,44,42,57,54)]) >= 1, 1, 0)
ds.leuk$orofacial.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(63,66,67)]) >= 1, 1, 0)
ds.leuk$gi.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(73,69,70,74)]) >= 1, 1, 0)
ds.leuk$genitourinary.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(78,80,77,166,167)]) >= 1, 1, 0)
ds.leuk$musculoskeletal.any.major.anomaly <- ifelse(rowSums(ds.leuk[c(84:86,87,88,90:92)] >= 1), 1, 0)
ds.leuk$any.major.comorbid.anomaly <- ifelse(rowSums(ds.leuk[162:169]) >= 1, 1, 0)

#' There are lots of columns kicking around that I don't see myself using.
#' The old organ system level variables, the chromosomal conditions other than DS 
#' and the cancer variables other than ALL come to mind.
cols <- c(which(grepl('conganomalies.', names(ds.leuk)) == TRUE), 65, 95:106, 113:151)

ds.leuk <- ds.leuk[, -cols]

save(ds.leuk, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180801.2.rdata')

rm(list = ls()); gc()



# Generate BD  (TX and NC) and cancer (all states) datasets ---------------

require(dplyr)

#' Need a clean environment for the references to ls() below to work.
rm(list = ls())

#' As of now, only TX and NC children in the ds.leuk dataset.
load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180801.2.rdata')
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/cancer.codes.v20180227.1.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.transpose.v20180614.rdata")
load("Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.txnc.v20180606.rdata")

data.frames.names <- paste0('ds.leuk.',ls()[1:3])

data.frames <- list(bd.codes.txnc, 
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
  
  path <- paste0('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/',data.frames.names[i],'.v20180801.rdata')
  
  save(list = data.frames.names[i], file = path)
  
}

rm(list = ls()); gc()



# Generate birth defects codes data frames: MI ----------------------------

require(dplyr)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180801.2.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.v20180712.rdata')
load('Z:/Jeremy/GOBACK/Datasets/Expanded datasets/bd.codes.mi.transpose.v20180712.rdata')

ds.leuk <- filter(ds.leuk, state == 'MI')

#' Results in 3102 of 3102 rows returned.
ds.leuk.bd.codes.mi <- filter(bd.codes.mi, studyid %in% ds.leuk$studyid)
ds.leuk.bd.codes.mi.transpose <- filter(bd.codes.mi.transpose, studyid %in% ds.leuk$studyid)

save(ds.leuk.bd.codes.mi, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.mi.v20180801.rdata')
save(ds.leuk.bd.codes.mi.transpose, file = 'Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.mi.transpose.v20180801.rdata')

rm(list = ls()); gc()


