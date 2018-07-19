
# Readme and TODO comments ------------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.10.
#' 
#' Generate tables for the DS-ALL birth defects project.
#' 
#' Table 1 gives summary statistics comparing non-DS ALL/AML, 
#' DS non-leukemia, and DS-leukemia groups with respect to sex birthweight,
#' gestational age, size for gestational age, and maternal and paternal 
#' ages.  
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#' TODO: Compute a size-for-gestational-age variable.

# Check distributions of table 1 variables --------------------------------

require(ggplot2)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')

varlist <- c('gest.age','m.age','f.age')

#' Histograms for continuous variables.
for (i in varlist){
  tmp <- floor(ds.leuk[,i])
  print(ggplot(data = ds.leuk) + 
          geom_histogram(aes(x = tmp), binwidth = 1) + 
          facet_wrap(~ds.leuk) +
          labs(x = i))
}

#' Birthweight with better binwidth.
print(ggplot(data = ds.leuk) + geom_histogram(aes(x = birth.wt)) + facet_wrap(~ds.leuk))

# Table 1 -----------------------------------------------------------------

require(gmodels); require(dplyr); require(xlsx)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180710.rdata')

#' For now, restrict to DS children.
ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML')

CrossTable(ds.leuk$sex, ds.leuk$ds.leuk, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE)

#' Report median and IQR for non-normally distributed variables.
varlist <- c('gest.age','m.age','f.age')

for (i in varlist){
  
  tmp <- data.frame(group = ds.leuk$ds.leuk, exposure = ds.leuk[,i])
  
  result1 <- rename(aggregate(exposure ~ group, data = tmp, median), median = exposure)
  result2 <- rename(aggregate(exposure ~ group, data = tmp, IQR), iqr = exposure)
  result <- left_join(result1, result2, by = 'group')
  
  write.xlsx(result, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/table1.v20180713.xlsx', 
             sheetName = i, row.names = FALSE, append = TRUE)
  
}

#' Significance tests.
for (i in varlist){
  
  tmp <- data.frame(group = ds.leuk$ds.leuk, exposure = ds.leuk[,i])
  
  print(i)
  
  print(kruskal.test(exposure ~ group, data = tmp))
  
}

#' Report mean and SD for normally distributed variables (actually variaBLE).
varlist <- 'birth.wt'

for (i in varlist){
  
  tmp <- data.frame(group = ds.leuk$ds.leuk, exposure = ds.leuk[,i])
  
  result1 <- rename(aggregate(exposure ~ group, data = tmp, mean), mean = exposure)
  result2 <- rename(aggregate(exposure ~ group, data = tmp, sd), sd = exposure)
  result <- left_join(result1, result2, by = 'group')
  
  write.xlsx(result, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/table1.v20180713.xlsx', 
             sheetName = i, row.names = FALSE, append = TRUE)
  
}

#' Significance tests.
for (i in varlist){
  
  tmp <- data.frame(group = ds.leuk$ds.leuk, exposure = ds.leuk[,i])
  
  print(i)
  
  print(kruskal.test(exposure ~ group, data = tmp))
  
}



# Table 1: only kids with DS ----------------------------------------------

require(gmodels); require(dplyr)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180710.rdata')

ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML')

CrossTable(ds.leuk$sex, ds.leuk$ds.leuk, prop.r = FALSE, prop.chisq = FALSE, chisq = TRUE)

#' Significance tests for non-normally distributed variables.
varlist <- c('gest.age','m.age','f.age')

for (i in varlist){
  
  tmp <- data.frame(group = ds.leuk$ds.leuk, exposure = ds.leuk[,i])
  
  print(i)
  
  print(kruskal.test(exposure ~ group, data = tmp))
  
}

#' Significance tests for normally distributed variables.
varlist <- 'birth.wt'

for (i in varlist){
  
  tmp <- data.frame(group = ds.leuk$ds.leuk, exposure = ds.leuk[,i])
  
  print(i)
  
  print(kruskal.test(exposure ~ group, data = tmp))
  
}




# Table 2: Birth defects burden in DS kids by leukemia status -------------

require(dplyr); require(gmodels)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')

#' This analysis will only involve kids with DS.
#' Given the problematic differences in birth defects prevalence between states, 
#' this analysis also temporarily only includes kids in Texas.
ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML' & state == 'TX')

CrossTable(ds.leuk$state, ds.leuk$ds.leuk, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE)

#' Compute a variable for whether any comorbid birth defect.
#' Will be a bad approach outside Texas where the minor defects are not being captured.
ds.leuk$any.comorbid.defect <- ifelse(ds.leuk$defect.total - 1 > 0, 1, 0)
CrossTable(ds.leuk$any.comorbid.defect, ds.leuk$ds.leuk, chisq = TRUE)

#' Existing organ system level variables - probably valid in Texas, probably not valid elsewhere.
table2.vars <- c(22,30,35,38,61,65,68,76,83,94)

#' New NBDPN-sourced organ system level variables for any MAJOR defect.  
#' Probably valid across all states with the exception of musculoskeletal.
#' Prevalence of craniosynostosis is too high in Texas, likely needs cleaning.
table2.vars.alt <- 159:167

for (i in table2.vars){
  print(names(ds.leuk[i]))
  print(CrossTable(ds.leuk[,i], ds.leuk$ds.leuk, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE, chisq = TRUE))
}

for (i in table2.vars.alt){
  print(names(ds.leuk[i]))
  print(CrossTable(ds.leuk[,i], ds.leuk$ds.leuk, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE, chisq = TRUE))
}

#' Parametric statistics and tests for continuous variables.
aggregate(defect.total -1 ~ ds.leuk, data = ds.leuk, mean)
aggregate(defect.total -1 ~ ds.leuk, data = ds.leuk, sd)
t.test(defect.total ~ ds.leuk, data = ds.leuk)

aggregate(majordefect.total -1 ~ ds.leuk, data = ds.leuk, mean)
aggregate(majordefect.total -1 ~ ds.leuk, data = ds.leuk, sd)
t.test(majordefect.total ~ ds.leuk, data = ds.leuk)



# Table 3: Specific major GI and CV defects by leukemia -------------------

require(dplyr); require(gmodels)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')

#' This analysis will only involve kids with DS.
#' Given the problematic differences in birth defects prevalence between states, 
#' this analysis also temporarily only includes kids in Texas.
ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML' & state == 'TX')

#' Overall frequency of digestive anomalies was higher in DS-leukemia group.
#' There is a suggestion of a greater frequency of CHD as well.
#' Do any specific major birth defects in these categories differ greatly in terms of prevalence?
gi <- 68:75
chd <- 38:60

tab.house <- data.frame(defect = as.character(),
                        count.ds.non.leu = as.numeric(),
                        prop.ds.non.leu = as.numeric(),
                        count.ds.leu = as.numeric(),
                        prop.ds.leu = as.numeric(),
                        p.chisq = as.numeric(),
                        p.fisher = as.numeric())

for (i in gi){
  
  tab <- CrossTable(ds.leuk[,i], ds.leuk$ds.leuk, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE, chisq = TRUE, fisher = TRUE)
  
  tab.out <- data.frame(defect = names(ds.leuk[i]),
                        count.ds.non.leu = tab$t[2,1],
                        prop.ds.non.leu = tab$prop.col[2,1],
                        count.ds.leu = tab$t[2,2],
                        prop.ds.leu = tab$prop.col[2,2],
                        p.chisq = tab$chisq$p.value,
                        p.fisher = tab$fisher.gt$p.value)
  
  tab.house <- rbind(tab.house, tab.out)
  
}

write.csv(tab.house, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/digestive.system.anomalies.by.leu.status.csv', row.names = FALSE)

tab.house <- data.frame(defect = as.character(),
                        count.ds.non.leu = as.numeric(),
                        prop.ds.non.leu = as.numeric(),
                        count.ds.leu = as.numeric(),
                        prop.ds.leu = as.numeric(),
                        p.chisq = as.numeric(),
                        p.fisher = as.numeric())

for (i in chd){
  
  tab <- CrossTable(ds.leuk[,i], ds.leuk$ds.leuk, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE, chisq = TRUE, fisher = TRUE)
  
  tab.out <- data.frame(defect = names(ds.leuk[i]),
                        count.ds.non.leu = tab$t[2,1],
                        prop.ds.non.leu = tab$prop.col[2,1],
                        count.ds.leu = tab$t[2,2],
                        prop.ds.leu = tab$prop.col[2,2],
                        p.chisq = tab$chisq$p.value,
                        p.fisher = tab$fisher.gt$p.value)
  
  tab.house <- rbind(tab.house, tab.out)
  
}

write.csv(tab.house, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/cardiovascular.anomalies.by.leu.status.csv', row.names = FALSE)



# Table 4: A closer look at other GI defects ------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.10.
#' 
#' Statistially speaking, the closest thing we have to a compelling result
#' right now is the higher likelihood of 'other digestive system anomalies' 
#' in DS children who develop leukemia.
#' 
#' Based on the prevalence of GI defects in this dataset (which can be 
#' obtained from the previously generated ds.leuk.comorbid.defect.counts.xlsx 
#' file), create new dummy variables for certain common GI defects.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr); require(xlsx); require(stringr); require(gmodels)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')
load('Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata')

#' Only dealing with DS cases diagnosed in TX for the moment.
ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML' & state == 'TX')

#' "Common" GI defects (say, >= 10 instances) not among those we've already computed variables for are macroglossia, high arched palate,
#' duodenal atresia or stenosis, other liver anomalies, other lip anomalies, anal atresia or stenosis, Hirschsprung's disease NOS, 
#' unspecified malrotation, annular pancreas, other palate anomalies.
#' One nice thing about the TX data is that I do not have to go through the trouble of dealing with the NC columns that have odd sixth digits.
#' To extend this code to NC, I'll need to map these additional columns to the diagnosis.
#' Read these into R.
codes <- read.xlsx(file = 'Z:/Jeremy/DS-ALL BD project/R outputs/ds.leuk.comorbid.defect.counts.xlsx', sheetIndex = 2, stringsAsFactors = FALSE,
                   colIndex = 1:2, rowIndex = 1:23, header = FALSE)

defect.codes <- codes$X1

#' Generate a character vector for column name assignment.
defect.names <- codes$X2
defect.names <- str_remove(defect.names, '[(].+[)]')
defect.names <- str_remove(defect.names, '[,]')
defect.names <- str_remove(defect.names, '[-]')
defect.names <- str_trim(defect.names, 'right')
defect.names <- tolower(defect.names)
defect.names <- str_replace_all(defect.names, ' ', '.')

#' Compute new columns.
for (i in 1:23){
  
  has.defect <- filter(ds.leuk.bd.codes.txnc.transpose, ds.leuk.bd.codes.txnc.transpose[, defect.codes[i]] == 1)
  
  ds.leuk[, defect.names[i]] <- ifelse(ds.leuk$studyid %in% has.defect$studyid, 1, 0)
  
}

#' Output table contents similar to table 3.
#' Have multipled prop.col by 100 to give percent with defect.
other.gi <- 168:190

tab.house <- data.frame(defect = as.character(),
                        count.ds.non.leu = as.numeric(),
                        pct.ds.non.leu = as.numeric(),
                        count.ds.leu = as.numeric(),
                        pct.ds.leu = as.numeric(),
                        p.chisq = as.numeric(),
                        p.fisher = as.numeric())

for (i in other.gi){
  
  tab <- CrossTable(ds.leuk[,i], ds.leuk$ds.leuk, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE, chisq = TRUE, fisher = TRUE)
  
  tab.out <- data.frame(defect = names(ds.leuk[i]),
                        count.ds.non.leu = tab$t[2,1],
                        pct.ds.non.leu = (tab$prop.col[2,1])*100,
                        count.ds.leu = tab$t[2,2],
                        pct.ds.leu = (tab$prop.col[2,2])*100,
                        p.chisq = tab$chisq$p.value,
                        p.fisher = tab$fisher.gt$p.value)
  
  tab.house <- rbind(tab.house, tab.out)
  
}

write.csv(tab.house, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/digestive.system.other.anomalies.by.leu.status.csv', row.names = FALSE)

rm(list = ls()); gc()


# Figure 1: Histogram of number of birth defects by leukemia status -------

require(ggplot2)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')

#' This analysis will only involve kids with DS.
#' Given the problematic differences in birth defects prevalence between states, 
#' this analysis also temporarily only includes kids in Texas.
ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML' & state == 'TX')

p <- ggplot(data = ds.leuk) + 
  geom_density(aes(x = number.of.comorbid.defects, fill = ds.leuk), alpha = 0.5) +
  labs(x = 'Number of Comorbid Birth Defects', y = 'Density') +
  guides(fill=guide_legend(title = NULL)) +
  geom_vline(aes(xintercept=9.5), linetype = 'dotdash') + 
  geom_vline(aes(xintercept = 10.4), linetype = 'solid') +
  theme_classic()
print(p)



# Figure 2: K-M curves by number of defects -------------------------------

require(survival); require(survminer); require(ggplot2)

load('Z:/Jeremy/DS-ALL BD project/Datasets/ds.leuk.v20180712.2.rdata')

#' This analysis will only involve kids with DS.
#' Given the problematic differences in birth defects prevalence between states, 
#' this analysis also temporarily only includes kids in Texas.
ds.leuk <- filter(ds.leuk, ds.leuk != 'Non-DS ALL/AML' & state == 'TX')

#' Identify BD cutpoints that could be used to break the sample in half or into tertiles.
quantile(ds.leuk$number.of.comorbid.defects, probs = c(0.33, 0.5, 0.67))

ds.leuk$defect.dichotomize <- factor(ifelse(ds.leuk$defect.total > 10, 1, 0),
                                     levels = c(0,1),
                                     labels = c('Median Number of Defects or Fewer', "Greater than Median Number of Defects"))

ds.leuk$defect.tertile <- factor(ifelse(ds.leuk$defect.total <= 7, 0, 
                                 ifelse(ds.leuk$defect.total >= 12, 2, 1)),
                          levels = c(0:2),
                          labels = c('First tertile', "Middle tertile", 'Upper tertile'))

#' Models for any leukemia, ALL and AML that can each be plotted using the code below.
fit <- survfit(Surv(person.yrs, leu.any) ~ defect.tertile, data = ds.leuk)
fit <- survfit(Surv(person.yrs, leu.any) ~ defect.dichotomize, data = ds.leuk)

new.plot <- ggsurvplot(fit, 
                       conf.int = TRUE, 
                       ylab = NULL,
                       ylim = c(0.9,1),
                       xlab = NULL,
                       font.tickslab = c(10, 'bold', 'black'),
                       linetype = 'strata', 
                       legend = "top",
                       legend.labs = c("< 10 Co-occurring defects",'10 or more co-occurring defects'))
print(new.plot)

#' A birthweight-adjusted Cox model for high vs. low number of defects.
dsleuk.surv <- data.frame(studyid = ds.leuk$studyid,<
                          time = ds.leuk$person.yrs, 
                          cancer = ds.leuk$leu.any, 
                          defect = ds.leuk$defect.dichotomize,
                          birth.wt = ds.leuk$birth.wt)

cox <- coxph(Surv(time, cancer) ~ defect + birth.wt, data = dsleuk.surv)
cox.coef <- as.data.frame(summary(cox)$coefficients)

estimates <- data.frame(var = rownames(cox.coef), 
                        hr = cox.coef$`exp(coef)`, 
                        ci.lower = exp(cox.coef$coef-(1.96*cox.coef$`se(coef)`)), 
                        ci.upper = exp(cox.coef$coef+(1.96*cox.coef$`se(coef)`)))

# Eye anomalies by state --------------------------------------------------

#' Including MI in this analysis has markedly reduced the overall prevalence of some organ system-level defects variables.
#' How different are the rates at which these are reported in each state?
prop.dx <- data.frame(defect = as.character(), pct.dx.mi = as.numeric(), pct.dx.nc = as.numeric(), pct.dx.tx = as.numeric())

for (i in table2.vars){
  
  prop <- as.numeric(round(CrossTable(ds.leuk[,i], ds.leuk$state, prop.r = FALSE, prop.chisq = FALSE, prop.t = FALSE, chisq = TRUE)$prop.col[2,] * 100, 2))
  prop <- data.frame(defect = names(ds.leuk[i]), pct.dx.mi = prop[1], pct.dx.nc = prop[2], pct.dx.tx = prop[3])
  prop.dx <- rbind(prop.dx, prop)
  
}

write.csv(prop.dx, file = 'Z:/Jeremy/DS-ALL BD project/R outputs/proportion.of.ds.cases.with.organ.system.defects.by.state.csv', row.names = FALSE)

#' Eye anomalies are an instructive example.
#' Prevalence is 4% in MI, 5% in NC and 71% in TX.
#' The prevalences of the named defects do not vary much.  However, TX kids are way more likely to have other eye defects == 1.
for (i in 31:34){
  ds.leuk[,i] <- ifelse(is.na(ds.leuk[,i]),0,ds.leuk[,i])
  print(names(ds.leuk[i]))
  print(CrossTable(ds.leuk[,i], ds.leuk$state, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE))
}

#' The prevalence of Brushfield spots of the iris in DS kids is something like 55% based on my reading.  The code for this is 743.800.
#' Pull the proportion of DS cases with this code by state.
load("Z:/Jeremy/DS-ALL BD project/Datasets/Expanded datasets/ds.leuk.bd.codes.txnc.transpose.v20180711.rdata")

which(names(ds.leuk.bd.codes.txnc.transpose) == '743.800')

tmp <- ds.leuk.bd.codes.txnc.transpose[, c(1,431:436)]

tmp <- left_join(tmp, select(ds.leuk, studyid, state), by = 'studyid')

for (i in 2:7){
  tmp[,i] <- ifelse(is.na(tmp[,i]), 0, tmp[,i])
  print(table(tmp[,i], tmp$state, useNA = 'ifany'))
}

tmp$`743.800` <- ifelse(is.na(tmp$`743.800`), 0, tmp$`743.800`)

#' Pretty much what I expected.
#' 68% of TX kids DX'd vs 3% of NC kids.
#' Even more infuriating, NC again has a number of codes that don't appear in the TX DHHS book that might be 
#' different forms of this 'other eye anomalies' diagnosis, but they're rarely if ever == 1.
table(tmp$`743.800`, tmp$state, useNA = 'ifany')



# Respiratory anomalies by state ------------------------------------------

for (i in 61:64){
  ds.leuk[,i] <- ifelse(is.na(ds.leuk[,i]),0,ds.leuk[,i])
  print(names(ds.leuk[i]))
  print(CrossTable(ds.leuk[,i], ds.leuk$state, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE))
}

which(names(ds.leuk.bd.codes.txnc.transpose) == '748.000')

tmp <- ds.leuk.bd.codes.txnc.transpose[, c(1,821:920)]

tmp <- left_join(tmp, select(ds.leuk, studyid, state), by = 'studyid')

for (i in 2:101){
  tmp[,i] <- ifelse(is.na(tmp[,i]), 0, tmp[,i])
  print(names(tmp[i]))
  print(table(tmp[,i], tmp$state, useNA = 'ifany'))
}
# Digestive anomalies by state --------------------------------------------

for (i in 68:75){
  ds.leuk[,i] <- ifelse(is.na(ds.leuk[,i]),0,ds.leuk[,i])
  print(names(ds.leuk[i]))
  print(CrossTable(ds.leuk[,i], ds.leuk$state, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE))
}

which(names(ds.leuk.bd.codes.txnc.transpose) == '750.000')
which(names(ds.leuk.bd.codes.txnc.transpose) == '751.900')

tmp <- ds.leuk.bd.codes.txnc.transpose[, c(1,969:1122)]

tmp <- left_join(tmp, select(ds.leuk, studyid, state), by = 'studyid')

differentially.reported.codes <- c()

for (i in 2:155){
  
  tmp[,i] <- ifelse(is.na(tmp[,i]), 0, tmp[,i])
  tab <- table(tmp[,i], tmp$state, useNA = 'ifany')
  
  if(dim(tab)[1] > 1){
      if(tab[2,2] > 10 & tab[2,2]/tab[2,1] > 3){
      new.diff.report.code <- names(tmp[i])
      differentially.reported.codes <- c(differentially.reported.codes, new.diff.report.code)
      }
      else{
        next
      }
  }
  
  else{
    next
  }
}
