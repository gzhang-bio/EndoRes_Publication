### library
# library(coin)
library(rstatix)

### read in the simple clinical data
all.anno.simple <- read.delim("all_anno_simple.txt", sep = "\t", check.names = FALSE)

all.anno.simple$NR_R <- factor(all.anno.simple$NR_R, levels = c("R", "NR"))
all.anno.simple$Histology_type_baseline <- factor(all.anno.simple$Histology_type_baseline, levels = c("NST", "ILBC", "Others"))
all.anno.simple$`Histology_type_post-pET` <- factor(all.anno.simple$`Histology_type_post-pET`, levels = c("NST", "ILBC", "Others"))

### simplify the values of pN
all.anno.simple$pN[all.anno.simple$pN == "2"] <- 1

p.dis <- NULL
p.dis.match <- NULL
p.dis.tam <- NULL
p.dis.tam.match <- NULL
p.dis.ai <- NULL
p.dis.ai.match <- NULL
p.dis.ai.bp <- NULL
p.val <- NULL
p.val.tam <- NULL
p.val.ai <- NULL


p.dis.tam.r.bp <- NULL
p.dis.tam.nr.bp <- NULL
p.dis.ai.r.bp <- NULL
p.dis.ai.nr.bp <- NULL
p.val.tam.r.bp <- NULL
p.val.tam.nr.bp <- NULL
p.val.ai.r.bp <- NULL
p.val.ai.nr.bp <- NULL

p.signif.dis <- NULL
p.signif.dis.match <- NULL
p.signif.dis.tam <- NULL
p.signif.dis.tam.match <- NULL
p.signif.dis.ai <- NULL
p.signif.dis.ai.match <- NULL
p.signif.val <- NULL
p.signif.val.tam <- NULL
p.signif.val.ai <- NULL

p.signif.dis.tam.r.bp <- NULL
p.signif.dis.tam.nr.bp <- NULL
p.signif.dis.ai.r.bp <- NULL
p.signif.dis.ai.nr.bp <- NULL
p.signif.val.tam.r.bp <- NULL
p.signif.val.tam.nr.bp <- NULL
p.signif.val.ai.nr.bp <- NULL
p.signif.val.ai.r.bp <- NULL

table.dis <- list()
table.dis.match <- list()
table.dis.tam <- list()
table.dis.tam.match <- list()
table.dis.ai <- list()
table.dis.ai.match <- list()
table.val <- list()
table.val.tam <- list()
table.val.ai <- list()

table.dis.tam.r.bp <- list()
table.dis.tam.nr.bp <- list()
table.dis.ai.r.bp <- list()
table.dis.ai.nr.bp <- list()
table.val.tam.r.bp <- list()
table.val.tam.nr.bp <- list()
table.val.ai.r.bp <- list()
table.val.ai.nr.bp <- list()

prop.table.dis <- list()
prop.table.dis.match <- list()
prop.table.dis.tam <- list()
prop.table.dis.tam.match <- list()
prop.table.dis.ai <- list()
prop.table.dis.ai.match <- list()
prop.table.val <- list()
prop.table.val.tam <- list()
prop.table.val.ai <- list()

prop.table.dis.tam.r.bp <- list()
prop.table.dis.tam.nr.bp <- list()
prop.table.dis.ai.r.bp <- list()
prop.table.dis.ai.nr.bp <- list()
prop.table.val.tam.r.bp <- list()
prop.table.val.tam.nr.bp <- list()
prop.table.val.ai.r.bp <- list()
prop.table.val.ai.nr.bp <- list()

features_all <- c("Age", 
                  "Histology_type_baseline", "Histology_type_post-pET", "Histology_type_baseline_post-pET", 
                  "pT", "pN",          
                  "Histology_grade_baseline", "Histology_grade_post-pET", "Histology_grade_baseline_post-pET",
                  "ER_percent_baseline", "ER_percent_post-pET", "ER_percent_baseline_post-pET",
                  "PR_percent_baseline", "PR_percent_post-pET", "PR_percent_baseline_post-pET",
                  "ER_PR_percent_baseline", "ER_PR_percent_post-pET", 
                  "HER2_IHC_Score_baseline", "HER2_IHC_Score_post-pET", "HER2_IHC_Score_baseline_post-pET",
                  "E_Cadherin_Status_baseline", "E_Cadherin_Status_post-pET", "E_Cadherin_Status_baseline_post-pET", 
                  "Ki67_baseline", "Ki67_post-pET", "Ki67_baseline_post-pET", 
                  "RS_Score_baseline")

features_bp <- c("Histology_type_baseline_post-pET",
                 "Histology_grade_baseline_post-pET",
                 "ER_percent_baseline_post-pET",
                 "PR_percent_baseline_post-pET",
                 "HER2_IHC_Score_baseline_post-pET",
                 "E_Cadherin_Status_baseline_post-pET",
                 "Ki67_baseline_post-pET")

features <- c()
### create the crorrect tables according to the feature, cohort and match
func_table <- function(table, match, feature) {
  condition <- paste(match, feature, sep = "_")
  switch(
    condition,
    "nomatch_HER2_IHC_Score_baseline" = {table <- table[c("IHC 0/1+", "IHC 2+, FISH-neg"), c("R", "NR")]},
    "nomatch_HER2_IHC_Score_post-pET" = {table <- table[c("IHC 0/1+", "IHC 2+, FISH-neg"), c("R", "NR")]},
    "nomatch_E_Cadherin_Status_baseline" = {table <- table[c("pos.", "neg."), c("R", "NR")]},
    "nomatch_E_Cadherin_Status_post-pET" = {table <- table[c("pos.", "neg."), c("R", "NR")]},
    "match_HER2_IHC_Score_baseline" = {table <- table[c("IHC 0/1+", "IHC 2+, FISH-neg"), c("IHC 0/1+", "IHC 2+, FISH-neg")]},
    "match_HER2_IHC_Score_post-pET" = {table <- table[c("IHC 0/1+", "IHC 2+, FISH-neg"), c("IHC 0/1+", "IHC 2+, FISH-neg")]},
    "match_HER2_IHC_Score_baseline_post-pET" = {table <- table[c("IHC 0/1+", "IHC 2+, FISH-neg"), c("IHC 0/1+", "IHC 2+, FISH-neg")]},
    "match_E_Cadherin_Status_baseline" = {table <- table[c("pos.", "neg."), c("pos.", "neg.")]},
    "match_E_Cadherin_Status_post-pET" = {table <- table[c("pos.", "neg."), c("pos.", "neg.")]},
    "match_E_Cadherin_Status_baseline_post-pET" = {table <- table[c("pos.", "neg."), c("pos.", "neg.")]}
  )
  return(table)
}

for(feature in features_all) {
  # feature <- "Histology_grade_baseline_post-pET"
  # feature <- "Ki67_post-pET"
  # feature <- "Ki67_baseline"
  # print(feature)
  
  for (cohort in c(1, 2)){
    # cohort <- 1
    switch (
      feature,
      "Age" = {
        data <- ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "Age"] < 60, "less than 60", "more equal 60")
        data <- factor(data, levels = c("less than 60", "more equal 60"))},
      "Histology_type_baseline" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, "Histology_type_baseline"]},
      "Histology_type_post-pET" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, "Histology_type_post-pET"]},
      "Histology_type_baseline_post-pET" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, c("Histology_type_baseline", "Histology_type_post-pET")]},
      "pT" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, "pT"]},
      "pN" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, "pN"]},
      "Histology_grade_baseline" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, "Histology_grade_baseline"]
        data <- factor(data, levels = c(1, 2, 3))},
      "Histology_grade_post-pET" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, "Histology_grade_post-pET"]
        data <- factor(data, levels = c(1, 2, 3))},
      "Histology_grade_baseline_post-pET" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, c("Histology_grade_baseline", "Histology_grade_post-pET")]
        colnames(data) <- c("Histology_grade_baseline", "Histology_grade_post-pET")
        data[, 1] <- factor(data[, 1], levels = c(1, 2, 3))
        data[, 2] <- factor(data[, 2], levels = c(1, 2, 3))},
      "ER_percent_baseline" = {
        data <- ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_baseline"] > 10, "more than 10", "less equal 10")
        data <- factor(data, levels = c("more than 10", "less equal 10"))}, 
      "ER_percent_post-pET" = {
        data <- ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_post-pET"] > 10, "more than 10", "less equal 10")
        data <- factor(data, levels = c("more than 10", "less equal 10"))},
      "ER_percent_baseline_post-pET" = {
        data <- data.frame(ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_baseline"] > 10, "more than 10", "less equal 10"),
                           ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_post-pET"] > 10, "more than 10", "less equal 10"))
        colnames(data) <- c("ER_percent_baseline", "ER_percent_post-pET") 
        data[, 1] <- factor(data[, 1], levels = c("more than 10", "less equal 10"))
        data[, 2] <- factor(data[, 2], levels = c("more than 10", "less equal 10"))},
      "PR_percent_baseline" = {
        data <- ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_baseline"] > 10, "more than 10", "less equal 10")
        data <- factor(data, levels = c("more than 10", "less equal 10"))},
      "PR_percent_post-pET" = {
        data <- ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_post-pET"] > 10, "more than 10", "less equal 10")
        data <- factor(data, levels = c("more than 10", "less equal 10"))},
      "PR_percent_baseline_post-pET" = {
        data <- data.frame(ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_baseline"] > 10, "more than 10", "less equal 10"),
                           ifelse(all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_post-pET"] > 10, "more than 10", "less equal 10"))
        colnames(data) <- c("PR_percent_baseline", "PR_percent_post-pET")
        data[, 1] <- factor(data[, 1], levels = c("more than 10", "less equal 10"))
        data[, 2] <- factor(data[, 2], levels = c("more than 10", "less equal 10"))},
      "ER_PR_percent_baseline" = { 
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_baseline"] > 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_baseline"] > 10] <- "ER-pos; PR-pos"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_baseline"] > 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_baseline"] <= 10] <- "ER-pos; PR-low"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_baseline"] <= 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_baseline"] > 10] <- "ER-low; PR-pos"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_baseline"]  <= 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_baseline"] <= 10] <- "ER-low; PR-low"
        data[is.na(all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_baseline"]) &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_baseline"] > 10] <- "ER-na; PR-pos"
        data <- factor(data, levels = c("ER-pos; PR-pos", "ER-pos; PR-low",
                                        "ER-low; PR-pos", "ER-low; PR-low",
                                        "ER-na; PR-pos"))},
      "ER_PR_percent_post-pET" = {
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_post-pET"] > 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_post-pET"] > 10] <- "ER-pos; PR-pos"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_post-pET"] > 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_post-pET"] <= 10] <- "ER-pos; PR-low"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_post-pET"] <= 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_post-pET"] > 10] <- "ER-low; PR-pos"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_post-pET"]  <= 10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_post-pET"] <= 10] <- "ER-low; PR-low"
        data[is.na(all.anno.simple[all.anno.simple$Cohort == cohort, "ER_percent_post-pET"]) &
               all.anno.simple[all.anno.simple$Cohort == cohort, "PR_percent_post-pET"] > 10] <- "ER-na; PR-pos"
        data <- factor(data, levels = c("ER-pos; PR-pos", "ER-pos; PR-low",
                                        "ER-low; PR-pos", "ER-low; PR-low",
                                        "ER-na; PR-pos"))},
      "HER2_IHC_Score_baseline" = {
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] %in% c(0, 1)] <- "IHC 0/1+"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_baseline"] == "0"] <- "IHC 2+, FISH-neg"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_baseline"] == "n.a."] <- "IHC 2+, FISH-n.a."
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_baseline"] == "1"] <- "IHC 2+, FISH-pos"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 3] <- "IHC 3+"
        data <- factor(data, levels = c("IHC 0/1+", "IHC 2+, FISH-neg",
                                        "IHC 2+, FISH-n.a.", "IHC 2+, FISH-pos", "IHC 3+"))},
      "HER2_IHC_Score_post-pET" = {
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] %in% c(0, 1)] <- "IHC 0/1+"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_post-pET"] == "0"] <- "IHC 2+, FISH-neg"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_post-pET"] == "n.a."] <- "IHC 2+, FISH-n.a."
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_post-pET"] == "1"] <- "IHC 2+, FISH-pos"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 3] <- "IHC 3+"
        data <- factor(data, levels = c("IHC 0/1+", "IHC 2+, FISH-neg",
                                        "IHC 2+, FISH-n.a.", "IHC 2+, FISH-pos", "IHC 3+"))},
      "HER2_IHC_Score_baseline_post-pET" = {
        data1 <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] %in% c(0, 1)] <- "IHC 0/1+"
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_baseline"] == "0"] <- "IHC 2+, FISH-neg"
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_baseline"] == "n.a."] <- "IHC 2+, FISH-n.a."
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 2 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_baseline"] == "1"] <- "IHC 2+, FISH-pos"
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_baseline"] == 3] <- "IHC 3+"
        
        data2 <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] %in% c(0, 1)] <- "IHC 0/1+"
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 2 &
                all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_post-pET"] == "0"] <- "IHC 2+, FISH-neg"
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 2 &
                all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_post-pET"] == "n.a."] <- "IHC 2+, FISH-n.a."
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 2 &
                all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_status_post-pET"] == "1"] <- "IHC 2+, FISH-pos"
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "HER2_IHC_Score_post-pET"] == 3] <- "IHC 3+"
        
        data <- data.frame(data1, data2)
        colnames(data) <- c("HER2_IHC_Score_baseline", "HER2_IHC_Score_post-pET")
        
        data[, 1] <- factor(data[, 1], levels = c("IHC 0/1+", "IHC 2+, FISH-neg",
                                                  "IHC 2+, FISH-n.a.", "IHC 2+, FISH-pos", "IHC 3+"))
        data[, 2] <- factor(data[, 2], levels = c("IHC 0/1+", "IHC 2+, FISH-neg",
                                                  "IHC 2+, FISH-n.a.", "IHC 2+, FISH-pos", "IHC 3+"))},
      "E_Cadherin_Status_baseline" = {
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_baseline"] %in% c("0.5", "1")] <- "pos."
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_baseline"] %in% c("0")] <- "neg."
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_baseline"] %in% c("n.a.")] <- "n.a."
        data <- factor(data, levels = c("pos.", "neg.", "n.a."))},
      "E_Cadherin_Status_post-pET" = {
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_post-pET"] %in% c("0.5", "1")] <- "pos."
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_post-pET"] %in% c("0")] <- "neg."
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_post-pET"] %in% c("n.a.")] <- "n.a."
        data <- factor(data, levels = c("pos.", "neg.", "n.a."))},
      "E_Cadherin_Status_baseline_post-pET" = {
        data1 <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_baseline"] %in% c("0.5", "1")] <- "pos."
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_baseline"] %in% c("0")] <- "neg."
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_baseline"] %in% c("n.a.")] <- "n.a."
        
        data2 <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_post-pET"] %in% c("0.5", "1")] <- "pos."
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_post-pET"] %in% c("0")] <- "neg."
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "E_Cadherin_Status_post-pET"] %in% c("n.a.")] <- "n.a."
        
        data <- data.frame(data1, data2)
        colnames(data) <- c("E_Cadherin_Status_baseline", "E_Cadherin_Status_post-pET")
        data[, 1] <- factor(data[, 1], levels = c("pos.", "neg.", "n.a."))
        data[, 2] <- factor(data[, 2], levels = c("pos.", "neg.", "n.a."))},
      "Ki67_baseline" = {
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >= 0 & 
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <= 9] <- "0-9"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >=10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <=19] <- "10-19"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >=20 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <=34] <- "20-34"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >=35 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <=100] <- "35-100"
        data <- factor(data, levels = c("0-9", "10-19", "20-34", "35-100"))},
      "Ki67_post-pET" = {
        data <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >= 0 & 
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <= 9] <- "0-9"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >=10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <=19] <- "10-19"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >=20 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <=34] <- "20-34"
        data[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >=35 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <=100] <- "35-100"
        data <- factor(data, levels = c("0-9", "10-19", "20-34", "35-100"))},
      "Ki67_baseline_post-pET" = {
        data1 <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >= 0 & 
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <= 9] <- "0-9"
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >=10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <=19] <- "10-19"
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >=20 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <=34] <- "20-34"
        data1[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] >=35 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_baseline"] <=100] <- "35-100"
        
        data2 <- rep(NA, nrow(all.anno.simple[all.anno.simple$Cohort == cohort,]))
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >= 0 & 
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <= 9] <- "0-9"
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >=10 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <=19] <- "10-19"
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >=20 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <=34] <- "20-34"
        data2[all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] >=35 &
               all.anno.simple[all.anno.simple$Cohort == cohort, "Ki67_post-pET"] <=100] <- "35-100"
        
        data <- data.frame(data1, data2)
        colnames(data) <- c("Ki67_baseline", "Ki67_post-pET")
        data[, 1] <- factor(data[, 1], levels = c("0-9", "10-19", "20-34", "35-100"))
        data[, 2] <- factor(data[, 2], levels = c("0-9", "10-19", "20-34", "35-100"))},
      "RS_Score_baseline" = {
        data <- all.anno.simple[all.anno.simple$Cohort == cohort, "RS_Score_baseline"]}
      ) 
  
    if (cohort == 1) {
      if (feature %in% features_bp) {
        ### only for baseline/post-pET matching comparing
        # table <- table(data[, 1], data[, 2])
        # table.dis.bp <- append(table.dis.bp, list(table))
        # prop.table.dis.bp <- append(prop.table.dis.bp, list(prop.table(table, 2)))
        # p.dis.bp <- append(p.dis.bp, mcnemar_test(func_table(table, "match", feature))$p)
        # p.signif.dis.bp <- append(p.signif.dis.bp, mcnemar_test(func_table(table, feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 2])
        table.dis.tam.r.bp <- append(table.dis.tam.r.bp, list(table))
        prop.table.dis.tam.r.bp <- append(prop.table.dis.tam.r.bp, list(prop.table(table, 2)))
        p.dis.tam.r.bp <- append(p.dis.tam.r.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.dis.tam.r.bp <- append(p.signif.dis.tam.r.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 2])
        table.dis.tam.nr.bp <- append(table.dis.tam.nr.bp, list(table))
        prop.table.dis.tam.nr.bp <- append(prop.table.dis.tam.nr.bp, list(prop.table(table, 2)))
        p.dis.tam.nr.bp <- append(p.dis.tam.nr.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.dis.tam.nr.bp <- append(p.signif.dis.tam.nr.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 2])
        table.dis.ai.r.bp <- append(table.dis.ai.r.bp, list(table))
        prop.table.dis.ai.r.bp <- append(prop.table.dis.ai.r.bp, list(prop.table(table, 2)))
        p.dis.ai.r.bp <- append(p.dis.ai.r.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.dis.ai.r.bp <- append(p.signif.dis.ai.r.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 2])
        table.dis.ai.nr.bp <- append(table.dis.ai.nr.bp, list(table))
        prop.table.dis.ai.nr.bp <- append(prop.table.dis.ai.nr.bp, list(prop.table(table, 2)))
        p.dis.ai.nr.bp <- append(p.dis.ai.nr.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.dis.ai.nr.bp <- append(p.signif.dis.ai.nr.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
      } else {
        features <- append(features, feature)
        ### matching
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR"],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R"])
        table.dis.match <- append(table.dis.match, list(table))
        prop.table.dis.match <- append(prop.table.dis.match, list(prop.table(table, 2)))
      
        p.dis.match <- append(p.dis.match, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.dis.match <- append(p.signif.dis.match, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR" & all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R" & all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"])
        table.dis.tam.match <- append(table.dis.tam.match, list(table))
        prop.table.dis.tam.match <- append(prop.table.dis.tam.match, list(prop.table(table, 2)))
        
        if (feature == "Ki67_post-pET") {
          p.dis.tam.match <- append(p.dis.tam.match, mcnemar_test(func_table(table[3:4, 1:2], "match", feature))$p)
          p.signif.dis.tam.match <- append(p.signif.dis.tam.match, mcnemar_test(func_table(table[3:4, 1:2], "match", feature))$p.signif) }
        else {
          p.dis.tam.match <- append(p.dis.tam.match, mcnemar_test(func_table(table, "match", feature))$p)
          p.signif.dis.tam.match <- append(p.signif.dis.tam.match, mcnemar_test(func_table(table, "match", feature))$p.signif)
        }
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR" & all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R" & all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"])
        table.dis.ai.match <- append(table.dis.ai.match, list(table))
        prop.table.dis.ai.match <- append(prop.table.dis.ai.match, list(prop.table(table, 2)))
        
        if (feature == "Ki67_post-pET") {
          p.dis.ai.match <- append(p.dis.ai.match, mcnemar_test(func_table(table[3:4, 1:2], "match", feature))$p)
          p.signif.dis.ai.match <- append(p.signif.dis.ai.match, mcnemar_test(func_table(table[3:4, 1:2], "match", feature))$p.signif) }
        else {
          p.dis.ai.match <- append(p.dis.ai.match, mcnemar_test(func_table(table, "match", feature))$p)
          p.signif.dis.ai.match <- append(p.signif.dis.ai.match, mcnemar_test(func_table(table, "match", feature))$p.signif)
        }
        
        ### non-matching
        table <- table(data, 
                       all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"])
        table.dis <- append(table.dis, list(table))
        prop.table.dis <- append(prop.table.dis, list(prop.table(table, 2)))
        
        if (feature %in% c("Ki67_baseline", "Ki67_post-pET")) {
          p.dis <- append(p.dis, chisq_test(table)$p)
          p.signif.dis <- append(p.signif.dis, chisq_test(func_table(table, "nomatch", feature))$p.signif)}
        else {
          p.dis <- append(p.dis, fisher_test(table)$p)
          p.signif.dis <- append(p.signif.dis, fisher_test(func_table(table, "nomatch", feature))$p.signif)}
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"],
                       all.anno.simple[all.anno.simple$Cohort == cohort & all.anno.simple$Treatment == "TAM", "NR_R"])
        table.dis.tam <- append(table.dis.tam, list(table))
        prop.table.dis.tam <- append(prop.table.dis.tam, list(prop.table(table, 2)))
        
        if (feature %in% c("Ki67_baseline", "Ki67_post-pET")) {
          p.dis.tam <- append(p.dis.tam, chisq_test(table)$p)
          p.signif.dis.tam <- append(p.signif.dis.tam, chisq_test(func_table(table, "nomatch", feature))$p.signif)}
        else {
          p.dis.tam <- append(p.dis.tam, fisher_test(table)$p)
          p.signif.dis.tam <- append(p.signif.dis.tam, fisher_test(func_table(table, "nomatch", feature))$p.signif)}
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"],
                       all.anno.simple[all.anno.simple$Cohort == cohort & all.anno.simple$Treatment == "AI", "NR_R"])
        table.dis.ai <- append(table.dis.ai, list(table))
        prop.table.dis.ai <- append(prop.table.dis.ai, list(prop.table(table, 2)))
        
        if (feature %in% c("Ki67_baseline", "Ki67_post-pET")) {
          p.dis.ai <- append(p.dis.ai, chisq_test(table)$p)
          p.signif.dis.ai <- append(p.signif.dis.ai, chisq_test(func_table(table, "nomatch", feature))$p.signif)}
        else {
          p.dis.ai <- append(p.dis.ai, fisher_test(table)$p)
          p.signif.dis.ai <- append(p.signif.dis.ai, fisher_test(func_table(table, "nomatch", feature))$p.signif)}
      }
    } 
    else if(cohort == 2) {
      if (feature %in% features_bp) {
        ### only for baseline/post-pET matching comparing
        # table <- table(data[, 1], data[, 2])
        # table.val.bp <- append(table.val.bp, list(table))
        # prop.table.val.bp <- append(prop.table.val.bp, list(prop.table(table, 2)))
        # p.val.bp <- append(p.val.bp, mcnemar_test(func_table(table, "match", feature))$p)
        # p.signif.val.bp <- append(p.signif.val.bp, mcnemar_test(func_table(table, feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 2])
        table.val.tam.r.bp <- append(table.val.tam.r.bp, list(table))
        prop.table.val.tam.r.bp <- append(prop.table.val.tam.r.bp, list(prop.table(table, 2)))
        p.val.tam.r.bp <- append(p.val.tam.r.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.val.tam.r.bp <- append(p.signif.val.tam.r.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 2])
        table.val.tam.nr.bp <- append(table.val.tam.nr.bp, list(table))
        prop.table.val.tam.nr.bp <- append(prop.table.val.tam.nr.bp, list(prop.table(table, 2)))
        p.val.tam.nr.bp <- append(p.val.tam.nr.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.val.tam.nr.bp <- append(p.signif.val.tam.nr.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "R", 2])
        table.val.ai.r.bp <- append(table.val.ai.r.bp, list(table))
        prop.table.val.ai.r.bp <- append(prop.table.val.ai.r.bp, list(prop.table(table, 2)))
        p.val.ai.r.bp <- append(p.val.ai.r.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.val.ai.r.bp <- append(p.signif.val.ai.r.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI" 
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 1],
                       data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"
                            & all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"] == "NR", 2])
        table.val.ai.nr.bp <- append(table.val.ai.nr.bp, list(table))
        prop.table.val.ai.nr.bp <- append(prop.table.val.ai.nr.bp, list(prop.table(table, 2)))
        p.val.ai.nr.bp <- append(p.val.ai.nr.bp, mcnemar_test(func_table(table, "match", feature))$p)
        p.signif.val.ai.nr.bp <- append(p.signif.val.ai.nr.bp, mcnemar_test(func_table(table, "match", feature))$p.signif)
        
      } else {
        table <- table(data, 
                       all.anno.simple[all.anno.simple$Cohort == cohort, "NR_R"])
        table.val <- append(table.val, list(table))
        prop.table.val <- append(prop.table.val, list(prop.table(table, 2)))
        
        if (feature %in% c("Ki67_baseline", "Ki67_post-pET")) {
          p.val <- append(p.val, chisq_test(table)$p)
          p.signif.val <- append(p.signif.val, chisq_test(func_table(table, "nomatch", feature))$p.signif)}
        else {
          p.val <- append(p.val, fisher_test(table)$p)
          p.signif.val <- append(p.signif.val, fisher_test(func_table(table, "nomatch", feature))$p.signif)}
        
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "TAM"],
                       all.anno.simple[all.anno.simple$Cohort == cohort & all.anno.simple$Treatment == "TAM", "NR_R"])
        table.val.tam <- append(table.val.tam, list(table))
        prop.table.val.tam <- append(prop.table.val.tam, list(prop.table(table, 2)))
        
        if (feature %in% c("Ki67_baseline", "Ki67_post-pET")) {
          p.val.tam <- append(p.val.tam, chisq_test(table)$p)
          p.signif.val.tam <- append(p.signif.val.tam, chisq_test(func_table(table, "nomatch", feature))$p.signif)}
        else {
          p.val.tam <- append(p.val.tam, fisher_test(table)$p)
          p.signif.val.tam <- append(p.signif.val.tam, fisher_test(func_table(table, "nomatch", feature))$p.signif)}
  
        table <- table(data[all.anno.simple[all.anno.simple$Cohort == cohort, "Treatment"] == "AI"],
                       all.anno.simple[all.anno.simple$Cohort == cohort & all.anno.simple$Treatment == "AI", "NR_R"])
        table.val.ai <- append(table.val.ai, list(table))
        prop.table.val.ai <- append(prop.table.val.ai, list(prop.table(table, 2)))
        
        if (feature %in% c("Ki67_baseline", "Ki67_post-pET")) {
          p.val.ai <- append(p.val.ai, chisq_test(table)$p)
          p.signif.val.ai <- append(p.signif.val.ai, chisq_test(func_table(table, "nomatch", feature))$p.signif)}
        else {
          p.val.ai <- append(p.val.ai, fisher_test(table)$p)
          p.signif.val.ai <- append(p.signif.val.ai, fisher_test(func_table(table, "nomatch", feature))$p.signif)}
      }
    }
  }
}

names(p.dis) <- features
names(p.dis.match) <- features
names(p.dis.tam) <- features
names(p.dis.tam.match) <- features
names(p.dis.ai) <- features
names(p.dis.ai.match) <- features
names(p.val) <- features
names(p.val.tam) <- features
names(p.val.ai) <- features

names(p.dis.tam.nr.bp) <- features_bp
names(p.dis.tam.r.bp) <- features_bp
names(p.dis.ai.nr.bp) <- features_bp
names(p.dis.ai.r.bp) <- features_bp
names(p.val.tam.nr.bp) <- features_bp
names(p.val.tam.r.bp) <- features_bp
names(p.val.ai.nr.bp) <- features_bp
names(p.val.ai.r.bp) <- features_bp

names(p.signif.dis.match) <- features
names(p.signif.dis.tam.match) <- features
names(p.signif.dis.ai.match) <- features
names(p.signif.dis) <- features
names(p.signif.dis.tam) <- features
names(p.signif.dis.ai) <- features
names(p.signif.val) <- features
names(p.signif.val.tam) <- features
names(p.signif.val.ai) <- features

names(p.signif.dis.tam.nr.bp) <- features_bp
names(p.signif.dis.tam.r.bp) <- features_bp
names(p.signif.dis.ai.nr.bp) <- features_bp
names(p.signif.dis.ai.r.bp) <- features_bp
names(p.signif.val.tam.nr.bp) <- features_bp
names(p.signif.val.tam.r.bp) <- features_bp
names(p.signif.val.ai.nr.bp) <- features_bp
names(p.signif.val.ai.r.bp) <- features_bp

names(table.dis) <- features
names(table.dis.match) <- features
names(table.dis.tam) <- features
names(table.dis.tam.match) <- features
names(table.dis.ai) <- features
names(table.dis.ai.match) <- features
names(table.val) <- features
names(table.val.tam) <- features
names(table.val.ai) <- features

names(table.dis.tam.nr.bp) <- features_bp
names(table.dis.tam.r.bp) <- features_bp
names(table.dis.ai.nr.bp) <- features_bp
names(table.dis.ai.r.bp) <- features_bp
names(table.val.tam.nr.bp) <- features_bp
names(table.val.tam.r.bp) <- features_bp
names(table.val.ai.nr.bp) <- features_bp
names(table.val.ai.r.bp) <- features_bp

names(prop.table.dis) <- features
names(prop.table.dis.match) <- features
names(prop.table.dis.tam) <- features
names(prop.table.dis.tam.match) <- features
names(prop.table.dis.ai) <- features
names(prop.table.dis.ai.match) <- features
names(prop.table.val) <- features
names(prop.table.val.tam) <- features
names(prop.table.val.ai) <- features

names(prop.table.dis.tam.nr.bp) <- features_bp
names(prop.table.dis.tam.r.bp) <- features_bp
names(prop.table.dis.ai.nr.bp) <- features_bp
names(prop.table.dis.ai.r.bp) <- features_bp
names(prop.table.val.tam.nr.bp) <- features_bp
names(prop.table.val.tam.r.bp) <- features_bp
names(prop.table.val.ai.nr.bp) <- features_bp
names(prop.table.val.ai.r.bp) <- features_bp

df.p.all <- data.frame(p.dis.match = p.dis.match,
                       p.signif.dis.match = p.signif.dis.match,
                       p.dis.tam.match = p.dis.tam.match,
                       p.signif.dis.tam.match = p.signif.dis.tam.match,
                       p.dis.ai.match = p.dis.ai.match,
                       p.signif.dis.ai.match = p.signif.dis.ai.match,
                       p.dis = p.dis,
                       p.signif.dis = p.signif.dis,
                       p.dis.tam = p.dis.tam,
                       p.signif.dis.tam = p.signif.dis.tam,
                       p.dis.ai = p.dis.ai,
                       p.signif.dis.ai = p.signif.dis.ai,
                       p.val = p.val,
                       p.signif.val = p.signif.val,
                       p.val.tam = p.val.tam,
                       p.signif.val.tam = p.signif.val.tam,
                       p.val.ai = p.val.ai,
                       p.signif.val.ai = p.signif.val.ai)

rownames(df.p.all) <- features

df.p.bp <- data.frame(p.dis.tam.r.bp = p.dis.tam.r.bp,
                      p.signif.dis.tam.r.bp = p.signif.dis.tam.r.bp,
                      p.dis.tam.nr.bp = p.dis.tam.nr.bp,
                      p.signif.dis.tam.nr.bp = p.signif.dis.tam.nr.bp,
                      p.dis.ai.r.bp = p.dis.ai.r.bp,
                      p.signif.dis.ai.r.bp = p.signif.dis.ai.r.bp,
                      p.dis.ai.nr.bp = p.dis.ai.nr.bp,
                      p.signif.dis.ai.nr.bp = p.signif.dis.ai.nr.bp,
                      p.val.tam.r.bp = p.val.tam.r.bp,
                      p.signif.val.tam.r.bp = p.signif.val.tam.r.bp,
                      p.val.tam.nr.bp = p.val.tam.nr.bp,
                      p.signif.val.tam.nr.bp = p.signif.val.tam.nr.bp,
                      p.val.ai.r.bp = p.val.ai.r.bp,
                      p.signif.val.ai.r.bp = p.signif.val.ai.r.bp,
                      p.val.ai.nr.bp = p.val.ai.nr.bp,
                      p.signif.val.ai.nr.bp = p.signif.val.ai.nr.bp)

rownames(df.p.bp) <- features_bp

### put everything into an object
setClass("clinSta", slots=list(p.dis = "vector",
                               p.dis.match = "vector",
                               p.dis.tam = "vector",
                               p.dis.tam.match = "vector",
                               p.dis.tam.r.bp = "vector",
                               p.dis.tam.nr.bp = "vector",
                               p.dis.ai = "vector",
                               p.dis.ai.match = "vector",
                               p.dis.ai.r.bp = "vector",
                               p.dis.ai.nr.bp = "vector",
                               p.val = "vector",
                               p.val.tam = "vector",
                               p.val.tam.r.bp = "vector",
                               p.val.tam.nr.bp = "vector",
                               p.val.ai = "vector",
                               p.val.ai.r.bp = "vector",
                               p.val.ai.nr.bp = "vector",
                               p.signif.dis = "vector",
                               p.signif.dis.match = "vector",
                               p.signif.dis.tam = "vector",
                               p.signif.dis.tam.match = "vector",
                               p.signif.dis.tam.r.bp = "vector",
                               p.signif.dis.tam.nr.bp = "vector",
                               p.signif.dis.ai = "vector",
                               p.signif.dis.ai.match = "vector",
                               p.signif.dis.ai.r.bp = "vector",
                               p.signif.dis.ai.nr.bp = "vector",
                               p.signif.val = "vector",
                               p.signif.val.tam = "vector",
                               p.signif.val.tam.r.bp = "vector",
                               p.signif.val.tam.nr.bp = "vector",
                               p.signif.val.ai = "vector",
                               p.signif.val.ai.r.bp = "vector",
                               p.signif.val.ai.nr.bp = "vector",
                               table.dis = "list",
                               table.dis.match = "list",
                               table.dis.tam = "list",
                               table.dis.tam.match = "list",
                               table.dis.tam.r.bp = "list",
                               table.dis.tam.nr.bp = "list",
                               table.dis.ai = "list",
                               table.dis.ai.match = "list",
                               table.dis.ai.r.bp = "list",
                               table.dis.ai.nr.bp = "list",
                               table.val = "list",
                               table.val.tam = "list",
                               table.val.tam.r.bp = "list",
                               table.val.tam.nr.bp = "list",
                               table.val.ai = "list",
                               table.val.ai.r.bp = "list",
                               table.val.ai.nr.bp = "list",
                               prop.table.dis = "list",
                               prop.table.dis.match = "list",
                               prop.table.dis.tam = "list",
                               prop.table.dis.tam.match = "list",
                               prop.table.dis.tam.r.bp = "list",
                               prop.table.dis.tam.nr.bp = "list",
                               prop.table.dis.ai = "list",
                               prop.table.dis.ai.match = "list",
                               prop.table.dis.ai.r.bp = "list",
                               prop.table.dis.ai.nr.bp = "list",
                               prop.table.val = "list",
                               prop.table.val.tam = "list",
                               prop.table.val.tam.r.bp = "list",
                               prop.table.val.tam.nr.bp = "list",
                               prop.table.val.ai = "list",
                               prop.table.val.ai.r.bp = "list",
                               prop.table.val.ai.nr.bp = "list"))

clins <- new("clinSta", 
             p.dis = p.dis,
             p.dis.match = p.dis.match,
             p.dis.tam = p.dis.tam,
             p.dis.tam.match = p.dis.tam.match,
             p.dis.tam.r.bp = p.dis.tam.r.bp,
             p.dis.tam.nr.bp = p.dis.tam.nr.bp,
             p.dis.ai = p.dis.ai,
             p.dis.ai.match = p.dis.ai.match,
             p.dis.ai.r.bp = p.dis.ai.r.bp,
             p.dis.ai.nr.bp = p.dis.ai.nr.bp,
             p.val = p.val,
             p.val.tam = p.val.tam,
             p.val.tam.r.bp = p.val.tam.r.bp,
             p.val.tam.nr.bp = p.val.tam.nr.bp,
             p.val.ai = p.val.ai,
             p.val.ai.r.bp = p.val.ai.r.bp,
             p.val.ai.nr.bp = p.val.ai.nr.bp,
             p.signif.dis = p.signif.dis,
             p.signif.dis.match = p.signif.dis.match,
             p.signif.dis.tam = p.signif.dis.tam,
             p.signif.dis.tam.match = p.signif.dis.tam.match,
             p.signif.dis.tam.r.bp = p.signif.dis.tam.r.bp,
             p.signif.dis.tam.nr.bp = p.signif.dis.tam.nr.bp,
             p.signif.dis.ai = p.signif.dis.ai,
             p.signif.dis.ai.match = p.signif.dis.ai.match,
             p.signif.dis.ai.r.bp = p.signif.dis.ai.r.bp,
             p.signif.dis.ai.nr.bp = p.signif.dis.ai.nr.bp,
             p.signif.val = p.signif.val,
             p.signif.val.tam = p.signif.val.tam,
             p.signif.val.tam.r.bp = p.signif.val.tam.r.bp,
             p.signif.val.tam.nr.bp = p.signif.val.tam.nr.bp,
             p.signif.val.ai = p.signif.val.ai,
             p.signif.val.ai.r.bp = p.signif.val.ai.r.bp,
             p.signif.val.ai.nr.bp = p.signif.val.ai.nr.bp,
             table.dis = table.dis,
             table.dis.match = table.dis.match,
             table.dis.tam = table.dis.tam,
             table.dis.tam.match = table.dis.tam.match,
             table.dis.tam.r.bp = table.dis.tam.r.bp,
             table.dis.tam.nr.bp = table.dis.tam.nr.bp,
             table.dis.ai = table.dis.ai,
             table.dis.ai.match = table.dis.ai.match,
             table.dis.ai.r.bp = table.dis.ai.r.bp,
             table.dis.ai.nr.bp = table.dis.ai.nr.bp,
             table.val = table.val,
             table.val.tam = table.val.tam,
             table.val.tam.r.bp = table.val.tam.r.bp,
             table.val.tam.nr.bp = table.val.tam.nr.bp,
             table.val.ai = table.val.ai,
             table.val.ai.r.bp = table.val.ai.r.bp,
             table.val.ai.nr.bp = table.val.ai.nr.bp,
             prop.table.dis = prop.table.dis,
             prop.table.dis.match = prop.table.dis.match,
             prop.table.dis.tam = prop.table.dis.tam,
             prop.table.dis.tam.match = prop.table.dis.tam.match,
             prop.table.dis.tam.r.bp = prop.table.dis.tam.r.bp,
             prop.table.dis.tam.nr.bp = prop.table.dis.tam.nr.bp,
             prop.table.dis.ai = prop.table.dis.ai,
             prop.table.dis.ai.match = prop.table.dis.ai.match,
             prop.table.dis.ai.r.bp = prop.table.dis.ai.r.bp,
             prop.table.dis.ai.nr.bp = prop.table.dis.ai.nr.bp,
             prop.table.val = prop.table.val,
             prop.table.val.tam = prop.table.val.tam,
             prop.table.val.tam.r.bp = prop.table.val.tam.r.bp,
             prop.table.val.tam.nr.bp = prop.table.val.tam.nr.bp,
             prop.table.val.ai = prop.table.val.ai,
             prop.table.val.ai.r.bp = prop.table.val.ai.r.bp,
             prop.table.val.ai.nr.bp = prop.table.val.ai.nr.bp)

saveRDS(clins, "clins.rds")

clins <- readRDS("clins.rds")
