# Load required packages
library(arules)
library(goeveg)
library(tidyr)

# Load goeveg dataset
schedenveg_bin <- goeveg::schedenveg

# Convert dataset to binary representation
schedenveg_bin[schedenveg_bin > 0] <- 1

# Check for missing values
sum(colSums(is.na(schedenveg_bin))) # No missing values

# Get frequency of species
species_frequency <- as.data.frame(sort(colSums((schedenveg_bin))))

# Reshape data to transaction table
pvt <- pivot_longer(schedenveg_bin, cols = colnames(schedenveg_bin))
pvt$fieldname <- rep(rownames(schedenveg_bin), each = 155)
pvt$group <- substr(pvt$fieldname, 1, 1)
pvt$value <- as.factor(pvt$value)
pvt <- pvt[pvt$value == 1, c("name", "fieldname")]
colnames(pvt) <- c("variable", "area")

# Write transaction table to file
write.table(pvt, "schedenveg_transactions.csv", row.names = FALSE)

# Read transaction data
tr <- read.transactions("schedenveg_transactions.csv", format = "single", header = TRUE, sep = " ", cols = c("area", "variable"), rm.duplicates = TRUE)

# Get species present in plot "A10_16"
items_onfield <- pvt[pvt$area == "A10_16", "variable"]$variable

# Generate association rules with apriori algorithm
rules <- apriori(tr, list(supp = 0.12, maxlen = 20, minlen = 1, confidence = 0.05, maxtime = 60, minval = 0.1))

# Filter out species already present in plot "A10_16"
useful_rules <- subset(rules, subset = !(rhs %in% items_onfield))

# Filter out species not present in plot "A10_16"
useful_rules <- subset(useful_rules, subset = (lhs %oin% items_onfield))


# add p-value of Chi-Square Test
useful_rules@quality$chiSquared <-interestMeasure(useful_rules,
                  tr,
                  measure = c("chiSquared"),
                  significance = TRUE)

# p value correction
apriori_df <- DATAFRAME(useful_rules)

correction_methods <- c("none", "holm", "BY", "fdr")


df <- data.frame(matrix(0, nrow(apriori_df), length(correction_methods)))
colnames(df) <- correction_methods
for (i in seq_along(correction_methods)) {
  df[,i] <- p.adjust(apriori_df$chiSquared, method = correction_methods[i])
}
df$rule_nr <- 1:nrow(df)
added_rules <- cbind(apriori_df, df)

# save ruleset table as csv file
write.table(added_rules, "./rules_schedenveg.csv", sep = ";", dec = ",", row.names = FALSE)

df <- df[order(df$none, decreasing = TRUE), ]

# Count the number of insignificant (p >= 0.05) rules
n_insignificant_rules <- nrow(df) - nrow(df[df$none >= 0.05, ])
print(paste("Number of insignificant (p >= 0.05) rules:", n_insignificant_rules))

##### # Creating Table 1 in Paper #####
# Plot rule correction methods
png(file = "rule_correction.png", width = 12.0, height = 12.0, units = "cm", res = 90)
plot(df$none, col = "black", type = "l", xlab = "Rule", ylab = "p-Value", ylim = c(0,1))
lines(df$holm, col = "yellow")
lines(df$BY, col = "green")
lines(df$fdr, col = "red")
abline(h = 0.05, lty = 1, lwd = 2)
abline(h = 0.01, lty = 2, lwd = 2)
abline(h = 0.001, lty = 3, lwd = 2)
legend(x = 34000, y = 0.91, fill = c("black", "yellow", "green", "red"), 
       legend = c("None", "Holm", "Benjamini-Yekutieli", "FDR"), 
       title = "Correction Method", bty = "n")
dev.off()

##### # Creating Table 1 in Paper #####

# 'best' rules, by 1. ChiSquare p-value (ASC), 2. Confidence (DESC), 3. Support (DESC)
edf <- apriori_df[order(apriori_df$chiSquared, -apriori_df$confidence, -apriori_df$support), 
                  c("LHS", "RHS", "support", "confidence", "chiSquared")]

# choose best rule per Species in RHS
#write.table(edf, "remaining_rules.csv", sep = ";", dec = ",", row.names = FALSE)
 


# Selecting the First (Best) Rule per Unique Consequent
t.first <- edf[match(unique(edf$RHS), edf$RHS), ]

# Get First Three Significant Figures
t.first[, 3:5] <- signif(t.first[, 3:5], 3)

# Sorting the Rules by Confidence and Support
t.first <- t.first[order(-t.first$confidence,-t.first$support), ]

# Writing the Top 10 Rules to a CSV File
write.table(
  t.first[1:10, ],
  file = "top_rules_schedenveg_chisquared_top10.csv",
  sep = ";",
  dec = ",",
  row.names = FALSE
)
