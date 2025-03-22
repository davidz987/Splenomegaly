# Individual-level data are not publicly available due to research participant privacy concerns; 
# however, requests from accredited researchers for access to individual-level data relevant 
# to this manuscript can be made by contacting the corresponding author. 

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

### Read in spleen features file
features <- data.table::fread('spleen_features.csv', sep=',', header=T)

calculatorOI = "pyradiomics"

# Read in splenectomy accessions
splenectomy <- data.table::fread("splenectomy_cases.txt")
splenectomy_vols <- features %>%
  filter(calculator == calculatorOI & measure == "shape" & metric == "MeshVolume") %>%
  rename(Label = label, rawval = value) %>%
  distinct() %>%
  mutate(rawval = rawval/1000) %>%   # mm^3 to mL 
  filter(accession %in% splenectomy$Accession)
splenectomy_threshold <- median(splenectomy_vols$rawval)

# Filter for volume
# - concert mm^3 to mL
# - remove cases of splenectomy
# - remove spleen volumes less than median of splenectomy cases
# - remove more than 4SDs from mean
df_hist = features %>%
  filter(calculator == calculatorOI & measure == "shape" & metric == "MeshVolume") %>%
  rename(Label = label, rawval = value) %>%
  distinct() %>%
  mutate(rawval = rawval/1000) # mm^3 to mL 
length(unique(df_hist$accession))

df_hist <- df_hist %>%
  filter(!(accession %in% splenectomy$Accession)) 
length(unique(df_hist$accession))

df_hist <- df_hist %>%
  filter(rawval > splenectomy_threshold) 
length(unique(df_hist$accession))

df_hist <- df_hist %>%
  filter((rawval > mean(rawval) - (4*sd(rawval))) & (rawval < mean(rawval) + (4*sd(rawval))))
length(unique(df_hist$accession))

### *BEGIN* Add correction factor for enhanced spleen volumes
# Read in contrast annotations
# - updated from classifier
new_contrast <- data.table::fread("contrast.txt")
df_hist <- merge(df_hist %>% select(id, accession, series_number, rawval),
                       new_contrast, by=c('id', 'accession', 'series_number'))
length(unique(df_hist$accession))

# identify accessions with both enhanced and unenhanced scans
duplicated_rows <- df_hist %>%
  group_by(id, accession) %>%
  filter(n_distinct(contrast) > 1) %>%
  ungroup() %>%
  group_by(id, accession, contrast) %>%
  summarise(value = mean(rawval))

mean(duplicated_rows$value[duplicated_rows$contrast == 0])
sd(duplicated_rows$value[duplicated_rows$contrast == 0])
mean(duplicated_rows$value[duplicated_rows$contrast == 1])
sd(duplicated_rows$value[duplicated_rows$contrast == 1])

duplicated_rows <- duplicated_rows %>%
  pivot_wider(
    names_from = contrast,
    values_from = value,
    names_prefix = "contrast_"
  )

ggscatter(duplicated_rows, x="contrast_1", y="contrast_0",
          color = "cornflowerblue", shape = 16, size = 3, alpha=0.5, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "orange", linetype = "dashed"), # Customize reg. line
          xlab = "Enhanced Splenic Volume (mL)",
          ylab = "Unnhanced Splenic Volume (mL)",
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) + 
  scale_x_continuous(breaks = seq(0, 800, by = 100)) +
  scale_y_continuous(breaks = seq(0, 800, by = 100))

res <- lm(contrast_0 ~ contrast_1, data=duplicated_rows)
summary(res)

# Get coefficient to correct enhanced to unenhanced (1/beta)
beta <- coef(summary(res))[2,1]

# Correct enhanced splenic volumes
df_hist$rawval[df_hist$contrast == 1] <- beta * df_hist$rawval[df_hist$contrast == 1] + coef(summary(res))[1,1]
### *END*

# Summarise series within an accession
df_hist <- df_hist %>%  
  group_by(id, accession) %>%
  summarise(value = mean(rawval))

# Read in other clinical characteristics
sex_age <- data.table::fread("age_sex.csv", sep=",")
sex_age <- sex_age %>% rename(SEX = Sex)
sex_age$Sex <- ifelse(sex_age$SEX == "Male", 1, 0)

race_ethnicity <- data.table::fread("race_ethnicity.csv", sep=",")

height <- data.table::fread("height.csv", sep=",")
weight <- data.table::fread("weight.csv", sep=",")

# Merge with age, sex, race, ethnicity, height, weight
# - filter height/weight for outliers --> to reduce disproportionate impact of extreme outliers
input_df <- merge(df_hist %>% select(id, accession, value), sex_age %>% select(id, birth_date_shift, SEX), by.x='id', by.y='id')
input_df <- merge(input_df, mapping %>% select(Accession, ShiftedDate), by.x="accession", by.y="Accession")
input_df$Age <- as.numeric(as.Date(input_df$ShiftedDate) - as.Date(input_df$birth_date_shift)) / 365.25
input_df <- input_df %>% select(-birth_date_shift, -ShiftedDate)

input_df <- merge(input_df, race_ethnicity %>% select(-id), by.x='id', by.y='id')

input_df <- merge(input_df, height %>% select(id, mean_value), by.x='id', by.y='id')
input_df$mean_value <- input_df$mean_value * 2.54 #change height from in to cm
colnames(input_df)[8] <- "height"
input_df <- input_df %>%
  filter((height > mean(height) - (4*sd(height))) & (height < mean(height) + (4*sd(height))))

input_df <- merge(input_df, weight %>% select(id, mean_value), by.x='id', by.y='id')
input_df$mean_value <- input_df$mean_value/35.274 #change weight oz to kg
colnames(input_df)[9] <- "weight"
input_df <- input_df %>%
  filter((weight > mean(weight) - (4*sd(weight))) & (weight < mean(weight) + (4*sd(weight))))

# Subset for IDs with phenotype info
pheno <- read.delim(file = "conditions_my_phecode_ro2.txt", 
                           sep = '\t', header = TRUE, check.names = FALSE)
pheno <- pheno %>% mutate_if(is.numeric,as.logical)

pheno_ids <- unique(mapping$id[mapping$id %in% pheno$id])
input_df <- input_df[input_df$id %in% pheno_ids,]
nrow(input_df)

# Compare chest v non-chest scans
chest <- data.table::fread("spleen_complete_chest.txt", sep="\t")
plot_df <- merge(input_df, chest %>% select(accession, chest), by='accession')

# plot of spleen volume stratified by chest scan
plot_df$chest <- factor(plot_df$chest, levels = c(1, 0))
ggplot(plot_df, aes(x=value, color=chest, fill=chest)) +
  geom_histogram(aes(y=after_stat(density)), alpha=0.5, position="identity", binwidth=20) +
  geom_density(alpha=0.6) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("Spleen Volume (mL)")

plot_df <- merge(input_df, chest %>% select(accession, chest) %>% distinct(), by='accession')
plot_df$chest[plot_df$chest == 0] <- "non-chest study"
plot_df$chest[plot_df$chest == 1] <- "chest study"
ggboxplot(plot_df, x = "chest", y = "value", alpha=0.05,
          xlab = "CT study type",
          ylab = "Spleen Volume (mL)")


# Get cases of splenomegaly - use radiology report
rad_splenomegaly <- data.table::fread("splenomegaly_cases.txt")
plot_df <- input_df
plot_df$Splenomegaly <- ifelse(plot_df$accession %in% rad_splenomegaly$Accession, "Yes", "No")
plot_df$Splenomegaly <- factor(plot_df$Splenomegaly, levels = c("Yes", "No"))
length(unique(plot_df$id[plot_df$Splenomegaly == "Yes"]))
nrow(plot_df[plot_df$Splenomegaly == "Yes",])

# Plot of spleen volume stratified by splenomegaly cases
ggdensity(plot_df, x = "value", color = "Splenomegaly", fill = "Splenomegaly",
          add = "mean", palette = c("#E7B800", "#00AFBB"),
          xlab = "Spleen Volume (mL)")

# Plot of spleen volume against sex
ggboxplot(plot_df, x="SEX", y="value",
          notch=TRUE, color="Splenomegaly", alpha=0.1,
          palette = c('#FF8B00','#999999'),
          ylab = "Spleen Volume (mL)")

# Plot of spleen volume after age/height/weight
ggscatter(plot_df, x="weight", y="value", color="Splenomegaly", alpha="Splenomegaly",
          palette = c('#FF8B00','#999999'), shape = 16, size = 2, # Points color, shape and size
          add = "loess", conf.int = TRUE,  # Add local regression fitting
          add.params = list(color = "cornflowerblue"), # Customize reg. line
          xlab = "Weight (kg)",
          ylab = "Spleen Volume (mL)") +
  scale_alpha_manual(values=c(0.8,0.1))

### PheWAS
library(PheWAS)

# My phecode annotations
my_phecode_info <- data.table::fread("pheinfo.txt", colClasses = list(character = "phecode"))

# Group to get max - for PheWAS
phewas_df <- input_df %>%
  group_by(id) %>%
  slice_max(value, n = 1, with_ties = FALSE) %>%
  ungroup()

phewas_df <- phewas_df %>%
  inner_join(mapping %>% select(-ShiftedDate) %>% rename(accession = Accession), by = c("id","accession"))
nrow(phewas_df)

RankNorm <- function(col) {
  k = 0.375
  n <- length(col)
  r <- rank(col, ties.method = "average")
  out <- qnorm((r - k) / (n - (2 * k) + 1))
  return(out)
}

geno_df <- phewas_df %>% select(id, value) %>%
  rename(max_vol = value)

# RINT transform
geno_df$norm_vol <- RankNorm(geno_df$max_vol)
geno_df <- geno_df %>% select(-max_vol)

cov_df <- phewas_df %>% select(id, SEX, Age, height, weight, Race)

results <- phewas(phenotypes = pheno, genotypes = geno_df, covariates = cov_df)
results$study<-'IDPs'

results2 <- merge(results, my_phecode_info %>% select(phecode, description, group), by.x="phenotype", by.y="phecode")
results2$adj_p <- p.adjust(results2$p, "bonferroni")
results2 <- results2[,c(1:7,19,8:18)]

myplot <- phewasManhattan(results, OR.direction=T, annotate.angle=0, annotate.size=4, point.size=2) +
  ggtitle(paste0("Spleen Volume - PheWAS")) +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1, size=8),
        plot.margin = margin(10,10,10,20))


### Patients w/ cirrhosis of liver phecode (571)
cirrhosis <- data.table::fread("phecode_571.txt")

# Get first date of diagnosis
first_cirrhosis <- cirrhosis %>%
  group_by(id) %>%
  slice_min(order_by = condition_start_date, n = 1) %>%
  ungroup()

cirrhosis_ids <- first_cirrhosis %>% select(id,condition_start_date) %>% distinct()

### Patients w/ lymphoma phecode (202)
lymphoma <- data.table::fread("phecode_202.txt")

# Get first date of diagnosis
first_lymphoma <- lymphoma %>%
  group_by(id) %>%
  slice_min(order_by = condition_start_date, n = 1) %>%
  ungroup()

lymphoma_ids <- first_lymphoma %>% select(id,condition_start_date) %>% distinct()

### Patients w/ hemolytic anemia phecode (282 / 283)
hemolytic <- data.table::fread("phecode_282.txt")
hemolytic2 <- data.table::fread("phecode_283.txt")
hemolytic <- rbind(hemolytic, hemolytic2)

# Get first date of diagnosis
first_hemolytic <- hemolytic %>%
  group_by(id) %>%
  slice_min(order_by = condition_start_date, n = 1) %>%
  ungroup()

hemolytic_ids <- first_hemolytic %>% select(id,condition_start_date) %>% distinct()

### Patients w/ thrombocytopenia phecode (287.3)
thrombocytopenia <- data.table::fread("phecode_2873.txt")

# Get first date of diagnosis
first_thrombocytopenia <- thrombocytopenia %>%
  group_by(id) %>%
  slice_min(order_by = condition_start_date, n = 1) %>%
  ungroup()

thrombocytopenia_ids <- first_thrombocytopenia %>% select(id,condition_start_date) %>% distinct()

### Patients w/ septicemia phecode (038)
septicemia <- data.table::fread("phecode_038.txt")

# Get first date of diagnosis
first_septicemia <- septicemia %>%
  group_by(id) %>%
  slice_min(order_by = condition_start_date, n = 1) %>%
  ungroup()

septicemia_ids <- first_septicemia %>% select(id,condition_start_date) %>% distinct()

### Patients w/ CHF phecode (428)
CHF <- data.table::fread("phecode_428.txt")

# Get first date of diagnosis
first_CHF <- CHF %>%
  group_by(id) %>%
  slice_min(order_by = condition_start_date, n = 1) %>%
  ungroup()

CHF_ids <- first_CHF %>% select(id,condition_start_date) %>% distinct()

# 'no conditions cohort'
# - if scan was done *before* first diagnoses of any condition
# - AND scan is not annotated as splenomegaly by rad report
all_ids <- rbind(cirrhosis_ids, lymphoma_ids, hemolytic_ids,
                 thrombocytopenia_ids, septicemia_ids, CHF_ids) %>%
              group_by(id) %>%
              slice_min(order_by = condition_start_date, n = 1) %>%
              ungroup() %>%
              distinct()
input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, all_ids, by="id", all.x=TRUE)
input_df_num$Other <- ifelse(is.na(input_df_num$condition_start_date) | (!is.na(is.na(input_df_num$condition_start_date)) & input_df_num$ShiftedDate < input_df_num$condition_start_date), 
                                    "Yes", "No")
input_df_num$Other <- ifelse((input_df_num$Other == "Yes") & !(input_df_num$accession %in% rad_splenomegaly$Accession), 
                             "Yes", "No")
length(unique(input_df_num$id[input_df_num$Other == "Yes"]))
nrow(input_df_num[input_df_num$Other == "Yes",])
input_df_num <- input_df_num %>%
  filter(Other == "Yes") %>%
  select(id, value) %>%
  mutate(phenotype = "All Others")
boxplot_df <- input_df_num

input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- input_df_num %>%
  filter(accession %in% rad_splenomegaly$Accession) %>%
  select(id, value) %>%
  mutate(phenotype = "Splenomegaly")
nrow(input_df_num)
length(unique((input_df_num$id)))
boxplot_df <- rbind(boxplot_df, input_df_num)

input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, cirrhosis_ids, by="id", all.x=TRUE)
input_df_num$Cirrhosis <- ifelse(!is.na(input_df_num$condition_start_date) & input_df_num$ShiftedDate >= input_df_num$condition_start_date, 
                                    "Yes", "No")
length(unique(input_df_num$id[input_df_num$Cirrhosis == "Yes"]))
nrow(input_df_num[input_df_num$Cirrhosis == "Yes",])
input_df_num <- input_df_num %>%
  filter(Cirrhosis == "Yes") %>%
  select(id, value) %>%
  mutate(phenotype = "Cirrhosis")
boxplot_df <- rbind(boxplot_df, input_df_num)

input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, lymphoma_ids, by="id", all.x=TRUE)
input_df_num$Lymphoma <- ifelse(!is.na(input_df_num$condition_start_date) & input_df_num$ShiftedDate >= input_df_num$condition_start_date, 
                                 "Yes", "No")
length(unique(input_df_num$id[input_df_num$Lymphoma == "Yes"]))
nrow(input_df_num[input_df_num$Lymphoma == "Yes",])
input_df_num <- input_df_num %>%
  filter(Lymphoma == "Yes") %>%
  select(id, value) %>%
  mutate(phenotype = "Lymphoma")
boxplot_df <- rbind(boxplot_df, input_df_num)

input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, hemolytic_ids, by="id", all.x=TRUE)
input_df_num$Hemolytic <- ifelse(!is.na(input_df_num$condition_start_date) & input_df_num$ShiftedDate >= input_df_num$condition_start_date, 
                                "Yes", "No")
length(unique(input_df_num$id[input_df_num$Hemolytic == "Yes"]))
nrow(input_df_num[input_df_num$Hemolytic == "Yes",])
input_df_num <- input_df_num %>%
  filter(Hemolytic == "Yes") %>%
  select(id, value) %>%
  mutate(phenotype = "Hemolytic anemia")
boxplot_df <- rbind(boxplot_df, input_df_num)

input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, thrombocytopenia_ids, by="id", all.x=TRUE)
input_df_num$Thrombocytopenia <- ifelse(!is.na(input_df_num$condition_start_date) & input_df_num$ShiftedDate >= input_df_num$condition_start_date, 
                                 "Yes", "No")
length(unique(input_df_num$id[input_df_num$Thrombocytopenia == "Yes"]))
nrow(input_df_num[input_df_num$Thrombocytopenia == "Yes",])
input_df_num <- input_df_num %>%
  filter(Thrombocytopenia == "Yes") %>%
  select(id, value) %>%
  mutate(phenotype = "Thrombocytopenia")
boxplot_df <- rbind(boxplot_df, input_df_num)

input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, septicemia_ids, by="id", all.x=TRUE)
input_df_num$Septicemia <- ifelse(!is.na(input_df_num$condition_start_date) & input_df_num$ShiftedDate >= input_df_num$condition_start_date, 
                                        "Yes", "No")
length(unique(input_df_num$id[input_df_num$Septicemia == "Yes"]))
nrow(input_df_num[input_df_num$Septicemia == "Yes",])
input_df_num <- input_df_num %>%
  filter(Septicemia == "Yes") %>%
  select(id, value) %>%
  mutate(phenotype = "Septicemia")
boxplot_df <- rbind(boxplot_df, input_df_num)

input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, CHF_ids, by="id", all.x=TRUE)
input_df_num$CHF <- ifelse(!is.na(input_df_num$condition_start_date) & input_df_num$ShiftedDate >= input_df_num$condition_start_date, 
                                  "Yes", "No")
length(unique(input_df_num$id[input_df_num$CHF == "Yes"]))
nrow(input_df_num[input_df_num$CHF == "Yes",])
input_df_num <- input_df_num %>%
  filter(CHF == "Yes") %>%
  select(id, value) %>%
  mutate(phenotype = "CHF")
boxplot_df <- rbind(boxplot_df, input_df_num)
boxplot_df$phenotype <- factor(boxplot_df$phenotype, levels = c("All Others", "CHF", 
                                                                "Cirrhosis", "Hemolytic anemia", 
                                                                "Lymphoma", "Septicemia", 
                                                                "Splenomegaly", "Thrombocytopenia"))

ggboxplot(boxplot_df, x="phenotype", y="value",
          notch=TRUE, alpha=0.05,
          xlab = "Phenotypes",
          ylab = "Spleen Volume (mL)")

# Use 'all other' cohort - correct for age, sex, height, weight
input_df_num <- merge(input_df, mapping, by.x="accession", by.y="Accession")
input_df_num <- merge(input_df_num, all_ids, by="id", all.x=TRUE)
input_df_num$Other <- ifelse(is.na(input_df_num$condition_start_date) | (!is.na(is.na(input_df_num$condition_start_date)) & input_df_num$ShiftedDate < input_df_num$condition_start_date), 
                             "Yes", "No")
input_df_num$Other <- ifelse((input_df_num$Other == "Yes") & !(input_df_num$accession %in% rad_splenomegaly$Accession), 
                             "Yes", "No")
length(unique(input_df_num$id[input_df_num$Other == "Yes"]))
nrow(input_df_num[input_df_num$Other == "Yes",])
all_others <- input_df_num[input_df_num$Other == "Yes",] %>%
  select(colnames(input_df))
sick_df <- anti_join(input_df, all_others)

# Fit a linear model to adjust the quantitative trait in the healthy cohort
model <- lm(value ~ SEX + Age + height + weight, data = all_others[all_others$accession %in% chest$accession[chest$chest == 0],])
summary(model)

# Splenomegaly threshold by bootstrapping for individual CIs
residuals <- model$residuals
sick_df$pred_value <- predict(model, newdata = sick_df)

# Bootstrapping:
bootstrap_ci <- function(residuals, predicted_volume){
  n_bootstrap <- 1000
  boot_predictions <- numeric(n_bootstrap)
  for (i in 1:n_bootstrap) {
    sampled_residual <- sample(residuals, size = 1, replace = TRUE)
    boot_predictions[i] <- predicted_volume + sampled_residual
  }
  ci <- quantile(boot_predictions, c(0.025, 0.975))
  return(ci[2])
}

set.seed(123)
sick_df$up_ci <- sapply(sick_df$pred_value, function(x) bootstrap_ci(residuals, x))

sick_df$pred_splenomegaly <- ifelse(sick_df$value > sick_df$up_ci, "Yes", "No")

### Compare with radiology reports
sick_df_num <- sick_df
sick_df_num$rad_splenomegaly <- ifelse(sick_df_num$accession %in% rad_splenomegaly$Accession, "Yes", "No")

TP <- nrow(sick_df_num[sick_df_num$pred_splenomegaly == "Yes" & sick_df_num$rad_splenomegaly == "Yes",])
FP <- nrow(sick_df_num[sick_df_num$pred_splenomegaly == "Yes" & sick_df_num$rad_splenomegaly == "No",])
PPV <- TP / (TP + FP)

TN <- nrow(sick_df_num[sick_df_num$pred_splenomegaly == "No" & sick_df_num$rad_splenomegaly == "No",])
FN <- nrow(sick_df_num[sick_df_num$pred_splenomegaly == "No" & sick_df_num$rad_splenomegaly == "Yes",])
NPV <- TN / (TN + FN)

sensitivity <- TP / (TP + FN)
specificity <- TN / (FP + TN)


# Compare spleen axes and correlations with volume
# - longest axis (mm) vs meshvolume
# - ellipsoid volume vs meshvolume
majoraxis = features %>%
  filter(calculator == calculatorOI & measure == "shape" & metric == "MajorAxisLength") %>%
  rename(Label = label, rawval = value) %>%
  distinct() %>%
  mutate(rawval = rawval/10) %>% # mm to cm 
  filter(accession %in% input_df$accession) %>%
  group_by(id, accession) %>%
  summarise(value = mean(rawval)) %>%
  rename(majoraxis = value)
majoraxis <- merge(majoraxis, input_df, by=c('id', 'accession'))

ggscatter(majoraxis, x="majoraxis", y="value",
          color = "cornflowerblue", shape = 16, size = 3, alpha=0.5, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "orange", linetype = "dashed"), # Customize reg. line
          xlab = "Major axis length (cm)",
          ylab = "Segmented splenic volume (mL)",
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) + 
  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
  ylim(0,850)


minoraxis = features %>%
  filter(calculator == calculatorOI & measure == "shape" & metric == "MinorAxisLength") %>%
  rename(Label = label, rawval = value) %>%
  distinct() %>%
  mutate(rawval = rawval/10) %>% # mm to cm 
  filter(accession %in% input_df$accession) %>%
  group_by(id, accession) %>%
  summarise(value = mean(rawval)) %>%
  rename(minoraxis = value)
leastaxis = features %>%
  filter(calculator == calculatorOI & measure == "shape" & metric == "LeastAxisLength") %>%
  rename(Label = label, rawval = value) %>%
  distinct() %>%
  mutate(rawval = rawval/10) %>% # mm to cm 
  filter(accession %in% input_df$accession) %>%
  group_by(id, accession) %>%
  summarise(value = mean(rawval)) %>%
  rename(leastaxis = value)
axes = merge(merge(leastaxis, minoraxis, by=c('id', 'accession')), majoraxis, by=c('id', 'accession'))

axes$ellipsoid <- (4/3)*pi*(axes$leastaxis/2)*(axes$minoraxis/2)*(axes$majoraxis/2)

ggscatter(axes, x="ellipsoid", y="value",
          color = "cornflowerblue", shape = 16, size = 3, alpha=0.5, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "orange", linetype = "dashed"), # Customize reg. line
          xlab = "Calculated ellipsoid splenic volume (mL)",
          ylab = "Segmented splenic volume (mL)",
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) + 
  scale_x_continuous(breaks = seq(0, 1500, by = 200)) +
  ylim(0,850)

res <- lm(value ~ ellipsoid, data=axes)
summary(res)
residuals <- resid(res)
residual_variance <- sum(residuals^2) / (length(residuals) - length(coef(res)))
rmse <- sqrt(mean(residuals^2))
residual_sd <- sd(residuals)


# Look at cases of splenosis and splenule
splenosis <- data.table::fread("splenosis_cases.txt")
splenule <- data.table::fread("splenule_cases.txt")
nrow(input_df[input_df$accession %in% splenosis$Accession,])
nrow(input_df[input_df$accession %in% splenule$Accession,])

