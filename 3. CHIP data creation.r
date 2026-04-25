library(data.table)
library(dplyr)

dat1 <- fread('data_participant.csv')
gene <- fread('GENE.csv')
dat2 <- fread('Death_register.csv')
dat2 <- dat2[,c(1:2)]

merged_data <- merge(dat1, dat2, by.x = "Participant ID", by.y = "Participant ID")
names(merged_data)
merged_data$Death_event <- ifelse(merged_data$`Date of death | Instance 0` == "", 0, 1)

merged_data$`Date of death | Instance 0` <- as.Date(merged_data$`Date of death | Instance 0`, format = "%Y/%m/%d")
merged_data$`Date of attending assessment centre | Instance 0` <- as.Date(merged_data$`Date of attending assessment centre | Instance 0`, format = "%Y-%m-%d")
fixed_date <- as.Date("2024-12-02", format = "%Y-%m-%d")

merged_data$Death_time <- ifelse(
  merged_data$Death_event == 1,
  merged_data$`Date of death | Instance 0` - merged_data$`Date of attending assessment centre | Instance 0`,
  fixed_date - merged_data$`Date of attending assessment centre | Instance 2`
)

valid_id <- gene %>%
  filter(!is.na(`Clonal haematopoiesis of indeterminate potential (CHIP) number of variants`),
         `Genomic ancestry` == "European ancestry (EUR)") %>%
  pull('Participant ID') %>%          
  unique() 
gene <- gene %>% filter(`Participant ID` %in% valid_id)
merged_data <- gene %>%
  left_join(merged_data, by = "Participant ID") %>%
  select(all_of(c(names(merged_data)[names(merged_data) != "Participant ID"], names(gene))))
merged_data <- merged_data %>% relocate(10, .before = 1)

write.csv(merged_data, "20251211_chip.csv", row.names = FALSE, na = "")
