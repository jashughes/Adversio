library(tidyverse)
setwd("~/Documents/Insight/CanVigAdRxn")

# There are some caveats about this approach that should be noted.
# Most importantly, numerical comparisons about frequency should be taken with
# a grain of salt since not all reactions are reported.
# Sometimes reactions are reported with slightly different terms, which
# means some common reactions could be removed from analysis.
# As a first approximation, this approach is reasonable for the
# intended application. A pharmaceutical company would likely have
# access to better records than I did.


#Load Drug Report Data & select most common (freq > 100) drugs
header <- read.delim("Report_Drug_header.csv", sep = ",") %>%
  select(Attribute.Physical.Name)
report_drug <- read.delim("report_drug.txt", sep = "$", header = FALSE,
                          col.names = t(header))
report_drug <- report_drug %>%
  select(-ends_with("_FR")) %>%
  mutate(DRUGNAME = gsub(" *", "",DRUGNAME))
freqDrug <- report_drug %>%
  group_by(DRUGNAME) %>%
  summarize(count = n()) %>%
  filter(count > 100) %>%
  dplyr::select(DRUGNAME)
report_drug <- report_drug %>%
  inner_join(freqDrug, by = "DRUGNAME") 

freqDrug_repID <- report_drug %>%
  dplyr::select(REPORT_ID) %>%
  distinct()
#####


#Load Reactions, limit to reports from most frequent drugs
header <- read.delim("Reactions_header.csv", sep = ",") %>%
  select(Attribute.Physical.Name)
reactions <- read.delim("reactions.txt", sep = "$", header = FALSE,
                          col.names = t(header))

rxn.bydrug <- reactions %>%
  select(-ends_with("_FR")) %>%
  inner_join(report_drug, by = "REPORT_ID") %>%
  select(REPORT_ID, DRUGNAME, PT_NAME_ENG, SOC_NAME_ENG) %>%
  distinct()
remove(report_drug)
remove(reactions)

#Calculate the most ferquent reactions, set conservatively here to 5% of reports
reps_per_drug <- rxn.bydrug %>%
  group_by(DRUGNAME) %>%
  select(DRUGNAME, REPORT_ID) %>%
  distinct()%>%
  summarize(reps_per_drug = n())

rxn.bydrug_withcounts <- rxn.bydrug %>%
  left_join(reps_per_drug, by = "DRUGNAME")

rxn.bydrug_withcounts <- rxn.bydrug_withcounts %>%
  group_by(DRUGNAME, PT_NAME_ENG, SOC_NAME_ENG, reps_per_drug) %>%
  summarize(rxn_count = n())

rxn.frequency <- rxn.bydrug_withcounts %>%
  mutate(Report_Frequency = rxn_count/reps_per_drug) %>%
  filter(Report_Frequency > .05)

#Remove unnecessary tables to free up space if necessary
remove(freqDrug, freqDrug_repID,reps_per_drug,rxn.bydrug,rxn.bydrug_withcounts)

#Graph, for rough informative purposes
rxn.frequency %>%
  group_by(DRUGNAME) %>%
  summarize(count = n()) %>%
  mutate(count = ifelse(count > 40, 40, count)) %>%
  ggplot() +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 40)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Number of common \nreported adverse effects") +
  ylab("Number of Drugs") +
  #annotate("text", x = 11.1, y = 500, label = label, size = 11, color = "navyblue") +
  scale_fill_manual(values = c('navyblue', 'navyblue')) +
  geom_bar(data = . %$% invisible(hist(count, breaks = 20)) %$% as_tibble(list(mids = mids, counts = counts, density = density)) 
           %>% mutate(fill = ifelse(mids <5,"removed","retained")), mapping = aes(x = mids, y = counts, fill = fill), stat = "identity")
  

AdverseReports <- rxn.frequency %>%
  group_by(DRUGNAME) %>%
  summarize(Toxicity = paste(PT_NAME_ENG, collapse = ","))

write_delim(AdverseReports, "Adverse Reports.txt", delim = "\t")

