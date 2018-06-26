library(tidyverse)
setwd("~/Documents/Insight")

#Import molecular targets database - download is available from DrugBank.ca
#However data is originally from UniProt, according to DrugBank.
targets <- read.csv("uniprot links.csv")

#Import toxicity information from DrugBank's database
#Requires first running xml_extract-drugbank-toxicity.R
df.t <- read.delim("toxicity-db.txt", sep = "\t")

#Because DrugBank does not have toxicity information for many targets (yet!)
#We will also combine this with toxocity information from Health Canada's
#Adverse Drug Event report system
#This requires first running Extract_Frequent_Adverse_Events.R
AdverseReports <- read.delim("Adverse Reports.txt", sep = "\t")

#Add toxicity information to molecular targets information
targets <- targets %>%
  left_join(df.t, by = "DrugBank.ID") %>%
  mutate(Name = toupper((Name))) %>%
  left_join(AdverseReports, by = c("Name" = "DRUGNAME")) %>%
  unite(AdverseEffects, c("ToxicityDescription", "Toxicity"), sep = ",") %>%
  filter(AdverseEffects != "NA,NA")

# targets <- targets %>%
#   mutate(UniProt.Name = gsub('Gamma-aminobutyric acid receptor subunit alpha-1', 
#                              'Gamma-aminobutyric acid receptor group',UniProt.Name)) %>%
#   mutate(UniProt.Name = gsub('Gamma-aminobutyric acid receptor subunit alpha-2', 
#                              'Gamma-aminobutyric acid receptor group',UniProt.Name))%>%
#   mutate(UniProt.Name = gsub('Gamma-aminobutyric acid receptor subunit alpha-3', 
#                              'Gamma-aminobutyric acid receptor group',UniProt.Name))%>%
#   mutate(UniProt.Name = gsub('Gamma-aminobutyric acid receptor subunit alpha-5', 
#                              'Gamma-aminobutyric acid receptor group',UniProt.Name)) %>%
#   mutate(UniProt.Name = gsub('Gamma-aminobutyric acid receptor subunit beta-1', 
#                              'Gamma-aminobutyric acid receptor group',UniProt.Name)) %>%
#   mutate(UniProt.Name = gsub('Gamma-aminobutyric acid receptor subunit beta-3', 
#                              'Gamma-aminobutyric acid receptor group',UniProt.Name))

#Data Exploration of drug counts and target counts
targets.bydrug <- targets %>%
  group_by(Name) %>%
  summarize(count = n())
targets.bydrug %>%
  mutate(count = ifelse(count>30,30,count)) %>%
  ggplot() +
  aes(x = count) +
  geom_histogram(bins = 30) +
  ylab("Number of Drugs") +
  xlab("Number of Molecular Targets") +
  theme_minimal()

drugs.bytarget <-targets %>%
  group_by(UniProt.Name) %>%
  summarize(count = n())

drugs.bytarget %>%
  mutate(count = ifelse(count>20,20,count)) %>%
  ggplot() +
  aes(x = count) +
  geom_histogram(bins = 20, aes(fill = ..)) +
  xlab("Number of Drugs") +
  ylab("Number of Molecular Targets") +
  theme_minimal()

targets.SE <- targets %>%
  mutate(Nausea = grepl("nausea", AdverseEffects, ignore.case = TRUE),
         Fever = grepl("fever|pyrexia", AdverseEffects, ignore.case = TRUE),
         Dyspnoea = grepl("dyspnoea|difficulty breathing|labored breathing", AdverseEffects, ignore.case = TRUE),
         Rash = grepl("rash", AdverseEffects, ignore.case = TRUE),
         Vomiting = grepl("vomiting", AdverseEffects, ignore.case = TRUE))

SE <- targets.SE %>%
  select(Nausea,Fever, Dyspnoea,Rash,Vomiting)

SE.breakdown <- as.data.frame(sapply(SE,table))
SE.breakdown$Total = c(0, SE.breakdown[1,1] + SE.breakdown[2,1])

SE.incidence <- SE.breakdown[2,] %>%
  gather(key = "SideEffect", value = "Incidence", 1:5) %>%
  mutate(Incidence = 100* Incidence/Total) %>%
  select(SideEffect, Incidence) 
SE.incidence %>%
  ggplot() +
  aes(x = reorder(SideEffect, Incidence), y = Incidence) +
  geom_bar(stat = 'identity', fill = 'navyblue') +
  coord_flip() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(text = element_text(size = 40)) +
  xlab('') +
  ylab('% Drugs causing effect')
  
targets.SE %>%
  select(-UniProt.ID, -UniProt.Name) %>%
  distinct() %>%
  group_by(Vomiting) %>%
  summarize(count = n())

freq.targets <- targets %>%
  group_by(UniProt.Name) %>%
  summarize(count = n()) %>%
  filter(count > 5) %>%
  select(UniProt.Name)

label = paste("Unique molecular targets\nretained:",
              as.character(dim(freq.targets)[1]))

drugs.bytarget %>%
  mutate(count = ifelse(count>20,20,count)) %>%
  ggplot() +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 40)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0,0)) +
  xlab("Drugs per target") +
  ylab("log(Frequency)") +
  annotate("text", x = 11.1, y = 500, label = label, size = 11, color = "navyblue") +
  scale_fill_manual(values = c('lightgrey', 'navyblue')) +
  geom_bar(data = . %$% invisible(hist(count, breaks = 20)) %$% as_tibble(list(mids = mids, counts = counts, density = density)) 
           %>% mutate(fill = ifelse(mids <5,"removed","retained")), mapping = aes(x = mids, y = counts, fill = fill), stat = "identity")

freq.targets.drugs <-targets %>%
  inner_join(freq.targets, by = "UniProt.Name") %>%
  group_by(Name) %>%
  summarize(n_targets = n())

label2 = paste("Unique drugs\n(observations):",
              as.character(dim(freq.targets.drugs)[1]))

targets.SE %>%
  inner_join(freq.targets, by = "UniProt.Name") %>%
  group_by(Name) %>%
  summarize(count = n()) %>%
  mutate(count = ifelse(count>20,20,count)) %>%
  ggplot() +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 40)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0,0)) +
  xlab("Targets per drug") +
  ylab("log(Frequency)") +
  annotate("text", x = 11.1, y = 200, label = label2, size = 11, color = "navyblue") +
  scale_fill_manual(values = c('navyblue', 'navyblue')) +
  geom_bar(data = . %$% invisible(hist(count, breaks = 20)) %$% as_tibble(list(mids = mids, counts = counts, density = density)) 
           %>% mutate(fill = ifelse(mids <5,"removed","retained")), mapping = aes(x = mids, y = counts, fill = fill), stat = "identity")
  

#Wide format array for learning, with a smaller number of features.
Nausea <- targets.SE %>%
  inner_join(freq.targets, by = "UniProt.Name") %>%
  select(Name, UniProt.Name, Nausea) %>%
  mutate(molecular_target_value = 1) %>%
  distinct() %>%
  group_by(Name, Nausea) %>%
  spread(key = UniProt.Name, value = molecular_target_value, fill = 0)

write_delim(Nausea, "Nausea.txt", delim = "\t")
