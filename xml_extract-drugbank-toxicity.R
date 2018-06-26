library(tidyverse)
library(xml2)
library(XML)
library(gdata)
setwd("~/Documents/Insight")


doc <- read_xml("full database.xml")
ns  <- c("db" = "http://www.drugbank.ca")

toxicity <- xml_find_all(doc, "//db:drug/db:toxicity/text()", ns = ns)

t.parent <- xml_parent(toxicity)
t.p.parent <- xml_parent(t.parent)

xml_find_first(t.p.parent[1],"//db:drug/db:drugbank-id", ns = ns)

df.t <- data.frame(matrix(ncol = 2, nrow = length(toxicity)))
colnames(df.t) <- c("DrugBank.ID",
                    "ToxicityDescription")
  
for (i in 1:length(toxicity)){
  df.t[i,1] <- xml_text(xml_find_first(t.p.parent[[i]], './db:drugbank-id', ns=ns))
  df.t[i,2] <- as.character(toxicity[i])
}

write_delim(df.t, "toxicity-db.txt", delim = "\t")
