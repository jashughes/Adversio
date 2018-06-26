# Adversio
This tool predicts drug side effects based on publicly available molecular information

# Compiling Data
1. Go to https://www.drugbank.ca/releases/latest and download their database (XML format). You will need to make an account, but the data is free to access.
2. Run xml_extract-drugbank-toxicity.R to extract toxicity (side effects) information.
3. Go to https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/canada-vigilance-online-database-data-extract.html to download the adverse reports database. This data is free to access.
4. Run Extract_Frequent_Adverse_Events.R to extract side effects information.
5. Go to https://www.drugbank.ca/releases/latest#external-links and download 'Target Drug-UniProt Links' for all drugs.
6. Run 'targets_by_drug.R' to join these datasets & output the data into a format ready for training the model. This script also produces some data exploration graphs (histograms, etc).
