library(tidyverse)
library(readxl)
ic50 <- read_excel("FinalProject/data/v17.3_fitted_dose_response.xlsx")  
head(ic50)


ic50 <- read_tsv("FinalProject/data/DREAM7_DrugSensitivity1_Drug_Response_Training.txt")
head(ic50)
ic50 %>% 
    unite(cell_drug, CELL_LINE_NAME, DRUG_NAME) %>%
    filter(duplicated(cell_drug, fromLast = TRUE) |
               duplicated(cell_drug, fromLast = FALSE)) %>%
    arrange(cell_drug) %>%
    data.frame() %>% head


ic50[duplicated(ic50$IC50_RESULTS_ID,fromLast = TRUE)|
         duplicated(ic50$IC50_RESULTS_ID,fromLast = FALSE),]

ic50[(duplicated(ic50$CELL_LINE_NAME,fromLast = TRUE)&
         duplicated(ic50$DRUG_ID,fromLast = TRUE)) |
         (duplicated(ic50$CELL_LINE_NAME,fromLast = FALSE)&
         duplicated(ic50$DRUG_ID,fromLast = FALSE)) , ] %>% 
    data.frame() %>% head 
    arrange(CELL_LINE_NAME, DRUG_ID)

    (duplicated(ic50$CELL_LINE_NAME,fromLast = TRUE)&
            duplicated(ic50$DRUG_ID,fromLast = TRUE)) |
        (duplicated(ic50$CELL_LINE_NAME,fromLast = FALSE)&
             duplicated(ic50$DRUG_ID,fromLast = FALSE)) %>% sum
    