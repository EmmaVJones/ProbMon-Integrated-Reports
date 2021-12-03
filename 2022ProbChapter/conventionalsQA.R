field <- read_excel('originalData/Lucy/FieldData_2019_2020_EVJ.xls')
analyte <- read_excel('originalData/Lucy/Analy_2019_2020_EVJ.xls')
metals <- read_excel('originalData/Lucy/Metals_2019_2020_EVJ.xls')

# all from 2020 so don't need to drop 2019 data

field <- field %>% 
  filter(!is.na(ID)) %>% 
  dplyr::select(FDT_STA_ID = ID, FDT_DATE_TIME = `Date Time`, FDT_FIELD_PH = `Field Ph`,
                FDT_TEMP_CELCIUS = `Temp Celsius`, DO_mg_L = `Do Optical`, 
                FDT_SPECIFIC_CONDUCTANCE = `Specific Conductance`) 

analyte <- analyte %>% 
  filter(!is.na(ID)) %>% 
  dplyr::rename(FDT_STA_ID = ID, FDT_DATE_TIME = `Date Time`, Ana_Sam_Mrs_Container_Id_Desc = `Container Description`) %>% 
  dplyr::select(FDT_STA_ID, FDT_DATE_TIME, Ana_Sam_Mrs_Container_Id_Desc, `Parm Name`, `Uncensored Value`) %>% 
  group_by(FDT_STA_ID, FDT_DATE_TIME, Ana_Sam_Mrs_Container_Id_Desc, `Parm Name` ) %>% 
  mutate(n = 1:n()) %>% 
  ungroup() %>% 
  group_by(FDT_STA_ID, FDT_DATE_TIME, Ana_Sam_Mrs_Container_Id_Desc) %>% 
  pivot_wider(names_from = `Parm Name`, values_from = `Uncensored Value`) %>% 
  rename('NITROGEN_mg_L' = 'NITROGEN, TOTAL (MG/L AS N)', 
         'PHOSPHORUS_mg_L'= 'PHOSPHORUS, TOTAL (MG/L AS P)',
         #'TDS RESIDUE,TOTAL FILTRABLE (DRIED AT 180C),MG/L TOTAL DISSOLVED SOLIDS',
         "AMMONIA_mg_L" = "NITROGEN, AMMONIA, TOTAL (MG/L AS N)",
         "NITRATE_mg_L" =  "NITRATE NITROGEN, TOTAL (MG/L AS N)" ,
         "NITROGEN_KJELDAHL_TOTAL_00625_mg_L" = "NITROGEN, KJELDAHL, TOTAL, (MG/L AS N)" , 
         'PHOSPHORUS_TOTAL_ORTHOPHOSPHATE_70507_mg_L' =  'PHOSPHORUS,IN TOTAL ORTHOPHOSPHATE (MG/L AS P)',
         'TSS_mg_L'=  "TSS RESIDUE, TOTAL NONFILTRABLE (MG/L) TOTAL SUSPENDED SOLIDS" ) %>% 
  dplyr::select(any_of(names(conventionals)))
metals <- metals %>% 
  filter(!is.na(ID)) %>% 
  dplyr::rename(FDT_STA_ID = ID, FDT_DATE_TIME = `Date Time`, Ana_Sam_Mrs_Container_Id_Desc = `Container Description`) %>% 
  dplyr::select(FDT_STA_ID, FDT_DATE_TIME, Ana_Sam_Mrs_Container_Id_Desc, `Parm Name`, `Uncensored Value`) %>% 
  group_by(FDT_STA_ID, FDT_DATE_TIME, Ana_Sam_Mrs_Container_Id_Desc) %>% 
  pivot_wider(names_from = `Parm Name`, values_from = `Uncensored Value`) 
  
logi <- full_join(field, analyte, by = c('FDT_STA_ID', 'FDT_DATE_TIME')) %>% 
  full_join(metals,  by = c('FDT_STA_ID', 'FDT_DATE_TIME')) %>% #, 'Ana_Sam_Mrs_Container_Id_Desc'))%>% 
  mutate(Source = 'Lucy') %>% 
  dplyr::select(Source, everything())


conventionalsSmash <- dplyr::select(conventionals, any_of(names(logi))) %>% 
  mutate(Source = 'Emma')%>% 
  dplyr::select(Source, everything())

allData <- bind_rows(logi, conventionalsSmash) %>% 
  dplyr::select(Source, FDT_STA_ID, FDT_DATE_TIME, Ana_Sam_Mrs_Container_Id_Desc, everything()) %>% 
  arrange(FDT_STA_ID, FDT_DATE_TIME, Source, Ana_Sam_Mrs_Container_Id_Desc)

write.csv(allData, 'processedData/QA/allData.csv', row.names = F, na = '')


dplyr::select(DO = DO_mg_L, 
              pH = FDT_FIELD_PH,
              SpCond = FDT_SPECIFIC_CONDUCTANCE, 
              TN = NITROGEN_mg_L, 
              TP = PHOSPHORUS_mg_L, 
              NH4 = AMMONIA_mg_L, 
              NO3 = NITRATE_mg_L, 
              TKN = NITROGEN_KJELDAHL_TOTAL_00625_mg_L,
              `Ortho-P` = PHOSPHORUS_TOTAL_ORTHOPHOSPHATE_70507_mg_L,
              Turb = `TURBIDITY,LAB NEPHELOMETRIC TURBIDITY UNITS, NTU`,
              TSS = TSS_mg_L, 
              Na = `SODIUM, DISSOLVED (MG/L AS NA)`, 
              K = `POTASSIUM, DISSOLVED (MG/L AS K)`,
              Cl = CHLORIDE_mg_L,
              Sf = SULFATE_mg_L,
              `70331VFine` = `SSC%Finer`,
              SSCCOARSE = SSC_COARSE,
              SSCFINE =  SSC_FINE, 
              SSCTOTAL = SSC_TOTAL,
              # sediment ppm data hasn't been collected since 2011
              ARSENICppm = as.numeric(NA),
              BERYLLIUMppm = as.numeric(NA),
              CADMIUMppm = as.numeric(NA),
              CHROMIUMppm = as.numeric(NA),
              COPPERppm = as.numeric(NA),
              LEADppm = as.numeric(NA),
              MANGppm = as.numeric(NA),
              NICKELppm = as.numeric(NA),
              SILVERppm = as.numeric(NA),
              ZINCppm = as.numeric(NA),
              ANTIMONYppm = as.numeric(NA),
              ALUMINUMppm = as.numeric(NA),
              SELENIUMppm = as.numeric(NA),
              IRONppm = as.numeric(NA),
              MERCURYppm = as.numeric(NA),
              THALLIUMppm = as.numeric(NA),
              CALCIUM = `CALCIUM, DISSOLVED (MG/L AS CA)`,
              MAGNESIUM = `MAGNESIUM, DISSOLVED (MG/L AS MG)`,
              ARSENIC = `ARSENIC, DISSOLVED  (UG/L AS AS)`,
              BARIUM = `BARIUM, DISSOLVED (UG/L AS BA)`,
              BERYLLIUM = `BERYLLIUM, DISSOLVED (UG/L AS BE)`,
              CADMIUM = `CADMIUM, DISSOLVED (UG/L AS CD)`,
              CHROMIUM = `CHROMIUM, DISSOLVED (UG/L AS CR)`,
              COPPER = `COPPER, DISSOLVED (UG/L AS CU)`,
              IRON = `IRON, DISSOLVED (UG/L AS FE)`,
              LEAD = `LEAD, DISSOLVED (UG/L AS PB)`,
              MANGANESE = `MANGANESE, DISSOLVED (UG/L AS MN)`,
              THALLIUM = `THALLIUM, DISSOLVED (UG/L AS TL)`,
              NICKEL = `NICKEL, DISSOLVED (UG/L AS NI)`,
              SILVER = `SILVER, DISSOLVED (UG/L AS AG)`,
              ZINC = `ZINC, DISSOLVED (UG/L AS ZN)`,
              ANTIMONY = `ANTIMONY, DISSOLVED (UG/L AS SB)`,
              ALUMINUM = `ALUMINUM, DISSOLVED (UG/L AS AL)`,
              SELENIUM = `SELENIUM, DISSOLVED (UG/L AS SE)`,
              HARDNESS = `HARDNESS, CA MG CALCULATED (MG/L AS CACO3) AS DISSOLVED`,
              MERCURY = `MERCURY-TL,FILTERED WATER,ULTRATRACE METHOD NG/L`,
              `Hg-C` = `RMK_50091`)
              
              
                        