# ProbMon-Integrated-Reports
This repository hold scripts for analysis and reporting on VDEQ Probabilistic Monitoring chapters for the integrated report cycles.

## Background
Each assessment cycle, the freshwater probabilistic monitoring program is responsible for reportin gon data collected during the current IR window (6 years). 

Generally, the process for producing any given chapter involves:
- interviewing regional biologists to understand stations sampled and why cetain stations were not sampled to appropriately adjust weights during analyses
- querying field, chemistry, and metals data from ODS for stations sampled
- cleaning and summarizing records
- scraping USGS StreamStats API to delineate appropriate watersheds for sampled sites and QAing all watersheds manually
- run landcover metrics for delineated watersheds using most appropriate to date NLCD, tiger road, and block population data
- combine all above information with probabilistic dataset into super wide dataset for publication
- adjust sample weights based on bio feedback, running CDF and relative risk across multiple subpopulations
- feed new prob dataset to report Rmd to generate latest chapter, updating any necessary text to best reflect results
- submit finalized chapter to CO for inclusion in IR report
