Files used for data preparation and analysis 

Data Preparation: 
- DrugProcessing:
  contains all steps to process immunosuppressive drug regimens
- Prophylax/Infection:
  contains the processing of prophylaxis, viral/bacterial/fungal infections and complications, and antibiotic treatments
- Basline Characteristics/Lab Values
  contains all extended baseline values such as BMI, blood type,... and the unit conversion of blood/urine values

Analysis: 
- diversity measurements: contains all alpha and beta diversity analysis for the different group comparison (Fig.1,2,3,5)
- prepMasterTable_paperRun:
   contains the generating of the meta input table (master table) adding baseline, lab value, antibiotic treatment, and viral infections based on the tables generated in the other files
   afterwards all association analysis using metadeconfoundR are listed together with different analysis (Fig.2,3,5)
- LongDat_completecohort: contains the long dat run of all samples post KT (Fig. 1)
- ReferenceCohortAnalysis: contains the analysis of the reference CKD study and comparison to the DZIF results (Fig. 6)
- heattree: visualization of differentially abundant bacteria on sev. taxa levels (Fig. 3)
- picrust: contains the analysis of functional capacity predicted by picrust, including ko preselection, and visualization of GOmixer module abundances (Fig.4)

- function.R contains several functions used to run the assocaition analysis with metadeconfoundR (preparation of meta input, feature input,...) + visualization functions (volcanos,heatmaps,...)
  
