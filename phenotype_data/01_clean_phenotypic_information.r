library(data.table)
library(dplyr)

df <- fread("BIP_phenotype_information.tsv")
# Looks to be a single sample that got an SM ID instead of external ID
df <- df %>% mutate(SAMPLE_ALIAS = SAMPLE_ALIAS.x) %>% select(-c(SAMPLE_ALIAS.x, SAMPLE_ALIAS.y))

# Now, need to clean things to get the Location, PI, and Phenotype information.

df <- df %>% mutate(LOCATION = COLLECTION, PI = COLLECTION,
	PHENOTYPE_FINE = PRIMARY_DISEASE)

# First, Phenotypes

# Here is the list of different strings present in the PRIMARY_DISEASE column

#  [1] ""                                                                                                
#  [2] "ADHD"                                                                                            
#  [3] "ADHD and personality disorder"                                                                   
#  [4] "Adjustment disorder, unspecified"                                                                
#  [5] "Agoraphobia without a history of panic disorder"                                                 
#  [6] "Alcohol abuse"                                                                                   
#  [7] "Alcohol dependence"                                                                              
#  [8] "BD after head trauma"                                                                            
#  [9] "Bipoalr 1 Disorder"                                                                              
# [10] "Bipoalr 2 Disorder"                                                                              
# [11] "Bipoalr Disorder"                                                                                
# [12] "Bipoalr Disorder 1"                                                                              
# [13] "Bipoalr Disorder 2"                                                                              
# [14] "Bipolar"                                                                                         
# [15] "Bipolar 1"                                                                                       
# [16] "Bipolar 1 Disorder"                                                                              
# [17] "Bipolar 2"                                                                                       
# [18] "Bipolar 2 Disorder"                                                                              
# [19] "Bipolar Affective Disorder"                                                                      
# [20] "Bipolar Disease"                                                                                 
# [21] "Bipolar Disease 1"                                                                               
# [22] "Bipolar disorder"                                                                                
# [23] "Bipolar Disorder"                                                                                
# [24] "Bipolar Disorder 1"                                                                              
# [25] "Bipolar Disorder 2"                                                                              
# [26] "Bipolar Disorder I"                                                                              
# [27] "Bipolar Disorder I Cykloid psykos"                                                               
# [28] "Bipolar Disorder II"                                                                             
# [29] "Bipolar disorder NOS"                                                                            
# [30] "Bipolar Disorder NOS"                                                                            
# [31] "Bipolar Disorder Type 1"                                                                         
# [32] "Bipolar Disorder UNS"                                                                            
# [33] "Bipolar Disorder Unspecified Type"                                                               
# [34] "Bipolar I"                                                                                       
# [35] "Bipolar I Disorder"                                                                              
# [36] "Bipolar II"                                                                                      
# [37] "Bipolar II Disorder"                                                                             
# [38] "Bipolar NOS"                                                                                     
# [39] "Bipolaraffective Disorder"                                                                       
# [40] "Borderline personality disorder"                                                                 
# [41] "BP"                                                                                              
# [42] "BP 3 Endast manier / hypomanier inducerade av antidepressiva"                                    
# [43] "BP 3&amp;frac12; Hypomani / depression endast associerade till substansmissbruk"                 
# [44] "BP 4 Depressioner och hypertymt temprament"                                                      
# [45] "BP 5 Unipol\x8ara (d/D) &amp;ge;5 recidiverande depressioner med hereditet f\x9ar bipol\x8ar sjd"
# [46] "BP UNS"                                                                                          
# [47] "BPD"                                                                                             
# [48] "Cannabis dependence"                                                                             
# [49] "CON"                                                                                             
# [50] "control"                                                                                         
# [51] "Control"                                                                                         
# [52] "CONTROL"                                                                                         
# [53] "Controls"                                                                                        
# [54] "Cyclothymic disorder"                                                                            
# [55] "Depressive disorder"                                                                             
# [56] "Depressive disorder and borderline personality disorder"                                         
# [57] "Depressive disorder NOS"                                                                         
# [58] "Depressive disorder; substance induced mood disorder"                                            
# [59] "Diagnosis unclear"                                                                               
# [60] "Diagnosis unclear (possibly dysthymic disorder)"                                                 
# [61] "Diagnosis unclear. Possibly depressive disorder. Dementia"                                       
# [62] "Dysthymic disorder"                                                                              
# [63] "Hypochondriasis"                                                                                 
# [64] "Induced BD"                                                                                      
# [65] "Induced mood disorder"                                                                           
# [66] "Major depressive disorder"                                                                       
# [67] "Major Depressive Disorder Recurrent"                                                             
# [68] "Major depressive disorder, recurrent, in full remission"                                         
# [69] "Major depressive disorder, single episode"                                                       
# [70] "Major depressive disorder, single episode, in full remission"                                    
# [71] "Other PsychosisManic Depressive Disorder"                                                        
# [72] "OtherBipolar Psychoses"                                                                          
# [73] "OtherHypomania"                                                                                  
# [74] "OtherMania/Hypomania"                                                                            
# [75] "Panic disorder without agoraphobia"                                                              
# [76] "Personality disorder NOS"                                                                        
# [77] "PSYCH NOS"                                                                                       
# [78] "Psychosis"                                                                                       
# [79] "psychosis at time of sample collection"                                                          
# [80] "Psychotic disorder NOS; Specific phobia; GAD"                                                    
# [81] "Schiz-Affective"                                                                                 
# [82] "Schizo-Affective Disorder"                                                                       
# [83] "Schizoaffective Bipolar Type"                                                                    
# [84] "Schizoaffective disorder"                                                                        
# [85] "Schizoaffective Disorder"                                                                        
# [86] "Schizoaffective Disorder Bipolar Type"                                                           
# [87] "Schizoaffective Disorder Depressive Type"                                                        
# [88] "Schizoaffective Disorder Manic"                                                                  
# [89] "Schizoaffective Schizophrenia"                                                                   
# [90] "Schizophrenia"                                                                                   
# [91] "Schizophreniform Disorder/NOS"                                                                   
# [92] "Schizophreniform Psychosis"                                                                      
# [93] "Schizotypal Personality Disorder"                                                                
# [94] "Social phobia"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in% 
	c("Bipolar",
	"Bipolar Affective Disorder",
	"Bipolaraffective Disorder",
	"BPD",
	"BP",
	"Schizoaffective Bipolar Type",
	"Schizoaffective Disorder Bipolar Type",
	"Schizoaffective Disorder Manic",
	"Bipoalr Disorder",
	"Bipolar Disease",
	"Bipolar disorder",
	"Bipolar Disorder")] <- "Bipolar Disorder"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in% 
	c("Bipolar Disorder NOS",
	"Bipolar NOS", 
	"Bipolar Disorder Unspecified Type",
	"BP UNS",
	"Bipolar disorder NOS",
	"Bipolar Disorder UNS")] <- "Bipolar Disorder NOS"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in%
	c("Schiz-Affective", 
	"Schizo-Affective Disorder",
	"Schizoaffective Disorder",
	"Schizoaffective disorder",
	"Schizoaffective Disorder Depressive Type",
	"Schizoaffective Schizophrenia")] <- "Schizoaffective"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in% 
	c("ADHD",
	"ADHD and personality disorder",
	"Adjustment disorder, unspecified",
	"Agoraphobia without a history of panic disorder",
	"Alcohol abuse",
	"Alcohol dependence",
	"BD after head trauma",
	"BP 3&amp;frac12; Hypomani / depression endast associerade till substansmissbruk", # This is bipolar only after substance abuse.
	"BP 4 Depressioner och hypertymt temprament", # Depression and hypertensive temprament
	"BP 5 Unipol\x8ara (d/D) &amp;ge;5 recidiverande depressioner med hereditet f\x9ar bipol\x8ar sjd", # Relapsing depressiokn with heredity for bipolar.
	"Borderline personality disorder",
	"BPD",
	"Cannabis dependence",
	"Cyclothymic disorder",
	"Depressive disorder",
	"Depressive disorder and borderline personality disorder",
	"Depressive disorder NOS",
	"Depressive disorder; substance induced mood disorder",
	"Dysthymic disorder",
	"Hypochondriasis",
	"Induced BD",
	"Induced mood disorder",
	"Major depressive disorder",
	"Major Depressive Disorder Recurrent",
	"Major depressive disorder, recurrent, in full remission",
	"Major depressive disorder, single episode",
	"Major depressive disorder, single episode, in full remission",
	"Other PsychosisManic Depressive Disorder",
	"OtherBipolar Psychoses",
	"OtherMania/Hypomania",
	"OtherHypomania",
	"Panic disorder without agoraphobia",
	"Personality disorder NOS",
	"PSYCH NOS",
	"Psychosis",
	"psychosis at time of sample collection",
	"Psychotic disorder NOS; Specific phobia; GAD",
	"Social phobia",
	"Schizophreniform Disorder/NOS",
	"Schizophreniform Psychosis",
	"Schizotypal Personality Disorder")] <- "Other"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in%
	c("CON",
	"control",
	"Control",
	"CONTROL",
	"Controls")]<- "Control"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in%
	c("",
	"Diagnosis unclear",
	"Diagnosis unclear (possibly dysthymic disorder)",
	"Diagnosis unclear. Possibly depressive disorder. Dementia",
	"CONTROL",
	"Controls")] <- "Unknown"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in%
	c("Bipoalr 2 Disorder",
	"Bipoalr Disorder 2",
	"Bipolar 2",
	"Bipolar 2 Disorder",
	"Bipolar Disorder 2",
	"Bipolar Disorder II",
	"Bipolar II",
	"Bipolar II Disorder")] <- "Bipolar Disorder 2"

df$PHENOTYPE_FINE[df$PHENOTYPE_FINE %in%
	c("BP 3 Endast manier / hypomanier inducerade av antidepressiva",
	"Bipoalr 1 Disorder",
	"Bipoalr Disorder 1",
	"Bipolar 1",
	"Bipolar 1 Disorder",
	"Bipolar Disease 1",
	"Bipolar Disorder 1",
	"Bipolar Disorder I",
	"Bipolar Disorder I Cykloid psykos",
	"Bipolar Disorder Type 1",
	"Bipolar I",
	"Bipolar I Disorder")] <- "Bipolar Disorder 1"

df <- df %>% mutate(PHENOTYPE_COARSE = PHENOTYPE_FINE)
df$PHENOTYPE_COARSE[df$PHENOTYPE_FINE %in%
	c("Bipolar Disorder 1", "Bipolar Disorder 2", "Bipolar Disorder NOS")] <- "Bipolar Disorder"

# Phenotype information complete!
# Next, let's determine the Location and PI from the Collection information.
df$COLLECTION <- gsub("^Neuropsychiatry / ", "", df$COLLECTION)
df$PI <- gsub(" \\(.*", "", df$COLLECTION)
df$PI[df$PI == "Michael O&apos;Donovan"] <- "Michael ODonovan"
df$PI[df$PI == "Fernando S. Goes, M.D."] <- "Fernando Goes"
df$PI[df$PI == "Adolfsson"] <- "Rolf Adolfsson"
df$PI[df$PI == "David St. Clair"] <- "David St Clair"
df$PI[df$PI == "McQuillin and Gurling"] <- "Andrew McQuillin"

df$LOCATION <- gsub(".*\\((.*)\\).*", "\\1", df$COLLECTION)
df$INSTITUTION <- df$LOCATION
df$LOCATION[df$LOCATION == "University of W\xfc\xbe\x99\xa3\xa4\xbcrzburg/Germany"] <- "Wurzburg, GER"
df$INSTITUTION[df$LOCATION == "Wurzburg, GER"] <- "University of Wurzburg"
df$LOCATION[df$LOCATION %in% c("Karolinska", "Karolinska/Sweden")] <- "Stockholm, SWE"
df$INSTITUTION[df$INSTITUTION %in% c("Karolinska", "Karolinska/Sweden")] <- "Karolinska Institute"
df$LOCATION[df$LOCATION == "VU University Medical Center- Amsterdam"] <- "Amsterdam, NED"
df$INSTITUTION[df$INSTITUTION == "VU University Medical Center- Amsterdam"] <- "VU University Medical Center"
df$LOCATION[df$LOCATION == "MGH"] <- "Boston, USA"
df$LOCATION[df$LOCATION == "Sweden"] <- "Umea, SWE"
df$INSTITUTION[df$INSTITUTION == "Sweden"] <- "Umea University"
df$LOCATION[df$LOCATION == "Johns Hopkins School of Medicine/USA"] <- "Baltimore, USA"
df$LOCATION[df$LOCATION %in% c("Edinburgh", "Aberdeen", "Cardiff")] <- paste0(df$LOCATION[df$LOCATION %in% c("Edinburgh", "Aberdeen", "Cardiff")], ", UK")
df$LOCATION[df$LOCATION %in% c("UCLA", "University of California Los Angeles")] <- "Amsterdam, NED"
df$LOCATION[df$LOCATION == "UK" & df$PI == "Willem Ouwehand"] <- "Cambridge, UK"
df$LOCATION[df$LOCATION == "UCL"] <- "London, UK"
df$LOCATION[df$LOCATION == "UCL"] <- "London, UK"
df$LOCATION[df$LOCATION == "Dublin"] <- "Dublin, IRE"
df$INSTITUTION[df$LOCATION == "Cambridge, UK" & df$PI == "Willem Ouwehand"] <- "University of Cambridge"
df$INSTITUTION[df$INSTITUTION %in% c("UCLA", "University of California Los Angeles")] <- "UCLA"
df$INSTITUTION[df$INSTITUTION == "Johns Hopkins School of Medicine/USA"] <- "Johns Hopkins School of Medicine"
df$INSTITUTION[df$INSTITUTION == "Edinburgh"] <- "University of Edinburgh"
df$INSTITUTION[df$INSTITUTION == "Aberdeen"] <- "University of Aberdeen"
df$INSTITUTION[df$INSTITUTION == "Cardiff"] <- "University of Cardiff"
df$INSTITUTION[df$INSTITUTION == "Dublin"] <- "National University of Ireland Galway"
df$INSTITUTION[df$INSTITUTION == "UCLA"] <- "University of California Los Angeles"
df$INSTITUTION[df$INSTITUTION == "UCL"] <- "University College London"
df$INSTITUTION[df$INSTITUTION == "MGH"] <- "Massachusetts General Hospital"

df$GENDER[df$GENDER %in% c("Not Reported", "Unknown", "")] <- "Unknown"

df <- df %>% select(-c("COLLECTION", "PRIMARY_DISEASE"))

# Final stage is to rename a small subset of the samples that got name changes in the .vcf

# df_pheno$SAMPLE_ALIAS[which(!(df_pheno$SAMPLE_ALIAS %in% df$s))]
#  [1] "1749.0"        "1968.0"        "1972.0"        "1974.0"       
#  [5] "1976.0"        "1979.0"        "1981.0"        "1982.0"       
#  [9] "1984.0"        "1985.0"        "1986.0"        "1987.0"       
# [13] "1988.0"        "1990.0"        "1991.0"        "1993.0"       
# [17] "1995.0"        "1996.0"        "1997.0"        "1998.0"       
# [21] "431-BG00852 D" "C0340_a"       "G1020,"        "zas8464_ 2"  

df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1749.0"] <- "1749_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1968.0"] <- "1968_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1972.0"] <- "1972_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1974.0"] <- "1974_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1976.0"] <- "1976_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1979.0"] <- "1979_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1981.0"] <- "1981_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1982.0"] <- "1982_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1984.0"] <- "1984_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1985.0"] <- "1985_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1986.0"] <- "1986_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1987.0"] <- "1987_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1988.0"] <- "1988_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1990.0"] <- "1990_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1991.0"] <- "1991_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1993.0"] <- "1993_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1995.0"] <- "1995_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1996.0"] <- "1996_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1997.0"] <- "1997_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "1998.0"] <- "1998_0"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "431-BG00852 D"] <- "431-BG00852_D"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "G1020,"] <- "G1020_"
df$SAMPLE_ALIAS[df$SAMPLE_ALIAS == "zas8464_ 2"] <-  "zas8464__2"

fwrite(df, file = "BIP_phenotype_information_cleaned.tsv", sep='\t')

