# Phenotype data curation steps

`00_merge_phenotypic_information.r`, `01_clean_phenotypic_information.r`, `02_merge_phenotypic_subtype_information.r`, `03_manual_curation_phenotypic_subtype_information.r`, `04_get_vcf_names.py`, `05_combine_pheno_file_and_vcf_names.r`, `06_checking_mcquillin_data.r` and `07_merge_psychosis_information.r` check the phenotype information from BSP and read in the latest subphenotype information from the collaborators. We carefully go through and ensure things match correctly, and that not much has changed as a sanity check.

This then generates BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv

Briefly, these are the steps that each script performs:

`00_merge_phenotypic_information.r`: Reads in two files from BSP (couldn't export all together as there are too many samples for the system), remove columns that are not required and outputs a smaller file `BIP_phenotype_information.tsv`. Compare this file to another set of files received from BSP (rather than Sam and Sinead) to ensure that everything is consistent. It is.

`01_clean_phenotypic_information.r`: This is an initial clean-up step to fix spelling errors and establish consistent naming conventions across the required columns. Where possible, we define subtype information and create a new column `PHENOTYPE_FINE`. The resultant phenotype columns in the file are `PHENOTYPE_COARSE` and `PHENOTYPE_FINE`, which have the following possible entries:

`PHENOTYPE_COARSE`
1. Bipolar Disorder
2. Control
3. Other
4. Schizoaffective
5. Schizophrenia
6. Unknown

`PHENOTYPE_FINE`
1. Bipolar Disorder
2. Bipolar Disorder 1
3. Bipolar Disorder 2
4. Bipolar Disorder NOS
5. Control
6. Other 
7. Schizoaffective
8. Schizophrenia
9. Unknown

`02_merge_phenotypic_subtype_information.r`: Here we read in raw files from the collaborators to attempt to increase the number of individuals for which we have subtype information. This again involved a number of data cleaning steps. The resultant file is called `BIP_phenotype_information_cleaned_new_subtype_information_added.tsv`.

From this file, we then noticed a few discrepancies between the coarse and fine phenotype definitions - these are manually curated in `03_manual_curation_phenotypic_subtype_information.r` and written to `BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised.tsv.`

Finally we checked the names of phenotypes in the vcf and compared them to those in this final phenotype file. 

To do this we just run `04_get_vcf_names.py` to read in the initial matrix table and export the sample IDs

We then read in the sample IDs and compare them to the cleaned phenotype IDs to identify any differences and obvious name changes.

The resultant cleaned names are then written to `BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv` in `05_combine_pheno_file_and_vcf_names.r`.

We also now have a further collection of files taken from the ICCBD project. We read these in and contrast the phenotype and subphenotype data with what is present here. We also incorporate psychosis information as a column in the final output phenotype file which we use for analysis. This is done in `07_merge_psychosis_information.r`. We run a further collection of checks and ensure that everything is consistent. Within `07_merge_psychosis_information.r`, `06_checking_mcquillin_data.r` is run to harmonise the McQuillin data.

The resultant final file is:
BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv

Note: Annabel noticed a discrepancy in the McQuillan data - there are only binary variable for much of the OPCRIT data which should have multiple categories. She sent me an older version of the data. I compared to see if this induces a difference in the definition of psychoses for this subset of the dataset. I check in `06_checking_mcquillin_data.r`, and correct the errors.

When updating the phenotype files, there are a number of scripts that need to be run to ensure that everything is updated both for the purposes of input into analyses, and for display of the QC pipeline. They are:

* `00_create_phenotype_keytable.py`: This creates the up to date hail table that is used to annotate with the latest phenotypes.
* `17_update_phenotypes.py`: This annotates the QC pipeline output with the latest phenotype information.
These two alterations are all that is required for downstream analyses.

To update the relevant numbers in the QC website, we additionally need to rerun:

* `03_initial_sample_qc_plot.r`
* `03_initial_sample_qc_filter.r`
* `15_final_sample_qc_plot.r`
* `15_filter_final_sample_qc.r`

# Raw data descriptions in the phenotypes files.

We begin with the subtype information that was originally provided to the Broad and pulled from the BSP database, this is then augmented with the data below in order to fill in subtype information which is not present in BSP.

Landen data:
Gender
Case status
BIP1 (TRUE/FALSE)
BIP2 (TRUE/FALSE)
BIP NOS (TRUE/FALSE)
SCZ BIP (TRUE/FALSE)
AOO and Psychosis information (TRUE/FALSE)

Reif data:
Sample name
Sample name (translated to remove umlauts etc)
Psychosis (TRUE/FALSE)
Bipolar classification (BIP1, BIP2, BIP NOS, UNKNOWN)

Yolken data: A lot more information
Primary disease (Bipolar Disorder/Control)
Age 
Gender
Ethnicity
Highest level of maternal education (momedu)
Level of proband education (educyrs)
Repeatable Battery for the Assessment of Neuropsychological Status Total Score (rbansatotal)
DSM IV Diagnosis (dsmivdiag)
Care Setting on Enrollment (care_setting_3)	
Young Mania Rating Scale (ymaniars)	
Brief Psychiatric Rating Scale (bprs)	
Hamilton Depression Rating Screen (hamtotscrn)
Cigarette Smoker (cigsmoker)	
Receiving Any Atypical Anti-psychotic Medication (anyatypical)
Receiving Lithium (lithium)
Receiving Antidepressant Medication	(antidepressant)
Age of first mood symptoms (moodonset)

Ophoff data:
Primary disease
Diagnosis (BIP1, BIP2, BIP NOS etc)

McQuillin data:
Primary disease
Diagnosis (BP_type: BIP1, BIP2, BIP NOS etc)
OPCRIT data - multiple files, the latest is incorrect.
Can generate psychosis TRUE/FALSE and AOO to match Sweden, MGH and Cardiff from this.

Smoller data: Split into lots of files, analysed by Chia-yen.
FID IID bp bp1 bp2 bp_nos bp_scz panic anxiety sub_dep alc_dep adhd asd dev_del int_dis epilepsy migraine psychosis suicidality age_fi_24 age_fi_40 age_fs_24 age_fs_40 age_d_24 age_d_40 fhx_bp fhx_mood fhx_psy li delivery_num postp_psychosis postp_dep peri_mood neurocog
