library(data.table)
library(dplyr)

dt_nov <- fread("McQuillin_data/annabel_121119_04pm.pheno")

# The old data is completely wrong.
dt_sept <- fread("McQuillin_data/correct_OPCRIT_data/OPCRIT_data_090919.pheno")
# There is one duplicate! Pick the first one...the OPCRIT numbers are differt but only slightly.
dt_sept <- dt_sept[-which(duplicated(dt_sept$id))[1],]
where <- which(dt_sept == 'na', arr.ind=TRUE)
dt_sept[[names(dt_sept)[where[2]]]][where[1]] <- NA
dt_sept[[names(dt_sept)[where[2]]]] <- as.integer(dt_sept[[names(dt_sept)[where[2]]]])
dt_july <- fread("McQuillin_data/correct_OPCRIT_data/annabel20jul19_200719_11pm.pheno")

check_differences <- function(younger_dt, older_dt, merge=FALSE, verbose=FALSE)
{
	# Dataset exclusive samples
	younger_dt_setdiff <- younger_dt %>% filter(id %in% setdiff(younger_dt$id, older_dt$id))
	older_dt_setdiff <- older_dt %>% filter(id %in% setdiff(older_dt$id, younger_dt$id))

	# Merge the remainder
	younger_dt <- younger_dt %>% filter(id %in% older_dt$id)
	older_dt <- older_dt %>% filter(id %in% younger_dt$id)
	dt <- merge(older_dt, younger_dt, by='id')

	percent_diff <- rep(0, 90)
	for(i in seq(1, 90))
	{
		if(i < 10) {
			i_str <- paste0('0', i)
		} else {
			i_str <- as.character(i)
		}

		not_equal <- which(dt[[paste0("opcrit.", i_str, ".x")]] != dt[[paste0("opcrit.", i_str, ".y")]])
		dt_tmp <- dt %>% 
			filter(
				!is.na(dt[[paste0("opcrit.", i_str, ".x")]]) & 
				!is.na(dt[[paste0("opcrit.", i_str, ".y")]])
			)
		
		if (merge) {
			# Replace the older value with the younger value
			dt[[paste0("opcrit.", i_str, ".y")]][not_equal] <- dt[[paste0("opcrit.", i_str, ".x")]][not_equal]
			# Replace missing in younger with non missing in older
			na_in_younger <- which(is.na(dt[[paste0("opcrit.", i_str, ".x")]]))
			dt[[paste0("opcrit.", i_str, ".x")]][na_in_younger] <- dt[[paste0("opcrit.", i_str, ".y")]][na_in_younger]
			# Replace missing in older with non missing in younger
			na_in_older <- which(is.na(dt[[paste0("opcrit.", i_str, ".y")]]))
			dt[[paste0("opcrit.", i_str, ".y")]][na_in_older] <- dt[[paste0("opcrit.", i_str, ".x")]][na_in_older]
			
			if (verbose) {
				# Check that everything now matches
				cat(all(dt[[paste0("opcrit.", i_str, ".y")]] == dt[[paste0("opcrit.", i_str, ".x")]], na.rm=TRUE), '\n')
				cat(all(is.na(dt[[paste0("opcrit.", i_str, ".y")]]) == is.na(dt[[paste0("opcrit.", i_str, ".x")]])), '\n')
			}
		}

		percent_diff[i] <- length(not_equal)/nrow(dt_tmp)

		if (verbose) {
			print(dt_tmp %>% 
				filter(dt_tmp[[paste0("opcrit.", i_str, ".x")]] != dt_tmp[[paste0("opcrit.", i_str, ".y")]]) %>% 
				select(paste0("opcrit.", i_str, ".x"), paste0("opcrit.", i_str, ".y"))
			)
		}
	}

	# Finally, check that the inferred BP subtype is the same in each of these three files.
	# Need to check if that colname exists, and if it does, check the difference.
	if("bp_type_inferred.x" %in% names(dt))
	{
		not_equal <- which(dt[["bp_type_inferred.x"]] != dt[["bp_type_inferred.y"]])
		dt_tmp <- dt %>% filter(!is.na(dt[["bp_type_inferred.x"]]) & !is.na(dt[["bp_type_inferred.y"]]))
		percent_bp_type_diff <- length(not_equal)/nrow(dt_tmp)
		cat("Percentage different between subtypes at the non-missing values...\n")
		cat(percent_bp_type_diff, '\n')
	}

	# Finally, add in the information that is only present in one of the datasets.
	if (merge) {

		dt <- dt[, -grep("\\.x", names(dt))]
		names(dt) <- gsub("\\.y", "", names(dt))

		# names(dt)[-grep("\\.x", names(dt))]
		nrow_merge <- nrow(dt)
		nrow_younger_dt_setdiff <- nrow(younger_dt_setdiff)
		nrow_older_dt_setdiff <- nrow(older_dt_setdiff)

		dt <- merge(dt, younger_dt_setdiff, all=TRUE)
		dt <- merge(dt, older_dt_setdiff, all=TRUE)
		
		# Check
		print(nrow(dt) == (nrow_merge + nrow_younger_dt_setdiff + nrow_older_dt_setdiff))
	}

	cat("Percentage different at the non-missing values...\n")
	print(percent_diff)

	cat("Size of the merged dataset...\n")
	cat(nrow(dt),"\n")

	return(dt)
}

dt_july_sept <- check_differences(dt_sept, dt_july, merge=TRUE)
dt_july_sept_nov <- check_differences(dt_nov, dt_july_sept, merge=TRUE)

fwrite(dt_july_sept_nov, "McQuillin_data/harmonised_opcrit_mcquillin_july_sept_nov.tsv", sep='\t')
