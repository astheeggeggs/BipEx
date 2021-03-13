library(dplyr)
library(data.table)

create_list <- function(dt, first_authors, last_authors)
{
	names(dt) <- gsub(" ", "_", names(dt))
	dt <- dt %>% mutate(name = gsub("  ", " ", paste(First_Name, Middle_Initial, Last_Name)))
	dt_first <- dt %>% filter(name %in% first_authors)
	dt_last <- dt %>% filter(name %in% last_authors)
	dt_middle <- dt %>% filter(!(name %in% first_authors | name %in% last_authors))
	dt_first <- dt_first %>% mutate(order = NA)
	dt_last <- dt_last %>% mutate(order = NA)
	i <- 1
	for (author in first_authors) { dt_first$order[grep(author, dt_first$name)] <- i; i <- i + 1}
	dt_middle <- dt_middle %>% arrange(Last_Name)
	dt_middle <- dt_middle %>% mutate(order = (seq(1, nrow(dt_middle)) + i - 1))
	i <- i + nrow(dt_middle)
	for (author in last_authors) { dt_last$order[grep(author, dt_last$name)] <- i; i <- i + 1}

	dt <- rbind(dt_first, dt_middle, dt_last)
	cols_for_affil <- c("order", "name", grep("Affiliation", names(dt), value=TRUE))
	id_vars <- c("order", "name")
	measure_vars <- setdiff(cols_for_affil, id_vars)
	dt_affil <- dt %>% arrange(order) %>% mutate(affiliation_used=FALSE) %>% select(one_of(cols_for_affil))
	long_dt_affil <- melt(dt_affil, id.vars=id_vars, measure.vars=measure_vars, value.name="Affiliation")
	setkeyv(long_dt_affil, c("order", "variable"))

	affiliations <- setdiff(unique(long_dt_affil$Affiliation), "")
	affil_list <- list()
	i <- 1
	for (affil in affiliations) {
		affil_list[[affil]] <- i
		i <- i+1
	}

	long_dt_affil <- long_dt_affil %>% filter(Affiliation != "") %>% mutate(id = unlist(affil_list[long_dt_affil$Affiliation]))
	dt <- merge(dt, dcast(long_dt_affil, name ~ ., fun.agg = function(x) paste(x, collapse=","), value.var="id"), by="name") %>% arrange(order)

	affiliations_to_print <- paste(paste(seq(1,length(affil_list)), names(affil_list), sep=". "), collapse="\n")
	authors_to_print <- paste(paste(dt$name, dt$`.`, sep=""), collapse=", ")
	return(list(names = authors_to_print, affils = affiliations_to_print))
}

author_list <- "Dalio_Author_List.tsv"
dt <- fread(author_list)
first_authors <- c("Duncan S Palmer", "Daniel P Howrigan", "Sinéad B Chapman")
last_authors <- "Benjamin M Neale"
authors <- create_list(dt, first_authors, last_authors)
cat(paste0(authors$names, "\n\n", authors$affils, "\n"))