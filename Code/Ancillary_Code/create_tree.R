library("rgbif")
devtools::install_github("taxize", "ropensci")
library("taxize")
library("ape")

setwd("H:/Projects/Stouffer_RDPFW/Data/Raw")
ids = as.matrix(read.csv("H:/Projects/Stouffer_RDPFW/Data/Raw/W_ids.csv"))

taxa.names = ids[,"Taxonomy"]
temp = gnr_resolve(taxa.names[1:10])
tsn = get_tsn(searchterm = taxa.names[1:100],searchtype = "scientific", accepted = TRUE)
tsn

out <- get_gbifid(taxa.names, verbose = FALSE)
cl <- classification(out, verbose = TRUE)
cl_naout <- cl[!is.na(cl)]

cl_col <- cbind(cl)

cl_col_fa <- cl_col

cl_col_fa$genus <- as.factor(cl_col$genus)
cl_col_fa$family <- as.factor(cl_col$family)
cl_col_fa$order <- as.factor(cl_col$order)
cl_col_fa$phylum <- as.factor(cl_col$phylum)
cl_col_fa$kingdom <- as.factor(cl_col$kingdom)
cl_col_mat <- data.frame(cl_col_fa)

cl_final <- cl_col_mat[complete.cases(cl_col_mat),]

tr <- as.phylo.formula(~kingdom/phylum/order/family/genus/species, data=cl_final)
tr_bl <- compute.brlen(tr,1)
plot(tr_bl)

classification(taxa.names, db = "itis")

get_uid(sciname = "Actinocyclus spiritus")
gnr_resolve("Actinocyclus spiritus")


tnrs(query = "Actinocyclus spiritus", source = "iPlant_TNRS")[ , -c(5:7)]

write.table(cl_col, file = "test.csv", sep = "/",row.names = FALSE, col.names = FALSE, quote = FALSE)
taxah <- read.csv("test.csv", header=FALSE)
