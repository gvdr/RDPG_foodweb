#Code use

The main function in `RDPG_function.R` is
```R
get_traits(Adj, dataframe=TRUE, names=FALSE)
```
which computes the full rank abstract functional traits of a food web.
It receives as input the adjacency matrix of the food web and returns
either a matrix (`dataframe=FALSE`) or a data.frame (`dataframe=TRUE`)
of abstract functional traits vectors of each species (both as
prey and as predators).
The function gets species names from the adjacency matrix (`names=FALSE`) or from
a vector of characters (`names=vec_of_names`).

The code to perform the various analysis presented in the paper is
in the Ancillary_Code folder, and I'll be adding readmes and howtos
shortly. If you want to know more details, please contact me!
