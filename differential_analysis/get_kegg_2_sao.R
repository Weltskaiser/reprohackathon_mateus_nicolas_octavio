# Generate useful Kegg data
# Has been done only once, the R script will now use our generated kegg_2_sao.json

library(KEGGREST)
library(jsonlite)
sao03010_json <- toJSON(keggGet("sao03010"), auto_unbox = TRUE)
sao00970_json <- toJSON(keggGet("sao00970"), auto_unbox = TRUE)
sao03010_list <- fromJSON(sao03010_json)
sao00970_list <- fromJSON(sao00970_json)
combined_list <- list(
  sao03010 = sao03010_list,
  sao00970 = sao00970_list
)
combined_json <- toJSON(combined_list, auto_unbox = TRUE, pretty = TRUE)
write(combined_json, "kegg_2_sao.json")