# Read a komodo2 list from a key-value file.
read_komodo2_file <- function(file.path){
  df <- utils::read.table(file = file.path, header = FALSE,
                          strip.white = TRUE,
                          comment.char = "#", sep = '=' ,
                          col.names = c('Key', 'Value'),
                          stringsAsFactors = FALSE)

  outlist <- as.list(df$Value)

  # Coerce numeric values to "numeric"
  outlist <- lapply(outlist,
                    function(v){
                      if(!is.na(suppressWarnings(as.numeric(v)))) {
                        return(as.numeric(v))
                      } else return(v)
                    })

  names(outlist) <- as.character(df$Key)
  class(outlist) <- c("KOMODO2", "list")
  return(outlist)
}
