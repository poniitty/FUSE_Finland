# resolve_higher_taxonomy_key(186341026)
# KEY <- 186341026
resolve_higher_taxonomy_key <- function(KEY = NULL, lib.loc = .libPaths()){
  
  if(is.null(KEY)){
    df <- data.frame(orig_key = NA,
                     taxon_name = NA,
                     genus = NA,
                     genus_key = NA,
                     family = NA,
                     family_key = NA,
                     note_higher = "Missing taxon key")
    return(df)
    stop()
  }
  
  KEY <- as.numeric(KEY)
  
  if(length(KEY) == 0 | is.na(KEY)){
    df <- data.frame(orig_key = NA,
                     taxon_name = NA,
                     genus = NA,
                     genus_key = NA,
                     family = NA,
                     family_key = NA,
                     note_higher = "Non-numeric taxon key")
    return(df)
    stop()
  }
  
  pckg <- require(rgbif, lib.loc = lib.loc)
  if(!pckg){stop("package ´rgbif´ not found")}
  
  ORIGNAME <- NA
  GENUS <- NA
  GENUSKEY <- NA
  FAMILY <- NA
  FAMILYKEY <- NA
  NOTE <- "OK"
  
  e <- try({
    nl <- name_usage(key = KEY)$data
    ORIGNAME <- nl$scientificName
    
    if(nrow(nl) == 1){
      if("genusKey" %in% names(nl)){
        nl <- name_usage(key = nl$genusKey)$data
        if(nrow(nl) == 1){
          if(nl$taxonomicStatus == "ACCEPTED"){
            GENUS <- nl$genus
            GENUSKEY <- nl$key
          } else {
            if("acceptedKey" %in% names(nl)){
              nl <- name_usage(key=nl$acceptedKey)$data
              if(nl$taxonomicStatus == "ACCEPTED"){
                GENUS <- nl$genus
                GENUSKEY <- nl$key
              } else {
                stop()
              }
            } else {
              if(nl$taxonomicStatus == "DOUBTFUL"){
                stop("DOUBTFUL")
              } else {
                stop()
              }
            }
          }
        } else {
          stop()
        }
      } else {
        if(nl$rank %in% c("SPECIES","FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
          GENUS <- strsplit(nl$scientificName, " ")[[1]][1]
          NOTE <- "GENUS EXTRACTED FROM SPECIES NAME"
        } else {
          stop()
        }
      }
    } else {
      stop()
    }
    
    if(!is.na(GENUSKEY) | !is.na(GENUS)){
      
      if(nrow(nl) == 1){
        if("familyKey" %in% names(nl)){
          nl <- name_usage(key = nl$familyKey)$data
          if(nrow(nl) == 1){
            if(nl$taxonomicStatus == "ACCEPTED"){
              FAMILY <- nl$family
              FAMILYKEY <- nl$key
            } else {
              if("acceptedKey" %in% names(nl)){
                nl <- name_usage(key=nl$acceptedKey)$data
                if(nl$taxonomicStatus == "ACCEPTED"){
                  FAMILY <- nl$family
                  FAMILYKEY <- nl$key
                } else {
                  stop()
                }
              } else {
                if(nl$taxonomicStatus == "DOUBTFUL"){
                  stop()
                } else {
                  stop()
                }
              }
            }
          } else {
            stop()
          }
        }
      } else {
        stop()
      }
    }
  })
  
  if(class(e) == "try-error"){
    df <- data.frame(orig_key = KEY,
                     taxon_name = NA,
                     genus = NA,
                     genus_key = NA,
                     family = NA,
                     family_key = NA,
                     note_higher = "FAILED")
    return(df)
  } else {
    df <- data.frame(orig_key = KEY,
                     taxon_name = ORIGNAME,
                     genus = GENUS,
                     genus_key = GENUSKEY,
                     family = FAMILY,
                     family_key = FAMILYKEY,
                     note_higher = NOTE)
    return(df)
  }
}
