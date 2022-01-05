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
                     dataset = NA,
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
                     dataset = NA,
                     note_higher = "Non-numeric taxon key")
    return(df)
    stop()
  }
  
  pckg <- require(rgbif, lib.loc = lib.loc)
  if(!pckg){stop("package ´rgbif´ not found")}
  
  tested_keys <- c(LCVP = "bae5856f-da10-4333-90a0-5a2135361b30", # The Leipzig catalogue of vascular plants
                   GBIF = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c", # GBIF BACKBONE
                   LIFE = "7ddf754f-d193-4cc9-b351-99906754a03b", # Catalogue of life
                   IPNI = "046bbc50-cae2-47ff-aa43-729fbf53f7c5", # International Plant Names Index
                   WCVP = "f382f0ce-323a-4091-bb9f-add557f3a9a2", # The World Checklist of Vascular Plants (WCVP)
                   TPL = "d9a4eedb-e985-4456-ad46-3df8472e00e8", # The plant list
                   ITIS = "9ca92552-f23a-41a8-a140-01abaa31c931") # Integrated Taxonomic Information System (ITIS)
  
  ORIGNAME <- NA
  GENUS <- NA
  GENUSKEY <- NA
  FAMILY <- NA
  FAMILYKEY <- NA
  NOTE <- "OK"
  DATASET <- NA
  
  e <- try({
    nl <- name_usage(key = KEY)$data
    ORIGNAME <- nl$scientificName
    if(nl$datasetKey %in% tested_keys){
      DATASET <- names(which(tested_keys == nl$datasetKey))
    } else {
      DATASET <- nl$datasetKey
    }
    
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
    } else {
      stop()
    }
  })
  
  if(class(e) == "try-error"){
    df <- data.frame(orig_key = KEY,
                     taxon_name = NA,
                     genus = NA,
                     genus_key = NA,
                     family = NA,
                     family_key = NA,
                     dataset = NA,
                     note_higher = "FAILED")
    return(df)
  } else {
    df <- data.frame(orig_key = KEY,
                     taxon_name = ORIGNAME,
                     genus = GENUS,
                     genus_key = GENUSKEY,
                     family = FAMILY,
                     family_key = FAMILYKEY,
                     dataset = DATASET,
                     note_higher = NOTE)
    return(df)
  }
}
