# orig_name <- "Sabulina rubella (Wahlenb.) Dillenb. & Kadereit"
# dataset_key <- "d9a4eedb-e985-4456-ad46-3df8472e00e8"
# resolve_taxon_name(orig_name)
# resolve_taxon_name(orig_name, dataset = "bae5856f-da10-4333-90a0-5a2135361b30")
resolve_taxon_name <- function(orig_name, dataset = NULL, lib.loc = .libPaths(), maxtry = 2){
  
  if(is.null(dataset)){
    dataset <- c(LCVP = "bae5856f-da10-4333-90a0-5a2135361b30", # The Leipzig catalogue of vascular plants
                 GBIF = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c", # GBIF BACKBONE
                 LIFE = "7ddf754f-d193-4cc9-b351-99906754a03b", # Catalogue of life
                 IPNI = "046bbc50-cae2-47ff-aa43-729fbf53f7c5", # International Plant Names Index
                 WCVP = "f382f0ce-323a-4091-bb9f-add557f3a9a2", # The World Checklist of Vascular Plants (WCVP)
                 TPL = "d9a4eedb-e985-4456-ad46-3df8472e00e8", # The plant list
                 ITIS = "9ca92552-f23a-41a8-a140-01abaa31c931") # Integrated Taxonomic Information System (ITIS)
  }
  
  tested_keys <- c(LCVP = "bae5856f-da10-4333-90a0-5a2135361b30", # The Leipzig catalogue of vascular plants
                   GBIF = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c", # GBIF BACKBONE
                   LIFE = "7ddf754f-d193-4cc9-b351-99906754a03b", # Catalogue of life
                   IPNI = "046bbc50-cae2-47ff-aa43-729fbf53f7c5", # International Plant Names Index
                   WCVP = "f382f0ce-323a-4091-bb9f-add557f3a9a2", # The World Checklist of Vascular Plants (WCVP)
                   TPL = "d9a4eedb-e985-4456-ad46-3df8472e00e8", # The plant list
                   ITIS = "9ca92552-f23a-41a8-a140-01abaa31c931") # Integrated Taxonomic Information System (ITIS)
  
  # Test if the given dataset id in the tested ones
  if(all(dataset %in% tested_keys)){
    
    
    pckg <- require(rgbif, lib.loc = lib.loc)
    if(!pckg){stop("package ´rgbif´ not found")}
    pckg <- require(taxize, lib.loc = lib.loc)
    if(!pckg){stop("package ´taxize´ not found")}
    pckg <- require(Taxonstand, lib.loc = lib.loc)
    if(!pckg){stop("package ´Taxonstand´ not found")}
    
    resolved_temp <- data.frame() # To store names
    for(dataset_key in dataset){
      # print(dataset_key)
      # Try N times if failed first
      NAME <- NULL
      NOTE <- "OK"
      for(tryn in 1:maxtry){
        if(is.null(NAME)){
          try({
            i <- orig_name
            
            try(if(Encoding(i) == "latin1"){
              ii <- i
              ii <- iconv(i, "latin1", "UTF-8")
              ii <- as.character(ii)
              if(nchar(i) > nchar(ii)){
                i <- ii
              }
            } else {
              ii <- i
              ii <- iconv(i, "UTF-8", "latin1")
              if(nchar(i) > nchar(ii)){
                i <- ii
              }
            }, silent = T)
            
            i <- gsub(" A-"," ×",i)
            i <- gsub(" x "," × ",i)
            i <- gsub(" X "," × ",i)
            i <- gsub("  "," ",i)
            i <- gsub("ë","e",i)
            i <- gsub("^x ","",i)
            i <- gsub("^X ","",i,ignore.case = T)
            i <- gsub("unknown ","",i,ignore.case = T)
            i <- gsub(" unknown","",i,ignore.case = T)
            i <- gsub("Uncertain ","",i,ignore.case = T)
            i <- gsub(" Uncertain","",i,ignore.case = T)
            i <- gsub("^ ","",i)
            i <- gsub("^ ","",i)
            i <- gsub("\\?","",i)
            i <- gsub("\\!","",i)
            i <- gsub("([0-9])","",i)
            i <- gsub("_", " ", i, perl = TRUE)
            
            parts <- strsplit(i, " ")[[1]]
            parts[1] <- gsub("[^\\p{L} ]", "", parts[1], perl = TRUE)
            parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), tolower(substr(parts[1], 2, 100)))
            if(length(parts) > 1){
              parts[2] <- tolower(parts[2])
              i <- paste(parts, collapse = " ")
              if(parts[2] %in% c("sp","sp.","spp","ssp","spp.","ssp.","spec.","species")){
                i <- parts[1]
              }
            }
            if(grepl(" = ",i)){
              parts <- strsplit(i, " = ")[[1]]
              p1 <- strsplit(parts[1], " ")[[1]]
              p2 <- strsplit(parts[2], " ")[[1]]
              if(p2[1] == "×"){
                p2 <- unlist(c(p1[1], p2))
              }
              if(nchar(p2[1]) < 3 & tolower(substr(p1[1],1,1)) == tolower(substr(p2[1],1,1))){
                p2[1] <- p1[1]
              }
              i <- paste(p2, collapse = " ")
            }
            if(grepl(" × ",i)){
              parts <- strsplit(i, " × ")[[1]]
              if(length(parts) == 2){
                p1 <- strsplit(parts[1], " ")[[1]]
                p2 <- strsplit(parts[2], " ")[[1]]
                
                if(nchar(p2[1]) < 3 & tolower(substr(p1[1],1,1)) == tolower(substr(p2[1],1,1))){
                  p2[1] <- p1[1]
                }
                i <- paste(paste(p1, collapse = " "), "×", paste(p2, collapse = " "))
              }
            }
            if(grepl(" × ",i)){
              parts <- strsplit(i, " ")[[1]]
              if(length(parts) > 2){
                parts[3] <- tolower(parts[3])
              }
              i <- paste(parts, collapse = " ")
            }
            if(grepl(" ×[[:alpha:]]+",i)){
              parts <- strsplit(i, " ")[[1]]
              if(length(parts) > 1){
                parts[2] <- tolower(parts[2])
              }
              i <- paste(parts, collapse = " ")
            }
            parts <- strsplit(i, " ")[[1]]
            parts[1] <- gsub("[^\\p{L} ]", "", parts[1], perl = TRUE)
            parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), tolower(substr(parts[1], 2, 100)))
            if(length(parts) > 1){
              parts[2] <- tolower(parts[2])
              i <- paste(parts, collapse = " ")
              if(parts[2] %in% c("sp","sp.","spp","ssp","spp.","ssp.","spec.","species")){
                i <- parts[1]
              }
            }
            
            nl <- name_suggest(q=i, datasetKey = dataset_key)$data
            
            if(nrow(nl) == 1){
              nl <- name_usage(key=nl$key)$data
              if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
              
              if(nrow(nl) == 1){
                if(!nl$synonym){
                  if(nl$taxonomicStatus == "ACCEPTED"){
                    NAME <- nl$scientificName
                    STATUS <- nl$taxonomicStatus
                    RANK <- nl$rank
                    if(nl$rank == "SPECIES"){
                      SPEC_NAME <- nl$scientificName
                    } else {
                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                        if("speciesKey" %in% names(nl)){
                          nl <- name_usage(key=nl$speciesKey)$data
                          SPEC_NAME <- nl$scientificName
                        } else {
                          SPEC_NAME <- NA
                        }
                      } else {
                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                          SPEC_NAME <- NA
                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                        } else {
                          stop()
                        }
                      }
                    }
                  } else {
                    if("acceptedKey" %in% names(nl)){
                      nl <- name_usage(key=nl$acceptedKey)$data
                      if(nl$taxonomicStatus == "ACCEPTED"){
                        NAME <- nl$scientificName
                        STATUS <- nl$taxonomicStatus
                        RANK <- nl$rank
                        if(nl$rank == "SPECIES"){
                          SPEC_NAME <- nl$scientificName
                        } else {
                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                            if("speciesKey" %in% names(nl)){
                              nl <- name_usage(key=nl$speciesKey)$data
                              SPEC_NAME <- nl$scientificName
                            } else {
                              SPEC_NAME <- NA
                            }
                          } else {
                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                              SPEC_NAME <- NA
                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
                            } else {
                              stop()
                            }
                          }
                        }
                      } else {
                        stop()
                      }
                    } else {
                      if(nl$taxonomicStatus == "DOUBTFUL"){
                        NAME <- nl$scientificName
                        STATUS <- nl$taxonomicStatus
                        RANK <- nl$rank
                        if(nl$rank == "SPECIES"){
                          SPEC_NAME <- nl$scientificName
                        } else {
                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                            if("speciesKey" %in% names(nl)){
                              nl <- name_usage(key=nl$speciesKey)$data
                              SPEC_NAME <- nl$scientificName
                            } else {
                              SPEC_NAME <- NA
                            }
                          } else {
                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                              SPEC_NAME <- NA
                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
                            } else {
                              stop()
                            }
                          }
                        }
                      } else {
                        stop()
                      }
                    }
                  }
                } else {
                  if("acceptedKey" %in% names(nl)){
                    nl <- name_usage(key=nl$acceptedKey)$data
                    if(nl$taxonomicStatus == "ACCEPTED"){
                      NAME <- nl$scientificName
                      STATUS <- nl$taxonomicStatus
                      RANK <- nl$rank
                      if(nl$rank == "SPECIES"){
                        SPEC_NAME <- nl$scientificName
                      } else {
                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                          if("speciesKey" %in% names(nl)){
                            nl <- name_usage(key=nl$speciesKey)$data
                            SPEC_NAME <- nl$scientificName
                          } else {
                            SPEC_NAME <- NA
                          }
                        } else {
                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                            SPEC_NAME <- NA
                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                          } else {
                            stop()
                          }
                        }
                      }
                    } else {
                      if(nl$taxonomicStatus == "DOUBTFUL"){
                        NAME <- nl$scientificName
                        STATUS <- nl$taxonomicStatus
                        RANK <- nl$rank
                        if(nl$rank == "SPECIES"){
                          SPEC_NAME <- nl$scientificName
                        } else {
                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                            if("speciesKey" %in% names(nl)){
                              nl <- name_usage(key=nl$speciesKey)$data
                              SPEC_NAME <- nl$scientificName
                            } else {
                              SPEC_NAME <- NA
                            }
                          } else {
                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                              SPEC_NAME <- NA
                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
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
              } else {
                if(nrow(nl) > 1){
                  
                  if(length(unique(nl$key)) == 1 & length(unique(nl$scientificName)) == 1){
                    nl <- nl[1,]
                    if(!nl$synonym){
                      if(nl$taxonomicStatus == "ACCEPTED"){
                        NAME <- nl$scientificName
                        STATUS <- nl$taxonomicStatus
                        RANK <- nl$rank
                        if(nl$rank == "SPECIES"){
                          SPEC_NAME <- nl$scientificName
                        } else {
                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                            if("speciesKey" %in% names(nl)){
                              nl <- name_usage(key=nl$speciesKey)$data
                              SPEC_NAME <- nl$scientificName
                            } else {
                              SPEC_NAME <- NA
                            }
                          } else {
                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                              SPEC_NAME <- NA
                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
                            } else {
                              stop()
                            }
                          }
                        }
                      } else {
                        if("acceptedKey" %in% names(nl)){
                          nl <- name_usage(key=nl$acceptedKey)$data
                          if(nl$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            if(nl$rank == "SPECIES"){
                              SPEC_NAME <- nl$scientificName
                            } else {
                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$speciesKey)$data
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            stop()
                          }
                        } else {
                          if(nl$taxonomicStatus == "DOUBTFUL"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            if(nl$rank == "SPECIES"){
                              SPEC_NAME <- nl$scientificName
                            } else {
                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$speciesKey)$data
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            stop()
                          }
                        }
                      }
                    } else {
                      if("acceptedKey" %in% names(nl)){
                        nl <- name_usage(key=nl$acceptedKey)$data
                        if(nl$taxonomicStatus == "ACCEPTED"){
                          NAME <- nl$scientificName
                          STATUS <- nl$taxonomicStatus
                          RANK <- nl$rank
                          if(nl$rank == "SPECIES"){
                            SPEC_NAME <- nl$scientificName
                          } else {
                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                              if("speciesKey" %in% names(nl)){
                                nl <- name_usage(key=nl$speciesKey)$data
                                SPEC_NAME <- nl$scientificName
                              } else {
                                SPEC_NAME <- NA
                              }
                            } else {
                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                SPEC_NAME <- NA
                                NOTE <- "LIKELY HIGHER LEVEL TAXON"
                              } else {
                                stop()
                              }
                            }
                          }
                        } else {
                          if("acceptedKey" %in% names(nl)){
                            nl <- name_usage(key=nl$acceptedKey)$data
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              stop()
                            }
                          } else {
                            if(nl$taxonomicStatus == "DOUBTFUL"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
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
              }
            } else {
              if(nrow(nl) > 1){
                if(length(gbif_parse(i)$canonicalname %in% nl$canonicalName) > 0){
                  if(gbif_parse(i)$canonicalname %in% nl$canonicalName){
                    nl <- nl[nl$canonicalName == gbif_parse(i)$canonicalname & (!is.na(nl$canonicalName)),]
                    if(nrow(nl) == 1){
                      nl <- name_usage(key=nl$key)$data
                      if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
                      
                      if(nrow(nl) == 1){
                        if(!nl$synonym){
                          if(nl$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            if(nl$rank == "SPECIES"){
                              SPEC_NAME <- nl$scientificName
                            } else {
                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$speciesKey)$data
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            if("acceptedKey" %in% names(nl)){
                              nl <- name_usage(key=nl$acceptedKey)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            } else {
                              if(nl$taxonomicStatus == "DOUBTFUL"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            }
                          }
                        } else {
                          if("acceptedKey" %in% names(nl)){
                            nl <- name_usage(key=nl$acceptedKey)$data
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              if(nl$taxonomicStatus == "DOUBTFUL"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
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
                      } else {
                        stop()
                      }
                    } else {
                      if(nrow(nl) > 1){
                        success <- FALSE
                        for(ikey in nl$key){
                          nl1 <- name_usage(key=ikey)$data
                          if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                            success <- TRUE
                            NAME <- nl1$scientificName
                            STATUS <- nl1$taxonomicStatus
                            RANK <- nl1$rank
                            if(nl1$rank == "SPECIES"){
                              SPEC_NAME <- nl1$scientificName
                            } else {
                              if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl1)){
                                  nl1 <- name_usage(key=nl1$speciesKey)$data
                                  SPEC_NAME <- nl1$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                            break
                          }
                        }
                        if(success == FALSE){
                          for(ikey in nl$key){
                            nl1 <- name_usage(key=ikey)$data
                            if("acceptedKey" %in% names(nl1)){
                              nl1 <- name_usage(key=nl1$acceptedKey)$data
                              
                              if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                success <- TRUE
                                NAME <- nl1$scientificName
                                STATUS <- nl1$taxonomicStatus
                                RANK <- nl1$rank
                                if(nl1$rank == "SPECIES"){
                                  SPEC_NAME <- nl1$scientificName
                                } else {
                                  if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl1)){
                                      nl1 <- name_usage(key=nl1$speciesKey)$data
                                      SPEC_NAME <- nl1$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                                break
                              }
                            }
                          }
                        }
                        if(success == FALSE){
                          NAME <- NA
                          STATUS <- NA
                          NOTE = "NOT FOUND"
                          RANK <- NA
                          SPEC_NAME <- NA
                        }
                      } else {
                        stop()
                      }
                    }
                  } else {
                    if(length(strsplit(i," ")[[1]]) == 1){
                      nl <- name_suggest(q=i, rank = "genus", datasetKey = dataset_key)$data
                      if(nrow(nl) > 0){
                        stop()
                      } else {
                        NAME <- NA
                        STATUS <- NA
                        RANK <- NA
                        NOTE = "LIKELY HIGHER LEVEL TAXON"
                        SPEC_NAME <- NA
                      }
                    } else {
                      nl <- nl[1,]
                      nl <- name_usage(key=nl$key)$data
                      if(nrow(nl) == 1){
                        if(!nl$synonym){
                          if(nl$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            if(nl$rank == "SPECIES"){
                              SPEC_NAME <- nl$scientificName
                            } else {
                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$speciesKey)$data
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            if("acceptedKey" %in% names(nl)){
                              nl <- name_usage(key=nl$acceptedKey)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            } else {
                              if(nl$taxonomicStatus == "DOUBTFUL"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            }
                          }
                        } else {
                          if("acceptedKey" %in% names(nl)){
                            nl <- name_usage(key=nl$acceptedKey)$data
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              stop()
                            }
                          } else {
                            stop()
                          }
                        }
                      } else {
                        if(nrow(nl) > 1){
                          
                          if(length(unique(nl$key)) == 1 & length(unique(nl$scientificName)) == 1){
                            nl <- nl[1,]
                            if(!nl$synonym){
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                if("acceptedKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$acceptedKey)$data
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
                                  }
                                } else {
                                  if(nl$taxonomicStatus == "DOUBTFUL"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              if("acceptedKey" %in% names(nl)){
                                nl <- name_usage(key=nl$acceptedKey)$data
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  stop()
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
                      }
                    }
                  }
                } else {
                  stop()
                }
              } else {
                ii <- gbif_parse(i)$canonicalnamecomplete
                if(!is.null(ii)){ i <- ii}
                nl <- name_suggest(q=i, datasetKey = dataset_key)$data
                if(nrow(nl) > 0 & nrow(nl) != 100){
                  if(nrow(nl) == 1){
                    nl <- name_usage(key=nl$key)$data
                    if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
                    
                    if(nrow(nl) == 1){
                      if(!nl$synonym){
                        if(nl$taxonomicStatus == "ACCEPTED"){
                          NAME <- nl$scientificName
                          STATUS <- nl$taxonomicStatus
                          RANK <- nl$rank
                          if(nl$rank == "SPECIES"){
                            SPEC_NAME <- nl$scientificName
                          } else {
                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                              if("speciesKey" %in% names(nl)){
                                nl <- name_usage(key=nl$speciesKey)$data
                                SPEC_NAME <- nl$scientificName
                              } else {
                                SPEC_NAME <- NA
                              }
                            } else {
                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                SPEC_NAME <- NA
                                NOTE <- "LIKELY HIGHER LEVEL TAXON"
                              } else {
                                stop()
                              }
                            }
                          }
                        } else {
                          if("acceptedKey" %in% names(nl)){
                            nl <- name_usage(key=nl$acceptedKey)$data
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              stop()
                            }
                          } else {
                            if(nl$taxonomicStatus == "DOUBTFUL"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              stop()
                            }
                          }
                        }
                      } else {
                        if("acceptedKey" %in% names(nl)){
                          nl <- name_usage(key=nl$acceptedKey)$data
                          if(nl$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            if(nl$rank == "SPECIES"){
                              SPEC_NAME <- nl$scientificName
                            } else {
                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$speciesKey)$data
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            if(nl$taxonomicStatus == "DOUBTFUL"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
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
                    } else {
                      stop()
                    }
                  } else {
                    if(nrow(nl) > 1){
                      success <- FALSE
                      for(ikey in nl$key){
                        nl1 <- name_usage(key=ikey)$data
                        if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                          success <- TRUE
                          NAME <- nl1$scientificName
                          STATUS <- nl1$taxonomicStatus
                          RANK <- nl1$rank
                          if(nl1$rank == "SPECIES"){
                            SPEC_NAME <- nl1$scientificName
                          } else {
                            if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                              if("speciesKey" %in% names(nl1)){
                                nl1 <- name_usage(key=nl1$speciesKey)$data
                                SPEC_NAME <- nl1$scientificName
                              } else {
                                SPEC_NAME <- NA
                              }
                            } else {
                              if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                SPEC_NAME <- NA
                                NOTE <- "LIKELY HIGHER LEVEL TAXON"
                              } else {
                                stop()
                              }
                            }
                          }
                          break
                        }
                      }
                      if(success == FALSE){
                        for(ikey in nl$key){
                          nl1 <- name_usage(key=ikey)$data
                          if("acceptedKey" %in% names(nl1)){
                            nl1 <- name_usage(key=nl1$acceptedKey)$data
                            
                            if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                              success <- TRUE
                              NAME <- nl1$scientificName
                              STATUS <- nl1$taxonomicStatus
                              RANK <- nl1$rank
                              if(nl1$rank == "SPECIES"){
                                SPEC_NAME <- nl1$scientificName
                              } else {
                                if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl1)){
                                    nl1 <- name_usage(key=nl1$speciesKey)$data
                                    SPEC_NAME <- nl1$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                              break
                            }
                          }
                        }
                      }
                      if(success == FALSE){
                        NAME <- NA
                        STATUS <- NA
                        NOTE = "NOT FOUND"
                        RANK <- NA
                        SPEC_NAME <- NA
                      }
                    } else {
                      stop()
                    }
                  }
                } else {
                  if(nrow(nl) == 100 & (grepl("×", orig_name)|grepl(" x ", orig_name)|grepl(" X ", orig_name))){
                    i <- orig_name
                    NAME <- NA
                    STATUS <- NA
                    RANK <- "HYBRID"
                    NOTE = "HYBRID TAXA! RESOLVING FAILED"
                    SPEC_NAME <- NA
                  } else {
                    ii <- gbif_parse(i)$canonicalname
                    if(!is.null(ii)){ i <- ii}
                    nl <- name_suggest(q=i, datasetKey = dataset_key)$data
                    if(nrow(nl) > 0){
                      if(nrow(nl) == 1){
                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                          nl <- name_usage(key=nl$key)$data
                          if(nl$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            SPEC_NAME <- NA
                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                          } else {
                            if("acceptedKey" %in% names(nl)){
                              nl <- name_usage(key=nl$acceptedKey)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      stop()
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            } else {
                              NAME <- NA
                              STATUS <- NA
                              SPEC_NAME <- NA
                              RANK <- NA
                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
                            }
                          }
                        } else {
                          nl <- name_usage(key=nl$key)$data
                          if(nrow(nl) == 1){
                            if(!nl$synonym){
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                if("acceptedKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$acceptedKey)$data
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
                                  }
                                } else {
                                  if(nl$taxonomicStatus == "DOUBTFUL"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              if("acceptedKey" %in% names(nl)){
                                nl <- name_usage(key=nl$acceptedKey)$data
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  if(nl$taxonomicStatus == "DOUBTFUL"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
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
                          } else {
                            if(nrow(nl) > 1){
                              
                              if(length(unique(nl$key)) == 1 & length(unique(nl$scientificName)) == 1){
                                nl <- nl[1,]
                                if(!nl$synonym){
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    if("acceptedKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$acceptedKey)$data
                                      if(nl$taxonomicStatus == "ACCEPTED"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        stop()
                                      }
                                    } else {
                                      if(nl$taxonomicStatus == "DOUBTFUL"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  if("acceptedKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$acceptedKey)$data
                                    if(nl$taxonomicStatus == "ACCEPTED"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      RANK <- nl$rank
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      stop()
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
                          }
                        }
                      } else {
                        success <- FALSE
                        for(ikey in nl$key){
                          nl1 <- name_usage(key=ikey)$data
                          if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                            success <- TRUE
                            NAME <- nl1$scientificName
                            STATUS <- nl1$taxonomicStatus
                            RANK <- nl1$rank
                            if(nl1$rank == "SPECIES"){
                              SPEC_NAME <- nl1$scientificName
                            } else {
                              if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl1)){
                                  nl1 <- name_usage(key=nl1$speciesKey)$data
                                  SPEC_NAME <- nl1$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                            break
                          }
                        }
                        if(success == FALSE){
                          for(ikey in nl$key){
                            nl1 <- name_usage(key=ikey)$data
                            if("acceptedKey" %in% names(nl1)){
                              nl1 <- name_usage(key=nl1$acceptedKey)$data
                              
                              if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                success <- TRUE
                                NAME <- nl1$scientificName
                                STATUS <- nl1$taxonomicStatus
                                RANK <- nl1$rank
                                if(nl1$rank == "SPECIES"){
                                  SPEC_NAME <- nl1$scientificName
                                } else {
                                  if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl1)){
                                      nl1 <- name_usage(key=nl1$speciesKey)$data
                                      SPEC_NAME <- nl1$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                                break
                              }
                            }
                          }
                        }
                        if(success == FALSE){
                          NAME <- NA
                          STATUS <- NA
                          RANK <- NA
                          NOTE = "NOT FOUND"
                          SPEC_NAME <- NA
                        }
                      }
                    } else {
                      ii <- gbif_parse(i)$canonicalname
                      if(!is.null(ii)){ i <- ii}
                      if(length(strsplit(i, " ")[[1]]) > 2){
                        i <- paste(strsplit(i, " ")[[1]][1:2], collapse = " ")
                        nl <- name_suggest(q=i, datasetKey = dataset_key)$data
                        if(nrow(nl) > 0){
                          if(nrow(nl) == 1){
                            nl <- name_usage(key=nl$key)$data
                            if(!nl$synonym){
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                NOTE = "PARSED TO SPECIES"
                                RANK <- nl$rank
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                if("acceptedKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$acceptedKey)$data
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    NOTE = "PARSED TO SPECIES"
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
                                  }
                                } else {
                                  if(nl$taxonomicStatus == "DOUBTFUL"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    NOTE = "PARSED TO SPECIES"
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              if("acceptedKey" %in% names(nl)){
                                nl <- name_usage(key=nl$acceptedKey)$data
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  NOTE = "PARSED TO SPECIES"
                                  RANK <- nl$rank
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  if("acceptedKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$acceptedKey)$data
                                    if(nl$taxonomicStatus == "ACCEPTED"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      RANK <- nl$rank
                                      NOTE = "PARSED TO SPECIES"
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      stop()
                                    }
                                  } else {
                                    if(nl$taxonomicStatus == "DOUBTFUL"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      RANK <- nl$rank
                                      NOTE = "PARSED TO SPECIES"
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
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
                            if(gbif_parse(i)$canonicalname %in% nl$canonicalName){
                              nl <- nl[nl$canonicalName == gbif_parse(i)$canonicalname & (!is.na(nl$canonicalName)),]
                              if(nrow(nl) == 1){
                                nl <- name_usage(key=nl$key)$data
                                if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
                                
                                if(nrow(nl) == 1){
                                  if(!nl$synonym){
                                    if(nl$taxonomicStatus == "ACCEPTED"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      NOTE = "PARSED TO SPECIES"
                                      RANK <- nl$rank
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      if("acceptedKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$acceptedKey)$data
                                        if(nl$taxonomicStatus == "ACCEPTED"){
                                          NAME <- nl$scientificName
                                          STATUS <- nl$taxonomicStatus
                                          RANK <- nl$rank
                                          NOTE = "PARSED TO SPECIES"
                                          if(nl$rank == "SPECIES"){
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                              if("speciesKey" %in% names(nl)){
                                                nl <- name_usage(key=nl$speciesKey)$data
                                                SPEC_NAME <- nl$scientificName
                                              } else {
                                                SPEC_NAME <- NA
                                              }
                                            } else {
                                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                SPEC_NAME <- NA
                                                NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                              } else {
                                                stop()
                                              }
                                            }
                                          }
                                        } else {
                                          stop()
                                        }
                                      } else {
                                        if(nl$taxonomicStatus == "DOUBTFUL"){
                                          NAME <- nl$scientificName
                                          STATUS <- nl$taxonomicStatus
                                          RANK <- nl$rank
                                          NOTE = "PARSED TO SPECIES"
                                          if(nl$rank == "SPECIES"){
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                              if("speciesKey" %in% names(nl)){
                                                nl <- name_usage(key=nl$speciesKey)$data
                                                SPEC_NAME <- nl$scientificName
                                              } else {
                                                SPEC_NAME <- NA
                                              }
                                            } else {
                                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                SPEC_NAME <- NA
                                                NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                              } else {
                                                stop()
                                              }
                                            }
                                          }
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    if("acceptedKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$acceptedKey)$data
                                      if(nl$taxonomicStatus == "ACCEPTED"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        NOTE = "PARSED TO SPECIES"
                                        RANK <- nl$rank
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        stop()
                                      }
                                    } else {
                                      stop()
                                    }
                                  }
                                } else {
                                  stop()
                                }
                              } else {
                                if(nrow(nl) > 1){
                                  success <- FALSE
                                  for(ikey in nl$key){
                                    nl1 <- name_usage(key=ikey)$data
                                    if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                      success <- TRUE
                                      NAME <- nl1$scientificName
                                      STATUS <- nl1$taxonomicStatus
                                      NOTE = "PARSED TO SPECIES"
                                      RANK <- nl1$rank
                                      if(nl1$rank == "SPECIES"){
                                        SPEC_NAME <- nl1$scientificName
                                      } else {
                                        if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl1)){
                                            nl1 <- name_usage(key=nl1$speciesKey)$data
                                            SPEC_NAME <- nl1$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                      break
                                    }
                                  }
                                  if(success == FALSE){
                                    for(ikey in nl$key){
                                      nl1 <- name_usage(key=ikey)$data
                                      if("acceptedKey" %in% names(nl1)){
                                        nl1 <- name_usage(key=nl1$acceptedKey)$data
                                        
                                        if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                          success <- TRUE
                                          NAME <- nl1$scientificName
                                          STATUS <- nl1$taxonomicStatus
                                          NOTE = "PARSED TO SPECIES"
                                          RANK <- nl1$rank
                                          if(nl1$rank == "SPECIES"){
                                            SPEC_NAME <- nl1$scientificName
                                          } else {
                                            if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                              if("speciesKey" %in% names(nl1)){
                                                nl1 <- name_usage(key=nl1$speciesKey)$data
                                                SPEC_NAME <- nl1$scientificName
                                              } else {
                                                SPEC_NAME <- NA
                                              }
                                            } else {
                                              if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                SPEC_NAME <- NA
                                                NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                              } else {
                                                stop()
                                              }
                                            }
                                          }
                                          break
                                        }
                                      }
                                    }
                                  }
                                  if(success == FALSE){
                                    NAME <- NA
                                    STATUS <- NA
                                    RANK <- NA
                                    NOTE = "NOT FOUND"
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  stop()
                                }
                              }
                            } else {
                              success <- FALSE
                              for(ikey in nl$key){
                                nl1 <- name_usage(key=ikey)$data
                                if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                  success <- TRUE
                                  NAME <- nl1$scientificName
                                  STATUS <- nl1$taxonomicStatus
                                  RANK <- nl1$rank
                                  NOTE = "PARSED TO SPECIES"
                                  if(nl1$rank == "SPECIES"){
                                    SPEC_NAME <- nl1$scientificName
                                  } else {
                                    if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl1)){
                                        nl1 <- name_usage(key=nl1$speciesKey)$data
                                        SPEC_NAME <- nl1$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                  break
                                }
                              }
                              if(success == FALSE){
                                for(ikey in nl$key){
                                  nl1 <- name_usage(key=ikey)$data
                                  if("acceptedKey" %in% names(nl1)){
                                    nl1 <- name_usage(key=nl1$acceptedKey)$data
                                    
                                    if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                      success <- TRUE
                                      NAME <- nl1$scientificName
                                      STATUS <- nl1$taxonomicStatus
                                      RANK <- nl1$rank
                                      NOTE = "PARSED TO SPECIES"
                                      if(nl1$rank == "SPECIES"){
                                        SPEC_NAME <- nl1$scientificName
                                      } else {
                                        if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl1)){
                                            nl1 <- name_usage(key=nl1$speciesKey)$data
                                            SPEC_NAME <- nl1$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                      break
                                    }
                                  }
                                }
                              }
                              if(success == FALSE){
                                NAME <- NA
                                STATUS <- NA
                                RANK <- NA
                                NOTE = "NOT FOUND"
                                SPEC_NAME <- NA
                              }
                            }
                          }
                        }
                      } else {
                        if(length(strsplit(i, " ")[[1]]) == 2){
                          if(strsplit(i, " ")[[1]][2] == "spec."){
                            NAME <- NA
                            STATUS <- NA
                            SPEC_NAME <- NA
                            RANK <- NA
                            NOTE <- "LIKELY HIGHER LEVEL TAXON"
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            if(is.null(NAME) | NOTE == "PARSED TO SPECIES"){
              i <- orig_name
              
              try(if(Encoding(i) == "latin1"){
                ii <- i
                ii <- iconv(i, "latin1", "UTF-8")
                ii <- as.character(ii)
                if(nchar(i) > nchar(ii)){
                  i <- ii
                }
              } else {
                ii <- i
                ii <- iconv(i, "UTF-8", "latin1")
                if(nchar(i) > nchar(ii)){
                  i <- ii
                }
              }, silent = T)
              
              i <- gsub(" A-"," ×",i)
              i <- gsub(" x "," × ",i)
              i <- gsub(" X "," × ",i)
              i <- gsub("  "," ",i)
              i <- gsub("ë","e",i)
              i <- gsub("^x ","",i)
              i <- gsub("^X ","",i,ignore.case = T)
              i <- gsub("unknown ","",i,ignore.case = T)
              i <- gsub(" unknown","",i,ignore.case = T)
              i <- gsub("Uncertain ","",i,ignore.case = T)
              i <- gsub(" Uncertain","",i,ignore.case = T)
              i <- gsub("^ ","",i)
              i <- gsub("^ ","",i)
              i <- gsub("\\?","",i)
              i <- gsub("\\!","",i)
              i <- gsub("([0-9])","",i)
              i <- gsub("_", " ", i, perl = TRUE)
              
              parts <- strsplit(i, " ")[[1]]
              parts[1] <- gsub("[^\\p{L} ]", "", parts[1], perl = TRUE)
              parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), tolower(substr(parts[1], 2, 100)))
              if(length(parts) > 1){
                parts[2] <- tolower(parts[2])
                i <- paste(parts, collapse = " ")
                if(parts[2] %in% c("sp","sp.","spp","ssp","spp.","ssp.","spec.","species")){
                  i <- parts[1]
                }
              }
              if(grepl(" = ",i)){
                parts <- strsplit(i, " = ")[[1]]
                p1 <- strsplit(parts[1], " ")[[1]]
                p2 <- strsplit(parts[2], " ")[[1]]
                if(p2[1] == "×"){
                  p2 <- unlist(c(p1[1], p2))
                }
                if(nchar(p2[1]) < 3 & tolower(substr(p1[1],1,1)) == tolower(substr(p2[1],1,1))){
                  p2[1] <- p1[1]
                }
                i <- paste(p2, collapse = " ")
              }
              if(grepl(" × ",i)){
                parts <- strsplit(i, " × ")[[1]]
                if(length(parts) == 2){
                  p1 <- strsplit(parts[1], " ")[[1]]
                  p2 <- strsplit(parts[2], " ")[[1]]
                  
                  if(nchar(p2[1]) < 3 & tolower(substr(p1[1],1,1)) == tolower(substr(p2[1],1,1))){
                    p2[1] <- p1[1]
                  }
                  i <- paste(paste(p1, collapse = " "), "×", paste(p2, collapse = " "))
                }
              }
              if(grepl(" × ",i)){
                parts <- strsplit(i, " ")[[1]]
                if(length(parts) > 2){
                  parts[3] <- tolower(parts[3])
                }
                i <- paste(parts, collapse = " ")
              }
              if(grepl(" ×[[:alpha:]]+",i)){
                parts <- strsplit(i, " ")[[1]]
                if(length(parts) > 1){
                  parts[2] <- tolower(parts[2])
                }
                i <- paste(parts, collapse = " ")
              }
              parts <- strsplit(i, " ")[[1]]
              parts[1] <- gsub("[^\\p{L} ]", "", parts[1], perl = TRUE)
              parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), tolower(substr(parts[1], 2, 100)))
              if(length(parts) > 1){
                parts[2] <- tolower(parts[2])
                i <- paste(parts, collapse = " ")
                if(parts[2] %in% c("sp","sp.","spp","ssp","spp.","ssp.","spec.","species")){
                  i <- parts[1]
                }
              }
              
              tp <- TPL(i, diffchar = 1, max.distance = 1)
              if(tp$Typo & tp$Plant.Name.Index){
                i <- gbif_parse(paste0(c(tp$New.Genus, tp$New.Hybrid.marker, tp$New.Species, tp$New.Infraspecific.rank, tp$New.Infraspecific, tp$New.Authority), collapse = " "))$canonicalnamecomplete
                
                nl <- name_suggest(q=i, datasetKey = dataset_key)$data
                
                if(nrow(nl) == 1){
                  nl <- name_usage(key=nl$key)$data
                  if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
                  
                  if(nrow(nl) == 1){
                    if(!nl$synonym){
                      if(nl$taxonomicStatus == "ACCEPTED"){
                        NAME <- nl$scientificName
                        STATUS <- nl$taxonomicStatus
                        RANK <- nl$rank
                        NOTE <- "FUZZY MATCHING"
                        if(nl$rank == "SPECIES"){
                          SPEC_NAME <- nl$scientificName
                        } else {
                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                            if("speciesKey" %in% names(nl)){
                              nl <- name_usage(key=nl$speciesKey)$data
                              SPEC_NAME <- nl$scientificName
                            } else {
                              SPEC_NAME <- NA
                            }
                          } else {
                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                              SPEC_NAME <- NA
                              NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                            } else {
                              stop()
                            }
                          }
                        }
                      } else {
                        if(nl$taxonomicStatus == "DOUBTFUL"){
                          NAME <- nl$scientificName
                          STATUS <- nl$taxonomicStatus
                          RANK <- nl$rank
                          NOTE <- "FUZZY MATCHING"
                          if(nl$rank == "SPECIES"){
                            SPEC_NAME <- nl$scientificName
                          } else {
                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                              if("speciesKey" %in% names(nl)){
                                nl <- name_usage(key=nl$speciesKey)$data
                                SPEC_NAME <- nl$scientificName
                              } else {
                                SPEC_NAME <- NA
                              }
                            } else {
                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                SPEC_NAME <- NA
                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
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
                      if("acceptedKey" %in% names(nl)){
                        nl <- name_usage(key=nl$acceptedKey)$data
                        if(nl$taxonomicStatus == "ACCEPTED"){
                          NAME <- nl$scientificName
                          STATUS <- nl$taxonomicStatus
                          RANK <- nl$rank
                          NOTE <- "FUZZY MATCHING"
                          if(nl$rank == "SPECIES"){
                            SPEC_NAME <- nl$scientificName
                          } else {
                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                              if("speciesKey" %in% names(nl)){
                                nl <- name_usage(key=nl$speciesKey)$data
                                SPEC_NAME <- nl$scientificName
                              } else {
                                SPEC_NAME <- NA
                              }
                            } else {
                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                SPEC_NAME <- NA
                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                              } else {
                                stop()
                              }
                            }
                          }
                        } else {
                          if(nl$taxonomicStatus == "DOUBTFUL"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            NOTE <- "FUZZY MATCHING"
                            if(nl$rank == "SPECIES"){
                              SPEC_NAME <- nl$scientificName
                            } else {
                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$speciesKey)$data
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
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
                  } else {
                    if(nrow(nl) > 1){
                      
                      if(length(unique(nl$key)) == 1 & length(unique(nl$scientificName)) == 1){
                        nl <- nl[1,]
                        if(!nl$synonym){
                          if(nl$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                            RANK <- nl$rank
                            NOTE <- "FUZZY MATCHING"
                            if(nl$rank == "SPECIES"){
                              SPEC_NAME <- nl$scientificName
                            } else {
                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                if("speciesKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$speciesKey)$data
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  SPEC_NAME <- NA
                                }
                              } else {
                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                  SPEC_NAME <- NA
                                  NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            if("acceptedKey" %in% names(nl)){
                              nl <- name_usage(key=nl$acceptedKey)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                NOTE <- "FUZZY MATCHING"
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            } else {
                              if(nl$taxonomicStatus == "DOUBTFUL"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                NOTE <- "FUZZY MATCHING"
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            }
                          }
                        } else {
                          if("acceptedKey" %in% names(nl)){
                            nl <- name_usage(key=nl$acceptedKey)$data
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              NOTE <- "FUZZY MATCHING"
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    stop()
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              stop()
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
                  }
                } else {
                  if(nrow(nl) > 1){
                    if(gbif_parse(i)$canonicalname %in% nl$canonicalName){
                      nl <- nl[nl$canonicalName == gbif_parse(i)$canonicalname & (!is.na(nl$canonicalName)),]
                      if(nrow(nl) == 1){
                        nl <- name_usage(key=nl$key)$data
                        if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
                        
                        if(nrow(nl) == 1){
                          if(!nl$synonym){
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              NOTE <- "FUZZY MATCHING"
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              if("acceptedKey" %in% names(nl)){
                                nl <- name_usage(key=nl$acceptedKey)$data
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  stop()
                                }
                              } else {
                                if(nl$taxonomicStatus == "DOUBTFUL"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            if("acceptedKey" %in% names(nl)){
                              nl <- name_usage(key=nl$acceptedKey)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                NOTE <- "FUZZY MATCHING"
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            } else {
                              stop()
                            }
                          }
                        } else {
                          stop()
                        }
                      } else {
                        if(nrow(nl) > 1){
                          success <- FALSE
                          for(ikey in nl$key){
                            nl1 <- name_usage(key=ikey)$data
                            if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                              success <- TRUE
                              NAME <- nl1$scientificName
                              STATUS <- nl1$taxonomicStatus
                              RANK <- nl1$rank
                              NOTE <- "FUZZY MATCHING"
                              if(nl1$rank == "SPECIES"){
                                SPEC_NAME <- nl1$scientificName
                              } else {
                                if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl1)){
                                    nl1 <- name_usage(key=nl1$speciesKey)$data
                                    SPEC_NAME <- nl1$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                              break
                            }
                          }
                          if(success == FALSE){
                            for(ikey in nl$key){
                              nl1 <- name_usage(key=ikey)$data
                              if("acceptedKey" %in% names(nl1)){
                                nl1 <- name_usage(key=nl1$acceptedKey)$data
                                
                                if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                  success <- TRUE
                                  NAME <- nl1$scientificName
                                  STATUS <- nl1$taxonomicStatus
                                  RANK <- nl1$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl1$rank == "SPECIES"){
                                    SPEC_NAME <- nl1$scientificName
                                  } else {
                                    if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl1)){
                                        nl1 <- name_usage(key=nl1$speciesKey)$data
                                        SPEC_NAME <- nl1$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                  break
                                }
                              }
                            }
                          }
                          if(success == FALSE){
                            NAME <- NA
                            STATUS <- NA
                            NOTE = "NOT FOUND"
                            RANK <- NA
                            SPEC_NAME <- NA
                          }
                        } else {
                          stop()
                        }
                      }
                    } else {
                      if(length(strsplit(i, " ")[[1]]) == 1){
                        nl <- name_suggest(q=i, rank = "genus", datasetKey = dataset_key)$data
                        if(nrow(nl) > 0){
                          stop()
                        } else {
                          NAME <- NA
                          STATUS <- NA
                          RANK <- NA
                          NOTE = "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                          SPEC_NAME <- NA
                        }
                      } else {
                        nl <- nl[1,]
                        nl <- name_usage(key=nl$key)$data
                        if(nrow(nl) == 1){
                          if(!nl$synonym){
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              NOTE <- "FUZZY MATCHING"
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              if("acceptedKey" %in% names(nl)){
                                nl <- name_usage(key=nl$acceptedKey)$data
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  stop()
                                }
                              } else {
                                if(nl$taxonomicStatus == "DOUBTFUL"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            if("acceptedKey" %in% names(nl)){
                              nl <- name_usage(key=nl$acceptedKey)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                NOTE <- "FUZZY MATCHING"
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            } else {
                              stop()
                            }
                          }
                        } else {
                          if(nrow(nl) > 1){
                            
                            if(length(unique(nl$key)) == 1 & length(unique(nl$scientificName)) == 1){
                              nl <- nl[1,]
                              if(!nl$synonym){
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  if("acceptedKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$acceptedKey)$data
                                    if(nl$taxonomicStatus == "ACCEPTED"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      RANK <- nl$rank
                                      NOTE <- "FUZZY MATCHING"
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      stop()
                                    }
                                  } else {
                                    if(nl$taxonomicStatus == "DOUBTFUL"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      RANK <- nl$rank
                                      NOTE <- "FUZZY MATCHING"
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                if("acceptedKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$acceptedKey)$data
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    NOTE <- "FUZZY MATCHING"
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
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
                        }
                      }
                    }
                  } else {
                    ii <- gbif_parse(i)$canonicalnamecomplete
                    if(!is.null(ii)){ i <- ii}
                    nl <- name_suggest(q=i, datasetKey = dataset_key)$data
                    if(nrow(nl) > 0 & nrow(nl) != 100){
                      if(nrow(nl) == 1){
                        nl <- name_usage(key=nl$key)$data
                        if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
                        
                        if(nrow(nl) == 1){
                          if(!nl$synonym){
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                              RANK <- nl$rank
                              NOTE <- "FUZZY MATCHING"
                              if(nl$rank == "SPECIES"){
                                SPEC_NAME <- nl$scientificName
                              } else {
                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$speciesKey)$data
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                            } else {
                              if("acceptedKey" %in% names(nl)){
                                nl <- name_usage(key=nl$acceptedKey)$data
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  stop()
                                }
                              } else {
                                if(nl$taxonomicStatus == "DOUBTFUL"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                  RANK <- nl$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl$rank == "SPECIES"){
                                    SPEC_NAME <- nl$scientificName
                                  } else {
                                    if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$speciesKey)$data
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  stop()
                                }
                              }
                            }
                          } else {
                            if("acceptedKey" %in% names(nl)){
                              nl <- name_usage(key=nl$acceptedKey)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                NOTE <- "FUZZY MATCHING"
                                if(nl$rank == "SPECIES"){
                                  SPEC_NAME <- nl$scientificName
                                } else {
                                  if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                    if("speciesKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$speciesKey)$data
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      SPEC_NAME <- NA
                                    }
                                  } else {
                                    if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                      SPEC_NAME <- NA
                                      NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                    } else {
                                      stop()
                                    }
                                  }
                                }
                              } else {
                                stop()
                              }
                            } else {
                              stop()
                            }
                          }
                        } else {
                          stop()
                        }
                      } else {
                        if(nrow(nl) > 1){
                          success <- FALSE
                          for(ikey in nl$key){
                            nl1 <- name_usage(key=ikey)$data
                            if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                              success <- TRUE
                              NAME <- nl1$scientificName
                              STATUS <- nl1$taxonomicStatus
                              RANK <- nl1$rank
                              NOTE <- "FUZZY MATCHING"
                              if(nl1$rank == "SPECIES"){
                                SPEC_NAME <- nl1$scientificName
                              } else {
                                if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                  if("speciesKey" %in% names(nl1)){
                                    nl1 <- name_usage(key=nl1$speciesKey)$data
                                    SPEC_NAME <- nl1$scientificName
                                  } else {
                                    SPEC_NAME <- NA
                                  }
                                } else {
                                  if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                    SPEC_NAME <- NA
                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                  } else {
                                    stop()
                                  }
                                }
                              }
                              break
                            }
                          }
                          if(success == FALSE){
                            for(ikey in nl$key){
                              nl1 <- name_usage(key=ikey)$data
                              if("acceptedKey" %in% names(nl1)){
                                nl1 <- name_usage(key=nl1$acceptedKey)$data
                                
                                if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                  success <- TRUE
                                  NAME <- nl1$scientificName
                                  STATUS <- nl1$taxonomicStatus
                                  RANK <- nl1$rank
                                  NOTE <- "FUZZY MATCHING"
                                  if(nl1$rank == "SPECIES"){
                                    SPEC_NAME <- nl1$scientificName
                                  } else {
                                    if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                      if("speciesKey" %in% names(nl1)){
                                        nl1 <- name_usage(key=nl1$speciesKey)$data
                                        SPEC_NAME <- nl1$scientificName
                                      } else {
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                        SPEC_NAME <- NA
                                        NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                  break
                                }
                              }
                            }
                          }
                          if(success == FALSE){
                            NAME <- NA
                            STATUS <- NA
                            NOTE = "NOT FOUND"
                            RANK <- NA
                            SPEC_NAME <- NA
                          }
                        } else {
                          stop()
                        }
                      }
                    } else {
                      if(nrow(nl) == 100 & grepl("×", orig_name)){
                        i <- orig_name
                        NAME <- NA
                        STATUS <- NA
                        RANK <- "HYBRID"
                        NOTE = "FUZZY MATCHING! HYBRID TAXA! RESOLVING FAILED"
                        SPEC_NAME <- NA
                      } else {
                        ii <- gbif_parse(i)$canonicalname
                        if(!is.null(ii)){ i <- ii}
                        nl <- name_suggest(q=i, datasetKey = dataset_key)$data
                        if(nrow(nl) > 0){
                          if(nrow(nl) == 1){
                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                              nl <- name_usage(key=nl$key)$data
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                                RANK <- nl$rank
                                SPEC_NAME <- NA
                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                              } else {
                                if("acceptedKey" %in% names(nl)){
                                  nl <- name_usage(key=nl$acceptedKey)$data
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    NOTE <- "FUZZY MATCHING"
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    stop()
                                  }
                                } else {
                                  NAME <- NA
                                  STATUS <- NA
                                  SPEC_NAME <- NA
                                  RANK <- NA
                                  NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                }
                              }
                            } else {
                              nl <- name_usage(key=nl$key)$data
                              if(nrow(nl) == 1){
                                if(!nl$synonym){
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    RANK <- nl$rank
                                    NOTE <- "FUZZY MATCHING"
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    if("acceptedKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$acceptedKey)$data
                                      if(nl$taxonomicStatus == "ACCEPTED"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        NOTE <- "FUZZY MATCHING"
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        stop()
                                      }
                                    } else {
                                      if(nl$taxonomicStatus == "DOUBTFUL"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        NOTE <- "FUZZY MATCHING"
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  if("acceptedKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$acceptedKey)$data
                                    if(nl$taxonomicStatus == "ACCEPTED"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      RANK <- nl$rank
                                      NOTE <- "FUZZY MATCHING"
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      if(nl$taxonomicStatus == "DOUBTFUL"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        NOTE <- "FUZZY MATCHING"
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
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
                              } else {
                                if(nrow(nl) > 1){
                                  
                                  if(length(unique(nl$key)) == 1 & length(unique(nl$scientificName)) == 1){
                                    nl <- nl[1,]
                                    if(!nl$synonym){
                                      if(nl$taxonomicStatus == "ACCEPTED"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        NOTE <- "FUZZY MATCHING"
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        if("acceptedKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$acceptedKey)$data
                                          if(nl$taxonomicStatus == "ACCEPTED"){
                                            NAME <- nl$scientificName
                                            STATUS <- nl$taxonomicStatus
                                            RANK <- nl$rank
                                            NOTE <- "FUZZY MATCHING"
                                            if(nl$rank == "SPECIES"){
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                                if("speciesKey" %in% names(nl)){
                                                  nl <- name_usage(key=nl$speciesKey)$data
                                                  SPEC_NAME <- nl$scientificName
                                                } else {
                                                  SPEC_NAME <- NA
                                                }
                                              } else {
                                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                  SPEC_NAME <- NA
                                                  NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                                } else {
                                                  stop()
                                                }
                                              }
                                            }
                                          } else {
                                            stop()
                                          }
                                        } else {
                                          if(nl$taxonomicStatus == "DOUBTFUL"){
                                            NAME <- nl$scientificName
                                            STATUS <- nl$taxonomicStatus
                                            RANK <- nl$rank
                                            NOTE <- "FUZZY MATCHING"
                                            if(nl$rank == "SPECIES"){
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                                if("speciesKey" %in% names(nl)){
                                                  nl <- name_usage(key=nl$speciesKey)$data
                                                  SPEC_NAME <- nl$scientificName
                                                } else {
                                                  SPEC_NAME <- NA
                                                }
                                              } else {
                                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                  SPEC_NAME <- NA
                                                  NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                                } else {
                                                  stop()
                                                }
                                              }
                                            }
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      if("acceptedKey" %in% names(nl)){
                                        nl <- name_usage(key=nl$acceptedKey)$data
                                        if(nl$taxonomicStatus == "ACCEPTED"){
                                          NAME <- nl$scientificName
                                          STATUS <- nl$taxonomicStatus
                                          RANK <- nl$rank
                                          NOTE <- "FUZZY MATCHING"
                                          if(nl$rank == "SPECIES"){
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                              if("speciesKey" %in% names(nl)){
                                                nl <- name_usage(key=nl$speciesKey)$data
                                                SPEC_NAME <- nl$scientificName
                                              } else {
                                                SPEC_NAME <- NA
                                              }
                                            } else {
                                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                SPEC_NAME <- NA
                                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                              } else {
                                                stop()
                                              }
                                            }
                                          }
                                        } else {
                                          stop()
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
                              }
                            }
                          } else {
                            success <- FALSE
                            for(ikey in nl$key){
                              nl1 <- name_usage(key=ikey)$data
                              if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                success <- TRUE
                                NAME <- nl1$scientificName
                                STATUS <- nl1$taxonomicStatus
                                SPEC_NAME <- NA
                                RANK <- nl1$rank
                                NOTE = "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                break
                              }
                            }
                            if(success == FALSE){
                              for(ikey in nl$key){
                                nl1 <- name_usage(key=ikey)$data
                                if("acceptedKey" %in% names(nl1)){
                                  nl1 <- name_usage(key=nl1$acceptedKey)$data
                                  
                                  if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                    success <- TRUE
                                    NAME <- nl1$scientificName
                                    STATUS <- nl1$taxonomicStatus
                                    SPEC_NAME <- NA
                                    RANK <- nl1$rank
                                    NOTE = "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                    break
                                  }
                                }
                              }
                            }
                            if(success == FALSE){
                              NAME <- NA
                              STATUS <- NA
                              RANK <- NA
                              NOTE = "NOT FOUND"
                              SPEC_NAME <- NA
                            }
                          }
                        } else {
                          ii <- gbif_parse(i)$canonicalname
                          if(!is.null(ii)){ i <- ii}
                          if(length(strsplit(i, " ")[[1]]) > 2){
                            i <- paste(strsplit(i, " ")[[1]][1:2], collapse = " ")
                            nl <- name_suggest(q=i, datasetKey = dataset_key)$data
                            if(nrow(nl) > 0){
                              if(nrow(nl) == 1){
                                nl <- name_usage(key=nl$key)$data
                                if(!nl$synonym){
                                  if(nl$taxonomicStatus == "ACCEPTED"){
                                    NAME <- nl$scientificName
                                    STATUS <- nl$taxonomicStatus
                                    NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                    RANK <- nl$rank
                                    if(nl$rank == "SPECIES"){
                                      SPEC_NAME <- nl$scientificName
                                    } else {
                                      if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                        if("speciesKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$speciesKey)$data
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          SPEC_NAME <- NA
                                        }
                                      } else {
                                        if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                          SPEC_NAME <- NA
                                          NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                        } else {
                                          stop()
                                        }
                                      }
                                    }
                                  } else {
                                    if("acceptedKey" %in% names(nl)){
                                      nl <- name_usage(key=nl$acceptedKey)$data
                                      if(nl$taxonomicStatus == "ACCEPTED"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        stop()
                                      }
                                    } else {
                                      if(nl$taxonomicStatus == "DOUBTFUL"){
                                        NAME <- nl$scientificName
                                        STATUS <- nl$taxonomicStatus
                                        RANK <- nl$rank
                                        NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                        if(nl$rank == "SPECIES"){
                                          SPEC_NAME <- nl$scientificName
                                        } else {
                                          if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                            if("speciesKey" %in% names(nl)){
                                              nl <- name_usage(key=nl$speciesKey)$data
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              SPEC_NAME <- NA
                                            }
                                          } else {
                                            if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                              SPEC_NAME <- NA
                                              NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        stop()
                                      }
                                    }
                                  }
                                } else {
                                  if("acceptedKey" %in% names(nl)){
                                    nl <- name_usage(key=nl$acceptedKey)$data
                                    if(nl$taxonomicStatus == "ACCEPTED"){
                                      NAME <- nl$scientificName
                                      STATUS <- nl$taxonomicStatus
                                      NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                      RANK <- nl$rank
                                      if(nl$rank == "SPECIES"){
                                        SPEC_NAME <- nl$scientificName
                                      } else {
                                        if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$speciesKey)$data
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                    } else {
                                      stop()
                                    }
                                  } else {
                                    stop()
                                  }
                                }
                              } else {
                                if(gbif_parse(i)$canonicalname %in% nl$canonicalName){
                                  nl <- nl[nl$canonicalName == gbif_parse(i)$canonicalname & (!is.na(nl$canonicalName)),]
                                  if(nrow(nl) == 1){
                                    nl <- name_usage(key=nl$key)$data
                                    if(nrow(nl) > 1){nl <- nl[!duplicated(nl),]}
                                    
                                    if(nrow(nl) == 1){
                                      if(!nl$synonym){
                                        if(nl$taxonomicStatus == "ACCEPTED"){
                                          NAME <- nl$scientificName
                                          STATUS <- nl$taxonomicStatus
                                          NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                          RANK <- nl$rank
                                          if(nl$rank == "SPECIES"){
                                            SPEC_NAME <- nl$scientificName
                                          } else {
                                            if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                              if("speciesKey" %in% names(nl)){
                                                nl <- name_usage(key=nl$speciesKey)$data
                                                SPEC_NAME <- nl$scientificName
                                              } else {
                                                SPEC_NAME <- NA
                                              }
                                            } else {
                                              if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                SPEC_NAME <- NA
                                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                              } else {
                                                stop()
                                              }
                                            }
                                          }
                                        } else {
                                          if("acceptedKey" %in% names(nl)){
                                            nl <- name_usage(key=nl$acceptedKey)$data
                                            if(nl$taxonomicStatus == "ACCEPTED"){
                                              NAME <- nl$scientificName
                                              STATUS <- nl$taxonomicStatus
                                              RANK <- nl$rank
                                              NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                              if(nl$rank == "SPECIES"){
                                                SPEC_NAME <- nl$scientificName
                                              } else {
                                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                                  if("speciesKey" %in% names(nl)){
                                                    nl <- name_usage(key=nl$speciesKey)$data
                                                    SPEC_NAME <- nl$scientificName
                                                  } else {
                                                    SPEC_NAME <- NA
                                                  }
                                                } else {
                                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                    SPEC_NAME <- NA
                                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                                  } else {
                                                    stop()
                                                  }
                                                }
                                              }
                                            } else {
                                              stop()
                                            }
                                          } else {
                                            if(nl$taxonomicStatus == "DOUBTFUL"){
                                              NAME <- nl$scientificName
                                              STATUS <- nl$taxonomicStatus
                                              RANK <- nl$rank
                                              NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                              if(nl$rank == "SPECIES"){
                                                SPEC_NAME <- nl$scientificName
                                              } else {
                                                if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                                  if("speciesKey" %in% names(nl)){
                                                    nl <- name_usage(key=nl$speciesKey)$data
                                                    SPEC_NAME <- nl$scientificName
                                                  } else {
                                                    SPEC_NAME <- NA
                                                  }
                                                } else {
                                                  if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                    SPEC_NAME <- NA
                                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                                  } else {
                                                    stop()
                                                  }
                                                }
                                              }
                                            } else {
                                              stop()
                                            }
                                          }
                                        }
                                      } else {
                                        if("acceptedKey" %in% names(nl)){
                                          nl <- name_usage(key=nl$acceptedKey)$data
                                          if(nl$taxonomicStatus == "ACCEPTED"){
                                            NAME <- nl$scientificName
                                            STATUS <- nl$taxonomicStatus
                                            NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                            RANK <- nl$rank
                                            if(nl$rank == "SPECIES"){
                                              SPEC_NAME <- nl$scientificName
                                            } else {
                                              if(nl$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                                if("speciesKey" %in% names(nl)){
                                                  nl <- name_usage(key=nl$speciesKey)$data
                                                  SPEC_NAME <- nl$scientificName
                                                } else {
                                                  SPEC_NAME <- NA
                                                }
                                              } else {
                                                if(nl$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                  SPEC_NAME <- NA
                                                  NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                                } else {
                                                  stop()
                                                }
                                              }
                                            }
                                          } else {
                                            stop()
                                          }
                                        } else {
                                          stop()
                                        }
                                      }
                                    } else {
                                      stop()
                                    }
                                  } else {
                                    if(nrow(nl) > 1){
                                      success <- FALSE
                                      for(ikey in nl$key){
                                        nl1 <- name_usage(key=ikey)$data
                                        if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                          success <- TRUE
                                          NAME <- nl1$scientificName
                                          STATUS <- nl1$taxonomicStatus
                                          NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                          RANK <- nl1$rank
                                          if(nl1$rank == "SPECIES"){
                                            SPEC_NAME <- nl1$scientificName
                                          } else {
                                            if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                              if("speciesKey" %in% names(nl1)){
                                                nl1 <- name_usage(key=nl1$speciesKey)$data
                                                SPEC_NAME <- nl1$scientificName
                                              } else {
                                                SPEC_NAME <- NA
                                              }
                                            } else {
                                              if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                SPEC_NAME <- NA
                                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                              } else {
                                                stop()
                                              }
                                            }
                                          }
                                          break
                                        }
                                      }
                                      if(success == FALSE){
                                        for(ikey in nl$key){
                                          nl1 <- name_usage(key=ikey)$data
                                          if("acceptedKey" %in% names(nl1)){
                                            nl1 <- name_usage(key=nl1$acceptedKey)$data
                                            
                                            if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                              success <- TRUE
                                              NAME <- nl1$scientificName
                                              STATUS <- nl1$taxonomicStatus
                                              NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                              RANK <- nl1$rank
                                              if(nl1$rank == "SPECIES"){
                                                SPEC_NAME <- nl1$scientificName
                                              } else {
                                                if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                                  if("speciesKey" %in% names(nl1)){
                                                    nl1 <- name_usage(key=nl1$speciesKey)$data
                                                    SPEC_NAME <- nl1$scientificName
                                                  } else {
                                                    SPEC_NAME <- NA
                                                  }
                                                } else {
                                                  if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                    SPEC_NAME <- NA
                                                    NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                                  } else {
                                                    stop()
                                                  }
                                                }
                                              }
                                              break
                                            }
                                          }
                                        }
                                      }
                                      if(success == FALSE){
                                        NAME <- NA
                                        STATUS <- NA
                                        RANK <- NA
                                        NOTE = "NOT FOUND"
                                        SPEC_NAME <- NA
                                      }
                                    } else {
                                      stop()
                                    }
                                  }
                                } else {
                                  success <- FALSE
                                  for(ikey in nl$key){
                                    nl1 <- name_usage(key=ikey)$data
                                    if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                      success <- TRUE
                                      NAME <- nl1$scientificName
                                      STATUS <- nl1$taxonomicStatus
                                      RANK <- nl1$rank
                                      NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                      if(nl1$rank == "SPECIES"){
                                        SPEC_NAME <- nl1$scientificName
                                      } else {
                                        if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                          if("speciesKey" %in% names(nl1)){
                                            nl1 <- name_usage(key=nl1$speciesKey)$data
                                            SPEC_NAME <- nl1$scientificName
                                          } else {
                                            SPEC_NAME <- NA
                                          }
                                        } else {
                                          if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                            SPEC_NAME <- NA
                                            NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                          } else {
                                            stop()
                                          }
                                        }
                                      }
                                      break
                                    }
                                  }
                                  if(success == FALSE){
                                    for(ikey in nl$key){
                                      nl1 <- name_usage(key=ikey)$data
                                      if("acceptedKey" %in% names(nl1)){
                                        nl1 <- name_usage(key=nl1$acceptedKey)$data
                                        
                                        if(nl1$taxonomicStatus == "ACCEPTED" & ifelse("kingdom" %in% names(nl1), nl1$kingdom == "Plantae", TRUE)){
                                          success <- TRUE
                                          NAME <- nl1$scientificName
                                          STATUS <- nl1$taxonomicStatus
                                          RANK <- nl1$rank
                                          NOTE = "FUZZY MATCHING! PARSED TO SPECIES"
                                          if(nl1$rank == "SPECIES"){
                                            SPEC_NAME <- nl1$scientificName
                                          } else {
                                            if(nl1$rank %in% c("FORM","SUBSPECIES","VARIETY","INFRASPECIFIC_NAME")){
                                              if("speciesKey" %in% names(nl1)){
                                                nl1 <- name_usage(key=nl1$speciesKey)$data
                                                SPEC_NAME <- nl1$scientificName
                                              } else {
                                                SPEC_NAME <- NA
                                              }
                                            } else {
                                              if(nl1$rank %in% c("GENUS","FAMILY","CLASS","KINGDOM","ORDER","PHYLUM","STRAIN","DOMAIN")){
                                                SPEC_NAME <- NA
                                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                                              } else {
                                                stop()
                                              }
                                            }
                                          }
                                          break
                                        }
                                      }
                                    }
                                  }
                                  if(success == FALSE){
                                    NAME <- NA
                                    STATUS <- NA
                                    RANK <- NA
                                    NOTE = "NOT FOUND"
                                    SPEC_NAME <- NA
                                  }
                                }
                              }
                            } else {
                              stop()
                            }
                          } else {
                            if(length(strsplit(i, " ")[[1]]) == 2){
                              if(strsplit(i, " ")[[1]][2] == "spec."){
                                NAME <- NA
                                STATUS <- NA
                                SPEC_NAME <- NA
                                RANK <- NA
                                NOTE <- "FUZZY MATCHING! LIKELY HIGHER LEVEL TAXON"
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
                
              } else {
                if(is.null(NAME)){
                  NAME <- NA
                  STATUS <- NA
                  NOTE = "NOT FOUND"
                  RANK <- NA
                  SPEC_NAME <- NA
                } else {
                  stop()
                }
              }
            }
          },silent=TRUE)
        }
      }
      
      if(is.null(NAME)){
        resolved_temp2 <- data.frame(orig_name = orig_name,
                   used_name = ifelse(is.null(i), orig_name, i),
                   found_name = "FAILED",
                   taxon_status = "FAILED",
                   taxon_rank = "FAILED",
                   species_name = "FAILED",
                   canonical_name = "FAILED",
                   canonical_species_name = "FAILED",
                   dataset = names(which(tested_keys == dataset_key)),
                   note = ifelse(grepl("×",orig_name) | grepl(" x ",orig_name) ,"HYBRID! RESOLVING FAILED","FAILED"))
      } else {
        
        if(grepl("Incertae sedis", NAME, ignore.case = T)){
          resolved_temp2 <- data.frame(orig_name = orig_name,
                                       used_name = ifelse(is.null(i), orig_name, i),
                                       found_name = NA,
                                       taxon_status = NA,
                                       taxon_rank = NA,
                                       species_name = NA,
                                       canonical_name = NA,
                                       canonical_species_name = NA,
                                       dataset = names(which(tested_keys == dataset_key)),
                                       note = "INCERTEA SEDIS! HIGHLY UNCERTAIN TAXON")
        } else {
          e <- try(if(is.na(NAME)){
            parsed <- NA
            spec_parsed <- NA
          } else {
            if(NOTE == "LIKELY HIGHER LEVEL TAXON"){
              parsed <- parsenames(NAME)$canonicalnamewithmarker
              spec_parsed <- NA
            } else {
              parsed <- parsenames(NAME)$canonicalnamewithmarker
              spec_parsed <- try(paste(strsplit(parsenames(NAME)$canonicalname, " ")[[1]][1:2], collapse = " "), silent = T)
              if(is.null(parsed) & parsenames(NAME)$type == "HYBRID"){
                parsed <- NA
                spec_parsed <- NA
                NOTE <- "HYBRID! RESOLVING MIGHT BE UNRELIABLE"
              }
              if(class(spec_parsed) == "try-error"){
                spec_parsed <- NA
              }
            }
          })
          
          if(class(e) == "try-error"){
            resolved_temp2 <- data.frame(orig_name = orig_name,
                                         used_name = ifelse(is.null(i), orig_name, i),
                                         found_name = "FAILED",
                                         taxon_status = "FAILED",
                                         taxon_rank = "FAILED",
                                         species_name = "FAILED",
                                         canonical_name = "FAILED",
                                         canonical_species_name = "FAILED",
                                         dataset = names(which(tested_keys == dataset_key)),
                                         note = "FAILED! PROPABLE NETWORK ERROR")
          } else {
            
            if(grepl("UNRANKED", RANK)){
              SPEC_NAME <- NA
            }
            
            resolved_temp2 <- data.frame(orig_name = orig_name,
                                         used_name = i,
                                         found_name = NAME,
                                         taxon_status = STATUS,
                                         taxon_rank = RANK,
                                         species_name = SPEC_NAME,
                                         canonical_name = parsed,
                                         canonical_species_name = spec_parsed,
                                         dataset = names(which(tested_keys == dataset_key)),
                                         note = NOTE)
          }
        }
      }
      
      if((grepl("×", orig_name)|grepl(" x ", orig_name)|grepl(" X ", orig_name)) & !grepl("HYBRID",resolved_temp2$taxon_rank)){
        resolved_temp2$taxon_rank <- ifelse(is.na(resolved_temp2$taxon_rank), "LIKELY HYBRID", paste0(resolved_temp2$taxon_rank, "! LIKELY HYBRID"))
      }
      
      resolved_temp <- rbind.data.frame(resolved_temp, resolved_temp2)
      
    }
    # resolved_temp
    return(resolved_temp)
    
  } else {
    stop("One or more of the given datasets not included in the tested datasets")
  }
}
