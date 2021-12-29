# orig_name <- "Pteridium pinetorum C. N. Page & R. R. Mill"
# dataset_key <- "d9a4eedb-e985-4456-ad46-3df8472e00e8"
# dataset <- c("f382f0ce-323a-4091-bb9f-add557f3a9a2","d9a4eedb-e985-4456-ad46-3df8472e00e8")
# resolve_taxon_name("Raphanus sativus var. sativus", dataset = dataset)
resolve_taxon_name <- function(orig_name, dataset = NULL, lib.loc = .libPaths()){
  
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
    
    resolved_temp <- data.frame() # To store names
    for(dataset_key in dataset){
      
      resolved_temp2 <- tryCatch({
        i <- orig_name
        nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
        
        if(nrow(nl) == 1){
          nl <- name_usage(key=nl$key)$data
          if(nrow(nl) > 1){
            nl <- nl[order(nl$taxonomicStatus),]
            nl <- nl[1,]
          }
          if("accepted" %in% names(nl)){
            if(nl$taxonomicStatus == "ACCEPTED"){
              NAME <- nl$scientificName
              STATUS <- nl$taxonomicStatus
            } else {
              nl <- name_usage(key=nl$acceptedKey)$data
              NAME <- nl$scientificName
              STATUS <- nl$taxonomicStatus
            }
          } else {
            
            if(nl$key == nl$speciesKey){
              if(nl$taxonomicStatus == "ACCEPTED"){
                NAME <- nl$scientificName
                STATUS <- nl$taxonomicStatus
              } else {
                if(nl$taxonomicStatus == "DOUBTFUL" & nl$synonym == FALSE){
                  NAME <- nl$scientificName
                  STATUS <- nl$taxonomicStatus
                } else {
                  NAME <- NA
                  STATUS <- NA
                }
              }
            } else {
              NAME <- NA
              STATUS <- NA
            }
          }
        } else {
          if(nrow(nl) == 0){
            i <- gbif_parse(i)$canonicalnamecomplete
            
            nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
            
            if(nrow(nl) == 1){
              
              nl <- name_usage(key=nl$key)$data
              if("acceptedKey" %in% names(nl)){
                nl <- name_usage(key=nl$acceptedKey)$data
              }
              if(nl$rank == "SUBSPECIES"){
                if(nl$taxonomicStatus == "ACCEPTED"){
                  NAME <- nl$scientificName
                  STATUS <- nl$taxonomicStatus
                } else {
                  NAME <- NA
                  STATUS <- NA
                }
              } else {
                if(nl$key == nl$speciesKey){
                  if(nl$taxonomicStatus == "ACCEPTED"){
                    NAME <- nl$scientificName
                    STATUS <- nl$taxonomicStatus
                  } else {
                    NAME <- NA
                    STATUS <- NA
                  }
                } else {
                  nl <- name_usage(key=nl$speciesKey)$data
                  if(nl$key == nl$speciesKey){
                    if(nl$taxonomicStatus == "ACCEPTED"){
                      NAME <- nl$scientificName
                      STATUS <- nl$taxonomicStatus
                    } else {
                      NAME <- NA
                      STATUS <- NA
                    }
                  } else {
                    NAME <- NA
                    STATUS <- NA
                  }
                }
              }
              
            } else {
              if(nrow(nl) > 1){
                success <- FALSE
                for(ikey in nl$key){
                  nl1 <- name_usage(key=ikey)$data
                  if(nl1$taxonomicStatus == "ACCEPTED"){
                    success <- TRUE
                    NAME <- nl1$scientificName
                    STATUS <- nl1$taxonomicStatus
                    break
                  }
                }
                if(success == FALSE){
                  for(ikey in nl$key){
                    nl1 <- name_usage(key=ikey)$data
                    if("acceptedKey" %in% names(nl1)){
                      nl1 <- name_usage(key=nl1$acceptedKey)$data
                      
                      if(nl1$taxonomicStatus == "ACCEPTED"){
                        success <- TRUE
                        NAME <- nl1$scientificName
                        STATUS <- nl1$taxonomicStatus
                        break
                      }
                    }
                  }
                }
                if(success == FALSE){
                  NAME <- NA
                  STATUS <- NA
                }
              } else {
                
                i <- gbif_parse(i)$canonicalname
                
                nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
                
                if(nrow(nl) == 1){
                  
                  nl <- name_usage(key=nl$key)$data
                  
                  if(nrow(nl) == 100){
                    NAME <- NA
                    STATUS <- NA
                  } else {
                    if(nl$taxonomicStatus == "ACCEPTED"){
                      NAME <- nl$scientificName
                      STATUS <- nl$taxonomicStatus
                    } else {
                      if("acceptedKey" %in% names(nl)){
                        nl1 <- name_usage(key=nl$acceptedKey)$data
                        
                        if(nl1$taxonomicStatus == "ACCEPTED"){
                          NAME <- nl1$scientificName
                          STATUS <- nl1$taxonomicStatus
                        } else {
                          NAME <- NA
                          STATUS <- NA
                        }
                      } else {
                        NAME <- NA
                        STATUS <- NA
                      }
                    }
                  }
                  
                } else {
                  if(nrow(nl) > 1){
                    
                    success <- FALSE
                    for(ikey in nl$key){
                      nl1 <- name_usage(key=ikey)$data
                      if(nl1$taxonomicStatus == "ACCEPTED"){
                        success <- TRUE
                        NAME <- nl1$scientificName
                        STATUS <- nl1$taxonomicStatus
                        break
                      }
                    }
                    if(success == FALSE){
                      for(ikey in nl$key){
                        nl1 <- name_usage(key=ikey)$data
                        if("acceptedKey" %in% names(nl1)){
                          nl1 <- name_usage(key=nl1$acceptedKey)$data
                          
                          if(nl1$taxonomicStatus == "ACCEPTED"){
                            success <- TRUE
                            NAME <- nl1$scientificName
                            STATUS <- nl1$taxonomicStatus
                            break
                          }
                        }
                      }
                    }
                    if(success == FALSE){
                      NAME <- NA
                      STATUS <- NA
                    }
                  } else {
                    
                    i <- gbif_parse(i)$canonicalname
                    
                    if(length(strsplit(i, " ")[[1]]) > 2){
                      
                      i <- paste(strsplit(i, " ")[[1]][1:2], collapse = " ")
                      
                      nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
                      
                      if(nrow(nl) == 1){
                        
                        nl <- name_usage(key=nl$key)$data
                        if("acceptedKey" %in% names(nl)){
                          nl <- name_usage(key=nl$acceptedKey)$data
                        }
                        if(nl$rank == "SUBSPECIES"){
                          if(nl$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl$scientificName
                            STATUS <- nl$taxonomicStatus
                          } else {
                            NAME <- NA
                            STATUS <- NA
                          }
                        } else {
                          if(nl$key == nl$speciesKey){
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                            } else {
                              NAME <- NA
                              STATUS <- NA
                            }
                          } else {
                            nl <- name_usage(key=nl$speciesKey)$data
                            if(nl$key == nl$speciesKey){
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                              } else {
                                NAME <- NA
                                STATUS <- NA
                              }
                            } else {
                              NAME <- NA
                              STATUS <- NA
                            }
                          }
                        }
                      } else {
                        if(nrow(nl) > 1){
                          success <- FALSE
                          for(ikey in nl$key){
                            nl1 <- name_usage(key=ikey)$data
                            if(nl1$taxonomicStatus == "ACCEPTED"){
                              success <- TRUE
                              NAME <- nl1$scientificName
                              STATUS <- nl1$taxonomicStatus
                              break
                            }
                          }
                          if(success == FALSE){
                            for(ikey in nl$key){
                              nl1 <- name_usage(key=ikey)$data
                              if("acceptedKey" %in% names(nl1)){
                                nl1 <- name_usage(key=nl1$acceptedKey)$data
                                
                                if(nl1$taxonomicStatus == "ACCEPTED"){
                                  success <- TRUE
                                  NAME <- nl1$scientificName
                                  STATUS <- nl1$taxonomicStatus
                                  break
                                }
                              }
                            }
                          }
                          if(success == FALSE){
                            NAME <- NA
                            STATUS <- NA
                          }
                        } else {
                          if(nrow(nl) == 0){
                            NAME <- NA
                            STATUS <- NA
                          } else {
                            NAME <- NA
                            STATUS <- NA
                          }
                        }
                        
                      }
                      
                    } else {
                      if(length(strsplit(i, " ")[[1]]) == 2){
                        NAME <- NA
                        STATUS <- NA
                      } else {
                        if(length(strsplit(i, " ")[[1]]) == 1){
                          NAME <- NA
                          STATUS <- NA
                        } else {
                          NAME <- NA
                          STATUS <- NA
                        }
                      }
                    }
                    
                  }
                  
                }
                
              }
              
            }
          } else {
            success <- FALSE
            for(ikey in nl$key){
              nl1 <- name_usage(key=ikey)$data
              if(nl1$taxonomicStatus == "ACCEPTED"){
                success <- TRUE
                NAME <- nl1$scientificName
                STATUS <- nl1$taxonomicStatus
                break
              }
            }
            if(success == FALSE){
              for(ikey in nl$key){
                nl1 <- name_usage(key=ikey)$data
                if("acceptedKey" %in% names(nl1)){
                  nl1 <- name_usage(key=nl1$acceptedKey)$data
                  
                  if(nl1$taxonomicStatus == "ACCEPTED"){
                    success <- TRUE
                    NAME <- nl1$scientificName
                    STATUS <- nl1$taxonomicStatus
                    break
                  }
                }
              }
            }
            if(success == FALSE){
              NAME <- NA
              STATUS <- NA
            }
          }
        }
        
        if(is.na(NAME)){
          parsed <- NA
        } else {
          parsed <- parsenames(NAME)$canonicalname
        }
        
        data.frame(orig_name = orig_name,
                   used_name = i,
                   found_name = NAME,
                   taxon_status = STATUS,
                   canonical_name = parsed,
                   dataset = names(which(tested_keys == dataset_key)))
        
      }, error = function(e) {
        tryCatch({
          Sys.sleep(0.5)
          
          i <- orig_name
          nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
          
          if(nrow(nl) == 1){
            nl <- name_usage(key=nl$key)$data
            if(nrow(nl) > 1){
              nl <- nl[order(nl$taxonomicStatus),]
              nl <- nl[1,]
            }
            if("accepted" %in% names(nl)){
              if(nl$taxonomicStatus == "ACCEPTED"){
                NAME <- nl$scientificName
                STATUS <- nl$taxonomicStatus
              } else {
                nl <- name_usage(key=nl$acceptedKey)$data
                NAME <- nl$scientificName
                STATUS <- nl$taxonomicStatus
              }
            } else {
              
              if(nl$key == nl$speciesKey){
                if(nl$taxonomicStatus == "ACCEPTED"){
                  NAME <- nl$scientificName
                  STATUS <- nl$taxonomicStatus
                } else {
                  if(nl$taxonomicStatus == "DOUBTFUL" & nl$synonym == FALSE){
                    NAME <- nl$scientificName
                    STATUS <- nl$taxonomicStatus
                  } else {
                    NAME <- NA
                    STATUS <- NA
                  }
                }
              } else {
                NAME <- NA
                STATUS <- NA
              }
            }
          } else {
            if(nrow(nl) == 0){
              i <- gbif_parse(i)$canonicalnamecomplete
              
              nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
              
              if(nrow(nl) == 1){
                
                nl <- name_usage(key=nl$key)$data
                if("acceptedKey" %in% names(nl)){
                  nl <- name_usage(key=nl$acceptedKey)$data
                }
                if(nl$rank == "SUBSPECIES"){
                  if(nl$taxonomicStatus == "ACCEPTED"){
                    NAME <- nl$scientificName
                    STATUS <- nl$taxonomicStatus
                  } else {
                    NAME <- NA
                    STATUS <- NA
                  }
                } else {
                  if(nl$key == nl$speciesKey){
                    if(nl$taxonomicStatus == "ACCEPTED"){
                      NAME <- nl$scientificName
                      STATUS <- nl$taxonomicStatus
                    } else {
                      NAME <- NA
                      STATUS <- NA
                    }
                  } else {
                    nl <- name_usage(key=nl$speciesKey)$data
                    if(nl$key == nl$speciesKey){
                      if(nl$taxonomicStatus == "ACCEPTED"){
                        NAME <- nl$scientificName
                        STATUS <- nl$taxonomicStatus
                      } else {
                        NAME <- NA
                        STATUS <- NA
                      }
                    } else {
                      NAME <- NA
                      STATUS <- NA
                    }
                  }
                }
                
              } else {
                if(nrow(nl) > 1){
                  success <- FALSE
                  for(ikey in nl$key){
                    nl1 <- name_usage(key=ikey)$data
                    if(nl1$taxonomicStatus == "ACCEPTED"){
                      success <- TRUE
                      NAME <- nl1$scientificName
                      STATUS <- nl1$taxonomicStatus
                      break
                    }
                  }
                  if(success == FALSE){
                    for(ikey in nl$key){
                      nl1 <- name_usage(key=ikey)$data
                      if("acceptedKey" %in% names(nl1)){
                        nl1 <- name_usage(key=nl1$acceptedKey)$data
                        
                        if(nl1$taxonomicStatus == "ACCEPTED"){
                          success <- TRUE
                          NAME <- nl1$scientificName
                          STATUS <- nl1$taxonomicStatus
                          break
                        }
                      }
                    }
                  }
                  if(success == FALSE){
                    NAME <- NA
                    STATUS <- NA
                  }
                } else {
                  
                  i <- gbif_parse(i)$canonicalname
                  
                  nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
                  
                  if(nrow(nl) == 1){
                    
                    nl <- name_usage(key=nl$key)$data
                    
                    if(nrow(nl) == 100){
                      NAME <- NA
                      STATUS <- NA
                    } else {
                      if(nl$taxonomicStatus == "ACCEPTED"){
                        NAME <- nl$scientificName
                        STATUS <- nl$taxonomicStatus
                      } else {
                        if("acceptedKey" %in% names(nl)){
                          nl1 <- name_usage(key=nl$acceptedKey)$data
                          
                          if(nl1$taxonomicStatus == "ACCEPTED"){
                            NAME <- nl1$scientificName
                            STATUS <- nl1$taxonomicStatus
                          } else {
                            NAME <- NA
                            STATUS <- NA
                          }
                        } else {
                          NAME <- NA
                          STATUS <- NA
                        }
                      }
                    }
                    
                  } else {
                    if(nrow(nl) > 1){
                      
                      success <- FALSE
                      for(ikey in nl$key){
                        nl1 <- name_usage(key=ikey)$data
                        if(nl1$taxonomicStatus == "ACCEPTED"){
                          success <- TRUE
                          NAME <- nl1$scientificName
                          STATUS <- nl1$taxonomicStatus
                          break
                        }
                      }
                      if(success == FALSE){
                        for(ikey in nl$key){
                          nl1 <- name_usage(key=ikey)$data
                          if("acceptedKey" %in% names(nl1)){
                            nl1 <- name_usage(key=nl1$acceptedKey)$data
                            
                            if(nl1$taxonomicStatus == "ACCEPTED"){
                              success <- TRUE
                              NAME <- nl1$scientificName
                              STATUS <- nl1$taxonomicStatus
                              break
                            }
                          }
                        }
                      }
                      if(success == FALSE){
                        NAME <- NA
                        STATUS <- NA
                      }
                    } else {
                      
                      i <- gbif_parse(i)$canonicalname
                      
                      if(length(strsplit(i, " ")[[1]]) > 2){
                        
                        i <- paste(strsplit(i, " ")[[1]][1:2], collapse = " ")
                        
                        nl <- name_suggest(q=i, rank = "species", datasetKey = dataset_key)$data
                        
                        if(nrow(nl) == 1){
                          
                          nl <- name_usage(key=nl$key)$data
                          if("acceptedKey" %in% names(nl)){
                            nl <- name_usage(key=nl$acceptedKey)$data
                          }
                          if(nl$rank == "SUBSPECIES"){
                            if(nl$taxonomicStatus == "ACCEPTED"){
                              NAME <- nl$scientificName
                              STATUS <- nl$taxonomicStatus
                            } else {
                              NAME <- NA
                              STATUS <- NA
                            }
                          } else {
                            if(nl$key == nl$speciesKey){
                              if(nl$taxonomicStatus == "ACCEPTED"){
                                NAME <- nl$scientificName
                                STATUS <- nl$taxonomicStatus
                              } else {
                                NAME <- NA
                                STATUS <- NA
                              }
                            } else {
                              nl <- name_usage(key=nl$speciesKey)$data
                              if(nl$key == nl$speciesKey){
                                if(nl$taxonomicStatus == "ACCEPTED"){
                                  NAME <- nl$scientificName
                                  STATUS <- nl$taxonomicStatus
                                } else {
                                  NAME <- NA
                                  STATUS <- NA
                                }
                              } else {
                                NAME <- NA
                                STATUS <- NA
                              }
                            }
                          }
                        } else {
                          if(nrow(nl) > 1){
                            success <- FALSE
                            for(ikey in nl$key){
                              nl1 <- name_usage(key=ikey)$data
                              if(nl1$taxonomicStatus == "ACCEPTED"){
                                success <- TRUE
                                NAME <- nl1$scientificName
                                STATUS <- nl1$taxonomicStatus
                                break
                              }
                            }
                            if(success == FALSE){
                              for(ikey in nl$key){
                                nl1 <- name_usage(key=ikey)$data
                                if("acceptedKey" %in% names(nl1)){
                                  nl1 <- name_usage(key=nl1$acceptedKey)$data
                                  
                                  if(nl1$taxonomicStatus == "ACCEPTED"){
                                    success <- TRUE
                                    NAME <- nl1$scientificName
                                    STATUS <- nl1$taxonomicStatus
                                    break
                                  }
                                }
                              }
                            }
                            if(success == FALSE){
                              NAME <- NA
                              STATUS <- NA
                            }
                          } else {
                            if(nrow(nl) == 0){
                              NAME <- NA
                              STATUS <- NA
                            } else {
                              NAME <- NA
                              STATUS <- NA
                            }
                          }
                          
                        }
                        
                      } else {
                        if(length(strsplit(i, " ")[[1]]) == 2){
                          NAME <- NA
                          STATUS <- NA
                        } else {
                          if(length(strsplit(i, " ")[[1]]) == 1){
                            NAME <- NA
                            STATUS <- NA
                          } else {
                            NAME <- NA
                            STATUS <- NA
                          }
                        }
                      }
                      
                    }
                    
                  }
                  
                }
                
              }
            } else {
              success <- FALSE
              for(ikey in nl$key){
                nl1 <- name_usage(key=ikey)$data
                if(nl1$taxonomicStatus == "ACCEPTED"){
                  success <- TRUE
                  NAME <- nl1$scientificName
                  STATUS <- nl1$taxonomicStatus
                  break
                }
              }
              if(success == FALSE){
                for(ikey in nl$key){
                  nl1 <- name_usage(key=ikey)$data
                  if("acceptedKey" %in% names(nl1)){
                    nl1 <- name_usage(key=nl1$acceptedKey)$data
                    
                    if(nl1$taxonomicStatus == "ACCEPTED"){
                      success <- TRUE
                      NAME <- nl1$scientificName
                      STATUS <- nl1$taxonomicStatus
                      break
                    }
                  }
                }
              }
              if(success == FALSE){
                NAME <- NA
                STATUS <- NA
              }
            }
          }
          
          if(is.na(NAME)){
            parsed <- NA
          } else {
            parsed <- parsenames(NAME)$canonicalname
          }
          
          data.frame(orig_name = orig_name,
                     used_name = i,
                     found_name = NAME,
                     taxon_status = STATUS,
                     canonical_name = parsed,
                     dataset = names(which(tested_keys == dataset_key)))
          
        }, error = function(e) {
          
          data.frame(orig_name = orig_name,
                     used_name = i,
                     found_name = "FAILED",
                     taxon_status = "FAILED",
                     canonical_name = "FAILED",
                     dataset = names(which(tested_keys == dataset_key)))
        })
      })
      
      resolved_temp <- rbind.data.frame(resolved_temp, resolved_temp2)
      
    }
    # resolved_temp
    return(resolved_temp)
    
  } else {
    stop("One or more of the given datasets not included in the tested datasets")
  }
}
