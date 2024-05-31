for(de in seq(65, 85, 1)){
  print(de)
  for(de2 in seq(65, 85, 1)){
    print(de2)
    if(de == 65 & de2 == 65){
      df <- data.frame(matrix(ncol = 3, nrow = 0))
      x <- c("deIn", "deOut", "difference")
      colnames(df) <- x
    }
    #### trial DE for enteric ####
    ## ------------ Notes --------------  ##
    ## you can comment out this section to test other values
    DEinside <- de # 'housed' value
    DEoutside <- de2 # outside value
    inDE <- paste0("housed_", DEinside)
    outDE <- paste0("sml_", DEoutside)
    ## ------------ ----- --------------  ##
    
    dairyCH4v2 <- dairyCH4 %>%
      dplyr::select(animal, category
                    , contains(c(inDE, outDE)))
    ## should be four columns. If so, rename
    if(ncol(dairyCH4v2) == 4){
      names(dairyCH4v2)[3:4] <- c(paste0("kgCH4GE.MJday", c("in", "out"))) 
      names(dairyCH4v2)
    } else {
      stop("wrong columns in 4a2")
    }
    dairyCH4v2 <- dairyCH4v2 %>%
      mutate(kgCH4GE.50pcScen = (kgCH4GE.MJdayout/2) + (kgCH4GE.MJdayin/2)) %>%
      rowwise() %>%
      mutate(litDifference = kgCH4GE.50pcScen - 125.44)
    head(dairyCH4v2)
    dairyCH4v3 <- dairyCH4v2 %>%
      filter(category != "heifer")
    mean(dairyCH4v3$litDifference)
    
    df <- bind_rows(df, 
                    bind_cols(
                      deIn = de
                      , deOut = de2
                      , difference = round(mean(dairyCH4v3$litDifference), 2)
                    )
    )
  }
}

# save output
x <- c("deIn", "deOut", "ent_ferm_difference")
colnames(df) <- x
fwrite(df
       , file.path("results", "trials", "ent_ferm_DEs.csv")
       , row.names = F)

#### 2 - VS check ####
for(de in seq(65, 85, 1)){
  print(de)
  for(de2 in seq(65, 85, 1)){
    print(de2)
    if(de == 65 & de2 == 65){
      df <- data.frame(matrix(ncol = 3, nrow = 0))
      x <- c("deIn", "deOut", "difference")
      colnames(df) <- x
    }
    ## ------------ Notes --------------  ##
    ## you can comment out this section to test other values
    DEinside <- de # 'housed' value
    DEoutside <- de2 # outside value
    inDE <- paste0("housed_", DEinside)
    outDE <- paste0("sml_", DEoutside)
    ## ------------ ----- --------------  ##
    
    vssTest <- VStable %>%
      filter(grepl("lactat", category))
    
    limitTest <- vssTest %>%
      dplyr::select(animal, category
                    , contains(c(inDE, outDE)))
    range(limitTest[, 3]); mean(limitTest[, 3])
    range(limitTest[, 4]); mean(limitTest[, 4])
    mean(mean(limitTest[, 3]), mean(limitTest[, 4]))
    diff <- mean(c(mean(limitTest[, 3])
                   , mean(limitTest[, 4]))
    ) - 5.73
    
    df <- bind_rows(df, 
                    bind_cols(
                      deIn = de
                      # , m1 = mean(limitTest[, 3])
                      # , m2 = mean(limitTest[, 4])
                      , deOut = de2
                      , difference = round(diff, 2)
                    )
    )
  }
}

# save output
x <- c("deIn", "deOut", "vs_difference")
colnames(df) <- x
fwrite(df
       , file.path("results", "trials", "VS_DEs.csv")
       , row.names = F)

#### 3 - check for N excretion ####
for(de in seq(65, 85, 1)){
  print(de)
  for(de2 in seq(65, 85, 1)){
    print(de2)
    if(de == 65 & de2 == 65){
      df <- data.frame(matrix(ncol = 3, nrow = 0))
      x <- c("deIn", "deOut", "N_exc_difference")
      colnames(df) <- x
    }
    ## ------------ Notes --------------  ##
    ## you can comment out this section to test other values
    DEinside <- 85 # 'housed' value
    DEoutside <- 79 # outside value
    inDE <- paste0("housed_", DEinside)
    outDE <- paste0("sml_", DEoutside)
    ## ------------ ----- --------------  ##
    
    # get yearly values based on scenario
    ## extract the two bits required
    NexcreteCows <- Nintake.retain.excreteCows %>%
      dplyr::select(c(animal, category, UK.Region
                      , kgDayNreten
                      , contains(c(
                        paste0("housed_", as.character(DEinside))
                        , paste0("sml_", as.character(DEoutside)))))) %>%
      # multiply each of the excretions by 6 months (i.e. half a year)
      ## half a year inside
      mutate(across(contains("kgNExcreteDay_housed"), ~ . * (365/2), .names = "kgNExcreteYr_housed")) %>%
      ## half a year outside
      mutate(across(contains("kgNExcreteDay_sml"), ~ . * (365/2), .names = "kgNExcreteYr_sml")) %>%
      ## sum together, to get the year total
      mutate(kgNExcreteYr = kgNExcreteYr_housed + kgNExcreteYr_sml)
    head(NexcreteCows)
    
    # should be near: 130.7 kg N/ dairy cow (Bougouin et al., 2022)
    ## get means for all of the same animal, category
    summariseN0 <- NexcreteCows %>%
      group_by(animal, category) %>%
      dplyr::summarise(across(where(is.numeric)
                              , list(mean = ~ mean(., na.rm = TRUE)))) %>%
      filter(grepl("dairy", category)) %>%
      dplyr::select(animal, category, contains("ExcreteYr"))
    head(summariseN0)
    
    fin <- mean(summariseN0$kgNExcreteYr_mean) - 130.8
    
    df <- bind_rows(df, 
                    bind_cols(
                      deIn = de
                      , deOut = de2
                      , N_exc_difference = round(fin, 2)
                    ))
                    
  }
}

# save output
x <- c("deIn", "deOut", "N_exc_difference")
colnames(df) <- x
fwrite(df
       , file.path("results", "trials", "N_excrete_DEs.csv")
       , row.names = F)

#### 4 - combine all together ####
pIn <- pblapply(list.files(file.path("results", "trials"), full.names = T, pattern = "DEs.csv"), fread)
pIn2 <- Reduce(function(x, y) merge(x, y, by = c("deIn", "deOut")
                            , all = TRUE),
       pIn) 
# save
fwrite(pIn2
       , file.path("results", "trials", "DE_error_comparison.csv")
       , row.names = F)

