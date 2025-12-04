## ----setup------------------------------------------------------------------------------------------------------
#Data handling
library(dplyr)
library(tidyr)
#API wrappers
library(rgbif)
library(rredlist)
#Handling spatial data
library(terra)
library(sf)
#Needed for manual GBIF API usage
library(httr)
library(jsonlite)
#Modelling and diagnostics
library(doParallel)
library(lme4)
library(broom.mixed)
#Plots
library(sjPlot)
library(kableExtra)

#Options for loading data

#Land use data: training years (match each gbif observation in these years to a land use category)
lu_years <- paste(seq(2010,2019))

#Land use categories:
# HILDA+ levels:
# 00 ocean
# 11 urban
# 22 cropland
# 33 pasture/rangeland
# 44 forest
# 55 unmanaged grass/shrubland
# 66 sparse/no vegetation
# 77 water
# 99 no data
lu_levels <- data.frame(
  ID = c(00,11,22,33,44,55,66,77,99),
  plum_category = c("Ocean","Urban","Cropland","Pasture/Rangeland","Forest", "Unmanaged grass/shrubland","Sparse/no vegetation","Water","No data")
)

#Habitat levels
hab_levels <- data.frame(
  ID = c("hab_1","hab_2","hab_3","hab_4","hab_5","hab_6","hab_8","hab_14_arable","hab_14_pasture","hab_14_urban","hab_14_degraded","hab_15"),
  iucn_habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetlands (inland)","Rocky Areas","Desert","Artificial (arable)","Artificial (pasture)","Artificial (urban and gardens)","Artificial (plantation and degraded)","Artificial (aquatic)")
)

#TEST - NOT RUN
# aves <- occ_download_prep(
#   pred("taxonKey", 3119195),
#   pred_gt("elevation", 5000),
#   user = "cmaney",
#   pwd = "RC)Y}GHr.i8$YE8",
#   email = "calummaney@gmail.com")

myredlistkey <- 'kYaDcjScQUgVR2EkcCQxU8ZTnY7v3wMs9LiF'


## ----gbifCall---------------------------------------------------------------------------------------------------
print("Checking downloads...")
if(length(list.files("raw_data/gbif/"))<2){
gbif_call <- occ_download(
  pred_or(
    pred("taxonKey", 212), #birds
    pred("taxonKey", 359), #mammals
    pred("taxonKey", 131), #amphibia
    pred("taxonKey", 358) #reptiles
  ),
  
  pred("hasGeospatialIssue", FALSE),
  
  pred_gte("year", 2010), #year start
  pred_lte("year", 2019), #year end
  
  #Points need coordinates
  pred_not(pred_isnull("decimalLatitude")),
  pred_not(pred_isnull("decimalLongitude")),
  
  #Keep good observations only
  pred_in("basisOfRecord",c("HUMAN_OBSERVATION","MACHINE_OBSERVATION","PRESERVED_SPECIMEN")),
  
  #Keep only low-uncertainty locations
  pred_lte("coordinateUncertaintyInMeters", 100),
  
  #Retrieve simplified data
  format = "SIMPLE_CSV",
  
  user = "cmaney",
  pwd = "RC)Y}GHr.i8$YE8",
  email = "calummaney@gmail.com")

occ_download_wait(gbif_call)
  
dir.create("raw_data/gbif/",showWarnings = F)

zip_file_path <- occ_download_get(key = gbif_call[1], path = "raw_data/gbif/", overwrite = TRUE)

saveRDS(zip_file_path,"raw_data/gbif/zip_file_path.rds")
}


## ----downloadHILDA----------------------------------------------------------------------------------------------
if(!file.exists("raw_data/hilda/hildap_vGLOB-1.0_geotiff_wgs84/hildap_GLOB-v1.0_lulc-states/")){
hilda_url <- "https://download.pangaea.de/dataset/921846/files/hildap_vGLOB-1.0_geotiff.zip"

hilda_zip <- download.file(hilda_url, destfile = "raw_data/hilda_plus.zip", method="curl")

unzip(zipfile = "raw_data/hilda_plus.zip",exdir = "raw_data/hilda/")

#Delete the .zip to tidy up
if (file.exists("raw_data/hilda_plus.zip")) {
  #Delete file if it exists
  file.remove("raw_data/hilda_plus.zip")
}
}


## ----gbifUnzip--------------------------------------------------------------------------------------------------
print("Loading GBIF data...")
zip_file_path <- readRDS("raw_data/gbif/zip_file_path.rds")

gbif_data <- occ_download_import(x = zip_file_path,path = ".") #est. 5 mins


## ----gbifGetIUCNIDfunction--------------------------------------------------------------------------------------
getIUCNtaxonID <- function(gbif_taxonKey = gbif_data[1,]$taxonKey){
  res <- GET(paste0("http://api.gbif.org/v1/species/",gbif_taxonKey,"/iucnRedListCategory"))
  
  if(length(res$content)>0){
  
  Sys.sleep(runif(n = 1,min=0,max=100)/1000)
  
  res <- try(fromJSON(rawToChar(res$content))$iucnTaxonID)
  
  return(res)
  } else{
    return(NA)
  }
}


## ----gbifGetAndSaveIUCNtags-------------------------------------------------------------------------------------
print("Checking IUCN Habitat data...")
if(!file.exists("wip/gbif_data_IUCN_codes.rds")){
#We only needed to run this once
gbif_taxonKeys <- unique(gbif_data$taxonKey) 

gbif_data_IUCN_codes <- lapply(gbif_taxonKeys,getIUCNtaxonID)

gbif_data_IUCN_codes[sapply(gbif_data_IUCN_codes, is.null)] <- NA

#Save the data
saveRDS(gbif_data_IUCN_codes,"wip/gbif_data_IUCN_codes.rds")
}


## ----gbifMargeIUCNKeys------------------------------------------------------------------------------------------
gbif_data_IUCN_codes <- readRDS("wip/gbif_data_IUCN_codes.rds")
gbif_taxonKeys <- unique(gbif_data$taxonKey) 

taxa.df <- data.frame(
  gbif_ID = gbif_taxonKeys,
  iucn_ID = gbif_data_IUCN_codes |> unlist()
)


## ----iucnGetTaxonCodes------------------------------------------------------------------------------------------
#iucn_hab <- read.csv("raw_data/Habitats/WCMC_Habitat_Info2_Aug2014_AA.csv")

iucnTaxonCodes <- gbif_data_IUCN_codes |> unlist()
iucnTaxonCodes <- iucnTaxonCodes[!is.na(iucnTaxonCodes)]
iucnTaxonCodes <- gsub("_.*","",iucnTaxonCodes)
iucnTaxonCodes <- as.numeric(iucnTaxonCodes)


## ----iucnGetHabsFromAPI-----------------------------------------------------------------------------------------
if(!file.exists("wip/iucnHabitats.rds")){
getHabitats <- function(i){
  #API call to retrieve the latest assessment of a given species, and 
  assessmentHabitats <- rl_sis_latest(
    i, 
    scope = "1",
    key = myredlistkey,
    parse = TRUE
    )$habitats
  
  habitats.out <- data.frame(
    taxonid = i,
    hab_code = assessmentHabitats$code,
    suitability = assessmentHabitats$suitability,
    season = assessmentHabitats$season,
    pause_base = 2,
    pause_cap = 3,
    pause_min = 1,
    times = 10
  )
  
  Sys.sleep((runif(1,0,75)/100)+0.5) #Need to wait for the API

  return(habitats.out)
}

iucnHabitats <- lapply(iucnTaxonCodes,getHabitats) |>
  bind_rows() |>
  unique.data.frame()

saveRDS(iucnHabitats,"wip/iucnHabitats.rds")
}


## ----joinDataPrepareIucnHabs------------------------------------------------------------------------------------
print("Joining habitat data...")
iucn_hab <- readRDS("wip/iucnHabitats.rds") #Sub in new API data

iucn_hab_suit <- iucn_hab[iucn_hab$suitability == "Suitable",]
iucn_hab_suit$taxonid <- paste(iucn_hab_suit$taxonid)

iucn_hab_suit <- iucn_hab_suit[,c("taxonid","hab_code","suitability")] |> unique.data.frame()



## ----joinRemoveGeneralists--------------------------------------------------------------------------------------
taxaToUse <- iucn_hab_suit |> 
  group_by(taxonid) |> 
  mutate(l1hab = gsub("_.*","",hab_code)) |>
  select(taxonid,l1hab,suitability) |>
  unique.data.frame() |>
  reframe(nHabs = sum(suitability=="Suitable")) |> 
  filter(nHabs <= 5) |> 
  pull(taxonid)

iucn_hab_suit <- iucn_hab_suit |> filter(taxonid %in% taxaToUse)


## ----joinIucnToGBIFdata-----------------------------------------------------------------------------------------
taxa_iucn.df <- left_join(taxa.df,iucn_hab_suit,by=c("iucn_ID"="taxonid")) |>
  filter(suitability == "Suitable")


## ----joinPivotHabitatsOut---------------------------------------------------------------------------------------
taxa_iucn_wide.df <- pivot_wider(taxa_iucn.df,names_from = hab_code,values_from = suitability)


## ----joinFilterNoHabitatData------------------------------------------------------------------------------------
gbif_withHabs <- gbif_data |> 
  filter(taxonKey %in% taxa_iucn.df$gbif_ID) |>
  transmute(
    gbifID,
    year,
    taxonKey,
    decimalLongitude,
    decimalLatitude,
    countryCode
  )


## ----joinHabitatsToObservations---------------------------------------------------------------------------------
obs_habs.df <- gbif_withHabs |>
  left_join(taxa_iucn_wide.df, by = c("taxonKey"="gbif_ID"))


## ----prepareHILDA-----------------------------------------------------------------------------------------------
print("Attaching land use data...")
hilda_files <- lapply(lu_years,function(x){list.files(path = "raw_data/hilda/hildap_vGLOB-1.0_geotiff_wgs84/hildap_GLOB-v1.0_lulc-states/",pattern = x,full.names = TRUE)}) |> unlist()

#Load and stack the path vector into one raster
hilda <- rast(hilda_files)

#Set the land use data to be categorical
hilda <- hilda |> as.factor() #this will take a while (est. 5-10 mins)

#Set names for the layers
names(hilda) <- paste0("X",lu_years)

#Set the correct levels for each layer in the data
for(thisYr in names(hilda)){
 levels(hilda[[thisYr]]) <- lu_levels
}

#Finished, prepared data.


## ----luAttachToObservations-------------------------------------------------------------------------------------
getLUfromObs <- function(obs, year = 2019){
  #Construct sf from the observation records for this year
  thisYear_obs <- obs[obs$year == year,]
  thisYear_obs.sf <- st_as_sf(thisYear_obs,coords=c("decimalLongitude","decimalLatitude"))
  
  thisLU <- extract(hilda[[paste0("X",year)]],thisYear_obs.sf)
  
  thisYear_obs.sf$HILDA_LU <- thisLU[,2]
  
  return(thisYear_obs.sf)
}

obs_habs_lu.sf <- lapply(lu_years,function(yr){getLUfromObs(obs_habs.df,year=yr)}) |> bind_rows() #est. 3 min


## ----mergeHabitats----------------------------------------------------------------------------------------------
print("Starting final modelling dataset...")
#Function to find column names matching a headline habitat code (level 1)
listMatchingHabs <- function(codeL1){
  allCols <- colnames(obs_habs_lu.sf)
  
  return(allCols[grepl(paste0("^",codeL1,"(_|$)"),allCols)])
}

#Reframe the dataset to merge the habitat categories
model.df <- obs_habs_lu.sf |>
  transmute(year,
            country = factor(countryCode),
            taxonKey,
            iucn_ID,
            HILDA_LU,
            
            #Habitats
            hab_1 = if_any(listMatchingHabs(1),~ !is.na(.x)),
            hab_2 = if_any(listMatchingHabs(2),~ !is.na(.x)),
            hab_3 = if_any(listMatchingHabs(3),~ !is.na(.x)),
            hab_4 = if_any(listMatchingHabs(4),~ !is.na(.x)),
            hab_5 = if_any(listMatchingHabs(5),~ !is.na(.x)),
            hab_6 = if_any(listMatchingHabs(6),~ !is.na(.x)),
            hab_7 = if_any(listMatchingHabs(7),~ !is.na(.x)),
            hab_8 = if_any(listMatchingHabs(8),~ !is.na(.x)),
            hab_9 = if_any(listMatchingHabs(9),~ !is.na(.x)),
            hab_10 = if_any(listMatchingHabs(10),~ !is.na(.x)),
            hab_11 = if_any(listMatchingHabs(11),~ !is.na(.x)),
            hab_12 = if_any(listMatchingHabs(12),~ !is.na(.x)),
            hab_13 = if_any(listMatchingHabs(13),~ !is.na(.x)),
            hab_14_arable = if_any(listMatchingHabs("14_(1)"),~ !is.na(.x)),
            hab_14_pasture = if_any(listMatchingHabs("14_(2)"),~ !is.na(.x)),
            hab_14_degraded = if_any(listMatchingHabs("14_(3|6)"),~ !is.na(.x)),
            hab_14_urban = if_any(listMatchingHabs("14_(4|5)"),~ !is.na(.x)),
            hab_15 = if_any(listMatchingHabs(15),~ !is.na(.x)),
            hab_16 = if_any(listMatchingHabs(16),~ !is.na(.x))
          ) |>
  data.frame() |>
  filter(HILDA_LU != "Ocean") |> #Remove points in the oceans/seas.
  filter(HILDA_LU != "No data") |> #Remove points with no land use data.
  droplevels() #Remove unneeded geographic levels


## ---------------------------------------------------------------------------------------------------------------
print("Removing duplicate observations...")

#Detecting duplicate observations
#model.df |> group_by(taxonKey,geometry,year) |> tally() |> group_by(n) |> tally()
#There definitely are some, so we need to delete them.

#Keeping year, lu, species, and habitat data, keeping unique records removes spatial duplicates:
model.df <- model.df |> unique.data.frame() # takes ~5 mins

model.df <- model.df |> select(-geometry) #We no longer need this now


## ---------------------------------------------------------------------------------------------------------------
print("Removing species with less than 10 observations...")

spCounts <- model.df |> group_by(taxonKey) |> tally()

lowNspecies <- spCounts$taxonKey[spCounts$n < 10]

model.df <- model.df |> filter(!taxonKey %in% lowNspecies)


## ---------------------------------------------------------------------------------------------------------------
print("Randomly subsampling species until no species makes up more than 2% of observations in its class...")

#Work out classes per species
classes <- gbif_data |> select(class,taxonKey) |> unique.data.frame() #takes a minute

#Join class data to modelling observations
model.df <- model.df |> left_join(classes,by="taxonKey")

subsampleByClass <- function(thisClass){
aves <- model.df |> filter(class == thisClass)

  spPercs <- aves |> 
  group_by(class,taxonKey) |> 
  tally() |>
  arrange(-n)
  
  S <- sum(spPercs$n)
  
  first50 <- spPercs$n[1:50]
  
  cSumFirst50 <- cumsum(first50)
  
  cSumDiff <- S - cSumFirst50
  
  x2 <- first50 - (cSumDiff / seq(49,0,-1))
  
  capIndex <- which(x2<0)[1]
  
  capFreq <- spPercs$n[capIndex]
  
  replacementFrequency <- floor((S - cSumFirst50[(capIndex-1)]) / (51 - capIndex))
  
  aves <- aves |>
    group_by(taxonKey) |>
    mutate(species_count = n()) |>
    mutate(toRemove = species_count - replacementFrequency) |>
    mutate(randomRank = runif(species_count,0,100000) |> rank()) |>
    filter(toRemove <= 0 | randomRank > toRemove) |>
    ungroup()
  
return(aves)
}

model.df <- lapply(unique(model.df$class),subsampleByClass) |> bind_rows()

#Finally, remove intermediate columns
model.df <- model.df |> select(-c(class,species_count,toRemove,randomRank))


## ----cleanEnvironment-------------------------------------------------------------------------------------------
#rm(list = setdiff(ls(), c("model.df","hilda")))
gc() #and tidy up memory


## ----testAssociationsTable--------------------------------------------------------------------------------------
table("LU forest" = model.df$HILDA_LU=="Forest", "Habitat forest" = model.df$hab_1)


## ----modelFunction----------------------------------------------------------------------------------------------
modelHabAssocs <- function(LandUse){

  #How many positive samples for "forest" land use class?
  lu_n_pos <- sum(model.df$HILDA_LU==LandUse)
  
  #Randomly subsample negatives to match the number of positives
  sampleRows_neg <- sample_n(
    model.df |> filter(HILDA_LU != LandUse),
    size = lu_n_pos
      )
  
  #Match this with all the positives
  sampleRows_pos <- model.df |> filter(HILDA_LU == LandUse)
  
  #Create the final modelling dataset
  lu.df <- rbind(sampleRows_pos, sampleRows_neg) |> 
    droplevels() |>
    mutate(idCol = row_number())
  
  #Variable to define if the land use matches the target
  lu.df$habitatMatch <- lu.df$HILDA_LU == LandUse
  
  mod.df <- sample_n(
    tbl = lu.df,
    size = floor(nrow(lu.df)*0.7)
    )
  
  eval.df <- lu.df |> filter(!idCol %in% mod.df$idCol)

  thisGlm <- glmer(
    formula = 
      habitatMatch ~ 
        hab_1 +
        hab_2 +
        hab_3 +
        hab_4 +
        hab_5 +
        hab_6 +
        #hab_7 +
        hab_8 +
        #hab_9 +
        #hab_10 +
        #hab_11 +
        #hab_12 +
        #hab_13 +
        hab_14_arable +
        hab_14_pasture +
        hab_14_degraded +
        hab_14_urban +
        hab_15 +
      (1|country),
    data = mod.df,
    family = "binomial"
  )
  
  return(
    list(
      LandUse = LandUse,
      mod.data = mod.df,
      eval.data = eval.df,
      model = thisGlm
      )
  )
}


## ----modelEvaluateAllParallel-----------------------------------------------------------------------------------
print("Starting model runs...")

lu_names <- unique(model.df$HILDA_LU)

modellers <- makeCluster(7)

clusterExport(modellers, "model.df")
clusterEvalQ(modellers,library(dplyr))
clusterEvalQ(modellers,library(lme4))

a <- parLapply(cl = modellers, X = lu_names, fun = modelHabAssocs)

stopCluster(modellers)

names(a) <- lu_names


## ----function_extractOddsRatios---------------------------------------------------------------------------------
getOddsRatios <- function(thisRegression){
  thisModel <- thisRegression$model
  
  theseOddsRatios <- tidy(thisModel,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
  
  
  oddsRatios_thisRow <- c(t(theseOddsRatios$estimate)) |> data.frame()
  rownames(oddsRatios_thisRow) <- theseOddsRatios$term
  colnames(oddsRatios_thisRow) <- thisRegression$LandUse

  return(oddsRatios_thisRow)
}


## ----extractOddsRatios------------------------------------------------------------------------------------------
obsFrame <- lapply(a,function(x){sum(model.df$HILDA_LU==x$LandUse)}) |> unlist() |> data.frame()
colnames(obsFrame) <- "Observations"


allOdds <- lapply(a,getOddsRatios) |> bind_cols() |>t() |> data.frame() |> select(-X.Intercept.)

#Add presentation names
colnames(allOdds) <- lapply(colnames(allOdds),function(x){hab_levels[hab_levels$ID==gsub("TRUE","",x),]$iucn_habitat}) |> unlist()

modelFigure <- cbind(obsFrame,allOdds) |> round(2) 

modelFigure |> kable()


## ----calculateTertiles------------------------------------------------------------------------------------------
posAssoc_tertiles <- modelFigure[,2:ncol(modelFigure)]
posAssoc_tertiles <- posAssoc_tertiles[posAssoc_tertiles>=1]
posAssoc_tertiles <- quantile(posAssoc_tertiles,probs = seq(0,1,1/3))

strongLimit <- as.numeric(posAssoc_tertiles[3])


## ----filterStrongAssocs_function--------------------------------------------------------------------------------
getStrongAssocs <- function(LandUse){
  theseCoefs <- allOdds[LandUse,]
  
  justCoefs <- as.numeric(theseCoefs)
  names(justCoefs) <- colnames(theseCoefs)
  
  return(names(justCoefs[justCoefs>=strongLimit]))
}

getStrongAssocs(LandUse = "Forest")


## ----finalHabMaps-----------------------------------------------------------------------------------------------
print("Generating habitat maps...")

createHabMap <- function(hab = "Forest",lu = hilda[[1]]){
  
  #Find all land uses associated with the habitat
  matches <- lapply(lapply(lu_names,getStrongAssocs),function(x){hab %in% x}) |> unlist()
  
  matching_lu <- lu_names[matches]
  habmap <- lu %in% matching_lu
  names(habmap) <- paste0(names(habmap),"_",hab)
  return(habmap)
}


mapHabs <- function(LU){
  return(lapply(hab_levels[,2],function(x){createHabMap(hab = x,lu = LU)}))
}

testYear <- mapHabs(hilda$X2019) |> rast()


## ---------------------------------------------------------------------------------------------------------------



## ---------------------------------------------------------------------------------------------------------------
knitr::purl("biofutures_life.qmd")

