library(here)
source(here("0-config.R"))
source(here("0-functions.R"))


#############################################################################################################
# AIM 1: Build PCA for each location, compare it to measure of multipathogen burden for each of 4 serological assays.
# WEALTH 
#############################################################################################################


#################################################################################
# Peru 
################################################################################
varlist_peru = c("ss.blender" , "ss.dresser", "ss.window", "ss.cell", "ss.floortype", "ss.cook",  "ss.handwash", 
                 "ss.water", "ss.soap", "ss.floor.f", "ss.cook.f", 
                 "ss.latrine.f", "td.drinksource.f", "td.defecate.f") # (-ss.hohjob.f) %>% # not directly related to WASH, -td.childfeces # UNCLEAR HOW THESE WERE CODED, -ss.latrine# UNCLEAR HOW THESE WERE CODED

peru_assets_environment = readr::read_csv(file = here("projects/6-multipathogen-burden/data/peru", "loretodata_swabs.csv")) %>%
  mutate(Sample = as.numeric(dbs.filtercode)) %>%
  dplyr::select(all_of(varlist_peru), district, community) 

####################################
# Check for variables with little variation 
check_dominant_vars <- function(peru_assets_environment, threshold = 0.9) {
  
  results <- data.frame( variable = character(),dominant_value = character(), prop = numeric() )
  
  for (v in names(peru_assets_environment)) {
    tab <- table(peru_assets_environment[[v]], useNA = "ifany")
    if (length(tab) > 0) {
      prop <- max(prop.table(tab, margin = NULL))
      if (prop >= threshold) {
        dom_val <- names(tab)[which.max(tab)]
        results <- rbind(
          results, data.frame(variable = v, dominant_value = dom_val, prop = round(prop, 3)))  }}}
  results }

dominant_vars <- check_dominant_vars(peru_assets_environment, threshold = 0.9)
dominant_vars_remove = unique(dominant_vars$variable)

####################################
# Check for lots of missingness 
check_na_props <- function(peru_assets_environment) {
  na_perc <- sapply(peru_assets_environment, function(x) mean(is.na(x)))
  data.frame(
    variable = names(na_perc),
    na_prop = round(na_perc, 3))}

# Run it
na_summary <- check_na_props(peru_assets_environment) %>%
  filter(na_prop > .5)
na_vars_remove = unique(na_summary$variable)

######################## Make cleaned dataset 

peru_assets_environment_final = peru_assets_environment %>%
  dplyr::select(-na_vars_remove, -dominant_vars_remove)

# Need to manually turn categorical variables into numerical 
sapply(peru_assets_environment_final, class)

# Check categorical variables 
table(peru_assets_environment_final$ss.floor.f)
table(peru_assets_environment_final$ss.latrine.f)
table(peru_assets_environment_final$td.drinksource.f)
table(peru_assets_environment_final$td.defecate.f)
table(peru_assets_environment_final$ss.floortype)

# We generally want higher numbers to correspond to better states of asset and environmental exposure.
peru_assets_environment_num = peru_assets_environment_final %>%
  mutate(floor_recode = ifelse(ss.floor.f == "dirt/sand", 1,
                               ifelse(ss.floor.f == "wood", 2, 
                                      ifelse(ss.floor.f == "cement/brick", 3, NA)))) %>%
  mutate(cook_recode = ifelse(ss.cook.f == "wood", 1,
                              ifelse(ss.cook.f == "natural gas", 2,
                                     ifelse(ss.cook.f == "liquid gas", 3, 
                                            ifelse(ss.cook.f == "elestricity", 4, NA))))) %>%
  mutate(latrine_recode = ifelse(ss.latrine.f == "no facilities/bush/field/surface water", 1,
                                 ifelse(ss.latrine.f == "pit latrine without slab/open pit" |
                                          ss.latrine.f == "pit latrine with slab", 2,
                                        ifelse(ss.latrine.f == "ventilation improved pit latrine", 3,
                                               ifelse(ss.latrine.f == "flush/pour to septic tank" |
                                                        ss.latrine.f == "flush/pour to open drains" | 
                                                        ss.latrine.f == "flush/pour to unkown place"|
                                                        ss.latrine.f == "flush/pour to pit latrine", 4, NA))))) %>%
  mutate(drink_recode = ifelse(td.drinksource.f == "surface water", 1,
                               ifelse(td.drinksource.f == "unprotected dug well" |
                                        td.drinksource.f == "unprotected spring", 2,
                                      ifelse(td.drinksource.f == "protected well" |
                                               td.drinksource.f == "protected spring", 3, 
                                             ifelse(td.drinksource.f == "public tap/stand pipe"|
                                                      td.drinksource.f == "piped to compound/plot", 4, NA))))) %>%
  mutate(defecate_recode = ifelse(td.defecate.f == "no structure/outside", 1,
                                  ifelse(td.defecate.f == "shared or public latrine", 2,
                                         ifelse(td.defecate.f == "private latrine", 3, NA )))) %>%
  dplyr::select(-ss.floor.f, -ss.latrine.f, -td.drinksource.f, -td.defecate.f, -ss.floortype , -ss.cook.f)
head(peru_assets_environment_num )


district = print(peru_assets_environment_num$district)
community = print(peru_assets_environment_num$community)
district_community_peru = data.frame(district, community)

# Make sure binary, where = NA, replae with mean values 
peru_clean = peru_assets_environment_num %>% 
  dplyr::select(-district, -community) %>%
  mutate(across(
    ss.blender:ss.handwash,           
    ~ ifelse(. %in% c(0, 1), ., NA)  )) %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ ifelse(. < 0 | . > 900, NA, .))) %>% 
  mutate(across(
    .cols = where(is.numeric),            
    .fns = ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
head(peru_clean)

# geo lebel 
# Make sure binary, where = NA, replae with mean values 
peru_geo_level = peru_assets_environment_num %>% 
  mutate(across(
    ss.blender:ss.handwash,           
    ~ ifelse(. %in% c(0, 1), ., NA)  )) %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ ifelse(. < 0 | . > 900, NA, .))) %>% 
  mutate(across(
    .cols = where(is.numeric),            
    .fns = ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  group_by(community) %>%
  summarize(
    across(
      c(ss.blender:ss.handwash, floor_recode:defecate_recode), 
      ~ mean(.x, na.rm = TRUE) ) ) %>%
  ungroup() %>%
  mutate(geo = community)
head(peru_geo_level)

# Convert data into matrix 
peru_matrix <-as.matrix(peru_clean)
pca_peru_matrix <- princomp(peru_matrix) 
summary(pca_peru_matrix)
screeplot(pca_peru_matrix)

## To get the first principal component in a variable ##
load <- loadings(pca_peru_matrix)[,1] 
print(load)

# Look at all loadings:
all_loadings <- loadings(pca_peru_matrix)
print(all_loadings[,1:3])

peru_pr.cp <- peru_matrix %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 
HHwealth_peru <- as.numeric(peru_pr.cp) ## Gives us the 1st PC in numeric form in pr.
summary(HHwealth_peru)

# Append wealth to district and community 
wealth_peru_place = district_community_peru %>%
  mutate(wealth = HHwealth_peru) %>%
  group_by(community) %>%
  summarize(across(wealth, ~ mean(.x, na.rm = TRUE))) %>%
  drop_na(community) %>% # %>%
 # left_join(method1, by = "community") %>%
  mutate(wealth = wealth - min(wealth)) %>%
  mutate(geo = community) %>% mutate(Location = "Peru") %>%
  dplyr::select(-community)
head(wealth_peru_place)


#################################################################################
# Bangladesh 
################################################################################
# Pearl's code
bangl_metadata <- read.csv(file = here("data/bangl/enrollment_wealth/untouched", "washb-bangladesh-enrol-public.csv"))
head(bangl_metadata)

df <- subset(bangl_metadata, select = -c(svyweek, svyyear))

# Vector of variables needed for the function asset_PCA
# removed "asset_phone" because not enough variation
# removed asset_tvbw and asset_tvcol, because they're redundant with asset_tv and cement which is the same with floor

# no WASH-related vars
varlist = c("landacre","roof","walls","floor",
            "elec","asset_radio","asset_refrig","asset_bike","asset_moto",
            "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_clock",
            "asset_khat","asset_chouki","asset_mobile")
varlist<-c("dataid","clusterid","hhid","block", varlist)

varlist_fac = c("roof","walls","floor",
                "elec","asset_radio","asset_refrig","asset_bike","asset_moto",
                "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_clock",
                "asset_khat","asset_chouki","asset_mobile")

stat_vars = c("landacre","walls","floor",
              "elec","asset_refrig","asset_bike","asset_moto",
              "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_khat",
              "asset_chouki","asset_mobile","wealthscore","wealth_tertile")

# Running the function assetPCA that computes the wealth score and divides them
# into tertiles. This function outputs a dataframe with the dataid variables +
# wealth scores and tertiles

assetPCA<-function(df, varlist, varlist_fac, stat_vars, reorder=F ){
  
  #Subset to only needed variables for subgroup analysis
  ret <- df %>%
    subset(select=c(varlist)) ## 22 vars including 4 ids
  
  #Select assets
  ret<-as.data.frame(ret) 
  id<-subset(df, select=c("dataid","clusterid","hhid","block")) ###only ID
  ret_assets_comp<-ret[,which(!(colnames(ret) %in% c("dataid","clusterid","hhid","block")))] ###only asset vars
  ret_assets_comp[,c(2:18)] <- lapply(ret_assets_comp[,c(2:18)], factor)
  
  print(head(ret_assets_comp))
  #Replace character blank with NA
  for(i in 2:ncol(ret_assets_comp)){
    ret_assets_comp[,i]<-ifelse(ret_assets_comp[,i]=="",NA,ret_assets_comp[,i])
  }
  
  #drop rows with no asset data
  id<-id[rowSums(is.na(ret_assets_comp[,5:ncol(ret_assets_comp)])) != ncol(ret_assets_comp)-4,]  
  ret_assets_comp<-ret_assets_comp[rowSums(is.na(ret_assets_comp[,5:ncol(ret_assets_comp)])) != ncol(ret_assets_comp)-4,]  
  
  
  #Drop assets with great missingness
  for(i in 1:ncol(ret_assets_comp)){
    cat(colnames(ret_assets_comp)[i],"\n")
    print(table(is.na(ret_assets_comp[,i])))
    print(class((ret_assets_comp[,i])))
  }
  
  cols.dont.want <- c("asset_clock")
  ret_assets_comp <- ret_assets_comp[, ! names(ret_assets_comp) %in% cols.dont.want, drop = F]
  
  #create level for missing factor levels
  table(is.na(ret_assets_comp))
  for(i in 2:ncol(ret_assets_comp)){
    ret_assets_comp[,i]<-as.character(ret_assets_comp[,i])
    ret_assets_comp[is.na(ret_assets_comp[,i]),i]<-"miss"
    ret_assets_comp[,i]<-as.factor(ret_assets_comp[,i])
    
  }
  
  ret_assets_comp$landacre[is.na(ret_assets_comp$landacre)] <- mean(ret_assets_comp$landacre, na.rm = T) ###repLace NAs with the mean
  table(is.na(ret_assets_comp))
  
  #Convert factors into indicators
  ret_assets_comp[,2:length(ret_assets_comp)]<-droplevels(ret_assets_comp[,2:length(ret_assets_comp)])

  Formula <- as.character(paste("~", "ret_assets_comp[,",2,"]", sep=''))
  for (i in 3:length(varlist_fac)) {
    Formula<-as.character(paste(Formula,"+","ret_assets_comp[,",i,"]",sep=""))
  }
  Formula<-as.formula(Formula)
  
  ###removing WASH-related vars
  ret_assets_compm <- data.frame(ret_assets_comp[, ! colnames(ret_assets_comp) %in% varlist_fac],
                                 model.matrix(Formula, ret_assets_comp))
  
  #Remove columns with almost no variance
  if(length(nearZeroVar(ret_assets_compm))>0){
    ret_assets_compm <-ret_assets_compm[,-nearZeroVar(ret_assets_compm)]
  } ## asset_radio and "roof" got removed
  
  ## Convert the data into matrix ##
  ret_assets_compm <-as.matrix(ret_assets_compm)
  
  for(y in 1:ncol(ret_assets_compm)) {
    for(x in 1:nrow(ret_assets_compm)) {
      if (is.na(ret_assets_compm[x,y])) ret_assets_compm[x,y] = 0
    }}
  
  ##Computing the principal component using eigenvalue decomposition ##
  princ.return <- princomp(ret_assets_compm) 
  print(summary(princ.return))
  screeplot(princ.return)
  
  ## To get the first principal component in a variable ##
  load <- loadings(princ.return)[,1] 
  
  pr.cp <- ret_assets_compm %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 
  HHwealth <- as.numeric(pr.cp) ## Gives us the 1st PC in numeric form in pr.
  
  #Create 3-level household weath index
  tertiles<-quantile(HHwealth, probs=seq(0, 1, 1/3))
  print(tertiles)
  ret_assets_compm<-as.data.frame(ret_assets_compm)
  ret_assets_compm$HHwealth<-HHwealth
  ret_assets_compm$HHwealth_ter<-rep(1, nrow(ret_assets_compm))
  ret_assets_compm$HHwealth_ter[HHwealth>=tertiles[1]]<-1
  ret_assets_compm$HHwealth_ter[HHwealth>=tertiles[2]]<-2
  ret_assets_compm$HHwealth_ter[HHwealth>=tertiles[3]]<-3
  table(ret_assets_compm$HHwealth_ter)
  
  if(reorder==T){
    ret_assets_compm$HHwealth_ter<-factor(ret_assets_compm$HHwealth_ter, levels=c("1", "2","3"))
  }else{
    levels(ret_assets_compm$HHwealth_ter)<-c("1", "2","3")
  }
  
  #Table assets by pca quintile to identify wealth/poverty levels
  d<-data.frame(id, ret_assets_compm)

  colnames(d) <- c("dataid","clusterid","hhid","block", stat_vars)
  
  #Save just the wealth data
  pca.wealth<-d %>% subset(select=c(dataid,clusterid,hhid,block,wealthscore,wealth_tertile))
  pca.wealth$dataid<-as.numeric(as.character(pca.wealth$dataid))
  
  d <-df %>% subset(., select=c("dataid","clusterid","hhid","block"))
  d$dataid<-as.numeric(as.character(d$dataid))
  d<-left_join(d, pca.wealth, by=c("dataid","clusterid","hhid","block"))
  
  return(as.data.frame(d))
  
}

d <- assetPCA(df, varlist = varlist, varlist_fac = varlist_fac, stat_vars = stat_vars, reorder = F)

wealth_bangl_place = d %>%
  group_by(block) %>%
  mutate(wealth = mean(wealthscore)) %>%
  distinct(block, wealth)  %>%
  ungroup() %>%
  mutate(block_r = 91 - block) %>%
  mutate(geo = block_r) %>% mutate(Location = "Bangladesh") %>% dplyr::select(-block_r, -block) 
head(wealth_bangl_place)
# geo is block r 

#################################################################################
# Kenya
################################################################################

kenya_meta = read_dta(file = here("data/kenya/master_files", "washk_household_vars_20190305.dta")) %>%
  mutate(improved_floor = ifelse(floor == 1, 1, 
                                 ifelse(floor == 9, NA, 0))) %>% # 1 = cement, 0 = earth/dung
  mutate(improved_food = ifelse(HHS == 1, 1, 
                                ifelse(HHS == 9, NA, 0))) %>% # 1 = little to no, 2 = moderate, 3 = severe
  mutate(improved_water = ifelse(prim_drink_ws_bl == 1 & water_time <= 30, 1, 0)) %>% # improved source and walking time < 30 min
  mutate(improved_hygiene = ifelse(hws_bl == 1, 1, 0)) %>% # soap and water present
  mutate(improved_sanitation = ifelse(imp_lat_bl == 1, 1,0 )) 
print(names(kenya_meta))

varlist_assets_kenya = c("roof", "walls", "elec", "radio", "tv", "mobilephone", "clock", "bicycle",
                         "motorcycle", "stove", "cooker", "car", "cow", "goat", "poultry")
varlist_wash_kenya = c( "improved_floor", "improved_food", "improved_water", "improved_hygiene", "improved_sanitation")

kenya_assets_environment = kenya_meta %>%
  dplyr::select(all_of(varlist_assets_kenya), all_of(varlist_wash_kenya), clusterid) 
  
# Check dominant variables 
check_dominant_vars <- function(kenya_assets_environment, threshold = 0.9) {
  
  results <- data.frame( variable = character(),dominant_value = character(), prop = numeric() )
  
  for (v in names(kenya_assets_environment)) {
    tab <- table(kenya_assets_environment[[v]], useNA = "ifany")
    if (length(tab) > 0) {
      prop <- max(prop.table(tab, margin = NULL))
      if (prop >= threshold) {
        dom_val <- names(tab)[which.max(tab)]
        results <- rbind(
          results, data.frame(variable = v, dominant_value = dom_val, prop = round(prop, 3)))  }}}
  results }

dominant_vars <- check_dominant_vars(kenya_assets_environment, threshold = 0.9)
dominant_vars_remove = unique(dominant_vars$variable)

# Check for lots of missingness 
check_na_props <- function(kenya_assets_environment) {
  na_perc <- sapply(kenya_assets_environment, function(x) mean(is.na(x)))
  data.frame(
    variable = names(na_perc),
    na_prop = round(na_perc, 3))}

# Run it
na_summary <- check_na_props(kenya_assets_environment) %>%
  filter(na_prop > .5)
na_vars_remove = unique(na_summary$variable)

######################## Make cleaned dataset 

kenya_assets_environment_final = kenya_assets_environment %>%
  dplyr::select(-na_vars_remove, -dominant_vars_remove)

# Need to manually turn categorical variables into numerical 
sapply(kenya_assets_environment_final, class)

# Check categorical variables 
table(kenya_assets_environment_final$roof)

# Location
clusterid = print(kenya_assets_environment_final$clusterid)

head(kenya_assets_environment_final)

# Make sure binary, where = NA, replae with mean values 
kenya_clean = kenya_assets_environment_final %>% 
  mutate(across(everything(), zap_labels)) %>%
  mutate(across(where(~ !is.numeric(.)), as.numeric)) %>%
  dplyr::select(-clusterid) %>%
  mutate(across(
    roof:stove,           
    ~ ifelse(. %in% c(0, 1), ., NA)  )) %>%
  mutate(across(
    cow:goat,           
    ~ ifelse(. > 1, 1, 0)  )) %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ ifelse(. > 8, NA, .))) %>% 
  mutate(across(
    .cols = where(is.numeric),            
    .fns = ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
head(kenya_clean)

# Convert data into matrix 
kenya_matrix <-as.matrix(kenya_clean)
pca_kenya_matrix <- princomp(kenya_matrix) 
summary(pca_kenya_matrix)
screeplot(pca_kenya_matrix)

## To get the first principal component in a variable ##
load <- loadings(pca_kenya_matrix)[,1] 
print(load)

# Look at all loadings:
all_loadings <- loadings(pca_kenya_matrix)
print(all_loadings[,1:3])

kenya_pr.cp <- kenya_matrix %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 
HHwealth_kenya <- as.numeric(kenya_pr.cp) ## Gives us the 1st PC in numeric form in pr.
summary(HHwealth_kenya)

# block assignment 
treatment_assignment_k = read.csv(file = here("data/kenya/public_ids", 
                                           "cluster_tx_masked.csv")) 
head(treatment_assignment_k)

wealth_kenya_place = data.frame(clusterid, HHwealth_kenya) %>%
  left_join(treatment_assignment_k, by = "clusterid") %>%
  group_by(block) %>%
  summarize(across(HHwealth_kenya, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(wealth = HHwealth_kenya - min(HHwealth_kenya)) %>% dplyr::select(-HHwealth_kenya ) %>%
  mutate(geo = block) %>% mutate(Location = "Kenya") %>% dplyr::select(-block)
head(wealth_kenya_place)

########## kenya geo level wash conditions 
household_kenya = read_dta(file = here("data/kenya/master_files", "washk_household_vars_20190305.dta")) %>%
  mutate(improved_floor = ifelse(floor == 1, 1, 
                                 ifelse(floor == 9, NA, 0))) %>% # 1 = cement, 0 = earth/dung
  mutate(improved_food = ifelse(HHS == 1, 1, 
                                ifelse(HHS == 9, NA, 0))) %>% # 1 = little to no, 2 = moderate, 3 = severe
  mutate(improved_water = ifelse(prim_drink_ws_bl == 1 & water_time <= 30, 1, 0)) %>% # improved source and walking time < 30 min
  mutate(improved_hygiene = ifelse(hws_bl == 1, 1, 0)) %>% # soap and water present
  mutate(improved_sanitation = ifelse(imp_lat_bl == 1, 1,0 )) %>% # improved latrine 
  mutate(improved_all = ifelse(improved_water == 1 & improved_hygiene == 1 &improved_sanitation == 1, 1, 0)) %>%
  dplyr::select(clusterid, improved_water, improved_floor, improved_hygiene, improved_sanitation, assetindex, improved_all) %>%
  left_join(treatment_assignment, by = "clusterid") %>%
  group_by(block) %>%
  summarize(across(improved_water:improved_all, ~ mean(., na.rm = TRUE))) %>%
  mutate(geo = block)
head(household_kenya )

################################################################################
# Cambodia 
################################################################################
# Wealth index has alredy been made 
# Load the .dta file
cam_wealth_data <- read_dta(file = here("projects/6-multipathogen-burden/data/cambodia", "KHPR61FL.DTA"))
head(cam_wealth_data)

cam_wealth_dat = cam_wealth_data %>%
  dplyr::select(hv001, hv270, hv271) %>%
  mutate(dhsclust = hv001) %>%
  group_by(dhsclust) %>%
  summarize(across(hv270:hv271, ~ mean(.x)))
table(cam_wealth_dat$dhsclust)

gps_cambodia =   readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_ea_dhs.csv")) 
head(gps_cambodia)

load(file = here("data/cambodia/cambodia_serology.Rdata"))
gps_dat = cambodia_serology %>%
  dplyr::select(dhsclust, psuid) %>% distinct()
head(gps_dat)

cam_gps_wealth = left_join(gps_cambodia, cam_wealth_dat ,  by = "dhsclust") %>%
  left_join(gps_dat, by = "dhsclust") %>%
  mutate(psuid = psuid + 1) %>%
  mutate(wealth = hv271) %>%
  dplyr::select(psuid, wealth) %>%
  mutate(geo = psuid) %>% mutate(Location = "Cambodia") %>% dplyr::select(-psuid)
head(cam_gps_wealth)

table(cam_gps_wealth$geo)
table(cambodia_multiplex$geo)

cam_wealth_place = cam_gps_wealth %>% 
  mutate(min = -min(wealth)) %>%
  mutate(wealth = wealth + min) %>% dplyr::select(-min)
head(cam_wealth_place)

##############################################################################
# Serology 
##############################################################################

##################### 
# PERU 
peru_serology = readr::read_csv(file = here("projects/6-multipathogen-burden/data/peru", "perusero.csv")) %>%
  mutate(Sample = as.numeric(Sample)) %>%
  distinct(Sample, .keep_all = TRUE)
length(unique(peru_serology$Sample))
head(peru_serology)

peru_meta = readr::read_csv(file = here("projects/6-multipathogen-burden/data/peru", "loretodata_swabs.csv")) %>%
  mutate(Sample = as.numeric(dbs.filtercode)) %>%
  dplyr::select(gps_long, gps_lat, Sample) %>%
  distinct() %>% drop_na(Sample)
length(unique(peru_meta$Sample))

# process 
peru_serology_antigens = left_join(peru_serology, peru_meta, by = c("Sample")) %>%
  pivot_longer(cols = pgp3:rbd591_pos, names_to = "antigen", values_to = "value") %>%
  filter(!str_detect(antigen, "_pos$")) %>%
  filter(!antigen %in% c("spike", "rbd591", "nucleo", "rbd541")) %>% # remove covid 
  filter(!antigen %in% c("Vero      (78)", "GST     (15)" )) %>% # remove controls 
  drop_na(Sample)
head(peru_serology_antigens)

antigen = unique(peru_serology_antigens$antigen)
print(antigen)
cutoff = c(212, 108, 235, 86, 830, 132, 99, 148, 
           339, 531, 62, 32, 2, 1759, 428, 20, 106,
           375, 107, 58, 422, 306, 155, 18, 140, 19, 14, 
           27, 102, 255, 172)
pathogen = c("Chlamydia trachomatis", "Chlamydia trachomatis", "Treponema palladium", "Treponema palladium", 
         "Cryptosporidium parvum", "Cryptosporidium parvum", "Giardia lamblia","Giardia lamblia", 
         "P. falciparum", "P. falciparum", "P. falciparum","P. falciparum","P. falciparum", "P. falciparum", 
         "P. malariae", "P. malariae", "P. vivax", "Wucheria bancrofti","Wucheria bancrofti", "Wucheria bancrofti", 
         "Onchocerca volvulus", "Onchocerca volvulus", "Strongyloides stercoralis", "T canis", "T gondii", "Taenia solium", "Taenia solium", 
         "Tetanus toxoid", "Diptheria toxoid", "wMEV", "RuV")
antigens_cutoff = data.frame(antigen, cutoff, pathogen)
antigen_path = antigen[1:27]
antigen_vax = antigen[28:31]

peru_serology_processed = left_join(peru_serology_antigens, antigens_cutoff, by = "antigen") %>%
  mutate(seropos = ifelse(value > cutoff & antigen %in% antigen_path , 1, 
                          ifelse(value < cutoff & antigen %in% antigen_vax, 1, 0))) %>%
  dplyr::select(Sample, community, antigen, pathogen, gps_lat, gps_long, seropos) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>% # only want antigens where at least 5% of the community is infected 
  filter(pop_seropos < .9499) %>%
   group_by(Sample, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  dplyr::distinct(Sample, community, pathogen, gps_lat, gps_long, path_pos) %>%
  drop_na()
head(peru_serology_processed)
table(peru_serology_processed$community)

# Community level seroprevalence
peru_multiplex = peru_serology_processed  %>%
  group_by(pathogen, community) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), 
         num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = community) %>%
  drop_na() %>%
  mutate(Location = "Peru") %>%
  distinct(community, pathogen, gps_lat, gps_long, fraction, geo, Location)
head(peru_multiplex)
 
##################### 
# Kenya 
treatment_assignment = read.csv(file = here("data/kenya/public_ids", 
                                            "cluster_tx_masked.csv")) 
head(treatment_assignment)
table(treatment_assignment$tr_masked)

treatment_assignment = read.csv(file = here("data/kenya/primary_outcomes", "endline-anthro.csv")) %>%
  distinct(block, clusterid, tr)
head(treatment_assignment)

antigen_vax = c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path = c("T. solium", "Cholera", "E. histolytica" , "Cryptosporidium", "P. falciparum", 
                 "Schistosomiasis", "P. ovale",  "Norovirus",  "P. malariae" , "Onchocerciasis" , "Dengue",         
                 "P. vivax" , "Campylobacter", "Zika", "Salmonella" , "Trachoma" ,"Giardia"  ,      
                 "Chikungunya" ,    "LT-ETEC" ,   "Shigella", "Strongyloides" )
  
luminex_kenya = read.csv(file = here("data/kenya/luminex/final", 
                                     "washb_kenya_luminex_igg_seropos_2025-09-21.csv")) %>%
  mutate(dataid = str_extract(childid, "(?<=-)\\d{5}(?=-)")) %>%
  left_join(treatment_assignment, by = "clusterid") %>%
  #filter(tr == c("Control", "Nutrition")) %>%
  filter(visit == 3) %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path , 1, 
                   ifelse(mfi < mficut & pathogen %in% antigen_vax, 1, 0))) %>%
  dplyr::select(clusterid, childid, antigen, pathogen, seropos, block) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>%
  filter(pop_seropos < .9499) %>%
  group_by(childid, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(clusterid, childid, block, pathogen, path_pos)
head(luminex_kenya)
print(length(unique(luminex_kenya$childid)))
print(length(unique(luminex_kenya$antigen)))
  
# Community level seroprevalence
kenya_multiplex = luminex_kenya  %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), 
         num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = block) %>%
  drop_na() %>%
  mutate(Location = "Kenya")
head(kenya_multiplex)

##################### 

# Bangladesh
public_ids = read.csv(file = here("data/bangl/public_ids", "public-ids.csv")) %>%
  distinct(dataid, clusterid, block, clusterid_r, block_r)
head(public_ids)

gps_dat = read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>% # had to add later? 
  left_join(public_ids, by = "dataid") %>% 
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat) 
head(gps_dat)

public_id_cluster = public_ids %>%
  distinct(clusterid, block, block_r, clusterid_r) %>%
  left_join(gps_dat, by = c("block", "block_r"))
head(public_id_cluster)

antigen_vax = c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path = c("T. solium", "Cholera", "E. histolytica" , "Cryptosporidium", "P. falciparum", 
                 "Schistosomiasis", "P. ovale",  "Norovirus",  "P. malariae" , "Onchocerciasis" , "Dengue",         
                 "P. vivax" , "Campylobacter", "Zika", "Salmonella" , "Trachoma" ,"Giardia"  ,      
                 "Chikungunya" ,    "LT-ETEC" ,   "Shigella", "Strongyloides" )

luminex_bangl <- read.csv(file = here("data/bangl/luminex/final", 
                                      "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = c("clusterid")) %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path , 1, 
                          ifelse(mfi < mficut & pathogen %in% antigen_vax, 1, 0))) %>%
 # filter(visit == 3) %>%
  dplyr::select(mfi, agemonth, clusterid, childid, antigen, pathogen, seropos, block, block_r, med_qgpslong, med_qgpslat) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>%
  filter(pop_seropos < .9499) %>%
  group_by(childid, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(mfi, agemonth, clusterid, childid, block, pathogen, path_pos, block_r, med_qgpslong, med_qgpslat)
head(luminex_bangl)
print(length(unique(luminex_bangl$childid)))
print(length(unique(luminex_bangl$antigen)))


# plot age and calendar time 
agetime_bangl <- read.csv(file = here("data/bangl/luminex/final", 
                                      "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = c("clusterid")) %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  distinct(childid, visit, agemonth, sampledate, eed) %>%
  mutate(eed = replace(eed, eed == "Endline sample", "Full study"))
head(agetime_bangl)
table(agetime_bangl$eed)

age = ggplot(data = agetime_bangl) +
  geom_histogram(aes(x = agemonth, fill = eed), alpha = .7,  position = "identity") +
  ylab("Number of samples") + xlab("Age in months") + theme_minimal() +
  scale_fill_manual(values = c( "#1F78B4", "#B2A07A"), name = "") +
#  theme(
#    strip.text = element_text(size = 15),
#    plot.title = element_text(size = 14),
#    axis.text.x = element_text(size = 12),
#    axis.text.y = element_text(size = 12),
#    axis.title.y = element_text(size = 12),
#    axis.title.x = element_text(size = 12)) +
  ggtitle("C") +
  theme(
    plot.title      = element_text(size = 20, face = "bold"),   # title
    plot.subtitle =   element_text(size = 20), 
    axis.title.x    = element_text(size = 20),   # x-axis title
    axis.title.y    = element_text(size = 20),   # y-axis title
    axis.text.x     = element_text(size = 14),   # x-axis tick labels
    axis.text.y     = element_text(size = 14) ,   # y-axis tick labels
    plot.tag         = element_text(face = "bold", size = 20),
    legend.position = "none",
    legend.text = element_text(size = 18), legend.title = element_text(size = 20))
age


time_plot = ggplot(data = agetime_bangl) +
  geom_histogram(aes(x = as.Date(sampledate), fill = eed), alpha = .7, position = "identity") +
  ylab("Number of samples") + xlab("Year") + theme_minimal() +
  scale_fill_manual(values = c( "#1F78B4", "#B2A07A"), name = "") +
  #theme(
  #  strip.text = element_text(size = 15),
  #  plot.title = element_text(size = 14),
  #  axis.text.x = element_text(size = 12),
  #  axis.text.y = element_text(size = 12),
   # axis.title.y = element_text(size = 12),
  #  axis.title.x = element_text(size = 12)) +
  ggtitle("D") +
  theme(
    plot.title      = element_text(size = 20, face = "bold"),   # title
    plot.subtitle =   element_text(size = 20), 
    axis.title.x    = element_text(size = 20),   # x-axis title
    axis.title.y    = element_text(size = 20),   # y-axis title
    axis.text.x     = element_text(size = 14),   # x-axis tick labels
    axis.text.y     = element_text(size = 14) ,   # y-axis tick labels
    plot.tag         = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    legend.text = element_text(size = 18), legend.title = element_text(size = 20))
time_plot

candd = plot_grid(age, time_plot, nrow = 2, rel_heights = c(.4, .5))
candd

library(lubridate)
month_plot = ggplot(
  data = agetime_bangl %>% 
    mutate(
      sampledate = as.Date(sampledate),
      month = month(sampledate)
    )) +
  geom_histogram(aes(x = month, fill = eed), alpha = .8, position = "dodge") +
  ylab("Density") + xlab("Month") + theme_minimal() +
  scale_fill_manual(values = c("darkslategray3", "tan")) +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12))
month_plot

plot_grid(age, time_plot, ncol = 1, labels = "C")


# Community level seroprevalence
bangl_multiplex = luminex_bangl  %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = block) %>%
  drop_na()  %>%
  mutate(Location = "Bangladesh")
head(bangl_multiplex)

lat_lon_bangl = bangl_multiplex  %>%
  distinct(med_qgpslong, med_qgpslat)
write.csv(lat_lon_bangl, "lat_lon_bangl.csv")

# Cambodia 

# Load disease data. 
cambodia_serology_public = readr::read_csv(file = here("data/cambodia", "cambodia_serology_public.csv")) %>%
  mutate(womanid = ...1)
head(cambodia_serology_public)
#ttmb     : Tetanous toxoid
#bm14     : Lymphatic filariasis bm14
#bm33     : Lymphatic filariasis bm33
#wb123    : Lymphatic filariasis wb123
#nie      : Strongyloides stercoralis NIE
#sag2a    : Toxoplasma gondii SAG2A
#t24      : Taenia solium T24
#pfmsp19  : Plasmodium falciparum MSP-1(19)
#pvmsp19  : Plasmodium vivax MSP-1(19)
# Load location data 
gps_cambodia =   readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_ea_dhs.csv")) 
head(gps_cambodia)
#load(file = here("projects/6-multipathogen-burden/data/cambodia/cambodia_serology.Rdata"))
gps_dat = gps_cambodia %>%
  mutate(psuid = ...1) %>%
  distinct(psuid, dhslat, dhslon)
head(gps_dat)


#######################################################
# Transform into seropositive vs seronegative using: https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004699#pntd-0004699-t001
# https://journals.asm.org/doi/10.1128/cvi.00052-16#T1 for tetanus

antigen = c("ttmb", "bm14", "wb123", "bm33", "nie", "sag2a", "t24", "pfmsp19", "pvmsp19")
cutoff = c(100, 65, 115, 966,792, 159, 486, 343, 196 )
pathogen = c("Tetanus", "Lymphatic filariasis", "Lymphatic filariasis", "Lymphatic filariasis", 
             "Strongyloides stercoralis", "Toxoplasma gondii", "Taenia solium", "Plasmodium falciparum",
             "Plasmodium vivax")
antigen_path_dat = data.frame(antigen, cutoff, pathogen)

antigen_vax = c( "Tetanus")
antigen_path = c("Lymphatic filariasis", "Strongyloides stercoralis", "Toxoplasma gondii", "Taenia solium", "Plasmodium falciparum",
                 "Plasmodium vivax")

cambodia_seropositivity = cambodia_serology_public %>%
  pivot_longer(cols = ttmb:gst, names_to = "antigen", values_to = "value") %>%
  filter(antigen != "gst") %>%
  filter(antigen != "bm33") %>% # we dont want this in our defintion
  mutate(psuid = psuid + 1) %>% 
  left_join(antigen_path_dat, by = "antigen") %>%
  mutate(seropos = ifelse(value > cutoff & pathogen %in% antigen_path , 1, 
                          ifelse(value < cutoff & pathogen %in% antigen_vax, 1, 0))) %>%
  dplyr::select(womanid, psuid, age, parity, antigen, pathogen, seropos) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>%
  filter(pop_seropos < .9499) %>%
  group_by(womanid, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  distinct(womanid, psuid, pathogen, path_pos)
head(cambodia_seropositivity)
table(cambodia_seropositivity$pathogen)
  
cambodia_multiplex = cambodia_seropositivity  %>%
  group_by(pathogen, psuid) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), 
         num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = psuid) %>%
  drop_na() %>%
  distinct(psuid, pathogen, num, non_na_count, fraction, geo) %>%
  left_join(gps_dat, by = "psuid") %>%
  mutate(Location = "Cambodia")
head(cambodia_multiplex)


############ geo analyiss
healthcare_dist = cambodia_multiplex %>%
  distinct(dhslat, dhslon, psuid)
write.csv(healthcare_dist, "healthcare_dist_cambodia.csv")


cam_distances = readr::read_csv(file = here("data/cambodia", "healthcare_dist_cambodia_with_distances.csv")) 
head(cam_distances)

#############################################################
#############################################################
# Calculate Rao across all of the pathogens 
peru_multiplex_dat = peru_multiplex %>% dplyr::select(pathogen, geo, Location, fraction) %>%
  distinct()
bangl_multiplex_dat = bangl_multiplex %>% dplyr::select(pathogen, geo, Location, fraction) %>%
  distinct()
kenya_multiplex_dat = kenya_multiplex %>% dplyr::select(pathogen, geo, Location, fraction) %>%
  distinct()
cambodia_multiplex_dat = cambodia_multiplex %>% dplyr::select(pathogen, geo, Location, fraction) %>%
  distinct()

multiplex_data = rbind(peru_multiplex_dat, bangl_multiplex_dat, kenya_multiplex_dat, cambodia_multiplex_dat )
head(multiplex_data)

product_results = list()
locations = unique(multiplex_data$Location)

for (l in seq_along(locations)) {
  
  i_location <- locations[l]
  print(paste("Processing location:", i_location))
  
  # Subset data for this location
  location_subset <- multiplex_data %>%
    filter(Location == i_location)
  
  geo_n <- unique(location_subset$geo)
  
  for (i in seq_along(geo_n)) {
    
    i_geo <- geo_n[i]
    print(paste("  Processing geo:", i_geo))
    
    geo_subset <- location_subset %>%
      filter(geo == i_geo) %>%
      ungroup() %>%
      dplyr::select(fraction, geo, Location)
    
    # Get unique combinations of pathogen prevalences 
    fractions <- geo_subset$fraction
    unique_combinations <- t(combn(fractions, 2, simplify = TRUE))
    
    # Compute products for each pair
    products <- apply(unique_combinations, 1, prod)
    
    products_df <- data.frame(
      Location = i_location,
      geo = i_geo, 
      fraction_mirrored = unique_combinations[, 1],
      fraction = unique_combinations[, 2],
      product = products )
    
    # Use combined key for list storage
    product_results[[paste(i_location, i_geo, sep = "_")]] <- products_df
  }
}

# Combine all results into a single data frame
product_results <- do.call(rbind, product_results)
head(product_results)

# Sum by geo and location
rao_location_final = product_results %>%
  group_by(Location, geo) %>%
  summarize(rao = sum(product), .groups = "drop") %>%
  group_by(Location) %>%
  mutate(max_rao = max(rao)) %>%
  ungroup() %>%
  mutate(rao = rao/max_rao) %>%
  dplyr::select(-max_rao) 
head(rao_location_final)

ggplot(data =  rao_location_final) +
  geom_histogram(aes(x = rao, fill = Location))+
  theme_minimal() +
  facet_wrap(vars(Location))


# cambodia helathcate seeking x multipathogen burden 
cam_distances = readr::read_csv(file = here("data/cambodia", "healthcare_dist_cambodia_with_distances.csv")) %>%
  mutate(Location = "Cambodia") %>%
  mutate(geo = psuid) %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("geo", "Location"))
head(cam_distances)

ggplot(data = cam_distances, aes(x = distance_km, y = rao)) +
  geom_point(alpha = .7, cex = 1.3) +
  geom_smooth(aes(x = distance_km, y = rao), method = "loess", se = TRUE, color = "darkslategray4", linetype = "solid") +
  theme_bw() +
  ylab("Multipathogen burden (Rao)") +
  xlab("Distance (km)") +
  theme(
    strip.text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    # axis.title.y = if(show_y) element_text() else element_blank(),
    legend.position = "none" ) 


#########################################################
#########################################################


wealth_indices = rbind(wealth_kenya_place, wealth_bangl_place, wealth_peru_place, cam_wealth_place) %>%
  group_by(Location) %>%
  mutate(wealth = (wealth - min(wealth)) / (max(wealth) - min(wealth))) %>%
  ungroup()

wealth_rao = left_join(rao_location_final, wealth_indices, by = c("geo", "Location"))
head(wealth_rao)

ggplot(data = wealth_rao) +
  geom_point(aes(x = wealth, y = rao), col = "darkslategray4") +
  geom_smooth(aes(x = wealth, y = rao), method = "lm", se = TRUE, color = "black", linetype = "solid") +
  facet_wrap(vars(Location), scales = "free")+ theme_minimal() +
  ylab("Multipathogen burden") + xlab("Wealth") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
     axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    # axis.title.y = if(show_y) element_text() else element_blank(),
    legend.position = "none" ) + # +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, .2, .4, .6, .8,1 ))


# check correlation 
locations = unique(wealth_rao$Location)
cor_results = list()

for (l in seq_along(locations)) {
  
  i_location <- locations[l]
  
  location_subset <- wealth_rao %>%
    filter(Location == i_location)
  
  cor_test <- cor.test(location_subset$rao, location_subset$wealth, 
                       method = "spearman")
  cor_results[[l]] <- data.frame(
    Location    = i_location,
    correlation = cor_test$estimate,
    p_value     = cor_test$p.value)
}

cor_results_final <- do.call(rbind, cor_results)
rownames(cor_results_final) <- NULL
head(cor_results_final)


###############################################
# try to predict continuous wealth 
library(xgboost)
library(caret)

# bangl 
head(wealth_bangl_place)
head(bangl_multiplex_dat)

# Bangl
predict_wealth_bangl = bangl_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(wealth_bangl_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final, by = c("geo", "Location"))
head(predict_wealth_bangl)

# Cambodia 
predict_wealth_cam = cambodia_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(cam_wealth_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("geo", "Location")) %>%
  mutate(wealth = wealth + 44001 ) %>% # add a bit so does not go to infintiy
  mutate(wealth = log(wealth))
head(predict_wealth_cam)
summary(predict_wealth_cam$wealth)

# Kenya 
predict_wealth_kenya = kenya_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(wealth_kenya_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("geo", "Location")) %>%
  drop_na()
head(predict_wealth_kenya)

# Peru
predict_wealth_peru = peru_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(wealth_peru_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final , by = c("geo", "Location")) %>%
  drop_na()
head(predict_wealth_peru)
head(predict_wealth_kenya)


set.seed(123)

# Predictor matrix
predictor_cols <- names(predict_wealth_cam)[c(3:8, 10)] # this is cam
predictor_cols <- names(predict_wealth_bangl)[c(3:12, 14)] # this is bangl
predictor_cols <- names(predict_wealth_kenya)[c(3:19, 21)] # this is kenya 
predictor_cols <- names(predict_wealth_peru)[c(3:16, 18)] # this is peru 
print(predictor_cols)

x_all <- as.matrix(predict_wealth_peru[, predictor_cols])

# Continuous outcome
y_all <- predict_wealth_peru$wealth

# Remove missing
complete_cases <- complete.cases(x_all, y_all)
x_all <- x_all[complete_cases, ]
y_all <- y_all[complete_cases]

cat(sprintf("Sample size: %d\n", length(y_all)))
cat(sprintf("Predictors: %d\n", ncol(x_all)))


####################################
#####################################
##### cross validate and out of sample 
set.seed(123)

# create randmomly sampled folds 
k <- 10
fold_ids <- sample(rep(1:k, length.out = nrow(x_all)))

# Standardize outcome value for all values - should be done at the fold level? 
y_mean <- mean(y_all)
y_sd   <- sd(y_all)
y_scaled <- (y_all - y_mean) / y_sd

# stores out-of-fold predictions
cv_preds   <- numeric(length(y_all)) 
cv_r2      <- numeric(k)
cv_rmse    <- numeric(k)
cv_mae     <- numeric(k)

# training on 90% of data, testing on 10% for each fold 
for (fold in 1:k) {
  
  fold_test_idx  <- which(fold_ids == fold)
  fold_train_idx <- which(fold_ids != fold)
  
  # Scale on training only
  y_mean_fold <- mean(y_all[fold_train_idx])
  y_sd_fold   <- sd(y_all[fold_train_idx])
  
  y_train_scaled <- (y_all[fold_train_idx] - y_mean_fold) / y_sd_fold
  y_test_scaled  <- (y_all[fold_test_idx]  - y_mean_fold) / y_sd_fold
  
  dtrain <- xgb.DMatrix(data = x_all[fold_train_idx, ], label = y_train_scaled)
  dval   <- xgb.DMatrix(data = x_all[fold_test_idx, ],  label = y_test_scaled)
  
  # best params 
  model <- xgb.train(
    params    = best_params, # hyperparmaters which are tuned within folds currently 
    data      = dtrain,
    nrounds   = 500, # max number of trees 
    early_stopping_rounds = 20, # stios training if validation loss hasnt impoved in 20 rounds (stops overfitting)
    watchlist = list(train = dtrain, eval = dval),
    verbose   = 0
  )
  
  # Predict outcoe fo wealth and rescale back to original units
  preds_scaled <- predict(model, dval)
  preds        <- preds_scaled * y_sd_fold + y_mean_fold
  actuals      <- y_all[fold_test_idx]
  
  # Store out-of-fold predictions
  cv_preds[fold_test_idx] <- preds
  
  # Per-fold metrics
  cv_rmse[fold] <- sqrt(mean((actuals - preds)^2))
  cv_mae[fold]  <- mean(abs(actuals - preds))
  cv_r2[fold]   <- cor(actuals, preds)^2
  
  cat(sprintf("Fold %d | R2: %.3f | RMSE: %.3f\n", fold, cv_r2[fold], cv_rmse[fold]))
}

# ── Overall metrics ───────────────────────────────────────────────────────────
cat("\n=== 5-Fold CV Performance ===\n")
cat(sprintf("R2:   %.3f (SD: %.3f)\n", mean(cv_r2),   sd(cv_r2))) # computed as squared Pearson correlation between observed and predicted, represents variance explained
cat(sprintf("RMSE: %.3f (SD: %.3f)\n", mean(cv_rmse), sd(cv_rmse)))
cat(sprintf("MAE:  %.3f (SD: %.3f)\n", mean(cv_mae),  sd(cv_mae)))

# ── Plot all out-of-fold predictions ─────────────────────────────────────────
plot_df <- data.frame(
  Observed  = y_all,
  Predicted = cv_preds,
  Fold      = as.factor(fold_ids))

peru_pred_plot = ggplot(plot_df, aes(x = Observed, y = Predicted)) +
  geom_point(aes(color = Fold), alpha = 0.9, cex = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
 # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
 # labs(
  #  title    = "Cambodia",
  #  subtitle = sprintf("R^2 = %.3f ± %.3f",
  #                     mean(cv_r2), sd(cv_r2)),
  #  x = "Observed Wealth",
  #  y = "Predicted Wealth") +
  labs(
    title    = "Peru",
    subtitle = sprintf("R^2 = %.3f",
                       mean(cv_r2)),
    x = "Observed Wealth",
    y = "Predicted Wealth" ) +
  theme_minimal() +
  scale_color_viridis_d(option = "D")
  

kenya_pred_plot
bangl_pred_plot
cam_pred_plot

plot_grid(bangl_pred_plot,  cam_pred_plot, kenya_pred_plot)

  ## What to report in your paper
#"Model performance was evaluated using 5-fold cross-validation. 
#The mean R² across folds was 0.XXX (SD = 0.XXX), RMSE = 0.XXX (SD = 0.XXX), 
#and MAE = 0.XXX (SD = 0.XXX)."


####################### aim 2: insights into specific living conditions 
#######################################################################
# Peru 
head(rao_location_final)
head(peru_geo_level)

peru_living_condition_analysis = rao_location_final %>%
  filter(Location == "Peru") %>%
  left_join(peru_geo_level, by = "geo")
head(peru_living_condition_analysis)

library(dplyr)
library(purrr)
library(broom)

# Select the columns of interest
cols_to_cor <- peru_living_condition_analysis %>%
  dplyr::select(ss.blender:defecate_recode) %>%
  names()

# Run Spearman correlations
spearman_results <- map_dfr(cols_to_cor, function(col_name) {
  test <- cor.test(
    peru_living_condition_analysis[[col_name]],
    peru_living_condition_analysis$rao,
    method = "spearman"
  )
  
  tibble(
    variable = col_name,
    spearman_rho = test$estimate,
    p_value = test$p.value
  )
})

# View results
spearman_results

################ Kenya 
head(household_kenya) 
kenya_living_condition_analysis = rao_location_final %>%
  filter(Location == "Kenya") %>%
  left_join(household_kenya %>% mutate(geo = as.character(geo)), by = "geo")
head(kenya_living_condition_analysis)

# Select the columns of interest
cols_to_cor <- kenya_living_condition_analysis %>%
  dplyr::select(improved_water:improved_sanitation) %>%
  names()

# Run Spearman correlations
spearman_results <- map_dfr(cols_to_cor, function(col_name) {
  test <- cor.test(
    kenya_living_condition_analysis[[col_name]],
    kenya_living_condition_analysis$rao,
    method = "spearman"
  )
  
  tibble(
    variable = col_name,
    spearman_rho = test$estimate,
    p_value = test$p.value
  )
})

# View results
spearman_results


######################################
######################################
######################################
# March 10 # best predicting model 


# Bangl
predict_wealth_bangl = bangl_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(wealth_bangl_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("geo", "Location"))
head(predict_wealth_bangl)

# Cambodia 
predict_wealth_cam = cambodia_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(cam_wealth_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("geo", "Location")) %>%
  mutate(wealth = wealth + 44001 ) %>% # add a bit so does not go to infintiy
  mutate(wealth = log(wealth))
head(predict_wealth_cam)
summary(predict_wealth_cam$wealth)

# Kenya 
predict_wealth_kenya = kenya_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(wealth_kenya_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("geo", "Location")) %>%
  drop_na()
head(predict_wealth_kenya)

# Peru
predict_wealth_peru = peru_multiplex_dat %>%
  distinct() %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  left_join(wealth_peru_place, by = c("geo", "Location")) %>%
  left_join(rao_location_final , by = c("geo", "Location")) %>%
  drop_na()
head(predict_wealth_peru)
head(predict_wealth_kenya)


set.seed(123)

# Predictor matrix
predictor_cols <- names(predict_wealth_cam)[c(3:8, 10)] # this is cam
predictor_cols <- names(predict_wealth_bangl)[c(3:12, 14)] # this is bangl
predictor_cols <- names(predict_wealth_kenya)[c(3:17, 19)] # this is kenya 
predictor_cols <- names(predict_wealth_peru)[c(3:16, 18)] # this is peru 

x_all <- as.matrix(predict_wealth_peru[, predictor_cols])

# Continuous outcome
y_all <- predict_wealth_peru$wealth

# Remove missing
complete_cases <- complete.cases(x_all, y_all)
x_all <- x_all[complete_cases, ]
y_all <- y_all[complete_cases]

cat(sprintf("Sample size: %d\n", length(y_all)))
cat(sprintf("Predictors: %d\n", ncol(x_all)))


# =============================================================================
# XGBoost Wealth Prediction — Nested Cross-Validation
# Hyperparameters are tuned *inside* each outer fold to prevent leakage.
# =============================================================================

library(xgboost)
library(ggplot2)
library(cowplot)
library(dplyr)

set.seed(123)

# -----------------------------------------------------------------------------
# 1.  COUNTRY-SPECIFIC PREDICTOR COLUMNS
#     Adjust the column indices to match each country's data frame.
# -----------------------------------------------------------------------------

get_predictor_cols <- function(df, country) {
  idx <- switch(country,
                cambodia   = c(3:8,  10),
                bangladesh = c(3:12, 14),
                kenya      = c(3:17, 19),
                peru       = c(3:16, 18),
                stop("Unknown country: ", country)
  )
  names(df)[idx]
}

# -----------------------------------------------------------------------------
# 2.  HYPERPARAMETER GRID
#     Searched inside every outer fold so no test-set information leaks into
#     model selection (proper nested / double cross-validation).
# -----------------------------------------------------------------------------

param_grid <- expand.grid(
  eta               = c(0.01, 0.05, 0.1),
  max_depth         = c(3L, 5L, 7L),
  subsample         = c(0.7, 0.9),
  colsample_bytree  = c(0.7, 0.9),
  min_child_weight  = c(1L, 5L),
  stringsAsFactors  = FALSE
)

# -----------------------------------------------------------------------------
# 3.  INNER CV — tune hyperparameters on the outer-fold training set
# -----------------------------------------------------------------------------

tune_xgb <- function(x_train, y_train,
                     param_grid,
                     inner_folds   = 3L,
                     max_rounds    = 500L,
                     early_stop    = 20L) {
  
  inner_fold_ids <- sample(rep(seq_len(inner_folds),
                               length.out = nrow(x_train)))
  best_score  <- Inf
  best_params <- param_grid[1L, ]
  best_nrounds <- max_rounds
  
  for (i in seq_len(nrow(param_grid))) {
    params_i <- as.list(param_grid[i, ])
    params_i$objective <- "reg:squarederror"
    params_i$eval_metric <- "rmse"
    params_i$nthread <- 1L
    
    inner_rmse <- numeric(inner_folds)
    
    for (j in seq_len(inner_folds)) {
      idx_tr <- which(inner_fold_ids != j)
      idx_vl <- which(inner_fold_ids == j)
      
      y_m  <- mean(y_train[idx_tr])
      y_s  <- sd(y_train[idx_tr])
      
      dtr  <- xgb.DMatrix(x_train[idx_tr, ],
                          label = (y_train[idx_tr] - y_m) / y_s)
      dvl  <- xgb.DMatrix(x_train[idx_vl, ],
                          label = (y_train[idx_vl] - y_m) / y_s)
      
      m <- xgb.train(
        params                = params_i,
        data                  = dtr,
        nrounds               = max_rounds,
        early_stopping_rounds = early_stop,
        watchlist             = list(train = dtr, eval = dvl),
        verbose               = 0
      )
      
      preds_scaled <- predict(m, dvl)
      preds        <- preds_scaled * y_s + y_m
      inner_rmse[j] <- sqrt(mean((y_train[idx_vl] - preds)^2))
    }
    
    mean_rmse <- mean(inner_rmse)
    if (mean_rmse < best_score) {
      best_score   <- mean_rmse
      best_params  <- param_grid[i, ]
      best_nrounds <- m$best_iteration
    }
  }
  
  list(params   = best_params,
       nrounds  = best_nrounds,
       inner_cv_rmse = best_score)
}

# -----------------------------------------------------------------------------
# 4.  OUTER CV — estimate generalisation performance
# -----------------------------------------------------------------------------

run_nested_cv <- function(df, country,
                          outer_folds = 10L,
                          inner_folds = 3L,
                          max_rounds  = 500L,
                          early_stop  = 20L,
                          seed        = 123L) {
  
  set.seed(seed)
  
  predictor_cols <- get_predictor_cols(df, country)
  x_all <- as.matrix(df[, predictor_cols])
  y_all <- df$wealth
  
  # Remove incomplete cases
  keep  <- complete.cases(x_all, y_all)
  x_all <- x_all[keep, ]
  y_all <- y_all[keep]
  
  cat(sprintf(
    "\n[%s]  n = %d  |  predictors = %d\n",
    toupper(country), length(y_all), ncol(x_all)
  ))
  
  fold_ids   <- sample(rep(seq_len(outer_folds), length.out = length(y_all)))
  cv_preds   <- numeric(length(y_all))
  cv_r2      <- numeric(outer_folds)
  cv_rmse    <- numeric(outer_folds)
  cv_mae     <- numeric(outer_folds)
  best_params_log <- vector("list", outer_folds)
  
  for (fold in seq_len(outer_folds)) {
    
    idx_te <- which(fold_ids == fold)
    idx_tr <- which(fold_ids != fold)
    
    # ── Tune inside outer training set ──────────────────────────────────────
    tuned <- tune_xgb(
      x_train     = x_all[idx_tr, ],
      y_train     = y_all[idx_tr],
      param_grid  = param_grid,
      inner_folds = inner_folds,
      max_rounds  = max_rounds,
      early_stop  = early_stop
    )
    
    best_params_log[[fold]] <- tuned$params
    
    # ── Fit final model on full outer training set ───────────────────────────
    y_m <- mean(y_all[idx_tr])
    y_s <- sd(y_all[idx_tr])
    
    final_params             <- as.list(tuned$params)
    final_params$objective   <- "reg:squarederror"
    final_params$eval_metric <- "rmse"
    final_params$nthread     <- 1L
    
    dtrain <- xgb.DMatrix(x_all[idx_tr, ],
                          label = (y_all[idx_tr] - y_m) / y_s)
    dtest  <- xgb.DMatrix(x_all[idx_te, ],
                          label = (y_all[idx_te] - y_m) / y_s)
    
    model <- xgb.train(
      params                = final_params,
      data                  = dtrain,
      nrounds               = max(tuned$nrounds, 10L),
      early_stopping_rounds = early_stop,
      watchlist             = list(train = dtrain, eval = dtest),
      verbose               = 0
    )
    
    preds   <- predict(model, dtest) * y_s + y_m
    actuals <- y_all[idx_te]
    
    cv_preds[idx_te] <- preds
    cv_rmse[fold]    <- sqrt(mean((actuals - preds)^2))
    cv_mae[fold]     <- mean(abs(actuals - preds))
    cv_r2[fold]      <- cor(actuals, preds)^2
    
    cat(sprintf("  Fold %2d | R² = %.3f | RMSE = %.3f | best_eta = %.3f | best_depth = %d\n",
                fold, cv_r2[fold], cv_rmse[fold],
                tuned$params$eta, tuned$params$max_depth))
  }
  
  # ── Summary ────────────────────────────────────────────────────────────────
  cat(sprintf("\n=== %s — %d-Fold Nested CV ===\n", toupper(country), outer_folds))
  cat(sprintf("  R²   : %.3f (SD = %.3f)\n", mean(cv_r2),   sd(cv_r2)))
  cat(sprintf("  RMSE : %.3f (SD = %.3f)\n", mean(cv_rmse), sd(cv_rmse)))
  cat(sprintf("  MAE  : %.3f (SD = %.3f)\n", mean(cv_mae),  sd(cv_mae)))
  
  list(
    country       = country,
    y_all         = y_all,
    cv_preds      = cv_preds,
    fold_ids      = fold_ids,
    cv_r2         = cv_r2,
    cv_rmse       = cv_rmse,
    cv_mae        = cv_mae,
    best_params   = best_params_log
  )
}

# -----------------------------------------------------------------------------
# 5.  RUN ALL COUNTRIES
# -----------------------------------------------------------------------------

results <- list(
  cambodia   = run_nested_cv(predict_wealth_cam,   "cambodia"),
  bangladesh = run_nested_cv(predict_wealth_bangl, "bangladesh"),
  kenya      = run_nested_cv(predict_wealth_kenya, "kenya"),
  peru       = run_nested_cv(predict_wealth_peru,  "peru")
)

# -----------------------------------------------------------------------------
# 6.  PLOTTING FUNCTION
# -----------------------------------------------------------------------------

country_labels <- c(
  cambodia   = "Cambodia",
  bangladesh = "Bangladesh",
  kenya      = "Kenya",
  peru       = "Peru"
)

make_pred_plot <- function(res) {
  label <- country_labels[[res$country]]
  r2    <- mean(res$cv_r2)
  
  plot_df <- data.frame(
    Observed  = res$y_all,
    Predicted = res$cv_preds,
    Fold      = factor(res$fold_ids)
  )
  
  ggplot(plot_df, aes(x = Observed, y = Predicted)) +
    geom_point(aes(color = Fold), alpha = 1, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    labs(
      title    = label,
    #subtitle = bquote(italic(R)^2 == .(sprintf("%.3f", r2))),
    subtitle = paste("R^2 =", round(r2, digits = 4)),
      x        = "Observed Wealth",
      y        = "Predicted Wealth"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "bottom"
    ) +
    scale_color_viridis_d(option = "D")
}

plots <- lapply(results, make_pred_plot)

# Individual plots
plots$cambodia
plots$bangladesh
plots$kenya
plots$peru

# Combined panel figure (3 countries shown as example; adjust as needed)
combined_plot <- plot_grid(
  plots$bangladesh,
  plots$cambodia,
  plots$kenya,
  plots$peru,
  nrow  = 2,
  ncol  = 2,
  align = "hv",
  labels = "AUTO"   # adds A, B, C, D panel labels
)

combined_plot

# Optionally save
# ggsave("wealth_prediction_cv.pdf", combined_plot,
#        width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 7.  SUMMARY TABLE (copy-paste into manuscript)
# -----------------------------------------------------------------------------

summary_table <- do.call(rbind, lapply(results, function(res) {
  data.frame(
    Country = country_labels[[res$country]],
    N       = length(res$y_all),
    R2_mean = round(mean(res$cv_r2),   3),
    R2_sd   = round(sd(res$cv_r2),     3),
    RMSE_mean = round(mean(res$cv_rmse), 3),
    RMSE_sd   = round(sd(res$cv_rmse),   3),
    MAE_mean  = round(mean(res$cv_mae),  3),
    MAE_sd    = round(sd(res$cv_mae),    3)
  )
}))

print(summary_table, row.names = FALSE)

# Methods text template:
# "Model performance was assessed using 10-fold nested cross-validation,
#  in which hyperparameters (learning rate, tree depth, subsampling rate,
#  column-sampling rate, and minimum child weight) were selected via 3-fold
#  inner cross-validation on each outer training set, preventing any leakage
#  of test-set information into model selection. Predictive performance was
#  summarised as the mean R² (squared Pearson correlation between observed and
#  out-of-fold predicted values), RMSE, and MAE across outer folds."


################################################

# =============================================================================
# XGBoost Wealth Prediction — Nested Cross-Validation
# Hyperparameters are tuned *inside* each outer fold to prevent leakage.
# =============================================================================

library(xgboost)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)   # rbindlist for importance aggregation

set.seed(000)

# -----------------------------------------------------------------------------
# 1.  COUNTRY-SPECIFIC PREDICTOR COLUMNS
#     Adjust the column indices to match each country's data frame.
# -----------------------------------------------------------------------------
head(predict_wealth_cam)
head(predict_wealth_bangl)

get_predictor_cols <- function(df, country) {
  idx <- switch(country,
                cambodia   = c(3:8,  10),
                bangladesh = c(3:12, 14),
                kenya      = c(3:17, 19),
                peru       = c(3:16, 18),
                stop("Unknown country: ", country)
  )
  names(df)[idx]
}

# -----------------------------------------------------------------------------
# 2.  HYPERPARAMETER GRID
#     Searched inside every outer fold so no test-set information leaks into
#     model selection (proper nested / double cross-validation).
# -----------------------------------------------------------------------------

param_grid <- expand.grid(
  eta               = c(0.01, 0.05, 0.1),
  max_depth         = c(3L, 5L, 7L),
  subsample         = c(0.7, 0.9),
  colsample_bytree  = c(0.7, 0.9),
  min_child_weight  = c(1L, 5L),
  stringsAsFactors  = FALSE
)

# -----------------------------------------------------------------------------
# 3.  INNER CV — tune hyperparameters on the outer-fold training set
# -----------------------------------------------------------------------------

tune_xgb <- function(x_train, y_train,
                     param_grid,
                     inner_folds   = 3L,
                     max_rounds    = 500L,
                     early_stop    = 20L) {
  
  inner_fold_ids <- sample(rep(seq_len(inner_folds),
                               length.out = nrow(x_train)))
  best_score  <- Inf
  best_params <- param_grid[1L, ]
  best_nrounds <- max_rounds
  
  for (i in seq_len(nrow(param_grid))) {
    params_i <- as.list(param_grid[i, ])
    params_i$objective <- "reg:squarederror"
    params_i$eval_metric <- "rmse"
    params_i$nthread <- 1L
    
    inner_rmse <- numeric(inner_folds)
    
    for (j in seq_len(inner_folds)) {
      idx_tr <- which(inner_fold_ids != j)
      idx_vl <- which(inner_fold_ids == j)
      
      y_m  <- mean(y_train[idx_tr])
      y_s  <- sd(y_train[idx_tr])
      
      dtr  <- xgb.DMatrix(x_train[idx_tr, ],
                          label = (y_train[idx_tr] - y_m) / y_s)
      dvl  <- xgb.DMatrix(x_train[idx_vl, ],
                          label = (y_train[idx_vl] - y_m) / y_s)
      
      m <- xgb.train(
        params                = params_i,
        data                  = dtr,
        nrounds               = max_rounds,
        early_stopping_rounds = early_stop,
        watchlist             = list(train = dtr, eval = dvl),
        verbose               = 0
      )
      
      preds_scaled <- predict(m, dvl)
      preds        <- preds_scaled * y_s + y_m
      inner_rmse[j] <- sqrt(mean((y_train[idx_vl] - preds)^2))
    }
    
    mean_rmse <- mean(inner_rmse)
    if (mean_rmse < best_score) {
      best_score   <- mean_rmse
      best_params  <- param_grid[i, ]
      best_nrounds <- m$best_iteration
    }
  }
  
  list(params   = best_params,
       nrounds  = best_nrounds,
       inner_cv_rmse = best_score)
}

# -----------------------------------------------------------------------------
# 4.  OUTER CV — estimate generalisation performance
# -----------------------------------------------------------------------------

run_nested_cv <- function(df, country,
                          outer_folds = 10L,
                          inner_folds = 3L,
                          max_rounds  = 500L,
                          early_stop  = 20L,
                          seed        = 123L) {
  
  set.seed(seed)
  
  predictor_cols <- get_predictor_cols(df, country)
  x_all <- as.matrix(df[, predictor_cols])
  y_all <- df$wealth
  
  # Remove incomplete cases
  keep  <- complete.cases(x_all, y_all)
  x_all <- x_all[keep, ]
  y_all <- y_all[keep]
  
  cat(sprintf(
    "\n[%s]  n = %d  |  predictors = %d\n",
    toupper(country), length(y_all), ncol(x_all)
  ))
  
  fold_ids        <- sample(rep(seq_len(outer_folds), length.out = length(y_all)))
  cv_preds        <- numeric(length(y_all))
  cv_r2           <- numeric(outer_folds)
  cv_rmse         <- numeric(outer_folds)
  cv_mae          <- numeric(outer_folds)
  best_params_log <- vector("list", outer_folds)
  importance_list <- vector("list", outer_folds)   # <-- collect per-fold importance
  
  for (fold in seq_len(outer_folds)) {
    
    idx_te <- which(fold_ids == fold)
    idx_tr <- which(fold_ids != fold)
    
    # ── Tune inside outer training set ──────────────────────────────────────
    tuned <- tune_xgb(
      x_train     = x_all[idx_tr, ],
      y_train     = y_all[idx_tr],
      param_grid  = param_grid,
      inner_folds = inner_folds,
      max_rounds  = max_rounds,
      early_stop  = early_stop
    )
    
    best_params_log[[fold]] <- tuned$params
    
    # ── Fit final model on full outer training set ───────────────────────────
    y_m <- mean(y_all[idx_tr])
    y_s <- sd(y_all[idx_tr])
    
    final_params             <- as.list(tuned$params)
    final_params$objective   <- "reg:squarederror"
    final_params$eval_metric <- "rmse"
    final_params$nthread     <- 1L
    
    dtrain <- xgb.DMatrix(x_all[idx_tr, ],
                          label = (y_all[idx_tr] - y_m) / y_s)
    dtest  <- xgb.DMatrix(x_all[idx_te, ],
                          label = (y_all[idx_te] - y_m) / y_s)
    
    model <- xgb.train(
      params                = final_params,
      data                  = dtrain,
      nrounds               = max(tuned$nrounds, 10L),
      early_stopping_rounds = early_stop,
      watchlist             = list(train = dtrain, eval = dtest),
      verbose               = 0
    )
    
    # ── Collect variable importance for this fold ────────────────────────────
    imp <- xgb.importance(model = model)
    if (!is.null(imp) && nrow(imp) > 0) {
      importance_list[[fold]] <- imp[, .(Feature, Gain, Cover, Frequency)]
    }
    
    preds   <- predict(model, dtest) * y_s + y_m
    actuals <- y_all[idx_te]
    
    cv_preds[idx_te] <- preds
    cv_rmse[fold]    <- sqrt(mean((actuals - preds)^2))
    cv_mae[fold]     <- mean(abs(actuals - preds))
    cv_r2[fold]      <- cor(actuals, preds)^2
    
    cat(sprintf("  Fold %2d | R² = %.3f | RMSE = %.3f | best_eta = %.3f | best_depth = %d\n",
                fold, cv_r2[fold], cv_rmse[fold],
                tuned$params$eta, tuned$params$max_depth))
  }
  
  # ── Summary ────────────────────────────────────────────────────────────────
  cat(sprintf("\n=== %s — %d-Fold Nested CV ===\n", toupper(country), outer_folds))
  cat(sprintf("  R²   : %.3f (SD = %.3f)\n", mean(cv_r2),   sd(cv_r2)))
  cat(sprintf("  RMSE : %.3f (SD = %.3f)\n", mean(cv_rmse), sd(cv_rmse)))
  cat(sprintf("  MAE  : %.3f (SD = %.3f)\n", mean(cv_mae),  sd(cv_mae)))
  
  # ── Aggregate importance across folds ─────────────────────────────────────
  # Average Gain/Cover/Frequency per feature over all folds that used it.
  # More stable than importance from any single model.
  imp_all <- rbindlist(importance_list, fill = TRUE)
  imp_agg <- imp_all[, .(
    Gain_mean      = mean(Gain,      na.rm = TRUE),
    Gain_sd        = sd(Gain,        na.rm = TRUE),
    Cover_mean     = mean(Cover,     na.rm = TRUE),
    Frequency_mean = mean(Frequency, na.rm = TRUE),
    n_folds        = .N
  ), by = Feature][order(-Gain_mean)]
  imp_agg[, Gain_norm := Gain_mean / sum(Gain_mean)]
  
  list(
    country       = country,
    y_all         = y_all,
    cv_preds      = cv_preds,
    fold_ids      = fold_ids,
    cv_r2         = cv_r2,
    cv_rmse       = cv_rmse,
    cv_mae        = cv_mae,
    best_params   = best_params_log,
    importance    = imp_agg
  )
}

# -----------------------------------------------------------------------------
# 5.  RUN ALL COUNTRIES
# -----------------------------------------------------------------------------

results <- list(
  cambodia   = run_nested_cv(predict_wealth_cam,   "cambodia"),
  bangladesh = run_nested_cv(predict_wealth_bangl, "bangladesh"),  #,
  kenya      = run_nested_cv(predict_wealth_kenya, "kenya"),
  peru       = run_nested_cv(predict_wealth_peru,  "peru")
)

# -----------------------------------------------------------------------------
# 6.  PLOTTING FUNCTION
# -----------------------------------------------------------------------------

country_labels <- c(
  cambodia   = "Cambodia",
  bangladesh = "Bangladesh",
  kenya      = "Kenya",
  peru       = "Peru"
)

make_pred_plot <- function(res) {
  label <- country_labels[[res$country]]
  r2    <- mean(res$cv_r2)
  
  plot_df <- data.frame(
    Observed  = res$y_all,
    Predicted = res$cv_preds,
    Fold      = factor(res$fold_ids)
  )
  
  ggplot(plot_df, aes(x = Observed, y = Predicted)) +
    geom_point(aes(color = Fold), alpha = 1, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    labs(
      title    = label,
      subtitle = paste("R^2 =", round(r2, digits = 4)),
   #   subtitle = bquote(italic(R)^2 == .(sprintf("%.3f", r2))),
      x        = "Observed Wealth",
      y        = "Predicted Wealth"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "none"
    ) +
    scale_color_viridis_d(option = "D")
}

plots <- lapply(results, make_pred_plot)

# Individual plots
plots$cambodia
plots$bangladesh
plots$kenya
plots$peru

# Combined panel figure (3 countries shown as example; adjust as needed)
combined_plot <- plot_grid(
  plots$bangladesh,
  plots$cambodia,
  plots$kenya,
  plots$peru,
  nrow  = 2,
  ncol  = 2,
  align = "hv",
  labels = "AUTO"   # adds A, B, C, D panel labels
)

combined_plot

# Optionally save
# ggsave("wealth_prediction_cv.pdf", combined_plot,
#        width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 7.  VARIABLE IMPORTANCE
#     Importance is averaged across all 10 outer folds for each country.
#     Gain  = average reduction in loss from splits on this feature (primary)
#     Cover = average number of training observations affected
#     Frequency = how often the feature is used across all trees
# -----------------------------------------------------------------------------

# ── Per-country importance plot ───────────────────────────────────────────────
make_importance_plot <- function(res, metric = "Gain_norm", top_n = 15L) {
  
  label  <- country_labels[[res$country]]
  imp    <- as.data.frame(res$importance)
  
  # Keep top N features by normalised Gain
  imp_top <- head(imp[order(-imp$Gain_norm), ], top_n)
  imp_top$Feature <- factor(imp_top$Feature,
                            levels = rev(imp_top$Feature))  # for coord_flip
  
  x_col  <- imp_top[[metric]]
  sd_col <- if (metric == "Gain_norm") imp_top$Gain_sd / sum(imp$Gain_mean) else NULL
  
  p <- ggplot(imp_top, aes(x = Feature, y = .data[[metric]])) +
    geom_col(fill = "#3B528B", width = 0.7) +
    coord_flip() +
    labs(
      title    = label,
    #  subtitle = sprintf("Top %d features by %s (fold-averaged)", top_n, metric),
      x        = NULL,
      y        = switch(metric,
                        Gain_norm      = "Normalized gain",
                        Cover_mean     = "Mean Cover",
                        Frequency_mean = "Mean Frequency",
                        metric)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40", size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # Add error bars for Gain (SD across folds)
  if (metric == "Gain_norm" && !is.null(sd_col)) {
    imp_top$ymin <- pmax(imp_top$Gain_norm - sd_col, 0)
    imp_top$ymax <- imp_top$Gain_norm + sd_col
    p <- p +
      geom_errorbar(data = imp_top,
                    aes(ymin = ymin, ymax = ymax),
                    width = 0.25, color = "gray30")
  }
  
  p
}

# Individual importance plots
imp_plots <- lapply(results, make_importance_plot)
imp_plots$cambodia
imp_plots$bangladesh
imp_plots$kenya
imp_plots$peru



# Combined importance panel
combined_imp <- plot_grid(
  imp_plots$bangladesh,
  imp_plots$cambodia,
  imp_plots$kenya,
  imp_plots$peru,
  nrow   = 2,
  #ncol   = 0,
  align  = "hv"
 # labels = "AUTO"
)
combined_imp

# ── Cross-country importance comparison ───────────────────────────────────────
# Shows which features are consistently important across countries.
# Only retains features appearing in >= 2 countries.

imp_combined <- rbindlist(lapply(results, function(res) {
  d <- as.data.frame(res$importance)
  d$country <- country_labels[[res$country]]
  d
}), fill = TRUE)

feature_counts <- imp_combined[, .(n_countries = uniqueN(country)), by = Feature]
shared_features <- feature_counts[n_countries >= 2, Feature]

imp_shared <- imp_combined[imp_combined$Feature %in% shared_features, ]
imp_shared$Feature <- factor(imp_shared$Feature,
                             levels = imp_combined[imp_combined$country == "Peru", ][order(-Gain_norm), Feature])

cross_country_plot <- ggplot(
  imp_shared,
  aes(x = Feature, y = Gain_norm, fill = country)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_fill_viridis_d(option = "D", name = "Country") +
  labs(
    title    = "Cross-Country Variable Importance",
    subtitle = "Normalised Gain (fold-averaged); features in ≥ 2 countries shown",
    x        = NULL,
    y        = "Normalised Gain"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "gray40", size = 9),
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom"
  )

cross_country_plot

# ── Print importance tables ───────────────────────────────────────────────────
for (nm in names(results)) {
  cat(sprintf("\n--- %s: Top 10 Features (Normalised Gain) ---\n",
              country_labels[[nm]]))
  top10 <- head(as.data.frame(results[[nm]]$importance), 10)
  top10 <- top10[, c("Feature", "Gain_norm", "Gain_sd", "Cover_mean",
                     "Frequency_mean", "n_folds")]
  top10[, 2:5] <- lapply(top10[, 2:5], round, 4)
  print(top10, row.names = FALSE)
}

# Optionally save
# ggsave("wealth_importance.pdf", combined_imp, width = 12, height = 9, dpi = 300)
# ggsave("wealth_importance_cross_country.pdf", cross_country_plot,
#        width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 8.  SUMMARY TABLE (copy-paste into manuscript)
# -----------------------------------------------------------------------------

summary_table <- do.call(rbind, lapply(results, function(res) {
  data.frame(
    Country = country_labels[[res$country]],
    N       = length(res$y_all),
    R2_mean = round(mean(res$cv_r2),   3),
    R2_sd   = round(sd(res$cv_r2),     3),
    RMSE_mean = round(mean(res$cv_rmse), 3),
    RMSE_sd   = round(sd(res$cv_rmse),   3),
    MAE_mean  = round(mean(res$cv_mae),  3),
    MAE_sd    = round(sd(res$cv_mae),    3)
  )
}))

print(summary_table, row.names = FALSE)

# Methods text template:
# "Model performance was assessed using 10-fold nested cross-validation,
#  in which hyperparameters (learning rate, tree depth, subsampling rate,
#  column-sampling rate, and minimum child weight) were selected via 3-fold
#  inner cross-validation on each outer training set, preventing any leakage
#  of test-set information into model selection. Predictive performance was
#  summarised as the mean R² (squared Pearson correlation between observed and
#  out-of-fold predicted values), RMSE, and MAE across outer folds."


## bangl performance 
loc = c("Bangladesh", "Cambodia", "Cambodia", "Bangladesh", "Cambodia", "Bangladesh", 
        "Cambodia", "Bangladesh", "Cambodia", "Bangladesh", "Cambodia", "Bangladesh", 
        "Bangladesh",  "Bangladesh",   "Bangladesh", "Bangladesh")
variables_included = c(11, 7, 6, 10, 5, 9, 4, 8, 3, 7, 2, 6, 5, 4, 3, 2)
r_2 = c(.55, .40, .36, .579, .40, .619, .40, .589, .36, .54, .33, .57, .58, .60, .61, .57)
dat_r2 = data.frame(loc, variables_included, r_2)

ggplot(data = dat_r2) +
  geom_line(aes(x = variables_included, y = r_2, col = loc)) 



##########################################
# age seroprevalence curves for bangladesh 

# serology data 
head(luminex_bangl)

serology_wash_bangl = luminex_bangl %>%
  separate(childid, into = c("prefix", "dataid", "suffix"), sep = "-", remove = FALSE) %>%
  mutate(dataid = as.numeric(dataid))
head(serology_wash_bangl)
print(unique(serology_wash_bangl$dataid))

bangl_wash <- read.csv(file = here("data/bangl/enrollment_wealth/untouched", "washb-bangladesh-enrol-public.csv")) %>%
  mutate(improved_water = ifelse(tubewell == 1 & watmin < 30, 1, 0 ),
         improved_hygiene = ifelse(hwsws == 1, 1, 0), 
         improved_sanitation = ifelse(latown == 1 & latseal == 1 & latfeces == 1, 1, 0), 
         improved_food = ifelse(hfiacat == "Food Secure", 1, 0)) %>%
  mutate( improved_floor = ifelse(cement == 1, 1, 0)) %>%
  dplyr::select(dataid, clusterid, block, improved_water, improved_sanitation, improved_hygiene,
                improved_food, improved_floor) 
head(bangl_wash)  
print(unique(bangl_wash$dataid))

### so 
serology_wash_bangl_model = left_join(bangl_wash ,  serology_wash_bangl, by = "dataid"  ) %>%
  drop_na() %>%
  filter(pathogen == "Giardia") %>%
  filter(agemonth  >= 6 & agemonth <= 28) %>%
  mutate(age_mid = floor(agemonth / 4) * 4 + 1) %>%
  group_by(age_mid, improved_water) %>%
  mutate(seropos = mean(path_pos, na.rm = TRUE)) %>%
  ungroup()
head(serology_wash_bangl_model)


water_bangl = ggplot(data = serology_wash_bangl_model)  +
  geom_point(aes(x = age_mid, y = seropos*100, col = as.character(improved_water))) +
  geom_smooth(aes(x = age_mid, y = seropos*100, col = as.character(improved_water)), 
              method = "loess",
              span = 0.7,
             # se = FALSE,
              linewidth = 1.4) +
  ylab("Seropositivity (%)") +
  theme_bw() +
  scale_color_manual(values = c( "tan", "darkslategray3"), name = "Improved water") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.position = "bottom") + 
 # ggtitle("Bangladesh: Giardia by WASH")
  xlab("Age (months)") +
  scale_x_continuous(
    breaks = seq(3, 28, by = 2))
water_bangl

serology_wash_bangl_model = left_join(bangl_wash ,  serology_wash_bangl, by = "dataid"  ) %>%
  drop_na() %>%
  filter(pathogen == "Giardia") %>%
  filter(agemonth  >= 6 & agemonth <= 28) %>%
  mutate(age_mid = floor(agemonth / 4) * 4 + 1) %>%
  group_by(age_mid, improved_sanitation) %>%
  mutate(seropos = mean(path_pos, na.rm = TRUE)) %>%
  ungroup()
head(serology_wash_bangl_model)

sanitation_bangl = ggplot(data = serology_wash_bangl_model)  +
  geom_point(aes(x = age_mid, y = seropos*100, col = as.character(improved_sanitation))) +
  geom_smooth(aes(x = age_mid, y = seropos*100, col = as.character(improved_sanitation)), 
              formula = y ~ s(x, k = 5), lwd = 1.4) +
  ylab("Seropositivity (%)") +
  theme_bw() +
  scale_color_manual(values = c( "tan", "darkslategray3"), name = "Improved sanitation") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.position = "bottom") +
  xlab("Age (months)") +
  scale_x_continuous(
    breaks = seq(3, 28, by = 2))
sanitation_bangl


serology_wash_bangl_model = left_join(bangl_wash ,  serology_wash_bangl, by = "dataid"  ) %>%
  drop_na() %>%
  filter(pathogen == "Giardia") %>%
  filter(agemonth  >= 6 & agemonth <= 28) %>%
  mutate(age_mid = floor(agemonth / 4) * 4 + 1) %>%
  group_by(age_mid, improved_hygiene) %>%
  mutate(seropos = mean(path_pos, na.rm = TRUE)) %>%
  ungroup()
head(serology_wash_bangl_model)

hygiene_bangl = ggplot(data = serology_wash_bangl_model)  +
  geom_point(aes(x = age_mid, y = seropos*100, col = as.character(improved_hygiene))) +
  geom_smooth(aes(x = age_mid, y = seropos*100, col = as.character(improved_hygiene)), 
              formula = y ~ s(x, k = 5), lwd = 1.4) +
  ylab("Seropositivity (%)") +
  theme_bw() +
  scale_color_manual(values = c( "tan", "darkslategray3"), name = "Improved hygiene") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.position = "bottom")  +
  xlab("Age (months)") +
  scale_x_continuous(
    breaks = seq(3, 28, by = 2))
hygiene_bangl


serology_wash_bangl_model = left_join(bangl_wash ,  serology_wash_bangl, by = "dataid"  ) %>%
  drop_na() %>%
  filter(pathogen == "Giardia") %>%
  filter(agemonth  >= 6 & agemonth <= 28) %>%
  mutate(age_mid = floor(agemonth / 4) * 4 + 1) %>%
  group_by(age_mid, improved_floor) %>%
  mutate(seropos = mean(path_pos, na.rm = TRUE)) %>%
  ungroup()
head(serology_wash_bangl_model)

flooring_bangl = ggplot(data = serology_wash_bangl_model)  +
  geom_point(aes(x = age_mid, y = seropos*100, col = as.character(improved_floor))) +
  geom_smooth(aes(x = age_mid, y = seropos*100, col = as.character(improved_floor)), 
              formula = y ~ s(x, k = 5), lwd = 1.4) +
  ylab("Seropositivity (%)") +
  theme_bw() +
  theme(legend_position = "bottom") +
  scale_color_manual(values = c( "tan", "darkslategray3"), name = "Improved floors") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.position = "bottom") +
  xlab("Age (months)") +
  scale_x_continuous(
    breaks = seq(3, 28, by = 2))

flooring_bangl


plot_grid(water_bangl,  sanitation_bangl, hygiene_bangl, nrow= 1)




##############################################
# Ben method 
serology_wash_bangl_model_ben = left_join(bangl_wash ,  serology_wash_bangl, by = "dataid"  ) %>%
  drop_na() %>%
  filter(pathogen == "Giardia") %>%
  filter(agemonth  >= 6 & agemonth <= 28) 
head(serology_wash_bangl_model_ben)


# Sanitation
fit_nosan = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                %>% filter(improved_sanitation == 0))
summary(fit_nosan)
fit_san = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                  %>% filter(improved_sanitation == 1))

di_pred <- data.frame(agemonth=seq(6,28,by=0.1))
di_pred$Unimproved_Sanitation <- predict(fit_nosan, newdata = di_pred)
di_pred$Improved_Sanitation <- predict(fit_san, newdata = di_pred)
head(di_pred)
di_pred = di_pred %>%
  pivot_longer(
    cols = c(Unimproved_Sanitation, Improved_Sanitation),
    names_to = "Sanitation",
    values_to = "value"
  )
head(di_pred)

san_ben = ggplot(data = di_pred) +
  geom_line(aes(x = agemonth, y = value*100, col = Sanitation), lwd =2) + 
  theme_bw() +
  scale_color_manual(values = c( "tan", "darkslategray3"))  +
  ylab("Seroprevalence (%)") +
  theme(legend_position = "bottom") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11), 
    legend.position = "bottom") +
  xlab("Age (months)")
san_ben

# Water 
fit_nowater = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                  %>% filter(improved_water == 0))
fit_water = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                %>% filter(improved_water == 1))

di_pred <- data.frame(agemonth=seq(6,28,by=0.1))
di_pred$Unimproved_Water <- predict(fit_nowater, newdata = di_pred)
di_pred$Improved_Water <- predict(fit_water, newdata = di_pred)
head(di_pred)
di_pred = di_pred %>%
  pivot_longer(
    cols = c(Unimproved_Water, Improved_Water),
    names_to = "Water",
    values_to = "value"
  )
head(di_pred)

h20_ben = ggplot(data = di_pred) +
  geom_line(aes(x = agemonth, y = value*100, col = Water), lwd =2) + 
  theme_bw() +
  scale_color_manual(values = c( "tan", "darkslategray3"))  +
  ylab("Seroprevalence (%)") +
  theme(legend_position = "bottom") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11), 
    legend.position = "bottom") +
  xlab("Age (months)")
h20_ben


# hygiene 
fit_nohygiene = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                    %>% filter(improved_hygiene == 0))
fit_hygiene = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                  %>% filter(improved_hygiene == 1))

di_pred <- data.frame(agemonth=seq(6,28,by=0.1))
di_pred$Unimproved_Hygiene <- predict(fit_nohygiene, newdata = di_pred)
di_pred$Improved_Hygiene <- predict(fit_hygiene, newdata = di_pred)
head(di_pred)
di_pred = di_pred %>%
  pivot_longer(
    cols = c(Unimproved_Hygiene, Improved_Hygiene),
    names_to = "Hygiene",
    values_to = "value"
  )
head(di_pred)

hygiene_ben = ggplot(data = di_pred) +
  geom_line(aes(x = agemonth, y = value*100, col = Hygiene), lwd =2) + 
  theme_bw() +
  scale_color_manual(values = c( "tan", "darkslategray3"))  +
  ylab("Seroprevalence (%)") +
  theme(legend_position = "bottom") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11), 
    legend.position = "bottom") +
  xlab("Age (months)")
hygiene_ben


# flooring
fit_nofloor = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                      %>% filter(improved_floor == 0))
fit_floor = loess(path_pos ~ agemonth , data = serology_wash_bangl_model_ben  
                    %>% filter(improved_floor == 1))

di_pred <- data.frame(agemonth=seq(6,28,by=0.1))
di_pred$Unimproved_Floor<- predict(fit_nofloor, newdata = di_pred)
di_pred$Improved_Floor <- predict(fit_floor, newdata = di_pred)
head(di_pred)
di_pred = di_pred %>%
  pivot_longer(
    cols = c(Unimproved_Floor, Improved_Floor),
    names_to = "Floor",
    values_to = "value"
  )
head(di_pred)

floor_ben = ggplot(data = di_pred) +
  geom_line(aes(x = agemonth, y = value*100, col = Floor), lwd =2) + 
  theme_bw() +
  scale_color_manual(values = c( "tan", "darkslategray3"))  +
  ylab("Seroprevalence (%)") +
  theme(legend_position = "bottom") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11), 
    legend.position = "bottom") +
  xlab("Age (months)")
floor_ben

plot_grid(h20_ben, san_ben, hygiene_ben, floor_ben, nrow =2 )



# =============================================================================
# Age-Seroprevalence Curves — Logistic GAM, stratified by improved sanitation
#
# Method: Generalised Additive Model with binomial(logit) family.
#   - Correct for binary outcomes (keeps predictions within [0, 1])
#   - Thin-plate regression spline on age (k controls wiggliness; mgcv
#     selects optimal smoothness via GCV/REML automatically)
#   - 95% confidence ribbon on the probability scale via delta method
#   - Cluster-robust standard errors via bam(..., cluster = clusterid)
#     to account for within-cluster correlation
# =============================================================================

library(mgcv)       # GAM fitting
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# 1.  PREPARE DATA
# -----------------------------------------------------------------------------

df <- serology_wash_bangl_model %>%
  mutate(
    sanitation = factor(
      improved_water,
      levels = c(0, 1),
      labels = c("Unimproved sanitation", "Improved sanitation")
    )
  ) %>%
  filter(!is.na(agemonth), !is.na(path_pos), !is.na(sanitation))

cat(sprintf("N = %d  |  Seropositive: %d (%.1f%%)\n",
            nrow(df), sum(df$path_pos), 100 * mean(df$path_pos)))
cat(sprintf("Age range: %.1f – %.1f months\n", min(df$agemonth), max(df$agemonth)))
print(table(df$sanitation, df$path_pos,
            dnn = c("Sanitation", "Seropositive")))

# -----------------------------------------------------------------------------
# 2.  FIT LOGISTIC GAMs — one per stratum
#     bam() is used over gam() for speed with larger datasets and supports
#     cluster= for robust SEs. k=10 gives enough flexibility for 0-60m range;
#     mgcv will penalise back toward a simpler curve if data don't support k.
# -----------------------------------------------------------------------------

fit_sero_gam <- function(data, k = 10) {
  mgcv::bam(
    path_pos ~ s(agemonth, bs = "tp", k = k),
    data    = data,
    family  = binomial(link = "logit"),
    cluster = data$clusterid,   # cluster-robust variance
    method  = "fREML"           # fast REML — appropriate for bam()
  )
}

gam_improved   <- fit_sero_gam(filter(df, sanitation == "Improved sanitation"))
gam_unimproved <- fit_sero_gam(filter(df, sanitation == "Unimproved sanitation"))

# Quick model summaries (check edf — if near k, increase k)
cat("\n--- Improved sanitation ---\n");    summary(gam_improved)
cat("\n--- Unimproved sanitation ---\n");  summary(gam_unimproved)

# -----------------------------------------------------------------------------
# 3.  PREDICTION GRID with 95% CI (delta method on probability scale)
# -----------------------------------------------------------------------------

age_grid <- data.frame(agemonth = seq(min(df$agemonth), max(df$agemonth),
                                      length.out = 200))

predict_ci <- function(model, newdata, label) {
  # predict() with se.fit returns fit + SE on the linear (logit) scale
  pred <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>%
    mutate(
      fit_logit = pred$fit,
      se_logit  = pred$se.fit,
      # Back-transform to probability scale
      fit  = plogis(fit_logit),
      lwr  = plogis(fit_logit - 1.96 * se_logit),
      upr  = plogis(fit_logit + 1.96 * se_logit),
      sanitation = label
    )
}

pred_improved   <- predict_ci(gam_improved,   age_grid, "Improved sanitation")
pred_unimproved <- predict_ci(gam_unimproved, age_grid, "Unimproved sanitation")

pred_all <- bind_rows(pred_improved, pred_unimproved) %>%
  mutate(sanitation = factor(sanitation,
                             levels = c("Unimproved sanitation",
                                        "Improved sanitation")))

# -----------------------------------------------------------------------------
# 4.  EMPIRICAL BINNED PROPORTIONS (overlaid as points for transparency)
#     Bin into ~3-month windows; show observed prevalence with 95% Wilson CI
# -----------------------------------------------------------------------------

bin_width <- 3   # months

empirical <- df %>%
  mutate(age_bin = floor(agemonth / bin_width) * bin_width + bin_width / 2) %>%
  group_by(sanitation, age_bin) %>%
  summarise(
    n_pos  = sum(path_pos),
    n      = n(),
    prev   = n_pos / n,
    # Wilson score interval (better than normal approx for small n / extreme p)
    lwr    = (n_pos + 1.96^2/2 - 1.96 * sqrt(n_pos*(n-n_pos)/n + 1.96^2/4)) /
      (n + 1.96^2),
    upr    = (n_pos + 1.96^2/2 + 1.96 * sqrt(n_pos*(n-n_pos)/n + 1.96^2/4)) /
      (n + 1.96^2),
    .groups = "drop"
  ) %>%
  filter(n >= 5)   # suppress bins with very few observations

# -----------------------------------------------------------------------------
# 5.  PLOT
# -----------------------------------------------------------------------------

palette <- c("Unimproved sanitation" = "#E76F51",
             "Improved sanitation"   = "#2A9D8F")

sero_plot <- ggplot() +
  
  # ── Confidence ribbons ───────────────────────────────────────────────────
 # geom_ribbon(
 #   data = pred_all,
 #   aes(x = agemonth, ymin = lwr, ymax = upr, fill = sanitation),
 #   alpha = 0.18
 # ) +
  
  # ── Smooth fitted curves ─────────────────────────────────────────────────
  geom_line(
    data = pred_all,
    aes(x = agemonth, y = fit, color = sanitation),
    linewidth = 1.1
  ) +
  
  # ── Empirical binned points (sized by n) ─────────────────────────────────
 # geom_pointrange(
  #  data = empirical,
  #  aes(x = age_bin, y = prev, ymin = lwr, ymax = upr,
  #      color = sanitation, size = n),
  #  alpha   = 0.7,
   # fatten  = 2.5,
  #  show.legend = TRUE
 # ) +
  scale_size_continuous(
    name   = "N in bin",
    range  = c(0.3, 1.5),
    breaks = c(10, 50, 100, 200)
  ) +
  
  # ── Scales & labels ──────────────────────────────────────────────────────
  scale_color_manual(values = palette, name = NULL) +
  scale_fill_manual( values = palette, name = NULL) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_x_continuous(
    breaks = seq(0, 60, by = 6),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title    = "Age-Seroprevalence Curve by Sanitation Status",
    subtitle = "Logistic GAM with 95% CI; points show 3-month binned observed prevalence (n ≥ 5)",
    x        = "Age (months)",
    y        = "Seroprevalence"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "gray40", size = 10),
    legend.position = "bottom",
    legend.box      = "vertical",
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.5)),
    fill  = guide_legend(order = 1),
    size  = guide_legend(order = 2)
  )

sero_plot




# -----------------------------------------------------------------------------
# 6.  DIFFERENCE CURVE — improved minus unimproved (on probability scale)
#     Useful for showing *where* in the age range sanitation matters most.
# -----------------------------------------------------------------------------

diff_df <- pred_improved %>%
  select(agemonth, fit_imp = fit, se_imp = se_logit, logit_imp = fit_logit) %>%
  left_join(
    pred_unimproved %>%
      select(agemonth, fit_unimp = fit, se_unimp = se_logit,
             logit_unimp = fit_logit),
    by = "agemonth"
  ) %>%
  mutate(
    # Difference in predicted probabilities + approximate SE via delta method
    diff     = fit_imp - fit_unimp,
    se_diff  = sqrt(
      (plogis(logit_imp) * (1 - plogis(logit_imp)) * se_imp)^2 +
        (plogis(logit_unimp) * (1 - plogis(logit_unimp)) * se_unimp)^2
    ),
    lwr_diff = diff - 1.96 * se_diff,
    upr_diff = diff + 1.96 * se_diff
  )

diff_plot <- ggplot(diff_df, aes(x = agemonth)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = lwr_diff, ymax = upr_diff), fill = "#457B9D", alpha = 0.2) +
  geom_line(aes(y = diff), color = "#457B9D", linewidth = 1.1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(0, 60, by = 6)) +
  labs(
    title    = "Difference in Seroprevalence: Improved − Unimproved Sanitation",
    subtitle = "Negative values = lower seroprevalence in improved sanitation group",
    x        = "Age (months)",
    y        = "Difference in seroprevalence"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title     = element_text(face = "bold"),
    plot.subtitle  = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank()
  )

diff_plot

# Combined two-panel figure
library(cowplot)
combined <- plot_grid(sero_plot, diff_plot,
                      nrow = 2, labels = "AUTO", rel_heights = c(1.4, 1))
combined

# Optionally save
# ggsave("seroprevalence_sanitation.pdf", combined, width = 9, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# 7.  FORMAL TEST — does sanitation modify the age-seroprevalence curve?
#     Fit a joint model with a sanitation × age interaction and test it.
# -----------------------------------------------------------------------------

gam_joint <- mgcv::bam(
  path_pos ~ sanitation +                          # main effect
    s(agemonth, bs = "tp", k = 10) +      # shared smooth
    s(agemonth, by = sanitation,           # difference smooth
      bs = "tp", k = 10),
  data    = df,
  family  = binomial(link = "logit"),
  cluster = df$clusterid,
  method  = "fREML"
)

cat("\n--- Joint model with sanitation × age interaction ---\n")
summary(gam_joint)
# The s(agemonth, by = sanitation) term tests whether the *shape* of the
# age curve differs between groups. A significant p-value means the curves
# are not parallel on the logit scale.







anthro = read.csv(file = here("data/bangl/anthro", "washb-bangladesh-anthro.csv"))
head(anthro)
print(unique(anthro$dataid))


bangl_wash <- read.csv(file = here("data/bangl/enrollment_wealth/untouched", "washb-bangladesh-enrol-public.csv")) %>%
       mutate(improved_water = ifelse(tubewell == 1 & watmin < 30, 1, 0 ),
       improved_hygiene = ifelse(hwsws == 1, 1, 0), 
       improved_sanitation = ifelse(latown == 1 & latseal == 1 & latfeces == 1, 1, 0), 
       improved_food = ifelse(hfiacat == "Food Secure", 1, 0)) %>%
       mutate( improved_floor = ifelse(cement == 1, 1, 0)) %>%
  dplyr::select(dataid, clusterid, block, improved_water, improved_sanitation, 
                improved_food, improved_floor) %>%
  left_join(anthro, by = c("dataid"))
head(bangl_wash)  
print(length(unique(bangl_wash$dataid)))

head(luminex_bangl )

# join this to serology 
serology_wash_bangl = luminex_bangl %>%
  separate(childid, into = c("prefix", "dataid", "suffix"), sep = "-", remove = FALSE) %>%
  mutate(dataid = as.numeric(dataid)) %>%
  left_join(anthro, by = c("dataid"))
print(length(unique(serology_wash_bangl$dataid)))

head(serology_wash_bangl)
  
  left_join(     )
luminex_bangl 


######################################################################################
#####################################################################################
######################################################################################
# March 14 

# Multidimensional Poverty Index 
######################################################################################


# Cambodia 
######################################################################################
cam_wealth_data <- read_dta(file = here("projects/6-multipathogen-burden/data/cambodia", "KHPR61FL.DTA"))
head(cam_wealth_data)


table(cam_wealth_data$hv213)
table(cam_wealth_data$hv214)
table(cam_wealth_data$hv215)

# make mpi 
mpi_cambodia = cam_wealth_data %>%
  mutate(nutrition_deprived = ifelse(is.na(ha40), NA,
                                    ifelse(ha40 < 1850, 1, 0))) %>% # BMI *100
  mutate(mortality = ifelse(is.na(hv111), NA, 
                            ifelse(hv111 == 0, 1, 0))) %>%
  mutate(mat_edu_deprived = ifelse(is.na(ha67), NA, 
                                   ifelse(ha67 == 99, NA, 
                                         ifelse(ha67  == 0 | ha67 == 1, 1, 0)))) %>% # 0 is no education, # 1 is pre-primary 
  mutate(cooking_deprived = ifelse(hv226 == 6 | hv226 == 7 | hv226 == 8 | hv226 == 9| hv226 == 10| hv226 == 11, 1, 0)) %>% 
           # coal, charcoal, wood, straw/grass/ ag crop residue, animal dung 
  mutate(sanitation_deprived = case_when(
      hv205 %in% c(15,23,31,41,42,43) ~ 1, # flush idk where, pit w/out slab, no facility, composting, bucket, hanging toilet 
      hv205 != 99 & hv225 == 1 ~ 1,
      hv205 != 99 & hv225 == 0 ~ 0,
      TRUE ~ NA_real_ ) ) %>%
  mutate(sh104 = replace(sh104, sh104 > 900, NA), sh104b = replace(sh104b, sh104b > 900, NA)) %>%
  mutate(water_deprived = case_when(
      hv237 == 0 ~ 1,
      hv237 == 1 & sh104  > 30 ~ 1,
      hv237 == 1 & sh104b > 30 ~ 1,
      hv237 == 1 ~ 0,
      TRUE ~ NA_real_)) %>%
  mutate( hv206 = replace(hv206, hv206 > 2, NA)) %>%
  mutate(electricity_deprived = ifelse(hv206 == 0, 1, 0 )) %>%
  mutate(housing_deprived = case_when(
      hv213 %in% c(11,21) |
      hv214 %in% c(11,12,13,21,22,23) |
      hv215 %in% c(11,12,21,22,23,24) ~ 1,
       hv213 == 96 | hv214 %in% c(96,99) | hv215 %in% c(96,99) ~ NA_real_,
    TRUE ~ 0)) %>%
  mutate(asset_count = hv208 + hv209 + hv210 + hv211 + hv212 + hv221 + hv243c, na.rm = TRUE ) %>%
  mutate(assets_deprived = ifelse(asset_count <= 1, 1, 0)) %>%
  mutate(poverty = rowSums(cbind(
    nutrition_deprived*(1/6),
    mortality*(1/6),
    mat_edu_deprived*(1/6),
    cooking_deprived*(1/18),
    sanitation_deprived*(1/18),
    water_deprived*(1/18),
    electricity_deprived*(1/18),
    housing_deprived*(1/18),
    assets_deprived*(1/18) ), na.rm = TRUE)) %>%
  mutate(deprived_status = ifelse(poverty > .2777, 1, 0)) %>%
  mutate(dhsclust = hv001) %>%
  group_by(dhsclust) %>%
  mutate(mean_dhs_deprived = mean(deprived_status, na.rm = TRUE)) %>%
  mutate(mean_dhs_intensity = mean(poverty[poverty > 0.2777], na.rm = TRUE)) %>% 
  distinct(mean_dhs_deprived, mean_dhs_intensity , dhsclust) %>%
  mutate(mpi = mean_dhs_deprived*mean_dhs_intensity)
  
ggplot(data = mpi_cambodia) +
  geom_histogram(aes(x = mpi))

# correlation with wealth 
cam_wealth_dat = cam_wealth_data %>%
  dplyr::select(hv001, hv270, hv271) %>%
  mutate(dhsclust = hv001) %>%
  group_by(dhsclust) %>%
  summarize(across(hv270:hv271, ~ mean(.x))) %>%
  left_join(mpi_cambodia, by = c("dhsclust")) %>% drop_na()
head(cam_wealth_dat )
cor(cam_wealth_dat$hv271, cam_wealth_dat$mpi, method = "spearman")

ggplot(data = cam_wealth_dat) +
  geom_point(aes(x = hv271, y = mpi)) + ylab("MPI") +
  xlab("DHS Wealth index") + theme_minimal()

# add gps data 
gps_cambodia =   readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_ea_dhs.csv")) 
head(gps_cambodia)

load(file = here("data/cambodia/cambodia_serology.Rdata"))
gps_dat = cambodia_serology %>%
  dplyr::select(dhsclust, psuid) %>% distinct()
head(gps_dat)

cam_gps_mpi = left_join(gps_cambodia, mpi_cambodia ,  by = "dhsclust") %>%
  left_join(gps_dat, by = "dhsclust") %>%
  mutate(psuid = psuid + 1) %>%
  dplyr::select(psuid, mpi) %>%
  mutate(geo = psuid) %>% mutate(Location = "Cambodia") %>% dplyr::select(-psuid)
head(cam_gps_mpi)

#################################
# serology 
# Cambodia 

# Load disease data. 
cambodia_serology_public = readr::read_csv(file = here("data/cambodia", "cambodia_serology_public.csv")) %>%
  mutate(womanid = ...1)
head(cambodia_serology_public)
#ttmb     : Tetanous toxoid
#bm14     : Lymphatic filariasis bm14
#bm33     : Lymphatic filariasis bm33
#wb123    : Lymphatic filariasis wb123
#nie      : Strongyloides stercoralis NIE
#sag2a    : Toxoplasma gondii SAG2A
#t24      : Taenia solium T24
#pfmsp19  : Plasmodium falciparum MSP-1(19)
#pvmsp19  : Plasmodium vivax MSP-1(19)
# Load location data 
gps_cambodia =   readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_ea_dhs.csv")) 
head(gps_cambodia)
#load(file = here("projects/6-multipathogen-burden/data/cambodia/cambodia_serology.Rdata"))
gps_dat = gps_cambodia %>%
  mutate(psuid = ...1) %>%
  distinct(psuid, dhslat, dhslon)
head(gps_dat)


#######################################################
# Transform into seropositive vs seronegative using: https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004699#pntd-0004699-t001
# https://journals.asm.org/doi/10.1128/cvi.00052-16#T1 for tetanus

antigen = c("ttmb", "bm14", "wb123", "bm33", "nie", "sag2a", "t24", "pfmsp19", "pvmsp19")
cutoff = c(100, 65, 115, 966,792, 159, 486, 343, 196 )
pathogen = c("Tetanus", "Lymphatic filariasis", "Lymphatic filariasis", "Lymphatic filariasis", 
             "Strongyloides stercoralis", "Toxoplasma gondii", "Taenia solium", "Plasmodium falciparum",
             "Plasmodium vivax")
antigen_path_dat = data.frame(antigen, cutoff, pathogen)

antigen_vax = c( "Tetanus")
antigen_path = c("Lymphatic filariasis", "Strongyloides stercoralis", "Toxoplasma gondii", "Taenia solium", "Plasmodium falciparum",
                 "Plasmodium vivax")

cambodia_seropositivity = cambodia_serology_public %>%
  pivot_longer(cols = ttmb:gst, names_to = "antigen", values_to = "value") %>%
  filter(antigen != "gst") %>%
  filter(antigen != "bm33") %>% # we dont want this in our defintion
  mutate(psuid = psuid + 1) %>% 
  left_join(antigen_path_dat, by = "antigen") %>%
  mutate(seropos = ifelse(value > cutoff & pathogen %in% antigen_path , 1, 
                          ifelse(value < cutoff & pathogen %in% antigen_vax, 1, 0))) %>%
  dplyr::select(womanid, psuid, age, parity, antigen, pathogen, seropos) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>%
  filter(pop_seropos < .9499) %>%
  group_by(womanid, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  distinct(womanid, psuid, pathogen, path_pos)
head(cambodia_seropositivity)
table(cambodia_seropositivity$pathogen)

cambodia_multiplex = cambodia_seropositivity  %>%
  group_by(pathogen, psuid) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), 
         num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = psuid) %>%
  drop_na() %>%
  distinct(psuid, pathogen, num, non_na_count, fraction, geo) %>%
  left_join(gps_dat, by = "psuid") %>%
  mutate(Location = "Cambodia")
head(cambodia_multiplex)

# diversity_metrics
cam_diversity_metrics = cambodia_multiplex %>%
  group_by(geo) %>%
  mutate(geo_cases = sum(num)) %>%
  ungroup() %>%
  mutate(p = num/geo_cases + .000001, plogp = p*log2(p), p_squared = p^2) %>%
  group_by(geo) %>%
  mutate(shannon = -sum(plogp), gini_simpson = 1 - sum(p_squared)) %>%
  group_by(geo) %>%
  mutate(alpha = sum(num > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(cam_gps_mpi, by = c("geo", "Location")) %>%
  distinct(geo, shannon, gini_simpson, alpha, mpi, Location) %>%
  drop_na() %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("Location", "geo"))
head(cam_diversity_metrics)


# 
head(rao_location_final)

cor.test(cam_diversity_metrics$shannon, cam_diversity_metrics$mpi, method = "pearson")
cor.test(cam_diversity_metrics$alpha, cam_diversity_metrics$mpi, method = "pearson")
cor.test(cam_diversity_metrics$gini_simpson, cam_diversity_metrics$mpi, method = "pearson")
cor.test(cam_diversity_metrics$rao, cam_diversity_metrics$mpi, method = "pearson")

diversity_met = c("alpha", "shannon", "gini_simpson", "rao")
mid = c(.295, .210, .144, .50)
low = c(.098, .0079, -.06, .33 )
up = c(.469, .39, .34, .64 )
plot_diversity_cam = data.frame(diversity_met, mid, low, up) %>%
  mutate(diversity_met = factor(diversity_met, levels = diversity_met[order(mid)]))

ggplot(plot_diversity_cam, aes(x = diversity_met, y = mid)) +
  geom_point(fill = "blue", width = 0.6) +            # bars
  ggtitle("Cambodia") +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0.1, color = "blue") +  # CI
  geom_hline(yintercept = 0, linetype = "dashed") +    # reference line at 0
  labs(
    x = "Diversity Metric",
    y = "Correlation (pearson) with MPI"
  ) + theme_minimal(base_size = 14)


############ geo analyiss
healthcare_dist = cambodia_multiplex %>%
  distinct(dhslat, dhslon, psuid)
write.csv(healthcare_dist, "healthcare_dist_cambodia.csv")

  

# Bangladesh 
######################################################################################
# treatment clusters
treatment_bangl = read.csv(file = here("data/bangl/parasites/untouched", "bangl_analysis_parasite.csv")) %>%
  dplyr::select(dataid, clusterid, block, personid, tr) %>%
  distinct(clusterid, block, tr)
head(treatment_bangl)
print(treatment_bangl$tr)

bangl_metadata <- read.csv(file = here("data/bangl/enrollment_wealth/untouched", "washb-bangladesh-enrol-public.csv"))
head(bangl_metadata)

table(bangl_metadata$hfiacat)
table(bangl_metadata$momedu)
table(bangl_metadata$latfeces)
table(bangl_metadata$elec)
table(bangl_metadata$roof)
table(bangl_metadata$floor)
table(bangl_metadata$walls)

mpi_bangladesh = bangl_metadata %>%
  left_join(treatment_bangl, by = c("clusterid", "block")) %>%
#  filter(tr %in% c("Control", "Nutrition")) %>%
  mutate(nutrition_deprived = ifelse(hfiacat == "Food Secure", 0, 1)) %>%
  # mortality not available 
  mutate(mat_edu_deprived  = ifelse(momedu == "Secondary (>5y)", 0, 1)) %>%
  # cooking not avaiable 
  mutate(sanitation_deprived = ifelse(latown == 1 & latseal == 1 & latfeces == 1 | tr %in% c("Sanitation") , 0, 1)) %>%
  mutate(water_deprived = ifelse(tubewell == 1 & watmin < 30 | tr %in% c("Water", "WSH") , 0, 1 )) %>%
  mutate(electricity_deprived = ifelse(elec == 1, 0, 1)) %>%
  mutate(housing_deprived = ifelse(roof == 1 | floor == 1 | walls == 1, 0, 1)) %>%
  mutate(asset_count = asset_radio + asset_tv + asset_phone + asset_bike + asset_refrig, na.rm = TRUE ) %>%
  mutate(assets_deprived = ifelse(asset_count <= 1, 1, 0)) %>%
  mutate(poverty = rowSums(cbind(
    nutrition_deprived*(1/6),
  #  mortality*(1/6),
    mat_edu_deprived*(1/6),
   # cooking_deprived*(1/18),
    sanitation_deprived*(1/18),
    water_deprived*(1/18),
    electricity_deprived*(1/18),
    housing_deprived*(1/18),
    assets_deprived*(1/18) ), na.rm = TRUE)) %>%
  mutate(deprived_status = ifelse(poverty > .2037, 1, 0)) %>%
  group_by(block) %>%
  mutate(mean_dhs_deprived = mean(deprived_status, na.rm = TRUE)) %>%
  mutate(mean_dhs_intensity = mean(poverty[poverty > 0.2037], na.rm = TRUE)) %>% 
  distinct(mean_dhs_deprived, mean_dhs_intensity , block) %>%
  mutate(mpi = mean_dhs_deprived*mean_dhs_intensity) %>%
  mutate(block_r = 91 - block) %>%
  mutate(geo = block_r) %>% mutate(Location = "Bangladesh") %>% dplyr::select(-block_r, -block) 

(1/6 + 1/6 + 5/18 )*.33333
head(mpi_bangladesh)

clusters_nc_dat =  bangl_metadata %>%
  left_join(treatment_bangl, by = c("clusterid", "block")) %>%
  filter(tr %in% c("Control", "Nutrition")) 
clusters_nc = unique(clusters_nc_dat$clusterid)

# serology 
# Bangladesh
public_ids = read.csv(file = here("data/bangl/public_ids", "public-ids.csv")) %>%
  distinct(dataid, clusterid, block, clusterid_r, block_r)
head(public_ids)

gps_dat = read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>% # had to add later? 
  left_join(public_ids, by = "dataid") %>% 
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat) 
head(gps_dat)

public_id_cluster = public_ids %>%
  distinct(clusterid, block, block_r, clusterid_r) %>%
  left_join(gps_dat, by = c("block", "block_r"))
head(public_id_cluster)

antigen_vax = c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path = c("T. solium", "Cholera", "E. histolytica" , "Cryptosporidium", "P. falciparum", 
                 "Schistosomiasis", "P. ovale",  "Norovirus",  "P. malariae" , "Onchocerciasis" , "Dengue",         
                 "P. vivax" , "Campylobacter", "Zika", "Salmonella" , "Trachoma" ,"Giardia"  ,      
                 "Chikungunya" ,    "LT-ETEC" ,   "Shigella", "Strongyloides" )

luminex_bangl <- read.csv(file = here("data/bangl/luminex/final", 
                                      "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = c("clusterid")) %>%
 # filter(clusterid %in% clusters_nc) %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path , 1, 
                          ifelse(mfi < mficut & pathogen %in% antigen_vax, 1, 0))) %>%
  filter(visit == 3) %>%
  dplyr::select(mfi, agemonth, clusterid, childid, antigen, pathogen, seropos, block, block_r, med_qgpslong, med_qgpslat) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>%
  filter(pop_seropos < .9499) %>%
  group_by(childid, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(mfi, agemonth, clusterid, childid, block, pathogen, path_pos, block_r, med_qgpslong, med_qgpslat)
head(luminex_bangl)
print(length(luminex_bangl$childid))

# Community level seroprevalence
bangl_multiplex = luminex_bangl  %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = block) %>%
  drop_na()  %>%
  mutate(Location = "Bangladesh")
head(bangl_multiplex)

# diversity 

bangl_diversity_metrics = bangl_multiplex %>%
  group_by(geo) %>%
  mutate(geo_cases = sum(num)) %>%
  ungroup() %>%
  mutate(p = num/geo_cases + .000001, plogp = p*log2(p), p_squared = p^2) %>%
  group_by(geo) %>%
  mutate(shannon = -sum(plogp), gini_simpson = 1- sum(p_squared)) %>%
  group_by(geo) %>%
  mutate(alpha = sum(num > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(mpi_bangladesh, by = c("geo", "Location")) %>%
  distinct(geo, shannon, gini_simpson, alpha, mpi, Location) %>%
  drop_na() %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("Location", "geo"))
head(bangl_diversity_metrics)

cor.test(bangl_diversity_metrics$shannon, bangl_diversity_metrics$mpi, method = "pearson")
cor.test(bangl_diversity_metrics$alpha, bangl_diversity_metrics$mpi, method = "pearson")
cor.test(bangl_diversity_metrics$gini_simpson, bangl_diversity_metrics$mpi, method = "pearson")
cor.test(bangl_diversity_metrics$rao, bangl_diversity_metrics$mpi, method = "pearson")

diversity_met = c("alpha", "shannon", "gini_simpson", "rao")
mid = c(-.035, -.0385, -.069, .374  )
low = c(-.24, -.24, -.27, .181   )
up = c( .17, .170, .14, .539  )
plot_diversity_bangl= data.frame(diversity_met, mid, low, up) %>%
  mutate(diversity_met = factor(diversity_met, levels = diversity_met[order(mid)]))

ggplot(plot_diversity_bangl, aes(x = diversity_met, y = mid)) +
  geom_point(fill = "blue", width = 0.6) +            # bars
  ggtitle("Bangladesh") + 
  geom_errorbar(aes(ymin = low, ymax = up), width = 0.1, color = "blue") +  # CI
  geom_hline(yintercept = 0, linetype = "dashed") +    # reference line at 0
  labs(
    x = "Diversity Metric",
    y = "Correlation (pearson) with MPI"
  ) + theme_minimal(base_size = 14)


#############################################################
# Kenya 
#washk_child_vars has whether they died 

treatment_assignment_k = read.csv(file = here("data/kenya/public_ids", 
                                              "cluster_tx_masked.csv")) 
head(treatment_assignment_k)

# household level 
kenya_meta = read_dta(file = here("data/kenya/master_files", "washk_household_vars_20190305.dta")) %>%
  left_join(treatment_assignment_k, by = "clusterid") 
head(kenya_meta)

table(kenya_meta$HHS_bi)
table(kenya_meta$mother_edu)
table(kenya_meta$tr_masked)

#child level 
child_death = read_dta(file = here("data/kenya/master_files", "washk_child_vars_20190305.dta")) %>%
  group_by(hhid) %>%
  mutate(mortality = ifelse(isdead_el %in% c(1, 3) | childdeath == 1, 1, 0)) %>%
  distinct(hhid, mortality)
head(child_death)

table(child_death$isdead_el)
table(child_death$childdeath)
table(child_death$mortality)

# blocks 
block_kenya = luminex_kenya = read.csv(file = here("data/kenya/luminex/final", 
                                     "washb_kenya_luminex_igg_seropos_pathogen_2026-01-29.csv")) %>%
  distinct(clusterid,block, tr)
head(block_kenya)
table(block_kenya$tr)

head(luminex_kenya)

mpi_kenya = kenya_meta %>% 
  left_join(child_death, by = "hhid") %>%
  left_join(block_kenya, by = "clusterid") %>%
  mutate(nutrition_deprived = ifelse(HHS_bi == 1, 1, 0)) %>% # moderate to severe hunger 
  mutate(mat_edu_deprived = ifelse(mother_edu == 0, 1, 0)) %>% # did not complete primary school
  mutate(sanitation_deprived = ifelse(imp_lat_el == 0 | tr == "WSH", 1, 0)) %>% # at endline improved latrine, maybe want to include treatment  
  mutate(water_deprived = ifelse(prim_drink_ws_bl == 1 & water_time < 31 | tr == "WSH" , 0, 1 )) %>%
  mutate(across(roof:car, ~replace(., . > 8, NA))) %>%
  mutate(electricity_deprived = ifelse(elec == 1, 0 ,1)) %>%
  mutate(cooking_deprived = ifelse(cooker == 1, 0, 1)) %>%
  mutate(housing_deprived = ifelse(roof == 1 | floor == 1 | walls == 1, 0, 1)) %>%
  mutate(asset_count = radio + tv + mobilephone + bicycle + car + motorcycle, na.rm = TRUE ) %>%
  mutate(assets_deprived = ifelse(asset_count <= 1, 1, 0)) %>%
  mutate(poverty = rowSums(cbind(
    nutrition_deprived*(1/6),
    mortality*(1/6),
    mat_edu_deprived*(1/6),
    cooking_deprived*(1/18),
    sanitation_deprived*(1/18),
    water_deprived*(1/18),
    electricity_deprived*(1/18),
    housing_deprived*(1/18),
    assets_deprived*(1/18) ), na.rm = TRUE)) %>%
  mutate(deprived_status = ifelse(poverty > .2778, 1, 0)) %>%
  group_by(block) %>%
  mutate(mean_dhs_deprived = mean(deprived_status, na.rm = TRUE)) %>%
  mutate(mean_dhs_intensity = mean(poverty[poverty > .2778], na.rm = TRUE)) %>% 
  distinct(mean_dhs_deprived, mean_dhs_intensity , block) %>%
  mutate(mpi = mean_dhs_deprived*mean_dhs_intensity)  %>%
  mutate(geo = block) %>%
  mutate(Location = "Kenya")
head(mpi_kenya)
  
# serology 
luminex_kenya = read.csv(file = here("data/kenya/luminex/final", 
                                     "washb_kenya_luminex_igg_seropos_2025-09-21.csv")) %>%
  mutate(dataid = str_extract(childid, "(?<=-)\\d{5}(?=-)")) %>%
 # left_join(treatment_assignment, by = "clusterid") %>%
 # filter(tr == c("Control", "Nutrition")) %>%
  filter(visit == 3) %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path , 1, 
                          ifelse(mfi < mficut & pathogen %in% antigen_vax, 1, 0))) %>%
  dplyr::select(clusterid, childid, antigen, pathogen, seropos) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>%
  filter(pop_seropos < .9499) %>%
  group_by(childid, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  left_join(block_kenya, by = "clusterid") %>%
  distinct(clusterid, childid, block, pathogen, path_pos)
head(luminex_kenya)
print(length(unique(luminex_kenya$childid)))

# Community level seroprevalence
kenya_multiplex = luminex_kenya  %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), 
         num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = block) %>%
  drop_na() %>%
  mutate(Location = "Kenya")
head(kenya_multiplex)

# diversity 
kenya_diversity_metrics = kenya_multiplex %>%
  group_by(geo) %>%
  mutate(geo_cases = sum(num)) %>%
  ungroup() %>%
  mutate(p = num/geo_cases + .000001, plogp = p*log2(p), p_squared = p^2) %>%
  group_by(geo) %>%
  mutate(shannon = -sum(plogp), gini_simpson = 1 - sum(p_squared)) %>%
  group_by(geo) %>%
  mutate(alpha = sum(num > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(mpi_kenya, by = c("geo", "Location")) %>%
  distinct(geo, shannon, gini_simpson, alpha, mpi, Location) %>%
  drop_na() %>%
  left_join(rao_location_final %>% mutate(geo = as.numeric(geo)), by = c("Location", "geo"))
head(kenya_diversity_metrics)

cor.test(kenya_diversity_metrics$shannon, kenya_diversity_metrics$mpi, method = "pearson")
cor.test(kenya_diversity_metrics$alpha, kenya_diversity_metrics$mpi, method = "pearson")
cor.test(kenya_diversity_metrics$gini_simpson, kenya_diversity_metrics$mpi, method = "pearson")
cor.test(kenya_diversity_metrics$rao, kenya_diversity_metrics$mpi, method = "pearson")

diversity_met = c("alpha", "shannon", "gini_simpson", "rao")
mid = c(.2696, .3144, .3131, .368 )
low = c(.0637, .1123, .11, .171  )
up = c(.4534, .4915, .49, .537  )
plot_diversity_kenya = data.frame(diversity_met, mid, low, up) %>%
  mutate(diversity_met = factor(diversity_met, levels = diversity_met[order(mid)]))

ggplot(plot_diversity_kenya, aes(x = diversity_met, y = mid)) +
  geom_point(fill = "blue", width = 0.6) +            # bars
  ggtitle("Kenya") +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0.1, color = "blue") +  # CI
  geom_hline(yintercept = 0, linetype = "dashed") +    # reference line at 0
  labs(
    x = "Diversity Metric",
    y = "Correlation (pearson) with MPI") + theme_minimal(base_size = 14)


############################################################################
# peru ####################################################################

# Peru 
################################################################################

peru_assets_environment = readr::read_csv(file = here("projects/6-multipathogen-burden/data/peru", "loretodata_swabs.csv")) %>%
  mutate(Sample = as.numeric(dbs.filtercode)) 
head(peru_assets_environment)

table(peru_assets_environment$ss.cook.f)
table(peru_assets_environment$ss.floor.f)

# We generally want higher numbers to correspond to better states of asset and environmental exposure.
mpi_peru = peru_assets_environment %>%
  # no nutrition 
  # no maternal education
  mutate(mortality = case_when(
      if_any(alive_1:alive_11, ~ . == 0) ~ 1,
      if_all(alive_1:alive_11, is.na) ~ NA_real_,
      TRUE ~ 0 )) %>%
  mutate(water_deprived = ifelse(td.drinksource.f %in% c("piped to compound/plot", "public tap/stand pipe") 
                                 & td.drinkmins < 31, 0, 1)) %>%
  mutate(sanitation_deprived = ifelse(ss.latrine.f %in% c("no facilities/bush/field/surface water",
                                                          "pit latrine without slab/open pit" ) |
                                        td.defecate.f %in% c("no structure/outside", "shared or public latrine"), 1, 0)) %>%
  mutate(cooking_deprived = ifelse(ss.cook.f == "wood", 1, 0) ) %>%
  # no electricity 
  mutate(housing_deprived = ifelse(ss.floor.f == "dirt/sand", 1, 0)) %>%
  mutate(across(ss.dresser:ss.moto, ~replace(., . > - 10, NA))) %>%
  mutate(asset_count = ss.moto + ss.cell + ss.sound, na.rm = TRUE ) %>%
  mutate(assets_deprived = ifelse(asset_count <= 1, 1, 0)) %>%
  mutate(poverty = rowSums(cbind(
   # nutrition_deprived*(1/6),
    mortality*(1/6),
   # mat_edu_deprived*(1/6),
    cooking_deprived*(1/18),
    sanitation_deprived*(1/18),
    water_deprived*(1/18),
   # electricity_deprived*(1/18),
    housing_deprived*(1/18),
    assets_deprived*(1/18) ), na.rm = TRUE)) %>%
  mutate(deprived_status = ifelse(poverty > .148, 1, 0)) %>%
  group_by(community) %>%
  mutate(mean_dhs_deprived = mean(deprived_status, na.rm = TRUE)) %>%
  mutate(mean_dhs_intensity = mean(poverty[poverty > .148], na.rm = TRUE)) %>% 
  distinct(mean_dhs_deprived, mean_dhs_intensity , community) %>%
  mutate(mpi = mean_dhs_deprived*mean_dhs_intensity)  %>%
  mutate(geo = community) %>%
  mutate(Location = "Peru")

# Serology 

peru_serology = readr::read_csv(file = here("projects/6-multipathogen-burden/data/peru", "perusero.csv")) %>%
  mutate(Sample = as.numeric(Sample)) %>%
  distinct(Sample, .keep_all = TRUE)
length(unique(peru_serology$Sample))
head(peru_serology)

peru_meta = readr::read_csv(file = here("projects/6-multipathogen-burden/data/peru", "loretodata_swabs.csv")) %>%
  mutate(Sample = as.numeric(dbs.filtercode)) %>%
  dplyr::select(gps_long, gps_lat, Sample) %>%
  distinct() %>% drop_na(Sample)
length(unique(peru_meta$Sample))

# process 
peru_serology_antigens = left_join(peru_serology, peru_meta, by = c("Sample")) %>%
  pivot_longer(cols = pgp3:rbd591_pos, names_to = "antigen", values_to = "value") %>%
  filter(!str_detect(antigen, "_pos$")) %>%
  filter(!antigen %in% c("spike", "rbd591", "nucleo", "rbd541")) %>% # remove covid 
  filter(!antigen %in% c("Vero      (78)", "GST     (15)" )) %>% # remove controls 
  drop_na(Sample)
head(peru_serology_antigens)

antigen = unique(peru_serology_antigens$antigen)
print(antigen)
cutoff = c(212, 108, 235, 86, 830, 132, 99, 148, 
           339, 531, 62, 32, 2, 1759, 428, 20, 106,
           375, 107, 58, 422, 306, 155, 18, 140, 19, 14, 
           27, 102, 255, 172)
pathogen = c("Chlamydia trachomatis", "Chlamydia trachomatis", "Treponema palladium", "Treponema palladium", 
             "Cryptosporidium parvum", "Cryptosporidium parvum", "Giardia lamblia","Giardia lamblia", 
             "P. falciparum", "P. falciparum", "P. falciparum","P. falciparum","P. falciparum", "P. falciparum", 
             "P. malariae", "P. malariae", "P. vivax", "Wucheria bancrofti","Wucheria bancrofti", "Wucheria bancrofti", 
             "Onchocerca volvulus", "Onchocerca volvulus", "Strongyloides stercoralis", "T canis", "T gondii", "Taenia solium", "Taenia solium", 
             "Tetanus toxoid", "Diptheria toxoid", "wMEV", "RuV")
antigens_cutoff = data.frame(antigen, cutoff, pathogen)
antigen_path = antigen[1:27]
antigen_vax = antigen[28:31]

peru_serology_processed = left_join(peru_serology_antigens, antigens_cutoff, by = "antigen") %>%
  mutate(seropos = ifelse(value > cutoff & antigen %in% antigen_path , 1, 
                          ifelse(value < cutoff & antigen %in% antigen_vax, 1, 0))) %>%
  dplyr::select(Sample, community, antigen, pathogen, gps_lat, gps_long, seropos) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>% # only want antigens where at least 5% of the community is infected 
  filter(pop_seropos < .9499) %>%
  group_by(Sample, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  dplyr::distinct(Sample, community, pathogen, gps_lat, gps_long, path_pos) %>%
  drop_na()
head(peru_serology_processed)
table(peru_serology_processed$community)

# Community level seroprevalence
peru_multiplex = peru_serology_processed  %>%
  group_by(pathogen, community) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), 
         num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = community) %>%
  drop_na() %>%
  mutate(Location = "Peru") # %>%
 # distinct(community, pathogen, gps_lat, gps_long, fraction, geo, Location)
head(peru_multiplex)

peru_diversity_metrics = peru_multiplex %>%
  group_by(geo) %>%
  mutate(geo_cases = sum(num)) %>%
  ungroup() %>%
  mutate(p = num/geo_cases + .000001, plogp = p*log2(p), p_squared = p^2) %>%
  group_by(geo) %>%
  mutate(shannon = -sum(plogp), gini_simpson = 1 - sum(p_squared)) %>%
  group_by(geo) %>%
  mutate(alpha = sum(num > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(mpi_peru, by = c("geo", "Location")) %>%
  distinct(geo, shannon, gini_simpson, alpha, mpi, Location) %>%
  drop_na() %>%
  left_join(rao_location_final, by = c("Location", "geo"))
head(peru_diversity_metrics)

cor.test(peru_diversity_metrics$shannon, peru_diversity_metrics$mpi, method = "pearson")
cor.test(peru_diversity_metrics$alpha, peru_diversity_metrics$mpi, method = "pearson")
cor.test(peru_diversity_metrics$gini_simpson, peru_diversity_metrics$mpi, method = "pearson")
cor.test(peru_diversity_metrics$rao, peru_diversity_metrics$mpi, method = "pearson")

diversity_met = c("alpha", "shannon", "gini_simpson", "rao")
mid = c( .52, .51, .45 , .012   )
low = c(.12 , .10, .03 , -.42    )
up = c( .78 , .77 , .74, .44   )
plot_diversity_peru= data.frame(diversity_met, mid, low, up) %>%
  mutate(diversity_met = factor(diversity_met, levels = diversity_met[order(mid)]))

ggplot(plot_diversity_peru, aes(x = diversity_met, y = mid)) +
  geom_point(fill = "blue", width = 0.6) +            # bars
  ggtitle("Peru") +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0.1, color = "blue") +  # CI
  geom_hline(yintercept = 0, linetype = "dashed") +    # reference line at 0
  labs(
    x = "Diversity Metric",
    y = "Correlation (pearson) with MPI") + theme_minimal(base_size = 14)


#############################################################################
# AIM2. multipathogen by rap 

# part 1: distance to healthcare facility - line plot 

healthcare_facility_bangl = read.csv(file = here("data/bangl", "lat_lon_bangl_with_distances.csv"))  %>%
  dplyr::select(-X)
head(healthcare_facility_bangl)
summary(healthcare_facility_bangl$distance)
table(healthcare_facility_bangl$nearest_facility_type)

bangl_rao_dist = bangl_multiplex %>%
  distinct(block, block_r, med_qgpslong,  med_qgpslat, geo) %>%
  left_join(healthcare_facility_bangl, by = c("med_qgpslong", "med_qgpslat") )  
head(bangl_rao_dist)

# spatial join 
library(sf)
fac_sf = st_as_sf(healthcare_facility_bangl,
                  coords = c("med_qgpslong","med_qgpslat"),
                  crs = 4326)

bangl_sf = bangl_multiplex %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat, geo) %>%
  st_as_sf(coords = c("med_qgpslong","med_qgpslat"), crs = 4326) %>%
  left_join(bangl_diversity_metrics, by = "geo")
bangl_rao_dist = st_join(bangl_sf, fac_sf, join = st_nearest_feature)
head(bangl_rao_dist)

ggplot(data = bangl_rao_dist) +
  geom_point(aes(x = distance, y = rao)) +
  geom_smooth(aes(x = distance, y = rao), method = "loess", se = TRUE, color = "darkslategray4", linetype = "solid") +
  theme_bw() +
  ylab("Multipathogen burden (Rao)") +
  xlab("Distance (km)") +
  theme(
    strip.text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    # axis.title.y = if(show_y) element_text() else element_blank(),
    legend.position = "none" ) 
  

# part 2: strata of relationship to development: stunting, wasting 


# part 3: strata of WASH conditions 

wash_characterize_bangl = bangl_metadata %>%
  left_join(treatment_bangl, by = c("clusterid", "block")) %>%
  #  filter(tr %in% c("Control", "Nutrition")) %>%
  mutate(nutrition_deprived = ifelse(hfiacat == "Food Secure", 0, 1)) %>%
  # mortality not available 
  mutate(mat_edu_deprived  = ifelse(momedu == "Secondary (>5y)", 0, 1)) %>%
  # cooking not avaiable 
  mutate(sanitation_deprived = ifelse(latown == 1 & latseal == 1 & latfeces == 1 | tr %in% c("Sanitation", "WSH") , 0, 1)) %>%
  mutate(water_deprived = ifelse(tubewell == 1 & watmin < 30 | tr %in% c("Water", "WSH") , 0, 1 )) %>%
  mutate(hygiene_deprived = ifelse(hwsws == 1 |  tr %in% c("Handwashing", "WSH") , 0, 1)) %>% # added
  mutate(floor_deprived = ifelse(cement == 0, 1, 0)) %>% # added
  mutate(electricity_deprived = ifelse(elec == 1, 0, 1)) %>%
  mutate(housing_deprived = ifelse(roof == 1 | floor == 1 | walls == 1, 0, 1)) %>%
  mutate(asset_count = rowSums(across(asset_radio:asset_refrig), na.rm = TRUE)) %>%
  mutate(assets_deprived = ifelse(asset_count <= 1, 1, 0)) %>%
  group_by(block) %>%
  summarize(across(nutrition_deprived:assets_deprived, ~mean(.x, na.rm = TRUE))) %>%
  ungroup() %>% dplyr::select(-asset_count, -nutrition_deprived, -mat_edu_deprived, 
                              -electricity_deprived, -assets_deprived ) %>%
  mutate(block_r = block) %>% dplyr::select(-block) %>%
  left_join(bangl_rao_dist, by = "block_r") %>%
  mutate(rao_quartile = ntile(rao, 4)) %>%
  dplyr::select(block, sanitation_deprived, water_deprived, hygiene_deprived, floor_deprived, block_r, geo,
                rao, rao_quartile) %>%
  pivot_longer(
    cols = c(sanitation_deprived, water_deprived, hygiene_deprived, floor_deprived),
    names_to = "condition",
    values_to = "deprivation"
  ) %>%
  mutate(condition = factor(condition, levels = c( "water_deprived", "hygiene_deprived",  "sanitation_deprived", "floor_deprived"))) %>%
  mutate(condition = gsub("_deprived", "", condition))
head(wash_characterize_bangl)


ggplot(data = wash_characterize_bangl) +
  geom_boxplot(aes(x = condition, y = deprivation*100, fill = as.character(rao_quartile))) +
  theme_minimal() +
  scale_fill_viridis_d(option = "G", name = "Rao", begin = 0.3) +
  ylab("Cluster deprivation (%)") +
  xlab("Environmental risk") +
  theme(
    strip.text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16) )  +
  ggtitle("Bangladesh")


### Child development 

anthro = read.csv(file = here("data/bangl/anthro", "washb-bangladesh-anthro.csv"))
head(anthro)
print(unique(anthro$dataid))

table(anthro$lazminus2) # stunt
table(anthro$whzminus2) # waste
table(anthro$wazminus2) # underweight

develop_bangl = anthro  %>%
  filter(svy == 2) %>%
  mutate(wasted_i = ifelse(whzminus2 == 1, 1, 0)) %>%
  mutate(stunted_i = ifelse(lazminus2 == 1, 1, 0)) %>%
  mutate(underweight_i = ifelse(wazminus2 == 1, 1, 0))  %>%
  group_by(block) %>%
  mutate(wasted = mean(wasted_i, na.rm = TRUE), 
         stunted = mean(stunted_i, na.rm = TRUE), 
         underweight = mean(underweight_i, na.rm = TRUE)) %>%
    ungroup() %>%
  distinct(block, wasted, stunted, underweight) %>%
 # mutate(block_r = block) %>%
  #dplyr::select(-block) %>%
  left_join(bangl_rao_dist, by = "block") %>%
  mutate(rao_quartile = ntile(rao, 2)) %>%
  dplyr::select(block, wasted, stunted, underweight, block_r, geo,
                rao, rao_quartile) %>%
  pivot_longer(
    cols = c(wasted, stunted, underweight),
    names_to = "condition",
    values_to = "deprivation" ) %>%
  mutate(condition = factor(condition, levels = c( "wasted", "underweight", "stunted"))) #%>%
 # mutate(condition = gsub("_deprived", "", condition))


ggplot(data = develop_bangl) +
  geom_boxplot(aes(x = condition, y = deprivation*100, fill = as.character(rao_quartile))) +
  theme_minimal() +
  scale_fill_viridis_d(option = "G", name = "Rao", begin = 0.3) +
  ylab("Cluster deprivation (%)") +
  xlab("Condition") +
  theme(
    strip.text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16) )  +
  ggtitle("Bangladesh")





#############################################################
# april 20 
# for a given location, loop through different combos of pathogens to be treated 
# calculate the effieciency of targeting at least 80% of all disease 
# as you move through new methods 

# start with Bangladesh 

# Bangladesh
public_ids = read.csv(file = here("data/bangl/public_ids", "public-ids.csv")) %>%
  distinct(dataid, clusterid, block, clusterid_r, block_r)
head(public_ids)

gps_dat = read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>% # had to add later? 
  left_join(public_ids, by = "dataid") %>% 
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat) 
head(gps_dat)

public_id_cluster = public_ids %>%
  distinct(clusterid, block, block_r, clusterid_r) %>%
  left_join(gps_dat, by = c("block", "block_r"))
head(public_id_cluster)

antigen_vax = c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path = c("T. solium", "Cholera", "E. histolytica" , "Cryptosporidium", "P. falciparum", 
                 "Schistosomiasis", "P. ovale",  "Norovirus",  "P. malariae" , "Onchocerciasis" , "Dengue",         
                 "P. vivax" , "Campylobacter", "Zika", "Salmonella" , "Trachoma" ,"Giardia"  ,      
                 "Chikungunya" ,    "LT-ETEC" ,   "Shigella", "Strongyloides" )

luminex_bangl <- read.csv(file = here("data/bangl/luminex/final", 
                                      "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = c("clusterid")) %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path , 1, 
                          ifelse(mfi < mficut & pathogen %in% antigen_vax, 1, 0))) %>%
  # filter(visit == 3) %>%
  dplyr::select(mfi, agemonth, clusterid, childid, antigen, pathogen, seropos, block, block_r, med_qgpslong, med_qgpslat) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > .0499) %>%
  filter(pop_seropos < .9499) %>%
  group_by(childid, pathogen ) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(mfi, agemonth, clusterid, childid, block, pathogen, path_pos, block_r, med_qgpslong, med_qgpslat)
head(luminex_bangl)
print(length(unique(luminex_bangl$childid)))

# Community level seroprevalence
bangl_multiplex = luminex_bangl  %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)), num = sum(path_pos)) %>% # number of obs per pathogen
  ungroup() %>%
  mutate(fraction = num/non_na_count) %>%
  mutate(geo = block) %>%
  drop_na()  %>%
  mutate(Location = "Bangladesh")
head(bangl_multiplex)

print(unique(bangl_multiplex$pathogen))

bangl_test = bangl_multiplex  %>%
  filter(pathogen %in% c("Rubella" , "Shigella" ,"Cholera", "Measles", "Salmonella"  ))
print(unique(bangl_test$pathogen))


# TRY IT OUT 
##############################################################
# Multipathogen targeting efficiency analysis
# Compares Shannon entropy, Gini-Simpson, and alpha diversity
# as ranking strategies for identifying minimal spatial coverage
# Goal: find minimum blocks needed to target >= 80% of each disease
##############################################################

library(tidyverse)
library(here)
library(vegan)   # for diversity measures
library(haven)
library(patchwork)

# ---------------------------------------------------------------
# 1. DATA PREP (your existing pipeline — included for completeness)
# ---------------------------------------------------------------

public_ids <- read.csv(file = here("data/bangl/public_ids", "public-ids.csv")) %>%
  distinct(dataid, clusterid, block, clusterid_r, block_r)

gps_dat <- read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>%
  left_join(public_ids, by = "dataid") %>%
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat)

public_id_cluster <- public_ids %>%
  distinct(clusterid, block, block_r, clusterid_r) %>%
  left_join(gps_dat, by = c("block", "block_r"))

antigen_vax  <- c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path <- c("T. solium", "Cholera", "E. histolytica", "Cryptosporidium",
                  "P. falciparum", "Schistosomiasis", "P. ovale", "Norovirus",
                  "P. malariae", "Onchocerciasis", "Dengue", "P. vivax",
                  "Campylobacter", "Zika", "Salmonella", "Trachoma", "Giardia",
                  "Chikungunya", "LT-ETEC", "Shigella", "Strongyloides")

luminex_bangl <- read.csv(
  file = here("data/bangl/luminex/final",
              "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = "clusterid") %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path, 1,
                          ifelse(mfi < mficut & pathogen %in% antigen_vax,  1, 0))) %>%
  dplyr::select(mfi, agemonth, clusterid, childid, antigen, pathogen,
                seropos, block, block_r, med_qgpslong, med_qgpslat) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > 0.0499, pop_seropos < 0.9499) %>%
  group_by(childid, pathogen) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(mfi, agemonth, clusterid, childid, block, pathogen,
           path_pos, block_r, med_qgpslong, med_qgpslat)

bangl_multiplex <- luminex_bangl %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)),
         num          = sum(path_pos)) %>%
  ungroup() %>%
  mutate(fraction = num / non_na_count,
         geo      = block) %>%
  drop_na() %>%
  mutate(Location = "Bangladesh")

# Subset to 5 pathogens of interest
bangl_test <- bangl_multiplex %>%
  filter(pathogen %in% c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella"))

# ---------------------------------------------------------------
# 2. BLOCK-LEVEL SEROPREVALENCE TABLE
#    One row per block_r × pathogen, with the mean seroprevalence
# ---------------------------------------------------------------

block_prev <- bangl_test %>%
  distinct(block_r, pathogen, fraction) %>%
  group_by(block_r, pathogen) %>%
  summarise(prev = mean(fraction, na.rm = TRUE), .groups = "drop")

measles_ordering = block_prev %>%
  filter(pathogen == "Measles") %>%
  arrange(-prev)
head(measles_ordering)
measles_ordering = measles_ordering$block_r

wealth_ordering = wealth_bangl_place %>%
  arrange(wealth)
head(wealth_ordering)
wealth_ordering = wealth_ordering$geo

# Wide format: rows = block_r, columns = pathogens
block_wide <- block_prev %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)

# ---------------------------------------------------------------
# 3. DIVERSITY MEASURE FUNCTIONS
#    All accept a numeric vector of prevalences for one block.
# ---------------------------------------------------------------

shannon_entropy <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)          # normalise so proportions sum to 1
  -sum(p * log(p))
}

gini_simpson <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)
  1 - sum(p^2)             # 1 - Simpson's dominance index
}

alpha_diversity <- function(p) {
  # Effective number of pathogens (Hill number q=1, equivalent to exp(Shannon))
  exp(shannon_entropy(p))
}

# ---------------------------------------------------------------
# 4. COMPUTE DIVERSITY SCORES PER BLOCK
# ---------------------------------------------------------------

pathogen_cols <- setdiff(names(block_wide), "block_r")

diversity_scores <- block_wide %>%
  rowwise() %>%
  mutate(
    shannon      = shannon_entropy(c_across(all_of(pathogen_cols))),
    gini_simpson = gini_simpson(c_across(all_of(pathogen_cols))),
    alpha_div    = alpha_diversity(c_across(all_of(pathogen_cols)))
  ) %>%
  ungroup()

# ---------------------------------------------------------------
# 5. HELPER: minimum blocks to reach >= threshold for all pathogens
#    in a given combination, using a specified ranking column
# ---------------------------------------------------------------

min_blocks_to_threshold <- function(block_prev_sub,   # long-format: block_r, pathogen, prev
                                    block_scores,      # diversity_scores filtered to same blocks
                                    rank_col,          # string: "shannon", "gini_simpson", "alpha_div"
                                    pathogens_in_combo,
                                    threshold = 0.80) {
  
  # Total disease burden per pathogen (denominator for cumulative %)
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  # Rank blocks high→low by the chosen diversity metric
  ranked_blocks <- block_scores %>%
    arrange(desc(.data[[rank_col]])) %>%
    pull(block_r)
  
  # Walk through ranked blocks, accumulating burden
  cumulative <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  for (i in seq_along(ranked_blocks)) {
    blk <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    # Check if all pathogens are at/above threshold
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total)
    
    if (all(pct_reached$pct >= threshold)) {
      return(tibble(
        n_blocks    = i,
        rank_method = rank_col,
        pct_details = list(pct_reached)
      ))
    }
  }
  # Threshold never reached
  return(tibble(
    n_blocks    = NA_integer_,
    rank_method = rank_col,
    pct_details = list(NULL)
  ))
}

# ---------------------------------------------------------------
# 6. LOOP OVER ALL COMBINATIONS OF 3 PATHOGENS (from 5)
# ---------------------------------------------------------------

all_pathogens <- c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella")
combos        <- combn(all_pathogens, 3, simplify = FALSE)
rank_methods  <- c("shannon", "gini_simpson", "alpha_div")

results <- map_dfr(combos, function(combo) {
  combo_label <- paste(sort(combo), collapse = " + ")
  
  # Filter block-level data to only the 3 pathogens in this combo
  block_prev_sub <- block_prev %>%
    filter(pathogen %in% combo)
  
  # Recompute diversity scores using only the 3 pathogens in this combo
  # (so entropy reflects only those pathogens' distribution)
  block_wide_sub <- block_prev_sub %>%
    pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
  
  sub_path_cols <- setdiff(names(block_wide_sub), "block_r")
  
  block_scores_sub <- block_wide_sub %>%
    rowwise() %>%
    mutate(
      shannon      = shannon_entropy(c_across(all_of(sub_path_cols))),
      gini_simpson = gini_simpson(c_across(all_of(sub_path_cols))),
      alpha_div    = alpha_diversity(c_across(all_of(sub_path_cols)))
    ) %>%
    ungroup()
  
  # Run for each diversity measure
  map_dfr(rank_methods, function(rm) {
    res <- min_blocks_to_threshold(
      block_prev_sub    = block_prev_sub,
      block_scores      = block_scores_sub,
      rank_col          = rm,
      pathogens_in_combo = combo
    )
    res %>% mutate(combo = combo_label, .before = 1)
  })
})

# ---------------------------------------------------------------
# 7. SUMMARY TABLE
# ---------------------------------------------------------------

summary_table <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  arrange(combo)

print(summary_table)

# ---------------------------------------------------------------
# 8. EFFICIENCY COMPARISON ACROSS DIVERSITY MEASURES
#    Compare mean / median blocks needed per measure
# ---------------------------------------------------------------

efficiency_summary <- results %>%
  group_by(rank_method) %>%
  summarise(
    mean_blocks   = mean(n_blocks, na.rm = TRUE),
    median_blocks = median(n_blocks, na.rm = TRUE),
    min_blocks    = min(n_blocks, na.rm = TRUE),
    max_blocks    = max(n_blocks, na.rm = TRUE),
    n_combos      = n(),
    n_reached     = sum(!is.na(n_blocks)),
    .groups = "drop"
  ) %>%
  arrange(mean_blocks)

print(efficiency_summary)

# ---------------------------------------------------------------
# 9. VISUALISATIONS
# ---------------------------------------------------------------

# --- 9a. Dot-and-range plot: blocks needed per combo × method ---

method_labels <- c(
  shannon      = "Shannon entropy",
  gini_simpson = "Gini-Simpson",
  alpha_div    = "Alpha diversity\n(Hill q=1)"
)

library(tidyverse)
p1 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  ggplot(aes(x = n_blocks, y = fct_reorder(combo, n_blocks, median),
             colour = rank_method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.85) +
  geom_vline(xintercept = length(unique(block_prev$block_r)) * 0.8,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  labs(
    title    = "Blocks needed to reach 80% coverage of each disease",
    subtitle = "Each row = one 3-pathogen combination; lower = more efficient",
    x        = "Number of blocks required",
    y        = NULL,
    colour   = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
p1

# --- 9b. Boxplot comparing measures across all combos ---

p2 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  ggplot(aes(x = rank_method, y = n_blocks, fill = rank_method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7, fill = "black") +
  labs(
    title    = "Distribution of blocks needed by ranking metric",
    subtitle = "Lower median = more efficient multipathogen targeting strategy",
    x        = NULL,
    y        = "Blocks required for ≥80% coverage"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
p2 

# --- 9c. Head-to-head: Shannon vs Gini-Simpson ---

p3 <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  ggplot(aes(x = shannon, y = gini_simpson)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8, colour = "#4472C4") +
  ggrepel::geom_text_repel(aes(label = combo), size = 2.5, max.overlaps = 10) +
  labs(
    title    = "Shannon vs Gini-Simpson: blocks required",
    subtitle = "Points above the line → Shannon needs fewer blocks",
    x        = "Shannon entropy (n blocks)",
    y        = "Gini-Simpson (n blocks)"
  ) +
  theme_bw(base_size = 11)
p3 

# Combine & save
combined_plot <- (p1 / (p2 | p3)) +
  plot_annotation(
    title   = "Multipathogen targeting efficiency — Bangladesh",
    caption = "80% coverage threshold; 3-pathogen combinations from 5 pathogens"
  )

combined_plot

# ---------------------------------------------------------------
# 10. OPTIONAL: CUMULATIVE COVERAGE CURVES
#     For a single combo & method, plot how coverage accrues
#     as you add blocks in ranked order
# ---------------------------------------------------------------

cumulative_coverage_curve <- function(block_prev_sub,
                                      block_scores,
                                      rank_col,
                                      pathogens_in_combo,
                                      threshold = 0.80) {
  
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  ranked_blocks <- block_scores %>%
    arrange(desc(.data[[rank_col]])) %>%
    pull(block_r)
  
  cumulative   <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  curve_rows   <- list()
  
  for (i in seq_along(ranked_blocks)) {
    blk       <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total, n_blocks = i)
    
    curve_rows[[i]] <- pct_reached
  }
  
  bind_rows(curve_rows)
}

# Example: first combo, Shannon ranking
example_combo  <- combos[[1]]
example_label  <- paste(sort(example_combo), collapse = " + ")

block_prev_ex  <- block_prev %>% filter(pathogen %in% example_combo)

block_wide_ex  <- block_prev_ex %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
sub_cols_ex    <- setdiff(names(block_wide_ex), "block_r")

block_scores_ex <- block_wide_ex %>%
  rowwise() %>%
  mutate(
    shannon      = shannon_entropy(c_across(all_of(sub_cols_ex))),
    gini_simpson = gini_simpson(c_across(all_of(sub_cols_ex))),
    alpha_div    = alpha_diversity(c_across(all_of(sub_cols_ex)))
  ) %>%
  ungroup()

curve_data <- map_dfr(rank_methods, function(rm) {
  cumulative_coverage_curve(
    block_prev_sub     = block_prev_ex,
    block_scores       = block_scores_ex,
    rank_col           = rm,
    pathogens_in_combo = example_combo
  ) %>% mutate(rank_method = recode(rm, !!!method_labels))
})

p_curve <- curve_data %>%
  ggplot(aes(x = n_blocks, y = pct, colour = pathogen, linetype = rank_method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(curve_data$n_blocks) * 0.02, y = 0.82,
           label = "80% threshold", hjust = 0, size = 3, colour = "grey40") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = paste("Cumulative coverage curve —", example_label),
    subtitle = "Blocks added in ranked order (high → low diversity)",
    x        = "Number of blocks included",
    y        = "Cumulative % of total disease burden captured",
    colour   = "Pathogen",
    linetype = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")
p_curve


#################################
# Add Rao

##############################################################
# Multipathogen targeting efficiency analysis
# Compares Shannon entropy, Gini-Simpson, and alpha diversity
# as ranking strategies for identifying minimal spatial coverage
# Goal: find minimum blocks needed to target >= 80% of each disease
##############################################################

library(tidyverse)
library(here)
library(vegan)   # for diversity measures
library(haven)
library(patchwork)

# ---------------------------------------------------------------
# 1. DATA PREP (your existing pipeline — included for completeness)
# ---------------------------------------------------------------

public_ids <- read.csv(file = here("data/bangl/public_ids", "public-ids.csv")) %>%
  distinct(dataid, clusterid, block, clusterid_r, block_r)

gps_dat <- read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>%
  left_join(public_ids, by = "dataid") %>%
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat)

public_id_cluster <- public_ids %>%
  distinct(clusterid, block, block_r, clusterid_r) %>%
  left_join(gps_dat, by = c("block", "block_r"))

antigen_vax  <- c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path <- c("T. solium", "Cholera", "E. histolytica", "Cryptosporidium",
                  "P. falciparum", "Schistosomiasis", "P. ovale", "Norovirus",
                  "P. malariae", "Onchocerciasis", "Dengue", "P. vivax",
                  "Campylobacter", "Zika", "Salmonella", "Trachoma", "Giardia",
                  "Chikungunya", "LT-ETEC", "Shigella", "Strongyloides")

luminex_bangl <- read.csv(
  file = here("data/bangl/luminex/final",
              "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = "clusterid") %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path, 1,
                          ifelse(mfi < mficut & pathogen %in% antigen_vax,  1, 0))) %>%
  dplyr::select(mfi, agemonth, clusterid, childid, antigen, pathogen,
                seropos, block, block_r, med_qgpslong, med_qgpslat) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > 0.0499, pop_seropos < 0.9499) %>%
  group_by(childid, pathogen) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(mfi, agemonth, clusterid, childid, block, pathogen,
           path_pos, block_r, med_qgpslong, med_qgpslat)

bangl_multiplex <- luminex_bangl %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)),
         num          = sum(path_pos)) %>%
  ungroup() %>%
  mutate(fraction = num / non_na_count,
         geo      = block) %>%
  drop_na() %>%
  mutate(Location = "Bangladesh")

# Subset to 5 pathogens of interest
bangl_test <- bangl_multiplex %>%
  filter(pathogen %in% c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella"))

# ---------------------------------------------------------------
# 2. BLOCK-LEVEL SEROPREVALENCE TABLE
#    One row per block_r × pathogen, with the mean seroprevalence
# ---------------------------------------------------------------

block_prev <- bangl_test %>%
  distinct(block_r, pathogen, fraction) %>%
  group_by(block_r, pathogen) %>%
  summarise(prev = mean(fraction, na.rm = TRUE), .groups = "drop")

measles_ordering = block_prev %>%
  filter(pathogen == "Measles") %>%
  arrange(-prev)
head(measles_ordering)
measles_ordering = measles_ordering$block_r

wealth_ordering = wealth_bangl_place %>%
  arrange(wealth)
head(wealth_ordering)
wealth_ordering = wealth_ordering$geo


# Wide format: rows = block_r, columns = pathogens
block_wide <- block_prev %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)

# ---------------------------------------------------------------
# 3. DIVERSITY MEASURE FUNCTIONS
#    All accept a numeric vector of prevalences for one block.
# ---------------------------------------------------------------

shannon_entropy <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)          # normalise so proportions sum to 1
  -sum(p * log(p))
}

gini_simpson <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)
  1 - sum(p^2)             # 1 - Simpson's dominance index
}

alpha_diversity <- function(p) {
  # Effective number of pathogens (Hill number q=1, equivalent to exp(Shannon))
  exp(shannon_entropy(p))
}

rao_quadratic <- function(p) {
  # Simplified Rao's quadratic entropy: sum of products of all pairwise
  # pathogen prevalences within a block (no phylogenetic distance term,
  # so d_ij = 1 for all i != j). Mirrors the products_df approach.
  if (length(p) < 2) return(0)
  pairs <- combn(p, 2, simplify = FALSE)
  sum(sapply(pairs, prod))
}

# ---------------------------------------------------------------
# 4. COMPUTE DIVERSITY SCORES PER BLOCK
# ---------------------------------------------------------------

pathogen_cols <- setdiff(names(block_wide), "block_r")

diversity_scores <- block_wide %>%
  rowwise() %>%
  mutate(
    shannon      = shannon_entropy(c_across(all_of(pathogen_cols))),
    gini_simpson = gini_simpson(c_across(all_of(pathogen_cols))),
    alpha_div    = alpha_diversity(c_across(all_of(pathogen_cols)))
  ) %>%
  ungroup()

# ---------------------------------------------------------------
# 5. HELPER: minimum blocks to reach >= threshold for all pathogens
#    in a given combination, using a specified ranking column
# ---------------------------------------------------------------

min_blocks_to_threshold <- function(block_prev_sub,   # long-format: block_r, pathogen, prev
                                    block_scores,      # diversity_scores filtered to same blocks
                                    rank_col,          # string: "shannon", "gini_simpson", "alpha_div"
                                    pathogens_in_combo,
                                    threshold = 0.80) {
  
  # Total disease burden per pathogen (denominator for cumulative %)
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  # Rank blocks high→low by the chosen diversity metric
  ranked_blocks <- block_scores %>%
    arrange(desc(.data[[rank_col]])) %>%
    pull(block_r)
  
  # Walk through ranked blocks, accumulating burden
  cumulative <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  for (i in seq_along(ranked_blocks)) {
    blk <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    # Check if all pathogens are at/above threshold
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total)
    
    if (all(pct_reached$pct >= threshold)) {
      return(tibble(
        n_blocks    = i,
        rank_method = rank_col,
        pct_details = list(pct_reached)
      ))
    }
  }
  # Threshold never reached
  return(tibble(
    n_blocks    = NA_integer_,
    rank_method = rank_col,
    pct_details = list(NULL)
  ))
}

# ---------------------------------------------------------------
# 6. LOOP OVER ALL COMBINATIONS OF 3 PATHOGENS (from 5)
# ---------------------------------------------------------------

all_pathogens <- c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella")
combos        <- combn(all_pathogens, 3, simplify = FALSE)
rank_methods  <- c("shannon", "gini_simpson", "alpha_div", "rao_quadratic")

results <- map_dfr(combos, function(combo) {
  combo_label <- paste(sort(combo), collapse = " + ")
  
  # Filter block-level data to only the 3 pathogens in this combo
  block_prev_sub <- block_prev %>%
    filter(pathogen %in% combo)
  
  # Recompute diversity scores using only the 3 pathogens in this combo
  # (so entropy reflects only those pathogens' distribution)
  block_wide_sub <- block_prev_sub %>%
    pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
  
  sub_path_cols <- setdiff(names(block_wide_sub), "block_r")
  
  block_scores_sub <- block_wide_sub %>%
    rowwise() %>%
    mutate(
      shannon      = shannon_entropy(c_across(all_of(sub_path_cols))),
      gini_simpson = gini_simpson(c_across(all_of(sub_path_cols))),
      alpha_div    = alpha_diversity(c_across(all_of(sub_path_cols))),
      rao_quadratic = rao_quadratic(c_across(all_of(sub_path_cols)))
    ) %>%
    ungroup()
  
  # Run for each diversity measure
  map_dfr(rank_methods, function(rm) {
    res <- min_blocks_to_threshold(
      block_prev_sub    = block_prev_sub,
      block_scores      = block_scores_sub,
      rank_col          = rm,
      pathogens_in_combo = combo
    )
    res %>% mutate(combo = combo_label, .before = 1)
  })
})

# ---------------------------------------------------------------
# 7. SUMMARY TABLE
# ---------------------------------------------------------------

summary_table <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  arrange(combo)

print(summary_table)

# ---------------------------------------------------------------
# 8. EFFICIENCY COMPARISON ACROSS DIVERSITY MEASURES
#    Compare mean / median blocks needed per measure
# ---------------------------------------------------------------

efficiency_summary <- results %>%
  group_by(rank_method) %>%
  summarise(
    mean_blocks   = mean(n_blocks, na.rm = TRUE),
    median_blocks = median(n_blocks, na.rm = TRUE),
    min_blocks    = min(n_blocks, na.rm = TRUE),
    max_blocks    = max(n_blocks, na.rm = TRUE),
    n_combos      = n(),
    n_reached     = sum(!is.na(n_blocks)),
    .groups = "drop"
  ) %>%
  arrange(mean_blocks)

print(efficiency_summary)

# ---------------------------------------------------------------
# 9. VISUALISATIONS
# ---------------------------------------------------------------

# --- 9a. Dot-and-range plot: blocks needed per combo × method ---

method_labels <- c(
  shannon       = "Shannon entropy",
  gini_simpson  = "Gini-Simpson",
  alpha_div     = "Alpha diversity",
  rao_quadratic = "Rao's quadratic\nentropy"
)

p1 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  ggplot(aes(x = n_blocks, y = fct_reorder(combo, n_blocks, median),
             colour = rank_method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.85) +
  geom_vline(xintercept = length(unique(block_prev$block_r)) * 0.8,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  labs(
    title    = "Blocks needed to reach 80% coverage of each disease",
    subtitle = "Each row = one 3-pathogen combination; lower = more efficient",
    x        = "Number of blocks required",
    y        = NULL,
    colour   = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
p1
# --- 9b. Boxplot comparing measures across all combos ---

p2 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  ggplot(aes(x = rank_method, y = n_blocks, fill = rank_method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7) +
  labs(
 #   title    = "Distribution of blocks needed by ranking metric",
  #  subtitle = "Lower median = more efficient multipathogen targeting strategy",
    x        = NULL,
    y        = "Blocks required for min 80% coverage of each pathogen"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
p2

# --- 9c. Head-to-head: Shannon vs Gini-Simpson ---

p3 <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  ggplot(aes(x = shannon, y = gini_simpson)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8, colour = "#4472C4") +
  ggrepel::geom_text_repel(aes(label = combo), size = 2.5, max.overlaps = 10) +
  labs(
    title    = "Shannon vs Gini-Simpson: blocks required",
    subtitle = "Points above the line → Shannon needs fewer blocks",
    x        = "Shannon entropy (n blocks)",
    y        = "Gini-Simpson (n blocks)"
  ) +
  theme_bw(base_size = 11)

# Combine & save
combined_plot <- (p1 / (p2 | p3)) +
  plot_annotation(
    title   = "Multipathogen targeting efficiency — Bangladesh",
    caption = "80% coverage threshold; 3-pathogen combinations from 5 pathogens"
  )

ggsave(
  filename = here("output", "multipathogen_efficiency.png"),
  plot     = combined_plot,
  width    = 14, height = 12, dpi = 300
)

# ---------------------------------------------------------------
# 10. OPTIONAL: CUMULATIVE COVERAGE CURVES
#     For a single combo & method, plot how coverage accrues
#     as you add blocks in ranked order
# ---------------------------------------------------------------

cumulative_coverage_curve <- function(block_prev_sub,
                                      block_scores,
                                      rank_col,
                                      pathogens_in_combo,
                                      threshold = 0.80) {
  
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  ranked_blocks <- block_scores %>%
    arrange(desc(.data[[rank_col]])) %>%
    pull(block_r)
  
  cumulative   <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  curve_rows   <- list()
  
  for (i in seq_along(ranked_blocks)) {
    blk       <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total, n_blocks = i)
    
    curve_rows[[i]] <- pct_reached
  }
  
  bind_rows(curve_rows)
}

# Example: first combo, Shannon ranking
example_combo  <- combos[[1]]
example_label  <- paste(sort(example_combo), collapse = " + ")

block_prev_ex  <- block_prev %>% filter(pathogen %in% example_combo)

block_wide_ex  <- block_prev_ex %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
sub_cols_ex    <- setdiff(names(block_wide_ex), "block_r")

block_scores_ex <- block_wide_ex %>%
  rowwise() %>%
  mutate(
    shannon       = shannon_entropy(c_across(all_of(sub_cols_ex))),
    gini_simpson  = gini_simpson(c_across(all_of(sub_cols_ex))),
    alpha_div     = alpha_diversity(c_across(all_of(sub_cols_ex))),
    rao_quadratic = rao_quadratic(c_across(all_of(sub_cols_ex)))
  ) %>%
  ungroup()

curve_data <- map_dfr(rank_methods, function(rm) {
  cumulative_coverage_curve(
    block_prev_sub     = block_prev_ex,
    block_scores       = block_scores_ex,
    rank_col           = rm,
    pathogens_in_combo = example_combo
  ) %>% mutate(rank_method = recode(rm, !!!method_labels))
})

p_curve <- curve_data %>%
  ggplot(aes(x = n_blocks, y = pct, colour = pathogen, linetype = rank_method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(curve_data$n_blocks) * 0.02, y = 0.82,
           label = "80% threshold", hjust = 0, size = 3, colour = "grey40") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = paste("Cumulative coverage curve —", example_label),
    subtitle = "Blocks added in ranked order (high → low diversity)",
    x        = "Number of blocks included",
    y        = "Cumulative % of total disease burden captured",
    colour   = "Pathogen",
    linetype = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

ggsave(
  filename = here("output", "coverage_curve_example.png"),
  plot     = p_curve,
  width    = 10, height = 6, dpi = 300
)

message("Done. Results saved to output/")








############ april 21 #################################################
########################################################################
# add comparator strategies 

##############################################################
# Multipathogen targeting efficiency analysis
# Compares Shannon entropy, Gini-Simpson, and alpha diversity
# as ranking strategies for identifying minimal spatial coverage
# Goal: find minimum blocks needed to target >= 80% of each disease
##############################################################

library(tidyverse)
library(here)
library(vegan)   # for diversity measures
library(haven)
library(patchwork)

# ---------------------------------------------------------------
# 1. DATA PREP (your existing pipeline — included for completeness)
# ---------------------------------------------------------------

public_ids <- read.csv(file = here("data/bangl/public_ids", "public-ids.csv")) %>%
  distinct(dataid, clusterid, block, clusterid_r, block_r)

gps_dat <- read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>%
  left_join(public_ids, by = "dataid") %>%
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat)

public_id_cluster <- public_ids %>%
  distinct(clusterid, block, block_r, clusterid_r) %>%
  left_join(gps_dat, by = c("block", "block_r"))

antigen_vax  <- c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path <- c("T. solium", "Cholera", "E. histolytica", "Cryptosporidium",
                  "P. falciparum", "Schistosomiasis", "P. ovale", "Norovirus",
                  "P. malariae", "Onchocerciasis", "Dengue", "P. vivax",
                  "Campylobacter", "Zika", "Salmonella", "Trachoma", "Giardia",
                  "Chikungunya", "LT-ETEC", "Shigella", "Strongyloides")

luminex_bangl <- read.csv(
  file = here("data/bangl/luminex/final",
              "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = "clusterid") %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path, 1,
                          ifelse(mfi < mficut & pathogen %in% antigen_vax,  1, 0))) %>%
  dplyr::select(mfi, agemonth, clusterid, childid, antigen, pathogen,
                seropos, block, block_r, med_qgpslong, med_qgpslat) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > 0.0499, pop_seropos < 0.9499) %>%
  group_by(childid, pathogen) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(mfi, agemonth, clusterid, childid, block, pathogen,
           path_pos, block_r, med_qgpslong, med_qgpslat)

bangl_multiplex <- luminex_bangl %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)),
         num          = sum(path_pos)) %>%
  ungroup() %>%
  mutate(fraction = num / non_na_count,
         geo      = block) %>%
  drop_na() %>%
  mutate(Location = "Bangladesh")

# Subset to 5 pathogens of interest
bangl_test <- bangl_multiplex %>%
  filter(pathogen %in% c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella"))
head(bangl_test)


# ---------------------------------------------------------------
# 2. BLOCK-LEVEL SEROPREVALENCE TABLE
#    One row per block_r × pathogen, with the mean seroprevalence
# ---------------------------------------------------------------

block_prev <- bangl_test %>%
  distinct(block_r, pathogen, fraction) %>%
  group_by(block_r, pathogen) %>%
  summarise(prev = mean(fraction, na.rm = TRUE), .groups = "drop")

# Wide format: rows = block_r, columns = pathogens
block_wide <- block_prev %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)

# ---------------------------------------------------------------
# 3. DIVERSITY MEASURE FUNCTIONS
#    All accept a numeric vector of prevalences for one block.
# ---------------------------------------------------------------

shannon_entropy <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)          # normalise so proportions sum to 1
  -sum(p * log(p))
}

gini_simpson <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)
  1 - sum(p^2)             # 1 - Simpson's dominance index
}

alpha_diversity <- function(p) {
  # Effective number of pathogens (Hill number q=1, equivalent to exp(Shannon))
  exp(shannon_entropy(p))
}

rao_quadratic <- function(p) {
  # Simplified Rao's quadratic entropy: sum of products of all pairwise
  # pathogen prevalences within a block (no phylogenetic distance term,
  # so d_ij = 1 for all i != j). Mirrors the products_df approach.
  if (length(p) < 2) return(0)
  pairs <- combn(p, 2, simplify = FALSE)
  sum(sapply(pairs, prod))
}

# ---------------------------------------------------------------
# 4. COMPUTE DIVERSITY SCORES PER BLOCK
# ---------------------------------------------------------------

pathogen_cols <- setdiff(names(block_wide), "block_r")

diversity_scores <- block_wide %>%
  rowwise() %>%
  mutate(
    shannon      = shannon_entropy(c_across(all_of(pathogen_cols))),
    gini_simpson = gini_simpson(c_across(all_of(pathogen_cols))),
    alpha_div    = alpha_diversity(c_across(all_of(pathogen_cols)))) %>% ungroup()

# ---------------------------------------------------------------
# 5. HELPER: minimum blocks to reach >= threshold for all pathogens
#    in a given combination, using a specified ranking column
# ---------------------------------------------------------------

min_blocks_to_threshold <- function(block_prev_sub,      # long-format: block_r, pathogen, prev
                                    block_scores,         # diversity_scores filtered to same blocks
                                    rank_col,             # string: column in block_scores, OR NULL
                                    pathogens_in_combo,
                                    threshold     = 0.80,
                                    preranked_blocks = NULL) { # optional external ordering vector
  
  # Total disease burden per pathogen (denominator for cumulative %)
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  # Use external ordering if supplied, otherwise rank by diversity column
  if (!is.null(preranked_blocks)) {
    # Keep only blocks present in block_prev_sub; preserve supplied order
    valid_blocks  <- unique(block_prev_sub$block_r)
    ranked_blocks <- preranked_blocks[preranked_blocks %in% valid_blocks]
  } else {
    ranked_blocks <- block_scores %>%
      arrange(desc(.data[[rank_col]])) %>%
      pull(block_r)
  }
  
  # Walk through ranked blocks, accumulating burden
  cumulative <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  for (i in seq_along(ranked_blocks)) {
    blk <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    # Check if all pathogens are at/above threshold
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total)
    
    if (all(pct_reached$pct >= threshold)) {
      return(tibble(
        n_blocks    = i,
        rank_method = rank_col,
        pct_details = list(pct_reached)
      ))
    }
  }
  # Threshold never reached
  return(tibble(
    n_blocks    = NA_integer_,
    rank_method = rank_col,
    pct_details = list(NULL)
  ))
}

# ---------------------------------------------------------------
# 6. LOOP OVER ALL COMBINATIONS OF 3 PATHOGENS (from 5)
# ---------------------------------------------------------------

all_pathogens <- c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella")
combos        <- combn(all_pathogens, 3, simplify = FALSE)
rank_methods  <- c("shannon", "alpha_div", "rao_quadratic")
comparator_methods <- c("measles_order", "wealth_order")

results <- map_dfr(combos, function(combo) {
  combo_label <- paste(sort(combo), collapse = " + ")
  
  # Filter block-level data to only the 3 pathogens in this combo
  block_prev_sub <- block_prev %>%
    filter(pathogen %in% combo)
  
  # Recompute diversity scores using only the 3 pathogens in this combo
  # (so entropy reflects only those pathogens' distribution)
  block_wide_sub <- block_prev_sub %>%
    pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
  
  sub_path_cols <- setdiff(names(block_wide_sub), "block_r")
  
  block_scores_sub <- block_wide_sub %>%
    rowwise() %>%
    mutate(
      shannon      = shannon_entropy(c_across(all_of(sub_path_cols))),
      gini_simpson = gini_simpson(c_across(all_of(sub_path_cols))),
      alpha_div    = alpha_diversity(c_across(all_of(sub_path_cols))),
      rao_quadratic = rao_quadratic(c_across(all_of(sub_path_cols)))
    ) %>%
    ungroup()
  
  # Run diversity-based methods
  diversity_results <- map_dfr(rank_methods, function(rm) {
    res <- min_blocks_to_threshold(
      block_prev_sub     = block_prev_sub,
      block_scores       = block_scores_sub,
      rank_col           = rm,
      pathogens_in_combo = combo
    )
    res %>% mutate(combo = combo_label, .before = 1)
  })
  
  # Run comparator strategies (pre-ranked external orderings)
  comparator_results <- map_dfr(
    list(measles_order = measles_ordering, wealth_order = wealth_ordering),
    function(ordering) {
      res <- min_blocks_to_threshold(
        block_prev_sub     = block_prev_sub,
        block_scores       = block_scores_sub,
        rank_col           = NULL,
        pathogens_in_combo = combo,
        preranked_blocks   = ordering
      )
      res
    },
    .id = "rank_method"
  ) %>%
    mutate(combo = combo_label, .before = 1)
  
  bind_rows(diversity_results, comparator_results)
})

# ---------------------------------------------------------------
# 7. SUMMARY TABLE
# ---------------------------------------------------------------

summary_table <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  arrange(combo)

print(summary_table)

# ---------------------------------------------------------------
# 8. EFFICIENCY COMPARISON ACROSS DIVERSITY MEASURES
#    Compare mean / median blocks needed per measure
# ---------------------------------------------------------------

efficiency_summary <- results %>%
  group_by(rank_method) %>%
  summarise(
    mean_blocks   = mean(n_blocks, na.rm = TRUE),
    median_blocks = median(n_blocks, na.rm = TRUE),
    min_blocks    = min(n_blocks, na.rm = TRUE),
    max_blocks    = max(n_blocks, na.rm = TRUE),
    n_combos      = n(),
    n_reached     = sum(!is.na(n_blocks)),
    .groups = "drop"
  ) %>%
  arrange(mean_blocks)

print(efficiency_summary)

# ---------------------------------------------------------------
# 9. VISUALISATIONS
# ---------------------------------------------------------------

# --- 9a. Dot-and-range plot: blocks needed per combo × method ---

method_labels <- c(
  shannon       = "Shannon entropy",
  alpha_div     = "Alpha diversity",
  rao_quadratic = "Rao's quadratic entropy",
  measles_order = "Measles-motivated",
  wealth_order  = "Wealth"
)

p1 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  ggplot(aes(x = n_blocks, y = fct_reorder(combo, n_blocks, median),
             colour = rank_method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.85) +
  geom_vline(xintercept = length(unique(block_prev$block_r)) * 0.8,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  labs(
    title    = "Blocks needed to reach 80% coverage of each disease",
    subtitle = "Each row = one 3-pathogen combination; lower = more efficient",
    x        = "Number of blocks required",
    y        = NULL,
    colour   = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
p1
# --- 9b. Boxplot comparing measures across all combos ---

p2 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  mutate(rank_method = fct_reorder(rank_method, n_blocks, median, na.rm = TRUE)) %>%
  ggplot(aes(x = rank_method, y = n_blocks, fill = rank_method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7) +
  labs(
  #  title    = "Distribution of blocks needed by ranking metric",
  #  subtitle = "Lower median = more efficient multipathogen targeting strategy",
    x        = NULL,
    y        = "Blocks required to target min 80% of all disease"
  ) +
  theme_bw(base_size = 11) +
  geom_hline(yintercept = 80, lty = "dashed", color = "gray") +
  theme(legend.position = "none", name = "strategy") +
  scale_fill_viridis_d(option = "D" , end = .7 , name = "strategy") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11), 
    legend.position = "bottom") 
p2

# --- 9c. Head-to-head: Shannon vs Rao ---

p3 <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  ggplot(aes(x = shannon, y = rao_quadratic)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8, colour = "#4472C4") +
  ggrepel::geom_text_repel(aes(label = combo), size = 2.5, max.overlaps = 10) +
  labs(
    title    = "Shannon vs Rao's quadratic entropy: blocks required",
    subtitle = "Points above the line → Shannon needs fewer blocks",
    x        = "Shannon entropy (n blocks)",
    y        = "Rao's quadratic entropy (n blocks)"
  ) +
  theme_bw(base_size = 11)
p3
# Combine & save
combined_plot <- (p1 / (p2 | p3)) +
  plot_annotation(
    title   = "Multipathogen targeting efficiency — Bangladesh",
    caption = "80% coverage threshold; 3-pathogen combinations from 5 pathogens"
  )

ggsave(
  filename = here("output", "multipathogen_efficiency.png"),
  plot     = combined_plot,
  width    = 14, height = 12, dpi = 300
)

# ---------------------------------------------------------------
# 10. OPTIONAL: CUMULATIVE COVERAGE CURVES
#     For a single combo & method, plot how coverage accrues
#     as you add blocks in ranked order
# ---------------------------------------------------------------

cumulative_coverage_curve <- function(block_prev_sub,
                                      block_scores,
                                      rank_col,
                                      pathogens_in_combo,
                                      threshold        = 0.80,
                                      preranked_blocks = NULL) {
  
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  if (!is.null(preranked_blocks)) {
    valid_blocks  <- unique(block_prev_sub$block_r)
    ranked_blocks <- preranked_blocks[preranked_blocks %in% valid_blocks]
  } else {
    ranked_blocks <- block_scores %>%
      arrange(desc(.data[[rank_col]])) %>%
      pull(block_r)
  }
  
  cumulative   <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  curve_rows   <- list()
  
  for (i in seq_along(ranked_blocks)) {
    blk       <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total, n_blocks = i)
    
    curve_rows[[i]] <- pct_reached
  }
  
  bind_rows(curve_rows)
}

# Example: first combo, Shannon ranking
example_combo  <- combos[[1]]
example_label  <- paste(sort(example_combo), collapse = " + ")

block_prev_ex  <- block_prev %>% filter(pathogen %in% example_combo)

block_wide_ex  <- block_prev_ex %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
sub_cols_ex    <- setdiff(names(block_wide_ex), "block_r")

block_scores_ex <- block_wide_ex %>%
  rowwise() %>%
  mutate(
    shannon       = shannon_entropy(c_across(all_of(sub_cols_ex))),
    alpha_div     = alpha_diversity(c_across(all_of(sub_cols_ex))),
    rao_quadratic = rao_quadratic(c_across(all_of(sub_cols_ex)))
  ) %>%
  ungroup()

# Diversity-based curves
curve_data_diversity <- map_dfr(rank_methods, function(rm) {
  cumulative_coverage_curve(
    block_prev_sub     = block_prev_ex,
    block_scores       = block_scores_ex,
    rank_col           = rm,
    pathogens_in_combo = example_combo
  ) %>% mutate(rank_method = recode(rm, !!!method_labels))
})

# Comparator curves (pre-ranked orderings)
curve_data_comparators <- map_dfr(
  list(measles_order = measles_ordering, wealth_order = wealth_ordering),
  function(ordering) {
    cumulative_coverage_curve(
      block_prev_sub     = block_prev_ex,
      block_scores       = block_scores_ex,
      rank_col           = NULL,
      pathogens_in_combo = example_combo,
      preranked_blocks   = ordering
    )
  },
  .id = "rank_method"
) %>%
  mutate(rank_method = recode(rank_method, !!!method_labels))

curve_data <- bind_rows(curve_data_diversity, curve_data_comparators)

p_curve <- curve_data %>%
  ggplot(aes(x = n_blocks, y = pct, colour = pathogen, linetype = rank_method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(curve_data$n_blocks) * 0.02, y = 0.82,
           label = "80% threshold", hjust = 0, size = 3, colour = "grey40") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = paste("Cumulative coverage curve —", example_label),
    subtitle = "Blocks added in ranked order (high → low diversity)",
    x        = "Number of blocks included",
    y        = "Cumulative % of total disease burden captured",
    colour   = "Pathogen",
    linetype = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")
p_curve


############################################################################################################
### april 22 
# now i want to assess morans i but only for the rao metric and compare it to intervention efficiency

##############################################################
# Multipathogen targeting efficiency analysis
# Compares Shannon entropy, Gini-Simpson, and alpha diversity
# as ranking strategies for identifying minimal spatial coverage
# Goal: find minimum blocks needed to target >= 80% of each disease
##############################################################

library(tidyverse)
library(here)
library(vegan)   # for diversity measures
library(haven)
library(patchwork)

# ---------------------------------------------------------------
# 1. DATA PREP (your existing pipeline — included for completeness)
# ---------------------------------------------------------------

public_ids <- read.csv(file = here("data/bangl/public_ids", "public-ids.csv")) %>%
  distinct(dataid, clusterid, block, clusterid_r, block_r)

gps_dat <- read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>%
  left_join(public_ids, by = "dataid") %>%
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat)

public_id_cluster <- public_ids %>%
  distinct(clusterid, block, block_r, clusterid_r) %>%
  left_join(gps_dat, by = c("block", "block_r"))

antigen_vax  <- c("Rubella", "Measles", "Tetanus", "Diptheria")
antigen_path <- c("T. solium", "Cholera", "E. histolytica", "Cryptosporidium",
                  "P. falciparum", "Schistosomiasis", "P. ovale", "Norovirus",
                  "P. malariae", "Onchocerciasis", "Dengue", "P. vivax",
                  "Campylobacter", "Zika", "Salmonella", "Trachoma", "Giardia",
                  "Chikungunya", "LT-ETEC", "Shigella", "Strongyloides")

luminex_bangl <- read.csv(
  file = here("data/bangl/luminex/final",
              "washb_bangl_luminex_igg_seropos_2025-09-21.csv")) %>%
  left_join(public_id_cluster, by = "clusterid") %>%
  filter(!pathogen %in% c("COVID19", "Schistosoma GST")) %>%
  mutate(seropos = ifelse(mfi > mficut & pathogen %in% antigen_path, 1,
                          ifelse(mfi < mficut & pathogen %in% antigen_vax,  1, 0))) %>%
  dplyr::select(mfi, agemonth, clusterid, childid, antigen, pathogen,
                seropos, block, block_r, med_qgpslong, med_qgpslat) %>%
  group_by(antigen) %>%
  mutate(pop_seropos = mean(seropos, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(pop_seropos > 0.0499, pop_seropos < 0.9499) %>%
  group_by(childid, pathogen) %>%
  mutate(path_pos = as.integer(all(seropos == 1))) %>%
  ungroup() %>%
  distinct(mfi, agemonth, clusterid, childid, block, pathogen,
           path_pos, block_r, med_qgpslong, med_qgpslat)

bangl_multiplex <- luminex_bangl %>%
  group_by(pathogen, block) %>%
  mutate(non_na_count = sum(!is.na(path_pos)),
         num          = sum(path_pos)) %>%
  ungroup() %>%
  mutate(fraction = num / non_na_count,
         geo      = block) %>%
  drop_na() %>%
  mutate(Location = "Bangladesh")

# Subset to 5 pathogens of interest
bangl_test <- bangl_multiplex %>%
  filter(pathogen %in% c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella"))

# ---------------------------------------------------------------
# 2. BLOCK-LEVEL SEROPREVALENCE TABLE
#    One row per block_r × pathogen, with the mean seroprevalence
# ---------------------------------------------------------------

block_prev <- bangl_test %>%
  distinct(block_r, pathogen, fraction) %>%
  group_by(block_r, pathogen) %>%
  summarise(prev = mean(fraction, na.rm = TRUE), .groups = "drop")

# Wide format: rows = block_r, columns = pathogens
block_wide <- block_prev %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)

# ---------------------------------------------------------------
# 3. DIVERSITY MEASURE FUNCTIONS
#    All accept a numeric vector of prevalences for one block.
# ---------------------------------------------------------------

shannon_entropy <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)          # normalise so proportions sum to 1
  -sum(p * log(p))
}

gini_simpson <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)
  1 - sum(p^2)             # 1 - Simpson's dominance index
}

alpha_diversity <- function(p) {
  # Effective number of pathogens (Hill number q=1, equivalent to exp(Shannon))
  exp(shannon_entropy(p))
}

rao_quadratic <- function(p) {
  # Simplified Rao's quadratic entropy: sum of products of all pairwise
  # pathogen prevalences within a block (no phylogenetic distance term,
  # so d_ij = 1 for all i != j). Mirrors the products_df approach.
  if (length(p) < 2) return(0)
  pairs <- combn(p, 2, simplify = FALSE)
  sum(sapply(pairs, prod))
}

# ---------------------------------------------------------------
# 4. COMPUTE DIVERSITY SCORES PER BLOCK
# ---------------------------------------------------------------

pathogen_cols <- setdiff(names(block_wide), "block_r")

diversity_scores <- block_wide %>%
  rowwise() %>%
  mutate(
    shannon      = shannon_entropy(c_across(all_of(pathogen_cols))),
    gini_simpson = gini_simpson(c_across(all_of(pathogen_cols))),
    alpha_div    = alpha_diversity(c_across(all_of(pathogen_cols)))
  ) %>%
  ungroup()

# ---------------------------------------------------------------
# 5. HELPER: minimum blocks to reach >= threshold for all pathogens
#    in a given combination, using a specified ranking column
# ---------------------------------------------------------------

min_blocks_to_threshold <- function(block_prev_sub,      # long-format: block_r, pathogen, prev
                                    block_scores,         # diversity_scores filtered to same blocks
                                    rank_col,             # string: column in block_scores, OR NULL
                                    pathogens_in_combo,
                                    threshold     = 0.80,
                                    preranked_blocks = NULL) { # optional external ordering vector
  
  # Total disease burden per pathogen (denominator for cumulative %)
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  # Use external ordering if supplied, otherwise rank by diversity column
  if (!is.null(preranked_blocks)) {
    # Keep only blocks present in block_prev_sub; preserve supplied order
    valid_blocks  <- unique(block_prev_sub$block_r)
    ranked_blocks <- preranked_blocks[preranked_blocks %in% valid_blocks]
  } else {
    ranked_blocks <- block_scores %>%
      arrange(desc(.data[[rank_col]])) %>%
      pull(block_r)
  }
  
  # Walk through ranked blocks, accumulating burden
  cumulative <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  for (i in seq_along(ranked_blocks)) {
    blk <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    # Check if all pathogens are at/above threshold
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total)
    
    if (all(pct_reached$pct >= threshold)) {
      return(tibble(
        n_blocks    = i,
        rank_method = rank_col,
        pct_details = list(pct_reached)
      ))
    }
  }
  # Threshold never reached
  return(tibble(
    n_blocks    = NA_integer_,
    rank_method = rank_col,
    pct_details = list(NULL)
  ))
}

# ---------------------------------------------------------------
# 6. LOOP OVER ALL COMBINATIONS OF 3 PATHOGENS (from 5)
# ---------------------------------------------------------------

all_pathogens <- c("Rubella", "Shigella", "Cholera", "Measles", "Salmonella")
combos        <- combn(all_pathogens, 3, simplify = FALSE)
rank_methods  <- c("shannon", "alpha_div", "rao_quadratic")
comparator_methods <- c("measles_order", "wealth_order")

results <- map_dfr(combos, function(combo) {
  combo_label <- paste(sort(combo), collapse = " + ")
  
  # Filter block-level data to only the 3 pathogens in this combo
  block_prev_sub <- block_prev %>%
    filter(pathogen %in% combo)
  
  # Recompute diversity scores using only the 3 pathogens in this combo
  # (so entropy reflects only those pathogens' distribution)
  block_wide_sub <- block_prev_sub %>%
    pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
  
  sub_path_cols <- setdiff(names(block_wide_sub), "block_r")
  
  block_scores_sub <- block_wide_sub %>%
    rowwise() %>%
    mutate(
      shannon       = shannon_entropy(c_across(all_of(sub_path_cols))),
      alpha_div     = alpha_diversity(c_across(all_of(sub_path_cols))),
      rao_quadratic = rao_quadratic(c_across(all_of(sub_path_cols)))
    ) %>%
    ungroup()
  
  # Run diversity-based methods
  diversity_results <- map_dfr(rank_methods, function(rm) {
    res <- min_blocks_to_threshold(
      block_prev_sub     = block_prev_sub,
      block_scores       = block_scores_sub,
      rank_col           = rm,
      pathogens_in_combo = combo
    )
    res %>% mutate(combo = combo_label, .before = 1)
  })
  
  # Run comparator strategies (pre-ranked external orderings)
  comparator_results <- map_dfr(
    list(measles_order = measles_ordering, wealth_order = wealth_ordering),
    function(ordering) {
      res <- min_blocks_to_threshold(
        block_prev_sub     = block_prev_sub,
        block_scores       = block_scores_sub,
        rank_col           = NULL,
        pathogens_in_combo = combo,
        preranked_blocks   = ordering
      )
      res
    },
    .id = "rank_method"
  ) %>%
    mutate(combo = combo_label, .before = 1)
  
  bind_rows(diversity_results, comparator_results)
})

# ---------------------------------------------------------------
# 7. SUMMARY TABLE
# ---------------------------------------------------------------

summary_table <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  arrange(combo)

print(summary_table)

# ---------------------------------------------------------------
# 8. EFFICIENCY COMPARISON ACROSS DIVERSITY MEASURES
#    Compare mean / median blocks needed per measure
# ---------------------------------------------------------------

efficiency_summary <- results %>%
  group_by(rank_method) %>%
  summarise(
    mean_blocks   = mean(n_blocks, na.rm = TRUE),
    median_blocks = median(n_blocks, na.rm = TRUE),
    min_blocks    = min(n_blocks, na.rm = TRUE),
    max_blocks    = max(n_blocks, na.rm = TRUE),
    n_combos      = n(),
    n_reached     = sum(!is.na(n_blocks)),
    .groups = "drop"
  ) %>%
  arrange(mean_blocks)

print(efficiency_summary)

# ---------------------------------------------------------------
# 9. VISUALISATIONS
# ---------------------------------------------------------------

# --- 9a. Dot-and-range plot: blocks needed per combo × method ---

method_labels <- c(
  shannon       = "Shannon entropy",
  alpha_div     = "Alpha diversity",
  rao_quadratic = "Rao's quadratic entropy",
  measles_order = "Measles-motivated",
  wealth_order  = "Wealth"
)

p1 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  ggplot(aes(x = n_blocks, y = fct_reorder(combo, n_blocks, median),
             colour = rank_method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.85) +
  geom_vline(xintercept = length(unique(block_prev$block_r)) * 0.8,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  labs(
    title    = "Blocks needed to reach 80% coverage of each disease",
    subtitle = "Each row = one 3-pathogen combination; lower = more efficient",
    x        = "Number of blocks required",
    y        = NULL,
    colour   = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
p1
# --- 9b. Boxplot comparing measures across all combos ---

p2 <- results %>%
  mutate(rank_method = recode(rank_method, !!!method_labels)) %>%
  mutate(rank_method = fct_reorder(rank_method, n_blocks, median, na.rm = TRUE)) %>%
  ggplot(aes(x = rank_method, y = n_blocks, fill = rank_method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7) +
  labs(
    title    = "Distribution of blocks needed by ranking metric",
    subtitle = "Lower median = more efficient multipathogen targeting strategy",
    x        = NULL,
    y        = "Blocks required for ≥80% coverage"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
p2
# --- 9c. Head-to-head: Shannon vs Rao ---

p3 <- results %>%
  dplyr::select(combo, rank_method, n_blocks) %>%
  pivot_wider(names_from = rank_method, values_from = n_blocks) %>%
  ggplot(aes(x = shannon, y = rao_quadratic)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8, colour = "#4472C4") +
  ggrepel::geom_text_repel(aes(label = combo), size = 2.5, max.overlaps = 10) +
  labs(
    title    = "Shannon vs Rao's quadratic entropy: blocks required",
    subtitle = "Points above the line → Shannon needs fewer blocks",
    x        = "Shannon entropy (n blocks)",
    y        = "Rao's quadratic entropy (n blocks)"
  ) +
  theme_bw(base_size = 11)
p3
# Combine & save
combined_plot <- (p1 / (p2 | p3)) +
  plot_annotation(
    title   = "Multipathogen targeting efficiency — Bangladesh",
    caption = "80% coverage threshold; 3-pathogen combinations from 5 pathogens" )


# ---------------------------------------------------------------
# 10. OPTIONAL: CUMULATIVE COVERAGE CURVES
#     For a single combo & method, plot how coverage accrues
#     as you add blocks in ranked order
# ---------------------------------------------------------------

cumulative_coverage_curve <- function(block_prev_sub,
                                      block_scores,
                                      rank_col,
                                      pathogens_in_combo,
                                      threshold        = 0.80,
                                      preranked_blocks = NULL) {
  
  total_burden <- block_prev_sub %>%
    filter(pathogen %in% pathogens_in_combo) %>%
    group_by(pathogen) %>%
    summarise(total = sum(prev, na.rm = TRUE), .groups = "drop")
  
  if (!is.null(preranked_blocks)) {
    valid_blocks  <- unique(block_prev_sub$block_r)
    ranked_blocks <- preranked_blocks[preranked_blocks %in% valid_blocks]
  } else {
    ranked_blocks <- block_scores %>%
      arrange(desc(.data[[rank_col]])) %>%
      pull(block_r)
  }
  
  cumulative   <- tibble(pathogen = pathogens_in_combo, cum_burden = 0)
  curve_rows   <- list()
  
  for (i in seq_along(ranked_blocks)) {
    blk       <- ranked_blocks[i]
    this_block <- block_prev_sub %>%
      filter(block_r == blk, pathogen %in% pathogens_in_combo) %>%
      dplyr::select(pathogen, prev)
    cumulative <- cumulative %>%
      left_join(this_block, by = "pathogen") %>%
      mutate(cum_burden = cum_burden + replace_na(prev, 0)) %>%
      dplyr::select(pathogen, cum_burden)
    
    pct_reached <- cumulative %>%
      left_join(total_burden, by = "pathogen") %>%
      mutate(pct = cum_burden / total, n_blocks = i)
    
    curve_rows[[i]] <- pct_reached
  }
  
  bind_rows(curve_rows)
}

# Example: first combo, Shannon ranking
example_combo  <- combos[[1]]
example_label  <- paste(sort(example_combo), collapse = " + ")

block_prev_ex  <- block_prev %>% filter(pathogen %in% example_combo)

block_wide_ex  <- block_prev_ex %>%
  pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0)
sub_cols_ex    <- setdiff(names(block_wide_ex), "block_r")

block_scores_ex <- block_wide_ex %>%
  rowwise() %>%
  mutate(
    shannon       = shannon_entropy(c_across(all_of(sub_cols_ex))),
    alpha_div     = alpha_diversity(c_across(all_of(sub_cols_ex))),
    rao_quadratic = rao_quadratic(c_across(all_of(sub_cols_ex)))
  ) %>%
  ungroup()

# Diversity-based curves
curve_data_diversity <- map_dfr(rank_methods, function(rm) {
  cumulative_coverage_curve(
    block_prev_sub     = block_prev_ex,
    block_scores       = block_scores_ex,
    rank_col           = rm,
    pathogens_in_combo = example_combo
  ) %>% mutate(rank_method = recode(rm, !!!method_labels))
})

# Comparator curves (pre-ranked orderings)
curve_data_comparators <- map_dfr(
  list(measles_order = measles_ordering, wealth_order = wealth_ordering),
  function(ordering) {
    cumulative_coverage_curve(
      block_prev_sub     = block_prev_ex,
      block_scores       = block_scores_ex,
      rank_col           = NULL,
      pathogens_in_combo = example_combo,
      preranked_blocks   = ordering
    )
  },
  .id = "rank_method"
) %>%
  mutate(rank_method = recode(rank_method, !!!method_labels))

curve_data <- bind_rows(curve_data_diversity, curve_data_comparators)

p_curve <- curve_data %>%
  ggplot(aes(x = n_blocks, y = pct, colour = pathogen, linetype = rank_method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(curve_data$n_blocks) * 0.02, y = 0.82,
           label = "80% threshold", hjust = 0, size = 3, colour = "grey40") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = paste("Cumulative coverage curve —", example_label),
    subtitle = "Blocks added in ranked order (high → low diversity)",
    x        = "Number of blocks included",
    y        = "Cumulative % of total disease burden captured",
    colour   = "Pathogen",
    linetype = "Ranking metric"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")
p_curve


# ---------------------------------------------------------------
# 11. SPATIAL AUTOCORRELATION (MORAN'S I) FOR RAO SCORES
#     One Moran's I per 3-pathogen combo, using block-level Rao
#     scores as the outcome and inverse-distance spatial weights.
# ---------------------------------------------------------------

library(spdep)

# Block-level GPS: one row per block_r with median coordinates
block_coords <- bangl_test %>%
  distinct(block_r, med_qgpslong, med_qgpslat) %>%
  drop_na()

# For each combo, compute Rao per block then run Moran's I
morans_results <- map_dfr(combos, function(combo) {
  combo_label <- paste(sort(combo), collapse = " + ")
  
  # Block-level Rao scores for this combo
  block_prev_sub <- block_prev %>% filter(pathogen %in% combo)
  
  block_rao <- block_prev_sub %>%
    pivot_wider(names_from = pathogen, values_from = prev, values_fill = 0) %>%
    rowwise() %>%
    mutate(rao = rao_quadratic(c_across(-block_r))) %>%
    ungroup() %>%
    dplyr::select(block_r, rao)
  
  # Merge with GPS, keep only blocks present in both
  spatial_df <- block_rao %>%
    inner_join(block_coords, by = "block_r") %>%
    drop_na()
  
  if (nrow(spatial_df) < 4) {
    return(tibble(combo = combo_label, moran_i = NA_real_,
                  moran_p = NA_real_, moran_z = NA_real_))
  }
  
  coords_mat <- as.matrix(spatial_df[, c("med_qgpslong", "med_qgpslat")])
  
  # Inverse-distance spatial weights (k=5 nearest neighbours as candidates,
  # then weight by 1/distance so nearby blocks contribute more)
  knn      <- knearneigh(coords_mat, k = min(5, nrow(spatial_df) - 1))
  nb       <- knn2nb(knn)
  dists    <- nbdists(nb, coords_mat)
  inv_dist <- lapply(dists, function(d) 1 / d)
  lw       <- nb2listw(nb, glist = inv_dist, style = "W", zero.policy = TRUE)
  
  mt <- moran.test(spatial_df$rao, lw, zero.policy = TRUE)
  
  tibble(
    combo    = combo_label,
    moran_i  = as.numeric(mt$estimate["Moran I statistic"]),
    moran_p  = mt$p.value,
    moran_z  = as.numeric(mt$statistic)
  )
})

print(morans_results)

# ---------------------------------------------------------------
# 12. EFFICIENCY VS MORAN'S I PLOT
#     Join Moran's I onto the Rao-only rows of results, then
#     scatter efficiency (n_blocks) against spatial clustering.
# ---------------------------------------------------------------

rao_efficiency <- results %>%
  filter(rank_method == "rao_quadratic") %>%
  dplyr::select(combo, n_blocks)

efficiency_morans <- rao_efficiency %>%
  left_join(morans_results, by = "combo") %>%
  drop_na(n_blocks, moran_i)

p_morans <- efficiency_morans %>%
  ggplot(aes(x = moran_i, y = n_blocks)) +
  geom_point(size = 3.5, alpha = 0.85, colour = "#4472C4") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linewidth = 0.7, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = combo), size = 2.8, max.overlaps = 10) +
  labs(
    title    = "Spatial clustering vs intervention efficiency (Rao ranking)",
    subtitle = "Higher Moran's I = more spatially clustered Rao burden;\nfewer blocks needed if clustering aligns with high-burden areas",
    x        = "Moran's I (Rao quadratic entropy)",
    y        = "Blocks required for ≥80% coverage"
  ) +
  theme_bw(base_size = 11)


p_morans <- efficiency_morans %>%
  ggplot(aes(x = moran_i, y = ((80 - n_blocks)/80)*100)) +
  geom_point(size = 3.5, alpha = 0.85, colour = "cadetblue") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linewidth = 0.7, linetype = "dashed", alpha = .2) +
  ggrepel::geom_text_repel(aes(label = combo), size = 3.8, max.overlaps = 10) +
  labs(
 #   title    = "Spatial clustering vs intervention efficiency (Rao ranking)",
  #  subtitle = "Higher Moran's I = more spatially clustered Rao burden;\nfewer blocks needed if clustering aligns with high-burden areas",
    x        = "Moran's I",
    y        = "Percent reduction in clusters needed\nto hit target (compared to random) (%)"
  ) +
  theme_bw() +
  scale_color_viridis_d(option = "D" , end = .7 , name = "strategy") +
  theme(
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11), 
    legend.position = "bottom") 

p_morans








