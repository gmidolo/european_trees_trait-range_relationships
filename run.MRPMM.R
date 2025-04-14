# Description: Run multi-response phylogenetic mixed models (MR–PMMs) 
# Tested on R version: 4.3.1
# Date: 11 March 2024
# Author: Gabriele Midolo 
# e-mail: midolo@fzp.czu.cz
# Reference: Midolo, G. (2024). Plant functional traits couple with range size and shape in European trees. Global Ecology and Biogeography.
#------------------------------------------------------------------------------------------------------------------------------------------#

# The following script is based on the work of Ben Halliwell, Luke A. Yates, and Barbara R. Holland.
# Check their TUTORIAL here: https://benjamin-halliwell.github.io/MR-PMM/MR-PMM_tutorial.html 
# and their preprint PAPER here: https://doi.org/10.1101/2022.12.13.520338 

#1. Load packages and prepare the data ####
set.seed(22)
library(tidyverse); library(ape); library(brms)

#Calulcate the correlation between plant height and range area
df <- data.frame(
  SPECIES = c('Abies alba', 'Abies borisii-regis', 'Abies cephalonica', 'Abies cilicica', 'Abies nordmanniana', 'Abies pinsapo', 'Acer campestre', 'Acer heldreichii', 'Acer monspessulanum', 'Acer platanoides', 'Acer pseudoplatanus', 'Aesculus hippocastanum', 'Alnus cordata', 'Alnus glutinosa', 'Alnus incana', 'Alnus alnobetula', 'Arbutus unedo', 'Betula pendula', 'Betula pubescens', 'Buxus balearica', 'Buxus sempervirens', 'Carpinus betulus', 'Carpinus orientalis', 'Castanea sativa', 'Cedrus libani', 'Celtis australis', 'Cornus mas', 'Cornus sanguinea', 'Corylus avellana', 'Euonymus europaeus', 'Fagus sylvatica', 'Frangula alnus', 'Fraxinus angustifolia', 'Fraxinus excelsior', 'Fraxinus ornus', 'Ilex aquifolium', 'Juglans regia', 'Juniperus communis', 'Juniperus drupacea', 'Juniperus excelsa', 'Juniperus foetidissima', 'Juniperus oxycedrus', 'Juniperus phoenicea', 'Juniperus thurifera', 'Larix decidua', 'Liquidambar orientalis', 'Olea europaea', 'Ostrya carpinifolia', 'Picea abies', 'Picea omorika', 'Picea orientalis', 'Pinus brutia', 'Pinus cembra', 'Pinus halepensis', 'Pinus heldreichii', 'Pinus mugo', 'Pinus nigra', 'Pinus peuce', 'Pinus pinaster', 'Pinus pinea', 'Pinus sylvestris', 'Pistacia atlantica', 'Pistacia lentiscus', 'Pistacia terebinthus', 'Platanus orientalis', 'Populus alba', 'Populus nigra', 'Populus tremula', 'Prunus avium', 'Prunus padus', 'Prunus spinosa', 'Quercus cerris', 'Quercus coccifera', 'Quercus faginea', 'Quercus frainetto', 'Quercus ilex', 'Quercus petraea', 'Quercus pubescens', 'Quercus pyrenaica', 'Quercus robur', 'Quercus suber', 'Quercus trojana', 'Salix alba', 'Salix caprea', 'Salix eleagnos', 'Sambucus nigra', 'Sorbus aria', 'Sorbus aucuparia', 'Sorbus domestica', 'Sorbus torminalis', 'Taxus baccata', 'Tetraclinis articulata', 'Tilia cordata', 'Tilia platyphyllos', 'Tilia tomentosa', 'Ulmus glabra', 'Ulmus laevis', 'Ulmus minor'),
  PLANT.HEIGHT = c(6.465,5.477,6.325,5.477,6.259,5.477,3.731,3.873,3.254,4.782,5.266,4.732,4.032,4.381,3.716,1.688,2.944,4.64,4.381,1.257,1.849,4.201,2.573,5.231,6.481,3.6,2.476,1.929,2.492,1.995,5.789,1.905,4.776,5.126,2.823,3.48,4.412,2.458,3.464,4.472,3.873,2.921,2.452,3.622,5.98,5.477,2.728,3.192,5.997,5.362,6.708,4.189,4.385,3.277,3.121,1.766,5.623,5,5.02,4.179,5.306,1.3,1.936,2,4.852,4.61,5.112,4.528,3.719,3.072,1.536,5.313,1.565,4.472,5.014,4.458,5.403,3.756,4.472,5.257,3.972,2.926,4.587,3.038,2.717,2.17,3.536,3.045,3.767,4.095,3.979,2.646,4.751,5.465,4.763,5.174,4.816,4.905),
  RANGE.AREA = c(58.609,19.206,19.396,27.555,38.932,14.251,82.144,43.094,68.2,86.193,71.636,20.345,26.835,91.61,100.216,101.997,64.611,112.959,104.989,33.151,49.425,77.608,54.926,48.84,33.512,64.68,71.724,84.369,88.933,82.225,75.277,97.49,81.475,87.753,61.964,69.003,84.791,111.544,29.588,66.075,43.14,72.279,62.965,41.492,41.693,20.847,69.018,57.519,80.593,11.533,33.434,43.89,38.024,50.975,20.381,46.85,52.579,15.987,48.655,40.825,102.269,78.063,62.349,69.644,54.872,100.962,91.811,114.642,78.551,108.543,89.709,63.645,61.917,45.999,54.088,62.754,77.281,70.533,48.224,87.163,49.997,41.121,97.16,100.37,60.889,85.877,85.195,110.686,66.53,72.476,66.336,45.731,90.209,72.794,52.633,88.311,84.59,85.689)
)
head(df)
#N.B.: plant height (= unit in meters) is square-root transformed to improve normality; range area was Box-Cox transformed (see SI online for the complete dataset).
#Plant height was retrieved from FloraVeg.EU (https://floraveg.eu/taxon/); 
#Check the manuscript for additional details on how trait data and range attributes were compiled

#Phylogenetic tree:
#N.B.: the phylogenetic tree was obtained from the ‘phylo.maker()’ function from the V.PhyloMaker R package (Jin & Qian, 2019; https://doi.org/10.1111/ecog.04434)
treetxt <- '((((((((Sambucus_nigra:102.692841,Ilex_aquifolium:102.692841)campanulids:4.04857,(Olea_europaea:17.739301,((Fraxinus_angustifolia:5.330423,Fraxinus_excelsior:5.330423):1.081153,Fraxinus_ornus:6.411576):11.327725):89.00211)mrcaott248ott320:5.599318,Arbutus_unedo:112.340729)mrcaott248ott650:2.225571,(Cornus_sanguinea:47.710427,Cornus_mas:47.710427):66.855873)mrcaott248ott27233:9.167937,((((((((Sorbus_domestica:2.706285,Sorbus_aucuparia:2.706285):5.893129,(Sorbus_torminalis:3.837911,Sorbus_aria:3.837911):4.761503):42.133854,((Prunus_spinosa:18.078495,Prunus_avium:18.078495):8.851155,Prunus_padus:26.92965):23.803618):48.228528,((Celtis_australis:79.170641,((Ulmus_minor:0.536514,Ulmus_glabra:0.536514):5.899536,Ulmus_laevis:6.43605):72.734592):6.322774,Frangula_alnus:85.493416):13.46838)Rosales.rn.d8s.tre:12.186212,(((((((Carpinus_betulus:8.226376,Carpinus_orientalis:8.226376):4.624557,Ostrya_carpinifolia:12.850933):19.198365,Corylus_avellana:32.049298):32.8096,(Betula_pendula:29.352098,Betula_pubescens:29.352098):35.5068):3.416551,(((Alnus_glutinosa:15.974943,Alnus_incana:15.974943):4.606542,Alnus_cordata:20.581485):18.979252,Alnus_alnobetula:39.560737):28.714712):22.373581,Juglans_regia:90.64903):7.230533,(((Quercus_robur:11.776698,Quercus_pubescens:11.776698,Quercus_ilex:11.776698,(((((Quercus_frainetto:3.966982,Quercus_faginea:3.966982):0.778652,Quercus_petraea:4.745634):0.424375,Quercus_pyrenaica:5.170009):3.507011,((Quercus_cerris:1.05467,Quercus_trojana:1.05467):0.52343,Quercus_suber:1.5781):7.09892):2.417941,Quercus_coccifera:11.094961):0.681738):0.532457,Castanea_sativa:12.309156):27.922252,Fagus_sylvatica:40.231407):57.648155):13.268445)mrcaott371ott2511:4.637557,((((Salix_caprea:0.308197,Salix_eleagnos:0.308197):7.034235,Salix_alba:7.342432):30.69623,((Populus_alba:1.574183,Populus_tremula:1.574183):1.537515,Populus_nigra:3.111698):34.926964):73.981924,Euonymus_europaeus:112.020586)mrcaott2ott1479:3.764979)mrcaott2ott371:2.793039,((Tilia_cordata:2.967151,(Tilia_platyphyllos:0.801068,Tilia_tomentosa:0.801068):2.166083):101.290908,(((((Acer_campestre:3.268696,Acer_platanoides:3.268696):9.197756,(Acer_heldreichii:6.709758,Acer_monspessulanum:6.709758):5.756694):1.629005,Acer_pseudoplatanus:14.095457):31.041674,Aesculus_hippocastanum:45.137131):34.78121,((Pistacia_atlantica:2.730903,Pistacia_terebinthus:2.730903):9.157292,Pistacia_lentiscus:11.888195):68.030146):24.339718)mrcaott96ott378:14.320545)mrcaott2ott96:3.826488,Liquidambar_orientalis:122.405092)mrcaott2ott2464:1.329145)Pentapetalae:4.81969,(Buxus_sempervirens:9.812122,Buxus_balearica:9.812122):118.741805)mrcaott2ott8379:1.770604,Platanus_orientalis:130.324531)mrcaott2ott969:194.725497,((((Juniperus_foetidissima:5.036479,((Juniperus_thurifera:1.964459,Juniperus_excelsa:1.964459):2.205496,Juniperus_phoenicea:4.169955):0.866525,((Juniperus_oxycedrus:3.777148,Juniperus_communis:3.777148):1.0901,Juniperus_drupacea:4.867248):0.169232):15.340157,Tetraclinis_articulata:20.376637):47.444059,Taxus_baccata:67.820696):59.002062,((((((((Pinus_mugo:4.876341,Pinus_sylvestris:4.876341):0.701058,Pinus_nigra:5.577399):3.316714,Pinus_heldreichii:8.894113):0.291606,((Pinus_pinea:4.318439,Pinus_pinaster:4.318439):2.865478,(Pinus_halepensis:2.227706,Pinus_brutia:2.227706):4.956211):2.001802):11.673147,(Pinus_cembra:5.599878,Pinus_peuce:5.599878):15.258988):14.157556,((Picea_abies:1.57471,Picea_omorika:1.57471):0.729596,Picea_orientalis:2.304306):32.712116):3.781157,Larix_decidua:38.797579):2.897999,((Abies_borisii-regis:4.548573,((((Abies_cilicica:0.692976,Abies_nordmanniana:0.692976):0.094111,Abies_cephalonica:0.787087):0.172634,Abies_pinsapo:0.959721):0.240809,Abies_alba:1.20053):3.348044):28.575012,Cedrus_libani:33.123586):8.571992):85.12718):198.22727)Spermatophyta;'
phy.tree <- read.tree(text = treetxt)

#Set names of phylogenetic tree and arrange the dataset accordingly
rownames(df) <- str_replace(df$SPECIES, ' ', '_') 
df <- df[match(phy.tree$tip.label, rownames(df)),]


#2. Run MR-PMMs ####

#2.1 Set the VCV matrix
C <- vcv.phylo(phy.tree, corr = T) 

#2.2 Run phylogenetic models

# rename traits:
names(df)[names(df) %in% 'PLANT.HEIGHT'] <- 'yi1'
names(df)[names(df) %in% 'RANGE.AREA'] <- 'yi2'
yi1 <- df[['yi1']]
yi2 <- df[['yi2']]

# fit MR-PMM:
#!!! N.B.: It may takes up to ~7 minutes to run !!!#
st <- Sys.time()

df$sp <- rownames(df)

model.formula <- bf(mvbind(yi2, yi1) ~ (1|a|gr(sp, cov = C))) + set_rescor(TRUE)

MRPMM <- brm(model.formula,
             family = gaussian(),
             data   = df, 
             data2  = list(C = C),
             control = list(adapt_delta = 0.99),
             cores  = 4)
et <- Sys.time()-st
message(paste0('DONE. Total elapsed time: ',paste(round(et,2)),' ',attr(et,'units')))


#3. Get model posterior stats####

#3.1 get 95% Credible intervals
median.ci95 <- as.data.frame(posterior_summary(MRPMM, probs = c(0.025, 0.975)))[c(5,8),]
median.ci95$dim <- ifelse(str_detect(rownames(median.ci95),'rescor'), 'ind', 'phy')

#3.2 get 50% Credible intervals
median.ci50 <- as.data.frame(posterior_summary(MRPMM, probs = c(0.25, 0.70)))[c(5,8),] %>% dplyr::select(-Estimate, -Est.Error)
median.ci50$dim <- ifelse(str_detect(rownames(median.ci50),'rescor'), 'ind', 'phy')

#3.3 inspect the results
res <- median.ci95 %>% left_join(median.ci50,'dim') %>%
  dplyr::rename(median=Estimate) %>%
  dplyr::select(median, dim, Q2.5, Q97.5, Q25, Q70)
res


#3.4 visualize the results
ggplot() +
  geom_vline(xintercept = 0, col = 'black', lty = 'longdash') +
  geom_linerange(data = res, aes(y=dim, xmin=Q2.5, xmax=Q97.5, col=dim),linewidth =.7) +
  geom_pointrange(data = res, aes(y= dim,  x=median, xmin=Q25, xmax=Q70, col=dim), 
                  size=.8, shape=21, linewidth =2, fill='black') +
  theme_bw() + 
  labs(y='', col='Corr.', x='Estimated correlation') + 
  lims(x=c(-1,1))
  

#4. Posterior Predictive Checks####
pp_check(MRPMM, resp= 'yi1', ndraws = 100) + labs(x='PLANT HEIGHT')
pp_check(MRPMM, resp= 'yi2', ndraws = 100) + labs(x='RANGE AREA')