#### 1. Get species of animals with occurrence in the Atlantic Forest ####
#Author: Weverton Trindade
#Last modification: May 30 2025

rm(list=ls())

#Load packages
#remotes::install_github("wevertonbio/faunabr") #If necessary, install faunabr
library(faunabr)
library(rgbif)
library(dplyr)
library(terra)
library(geobr)
library(sf)
library(data.table)

#Import Atlantic forest shapefile
af <- vect("Vectors/AF.gpkg")
#Import br shapefile
br <- read_state() %>% vect()
#Get states covered by AF
br_af <- terra::crop(br, af)
plot(br_af)
af_states <- br_af$abbrev_state %>% unique()

# #### Species list based on Fauna do Brasil ####
# #Get last version of faunabr
my_dir <- "Data/faunabr_data/"
# dir.create(my_dir, recursive = T)
# get_faunabr(output_dir = my_dir, data_version = "1.37") #Run this only once)
# #Using version 1.37


#Load data
bf <- load_faunabr(data_dir = my_dir, data_version = "1.37")

#Select only chordata with accepted names
# Only "Amphibia", "Aves", "Mammalia", and "Reptilia"
#In states covered by Atlantic Forest
faunabr::fauna_attributes(bf, attribute = "taxonomicStatus")

chord <- select_fauna(data = bf, phylum = "Chordata",
                      class = c("Amphibia", "Aves", "Mammalia", "Reptilia"),
                      states = af_states, filter_states = "in",
                      taxonomicStatus = "accepted",
                      origin = c("Native", "cryptogenic"))

#Save as a csv file
write.csv(chord, file = "Data/species_list_faunabr_AF.csv", row.names = TRUE)


#### Search species based on GBIF occurrences ####
#Species with occurrences inside Atlantic Forest

#Add a buffer of 50km in af
af50 <- buffer(af, width = 50 * 1000)

#Simplify geometry
af_simp <- simplifyGeom(af) %>% makeValid()
plot(af_simp)

#Get geometry
af_wkt <- af_simp %>% sf::st_as_sf() %>% sf::st_geometry(af_simp) %>% sf::st_as_text()
#Fix geometry
af_wkt <- af_wkt %>% wk::wkt() %>%
  wk::wk_orient()

#Extract species from Atlantic Forest
spp_gbif <- occ_download(pred_within(af_wkt), #Atlantic Forest
                         pred("taxonKey", 44), #Chordata
                         user = user,
                         pwd = pwd,
                         email = email,
                         format = "SPECIES_LIST")
#Save


occ_download_wait(spp_gbif) #Wait... It will take some minutes (maybe hours)...
#https://www.gbif.org/occurrence/download/0012627-250515123054153


spp_l <- occ_download_get(key = '0012627-250515123054153', #Paste the key here
                          path = "Data/GBIF_Species", overwrite = TRUE)
#Import list of species
d_gbif <- occ_download_import(spp_l, na.strings = c("", NA))

#Get only species with accepted names
head(d_gbif)
# seleciona apenas as entradas que estão classificadas no nível de espécie
spp_gbif <- d_gbif %>% filter(taxonRank == "SPECIES",
                              taxonomicStatus == "ACCEPTED")

#Subset classes
spp_gbif$class %>% table()

spp_gbif <- spp_gbif %>% filter(class %in% c("Amphibia", "Aves", "Mammalia",
                                             "Crocodylia",
                                             "Squamata", "Testudines")) #Reptillia

#Get unique species
spp_gbif <- spp_gbif %>% select(scientificName, numberOfOccurrences, class,
                                order, family, speciesKey,
                                iucnRedListCategory) %>% distinct()

#Get binomial name
spp_gbif <- spp_gbif %>% mutate(species = faunabr::extract_binomial(spp_gbif$scientificName),
                                .before = 1)

#Save
write.csv(spp_gbif, "Data/species_list_gbif_AF.csv")

####Matching dataset with faunabr####
#Impor data again
spp_gbif <- read.csv("Data/species_list_gbif_AF.csv")
spp_fauna <- read.csv("Data/species_list_faunabr_AF.csv")
my_dir <- "Data/faunabr_data/"
#Load data
bf <- load_faunabr(data_dir = my_dir, data_version = "1.37")

#Get species from gbif absent from faunapr
spp_out <- setdiff(spp_gbif$species, spp_fauna$species)
#Identifica as espécies que estão no spp_gbif mas não em chord

#From the species absent from faunapr, which one are present in faunabr?
sp_gbif_in_fauna <- subset_fauna(data = bf, species = spp_out)
table(sp_gbif_in_fauna$origin)

#Fix synonyms
sp_gbif_in_fauna$species[which(sp_gbif_in_fauna$taxonomicStatus == "synonym")] <- sp_gbif_in_fauna$acceptedName[which(sp_gbif_in_fauna$taxonomicStatus == "synonym")] %>% extract_binomial()

sp_gbif_in_fauna <- bf %>% filter(species %in% sp_gbif_in_fauna$species,
                                  taxonRank == "species") %>%
  filter(!(species %in% spp_fauna$species))

#Remove domesticaded and introduced species
sp_gbif_in_fauna$origin %>% table()
sp_gbif_in_fauna <- sp_gbif_in_fauna %>% filter(origin != "domesticated",
                                                origin != "introduced",
                                                origin != "exotic")
#Remove species already in faunapr
sp_gbif_in_fauna <- sp_gbif_in_fauna %>%
  filter(!(species %in% spp_fauna$species))

#Get species absent from faunabr
spp_out_br <- setdiff(spp_out, bf$species)

#Try to fix names
check_out <- florabr::match_names(species = spp_out_br,
                                  species_to_match = bf$species, parallel = TRUE,
                                  ncores = 8, progress_bar = TRUE)

#Get unique values
check_out3 <- check_out %>% distinct()

#Get only species with distance <= 3 a
check_out3 <- check_out %>% filter(Distance <= 3, !is.na(Suggested_name)) %>%
  distinct()

#Identify source of species in spp_fauna
spp_fauna <- spp_fauna %>% mutate(source = "faunabr")
#Add species in spp_fauna from gbif (species in faunabr as well)
spp_gbif_in_fauna <- subset_fauna(data = bf,
                                  species = check_out3$Suggested_name)
#Remove species already in spp_fauna
spp_gbif_in_fauna <- spp_gbif_in_fauna %>%
  filter(!(species %in% spp_fauna$species)) %>% filter(phylum == "Chordata")
#Identify source of species in spp_fauna
spp_gbif_in_fauna <- spp_gbif_in_fauna %>%
  mutate(source = "gbif_faunabr")

#Get species only from gbif
spp_only_gbif <- check_out %>%
  filter(!(input_name %in% check_out3$input_name)) %>% distinct()
spp_only_gbif <- spp_gbif %>% filter(species %in% spp_only_gbif$input_name) %>%
  distinct()
#Identify source of species
spp_only_gbif <- spp_only_gbif %>%
  mutate(source = "gbif_only")

#Create dataframe with all species
## First, get and rename columns
df_faunabr <- spp_fauna %>% select(species, class, order, family, source)
df_gbif_faunabr <- spp_gbif_in_fauna %>% select(species, class, order, family, source)
df_only_gbif <- spp_only_gbif %>% select(species, class, order, family, source)
df <- bind_rows(df_faunabr, df_gbif_faunabr, df_only_gbif) %>% distinct()

#Save data
write.csv(df, "Data/All_species_AF_unfiltered.csv", row.names = FALSE)

#Read data again
df <- read.csv("Data/All_species_AF_unfiltered.csv")

#Now, get all species and ask for some IA (here, gemini):
# - which specie are not native from south america
# - which species are fossils
spp <- df$species %>% unique()
spp %>% dput() #paste it in IA
#Extinct species accordind to Gemini
fossil_species <- c(
  "Pampatherium paulacoutoi", "Dusicyon avus", "Holmesina paulacoutoi",
  "Leontinia gaudryi", "Mirandatherium alipioi", "Asmithwoodwardia scotti",
  "Guggenheimia brasiliensis", "Nothrotherium maquinense", "Eutreptodactylus itaboraiensis",
  "Glyptodon clavipes", "Pristiguana brasilensis", "Smilodon populator",
  "Diogenornis fragilis", "Ahytherium aureum", "Rhynchippus brasiliensis",
  "Taubatherium paulacoutoi", "Derorhynchus singularis", "Australonyx aquae",
  "Peiropemys mezzalirai", "Monodelphopsis travassosi", "Guggenheimia crocheti",
  "Minusculodelphis minimus", "Minusculodelphis modicum", "Procaroloameghinia pricei",
  "Patene simpsoni", "Epidolops ameghinoi", "Riostegotherium yanei",
  "Catonyx cuvieri", "Carodnia vieirai", "Zeusdelphys complicatus",
  "Palaeocladosictis mosesi", "Megatherium americanum", "Victorlemoineia prototypica",
  "Itaboraidelphys camposi", "Protodidelphis vanzolinii", "Apodops pricei",
  "Eobrasilia coutoi", "Itaboravis elaphrocnemoides", "Pricemys caiera",
  "Gashternia carioca", "Protolipterna ellipsodontoides", "Noronhomys vespuccii",
  "Taubatherium major", "Paleopsilopterus itaboraiensis",
  "Dicolpomys fossor", "Eurycephalella alcinae",
  "Sarcosuchus hartti", "Lamegoia conodonta", "Eocaiman itaboraiensis",
  "Homalostylops atavus", "Glyptodon reticulatus", "Gaylordia mater",
  "Gaylordia macrocynodonta", "Ocnopus gracilis", "Carolocoutoia ferigoloi",
  "Arariphrynus placidoi"
)
#Get info and check manually
df_fossil <- df %>% filter(species %in% fossil_species)
#Castoria angustidens is not extinct
df_fossil <- df_fossil %>% filter(species != "Castoria angustidens")

#Check species non-native from south america
setdiff(spp, df_fossil$species) %>% dput()
#After searching in gemini...
non_native_species <- c(
  "Lobodon carcinophaga", # Foca-caranguejeira: Antártida e ilhas subantárticas
  "Mirounga leonina", # Elefante-marinho-do-sul: Oceanos do sul (Antártida, Argentina, Chile, África do Sul, Nova Zelândia, Austrália)
  "Arctocephalus tropicalis", # Lobo-marinho-subtropical: Ilhas do Atlântico Sul e Oceano Índico Sul
  "Arctocephalus gazella", # Lobo-marinho-antártico: Ilhas subantárticas e antárticas
  "Morus bassanus", # Ganso-patola-do-norte: Atlântico Norte
  "Morus capensis", # Ganso-patola-do-cabo: África do Sul, Namíbia, Moçambique
  "Morus serrator", # Ganso-patola-australiano: Austrália, Nova Zelândia
  "Chionis albus", # Pomba-das-neves: Península Antártica e ilhas subantárticas do Atlântico Sul
  "Glareola pratincola", # Corredor-de-asa-vermelha: Sul da Europa, África, sul da Ásia
  "Chroicocephalus ridibundos", # Gaivota-de-cabeça-preta: Europa, Ásia
  "Xema sabini", # Gaivota-de-sabine: Ártico
  "Calidris ferruginea", # Seixoeira-de-bico-curvo: Sibéria Ártica
  "Calidris minuta", # Maçarico-anão: Norte da Sibéria, Escandinávia
  "Calidris pugnax", # Combatente: Tundra do norte da Eurásia
  "Calidris subruficollis", # Maçarico-de-bico-grosso: Tundra do Alasca e Canadá
  "Limosa lapponica", # Fuselo: Tundra ártica da Eurásia e América do Norte
  "Numenius hudsonicus", # Maçarico-de-bico-curvo-americano: América do Norte
  "Phalaropus fulicarius", # Falaropo-de-bico-grosso: Tundra ártica
  "Phalaropus lobatus", # Falaropo-de-bico-fino: Ártico e subártico da América do Norte e Eurásia
  "Tringa glareola", # Maçarico-bique-bique: Norte da Eurásia
  "Tringa totanus", # Perna-vermelha-comum: Europa, Ásia
  "Xenus cinereus", # Maçarico-malhado: Sibéria, Ásia
  "Stercorarius maccormicki", # Stercorarius-de-mccormick: Antártida
  "Stercorarius parasiticus", # Moleiro-parasítico: Regiões árticas e subárticas do Hemisfério Norte
  "Stercorarius pomarinus", # Moleiro-pomarinho: Tundra ártica
  "Stercorarius skua", # Moleiro-grande: Atlântico Norte
  "Chlidonias leucopterus", # Gaivina-de-asa-branca: Europa Oriental, Ásia, África
  "Gygis alba", # Andorinha-do-mar-branca: Ilhas tropicais e subtropicais dos oceanos Pacífico, Índico e Atlântico
  "Sterna paradisaea", # Andorinha-do-mar-ártica: Regiões árticas e subárticas da Eurásia e América do Norte
  "Sterna vittata", # Andorinha-do-mar-antártica: Ilhas e costas subantárticas
  "Phaethon rubricauda", # Rabo-de-palha-de-cauda-vermelha: Oceanos tropicais e subtropicais
  "Oceanodroma castro", # Painho-de-castro: Atlântico tropical, Pacífico
  "Lugensa brevirostris", # Petrel-de-kermadec: Oceano Índico e Pacífico Sul
  "Pterodroma macroptera", # Petrel-de-cara-branca: Oceanos do sul (Austrália, Nova Zelândia)
  "Crex crex", # Codornizão: Europa, Ásia
  "Porphyrio alleni", # Frango-d'água-de-allen: África subsaariana
  "Porzana carolina", # Frango-d'água-americano: América do Norte
  "Pavo cristatus", # Pavão: Índia, Sri Lanka
  "Elephantulus brachyrhynchus", # Rato-elefante-de-nariz-curto: África Austral e Oriental
  "Lophura nycthemera", # Faisão-prateado: Sudeste Asiático (China, Laos, Vietnã)
  "Ovis dalli", # Ovelha-de-dall: Alasca, Yukon, Colúmbia Britânica
  "Lepus oiostolus", # Lebre-tibetan: Planalto Tibetano
  "Dromaius novaehollandiae", # Emu: Austrália
  "Acrochordus javanicus", # Cobra-javanense-de-água: Sudeste Asiático (Indonésia, Malásia, Tailândia)
  "Cygnus atratus", # Cisne-negro: Austrália
  "Anas platyrhynchos", # Pato-real: Hemisfério Norte (América do Norte, Europa, Ásia)
  "Troglodytes aedon", # Corruíra-comum: Américas (do Canadá à Patagônia)
  "Thalasseus sandvicensis", # Andorinha-do-mar-de-bico-preto: Costas da América do Norte, Europa, África
  "Dorcopsulus vanheurni", # Canguru-arborícola-de-van-heurn: Nova Guiné
  "Herpestes javanicus", # Mangusto-javanes: Sudeste Asiático
  "Taeniopygia guttata", # Diamante-mandarim: Austrália
  "Numida meleagris", # Galinha-d'angola: África
  "Psittacula krameri", # Periquito-de-colar: África, subcontinente indiano
  "Branta canadensis", # Ganso-canadense: América do Norte
  "Tragelaphus spekii", # Sitatunga: África Central
  "Eumetopias jubatus", # Leão-marinho-de-steller: Pacífico Norte
  "Anser albifrons", # Ganso-bravo-de-testa-branca: Norte da Eurásia e América do Norte
  "Xenopus laevis", # Rã-de-unhas-africana: África subsaariana
  "Alopochen aegyptiaca", # Ganso-do-nilo: África subsaariana, Vale do Nilo
  "Panthera pardus", # Leopardo: África, Ásia
  "Fregata ariel", # Fragata-pequena: Oceanos Índico e Pacífico
  "Georychus capensis", # Toupeira-rato-do-cabo: África do Sul
  "Myosorex varius", # Musaranho-variado: África do Sul
  "Oryctolagus cuniculus", # Coelho-europeu: Península Ibérica
  "Netta rufina", # Pato-vermelho: Sul da Europa, Ásia Central, Norte da África
  "Agapornis personatus", # Periquito-mascarado: Tanzânia
  "Sarkidiornis melanotos", # Pato-de-crista: África subsaariana, Madagascar, sul da Ásia
  "Quiscalus mexicanus", # Graúna-mexicana: América do Norte e Central
  "Dumetella carolinensis", # Tordo-cinzento: América do Norte
  "Cuon alpinus", # Dhole: Ásia Central e Sudeste Asiático
  "Microtus longicaudus", # Rato-do-campo-de-cauda-longa: América do Norte
  "Marmota bobak", # Marmota-das-estepes: Estepe da Europa Oriental e Ásia Central
  "Puffinus gravis", # Pardela-preta-de-bico-amarelo: Atlântico Sul
  "Sciurus ignitus", # Nome não encontrado/válido
  "Anomalurus derbianus", # Esquilo-voador-de-derby: África Ocidental e Central
  "Epomops buettikoferi", # Morcego-fruteiro-de-buettikofer: África Ocidental
  "Micrurus fulvius", # Cobra-coral-oriental: Sudeste dos Estados Unidos
  "Catharus ustulatus", # Tordo-de-swainson: América do Norte
  "Lepus capensis", # Lebre-do-cabo: África, Arábia, Índia
  "Scinax pixinguinha", # Perereca-pixinguinha: Brasil
  "Basileuterus vermivorus", # Mariquita-de-peito-amarelo: Leste da América do Norte
  "Tyto alba", # Coruja-das-torres: Cosmopolita
  "Eptesicus fuscus", # Morcego-pardo-grande: América do Norte e Central
  "Cynanthus latirostris", # Beija-flor-de-bico-largo: Sudoeste dos Estados Unidos, México
  "Aix sponsa", # Pato-mandarim-americano: América do Norte
  "Mertensophryne lindneri", # Sapo-de-lindner: África Central e Oriental
  "Hemidactylus turcicus", # Lagartixa-doméstica-mediterrânica: Bacia do Mediterrâneo
  "Cardellina canadensis", # Mariquita-canadense: Leste da América do Norte
  "Struthio camelus", # Avestruz: África
  "Anas aucklandica", # Marreca-de-auckland: Nova Zelândia
  "Chamaeleo dilepis", # Camaleão-comum-africano: África Central e Oriental
  "Platycercus adscitus", # Rosela-de-cabeça-azul: Nordeste da Austrália
  "Meleagris gallopavo", # Peru: América do Norte
  "Xenopus romeri", # Rã-de-garras-de-romer: África do Sul
  "Tadorna ferruginea", # Ganso-ferrugem: Eurásia, África
  "Aix galericulata", # Pato-mandarim: Leste da Ásia
  "Kobus megaceros", # Kob-do-nilo: Etiópia (pantanais de Gambela)
  "Myonycteris torquata", # Morcego-fruteiro-de-colar: África Ocidental e Central
  "Sigmodon hispidus", # Rato-algodoeiro: Sul dos Estados Unidos, México, América Central
  "Apus apus", # Andorinhão-preto: Europa, Ásia, África
  "Lepidodactylus lugubris", # Lagartixa-de-luto: Sudeste Asiático, Oceania
  "Setophaga adelaidae", # Mariquita-de-porto-rico: Porto Rico
  "Hemidactylus garnotii", # Lagartixa-de-garnot: Sudeste Asiático, Oceania
  "Lithobates catesbeianus", # Rã-touro-americana: Leste dos Estados Unidos e Canadá
  "Corvus brachyrhynchos", # Corvo-americano: América do Norte
  "Corvus macrorhynchos", # Corvo-de-bico-grosso: Ásia
  "Lasiurus borealis", # Morcego-vermelho-oriental: América do Norte
  "Puffinus griseus", # Pardela-escura: Oceanos do sul (Nova Zelândia, Chile, Ilhas Malvinas)
  "Anser anser", # Ganso-comum: Europa, Ásia
  "Melopsittacus undulatus", # Periquito-australiano: Austrália
  "Callorhinus ursinus", # Lobo-marinho-do-norte: Pacífico Norte
  "Synaptomys borealis", # Rato-do-campo-de-barro-boreal: América do Norte (Canadá e Alasca)
  "Pagodroma nivea", # Petrel-das-neves: Antártida e ilhas subantárticas
  "Miniopterus pusillus", # Morcego-de-asa-longa-pequeno: Sudeste Asiático, Indonésia, Oceania
  "Cygnus olor", # Cisne-branco: Europa, Ásia
  "Oxyura australis", # Pato-de-rabo-espinho-australiano: Austrália
  "Anser canagicus", # Ganso-imperador: Alasca, Sibéria
  "Trachylepis maculata", # Nome não encontrado/válido
  "Sorex vagrans", # Musaranho-errante: Oeste da América do Norte
  "Tyrannus couchii", # Suiriri-de-couch: Sul dos Estados Unidos, México, América Central
  "Gallinula chloropus", # Galinha-d'água-comum: Cosmopolita
  "Pseudomys gouldii", # Rato-de-gould: Austrália (extinto na natureza)
  "Turdus fulviventris", # Sabiá-de-barriga-amarela: América Central e do Sul
  "Hyla meridionalis", # Rã-arbórea-meridional: Sul da Europa, Norte da África
  "Charadrius asiaticus", # Borrelho-asiático: Ásia Central
  "Streptopelia decaocto", # Rola-turca: Europa, Ásia
  "Serinus canaria", # Canário: Ilhas Canárias, Açores, Madeira
  "Ardea herodias", # Garça-azul-grande: América do Norte e Central
  "Crotalus horridus", # Cascavel-do-leste: Leste dos Estados Unidos
  "Tringa nebularia", # Perna-verde-comum: Norte da Eurásia
  "Tiaris olivaceus", # Tiaris-oliváceo: Caribe, América Central
  "Euphagus carolinus", # Melro-de-rusty: América do Norte
  "Chrysolophus pictus", # Faisão-dourado: China
  "Estrilda troglodytes", # Estrilda-de-bochechas-vermelhas: África Ocidental
  "Mustela putorius", # Doninha-europeia: Europa, Norte da África
  "Camelus bactrianus", # Camelo-bactriano: Ásia Central
  "Procyon lotor", # Guaxinim: América do Norte
  "Colinus virginianus", # Codorna-virgínia: Leste da América do Norte, México
  "Nymphicus hollandicus", # Calopsita: Austrália
  "Microtus oeconomus", # Rato-do-campo-de-turfeira: Norte da Eurásia e América do Norte
  "Bandicota bengalensis", # Rato-toupeira-indiano: Sul e Sudeste Asiático
  "Peromyscus boylii", # Rato-de-boyle: América do Norte (México, sudoeste dos EUA)
  "Equus asinus", # Asno: África, Oriente Médio
  "Rhabdomys pumilio", # Rato-listrado-de-quatro-listras: África Austral
  "Garrulus glandarius", # Gaio-comum: Europa, Ásia, Norte da África
  "Phoenicopterus roseus", # Flamingo-comum: África, Ásia, sul da Europa
  "Charadrius hiaticula", # Borrelho-grande-de-coleira: Europa, Ásia, América do Norte
  "Chalcophaps indica", # Pomba-esmeralda-comum: Sul e Sudeste Asiático, Austrália
  "Pantherophis guttatus", # Cobra-do-milho: Sudeste dos Estados Unidos
  "Geopelia cuneata", # Diamante-pomba: Austrália
  "Vanellus vanellus" # Abibe-comum: Europa, Ásia
)
#Subset
df_non_native <- df %>% filter(species %in% non_native_species)

#Get species to remove
to_remove <- c(df_non_native$species,
               df_fossil$species) %>% unique()

#Remove from df
df_final <- df %>% filter(!(species %in% to_remove))

#Save
write.csv(df_final, "Data/All_AF_species.csv", row.names = FALSE)

#### Compare species with species in napibio ####
#rm(list=ls())
#Import data again
df_final <- read.csv("Data/All_AF_species.csv")
#Import species list from napibio
napibio <- fread("../napibio/Data/Animals/n_records_in_gbif.gz")
#Remove fossils and non-native species
to_remove <- c(fossil_species, non_native_species)
napibio_to_remove <- napibio %>% filter(species %in% to_remove)
#Remove
napibio <- napibio %>% filter(!(species %in% napibio_to_remove$species))

#See species to search (not in napibio)
to_search <- df_final %>% filter(!(species %in% napibio$species))
#Save
write.csv(to_search, "Data/species_to_search_AF.csv", row.names = FALSE)

#### Standardize names and get synonyms using gbif ####
library(rgbif)
library(faunabr)
library(pbapply)
library(dplyr)
library(data.table)

#Get species list
d <- read.csv("Data/species_to_search_AF.csv")
spp <- d$species %>% unique()

#Get synonyms from faunabr
#Load data
my_dir <- "Data/faunabr_data/"
bf <- load_faunabr(my_dir)

all_syn <- fauna_synonym(data = bf, species = spp,
                         include_subspecies = FALSE) %>%
  filter(!is.na(synonym))

#Ranks to consider
r <- c("VARIETY", "SPECIES", "FORM", "SUBSPECIES")

#Create cluster
cl <- parallel::makeCluster(8)
parallel::clusterEvalQ(cl, {library("rgbif");
  library("dplyr");
  library(data.table);
  library(faunabr)})
parallel::clusterExport(cl, varlist = c("spp", "r", "all_syn"))
d_gbif <- pblapply(spp, function(sp){

  #For test
  #sp = spp[1]
  #print(sp)

  #Get synonyms
  sp_fauna_syn <- all_syn %>% filter(acceptedName %in% sp)
  sp_fauna_syn <- sp_fauna_syn$synonym
  all_sp <- unique(c(sp, sp_fauna_syn))

  #Get name backbone
  sp_backbone <- rbindlist(lapply(all_sp, rgbif::name_backbone,
                                  kingdom = "Animalia"), fill = TRUE)
  #Subset ranks
  is_rank <- "rank" %in% colnames(sp_backbone)
  is_species <- if(is_rank) {
    any(r %in% sp_backbone$rank)} else {FALSE}

  #Try with name_suggest
  if(!is_rank | !is_species){
    #Get name backbone
    sp_backbone <- lapply(all_sp, rgbif::name_suggest,
                          field = c('key','canonicalName', "kingdom", "rank", "class"))
    sp_backbone <- lapply(sp_backbone, function(x) x$data)
    sp_backbone <- unique(rbindlist(sp_backbone))
    if(nrow(sp_backbone)>0){
      sp_backbone <- sp_backbone %>% dplyr::rename(speciesKey = key,
                                                   species = canonicalName)
      #Subset ranks
      is_rank <- "rank" %in% colnames(sp_backbone)
      is_species <- if(is_rank) {
        any(r %in% sp_backbone$rank)} else {FALSE}
    }
  }

  if(is_rank & is_species) {
    sp_backbone <- sp_backbone %>% filter(rank %in% r)

    sp_synonyms <- lapply(sp_backbone$speciesKey, name_usage)
    syn_gbif <- rbindlist(lapply(sp_synonyms, function(x) x$data), fill = TRUE)
    if("basionymKey" %in% colnames(syn_gbif)){
      sp_syns <- extract_binomial(syn_gbif$basionym)
    } else {sp_syns <- NA}
    #Remove wrong syns with "?" character
    sp_syns[grepl("\\?", sp_syns)] <- NA

    if(length(sp_syns) > 1){
      sp_syns <- sp_syns[!is.na(sp_syns)]
      sp_syns <- paste(unique(sp_syns), collapse = ";")
    }

    #Check class
    sp_class <- unique(sp_backbone$class)
    if(is.null(sp_class)){
      sp_class <- NA
    }

    #Create dataframe with old and new name and synonyms
    df <- data.frame(species = sp,
                     accepted_gbif = paste(unique(sp_backbone$species), collapse = ";"),
                     synonyms = sp_syns,
                     class = sp_class)} else {
                       df <- data.frame(species = sp,
                                        accepted_gbif = NA,
                                        synonyms = NA,
                                        class = NA)
                     }
  return(df)
}, cl = cl)
parallel::stopCluster(cl)
d_gbif <- rbindlist(d_gbif)
#Bind column with source
df <- merge(d_gbif, d, by = "species")

#Identify synonyms gbif and synomyms faunabr
all_syn2 <- rbindlist(lapply(unique(all_syn$acceptedName), function(x){
  syn_x <- all_syn %>% filter(acceptedName == x)
  data.frame(species = x,
             synonym_faunabr = paste(syn_x$synonym, collapse = ";"))
}))

#Merge data
df2 <- df %>% rename(synonyms_gbif = synonyms) %>%
  left_join(., all_syn2, by = "species")

#Remove duplicates species
spp_dup <- df2$accepted_gbif[duplicated(df2$accepted_gbif)] %>%
  na.omit()
df2 %>% filter(accepted_gbif %in% spp_dup) %>% View()
unique(df2$source)
preferred_order <- c("faunabr", "gbif_faunabr", "gbif_only")
#Keep only unique values
df3 <- df2 %>%
  group_by(accepted_gbif) %>%
  arrange(match(source, preferred_order)) %>%
  filter(row_number() == 1) %>%
  ungroup()
#Check
df3$accepted_gbif[duplicated(df3$accepted_gbif)] %>%
  na.omit() #Should be 0
#Check if there are synonums in gbif as accepted gbif
sp_dup2 <- intersect(df3$synonyms_gbif, df3$accepted_gbif) %>% na.omit()
#Fix manually: they are not synonyms
df3$synonyms_gbif[df3$synonyms_gbif == "Micrurus pyrrhocryptus"] <- NA
df3$synonyms_gbif[df3$synonyms_gbif == "Callicebus personatus"] <- NA
unique(df3$accepted_gbif) %>% length()

#Fix duplicates species: synonyms
all_syn_gbif <- strsplit(df3$synonyms_gbif, ";") %>% unlist() %>%
  na.omit()
spp_dup <- all_syn_gbif[duplicated(all_syn_gbif)]
df3 %>% filter(grepl(paste(spp_dup, collapse = "|"), synonyms_gbif)) %>%
  View()
#Penelope araucuan is synonyms only of Ortalis araucuan, not of Ortalis guttata
df3$synonyms_gbif[df3$species == "Ortalis guttata"] <- "Penelope guttata"

#Fill and fix NA classes manually
df4 <- df3 %>% select(-class.x) %>% rename(class = class.y) #Fix column names
df4$species[is.na(df4$class)]
#Check classes
df4 %>% count(class)
#Remove fish
df4 <- df4 %>% filter(class != "Actinopterygii")

#Get key number and check number of records in GBIF
spp <- df4$species
spp %>% as.data.frame() %>% View()
#Ranks to consider
r <- c("VARIETY", "SPECIES", "FORM", "SUBSPECIES")

#Check keys and number of records in GBIF
cl <- parallel::makeCluster(8)
parallel::clusterEvalQ(cl, {library("rgbif");
  library("dplyr");
  library("faunabr");
  library("data.table")})
parallel::clusterExport(cl, varlist = c("spp", "r", "df4"))

sp_gbif <- pblapply(spp, function(i){
  #try({
  #print(i)
  #Get all synonyms
  spp_i <- df4 %>% filter(species == i)
  all_sp <- spp_i %>% dplyr::select(-source) %>%
    as.character() %>% na.omit() %>% unique()

  #Get name backbone
  sp_backbone <- rbindlist(lapply(all_sp, rgbif::name_backbone,
                                  kingdom = "Animalia"), fill = TRUE)
  #Subset ranks
  is_rank <- "rank" %in% colnames(sp_backbone)
  is_species <- if(is_rank) {
    any(r %in% sp_backbone$rank)} else {FALSE}

  if(is_rank & is_species) {
    sp_backbone <- sp_backbone %>% filter(rank %in% r)}

  #Try with name_suggest
  if(!is_rank | !is_species){
    #Get name backbone
    sp_backbone <- lapply(all_sp, rgbif::name_suggest,
                          field = c('key','canonicalName', "kingdom", "rank"))
    sp_backbone <- lapply(sp_backbone, function(x) x$data)
    sp_backbone <- unique(rbindlist(sp_backbone))
    if(nrow(sp_backbone)>0){
      sp_backbone <- sp_backbone %>% dplyr::rename(speciesKey = key,
                                                   species = canonicalName)
      #Subset ranks
      is_rank <- "rank" %in% colnames(sp_backbone)
      is_species <- if(is_rank) {
        any(r %in% sp_backbone$rank)} else {FALSE}
      #Check if species is in the backbone
      if(is_species){
        is_species <- i %in% sp_backbone$species
      }
    }
  }

  #Get occurrences
  #First, check which columns with key is present
  k <- intersect(c("usageKey", "key"), colnames(sp_backbone))
  #Look for more synonyms
  if(is_rank & is_species & nrow(sp_backbone) > 0) {
    keys <- unique(c(sp_backbone[[k]]))
    #Get number of occurrence
    n_occ <- rgbif::occ_count(taxonKey = paste(keys, collapse = ";"))
    n_occ_with_coordinates <- rgbif::occ_count(taxonKey = paste(keys,
                                                                collapse = ";"),
                                               hasCoordinate = TRUE)
    #Get dataframe
    keys_species <- rbindlist(lapply(keys, function(x) name_usage(x)$data),
                              fill = TRUE)$canonicalName

    df <- data.frame(species = i,
                     occ_total = n_occ,
                     occ_with_coordinates = n_occ_with_coordinates,
                     keys = paste(keys, collapse = ";"),
                     species_keys = paste(keys_species, collapse = ";"))
  } else {
    df <- data.frame(species = i,
                     occ_total = NA,
                     occ_with_coordinates = NA,
                     keys = NA)
  }
  return(df)
  #})
}, cl = cl)
parallel::stopCluster(cl)
names(sp_gbif) <- spp

#Remove errors
sp_gbif <- sp_gbif[sapply(sp_gbif, function(x) is.data.frame(x))]
df_gbif <- rbindlist(sp_gbif, fill = TRUE)
#Absent species
spp_out <- setdiff(spp, unique(df_gbif$species)) #Must be 0
#Try this one
#i = spp_out
#df_gbif <- bind_rows(df_gbif, df)

#Check species with a lot of records again with IA
df_gbif$species %>% unique() %>% dput()


#Remove species that we know is introduced
introduced <- c(
    "Agapornis fischeri", # Periquito-de-fischer, África Oriental
    "Agapornis roseicollis", # Periquito-de-colar-rosa, África Austral
    "Anser cygnoides", # Ganso-cisne, Ásia Oriental
    "Aptenodytes patagonicus", # Pinguim-rei, Ilhas Subantárticas
    "Basiliscus vittatus", # Basilisco-comum, México e América Central
    "Bubulcus ibis", # Garça-vaqueira, Europa, África e Ásia
    "Cavia porcellus", # Porquinho-da-índia (doméstico), Andes (domesticado a partir de espécies selvagens)
    "Chenonetta jubata", # Pato-jubado, Austrália
    "Corvus splendens", # Corvo-doméstico, Ásia
    "Cuculus canorus", # Cuco-comum, Europa, Ásia e África
    "Egretta garzetta", # Garça-pequena, Europa, África, Ásia e Oceania
    "Egretta gularis", # Garça-dos-recifes, África Ocidental e Ásia
    "Equus caballus", # Cavalo (doméstico), Eurásia (domesticado)
    "Falco aesalon", # Esmerilhão, Europa, Ásia e África
    "Falco tinnunculus", # Peneireiro-vulgar, Europa, Ásia e África
    "Felis catus", # Gato-doméstico, Oriente Médio (domesticado)
    "Gallinago gallinago", # Narceja-comum, Europa, Ásia e África
    "Gallus gallus", # Galinha (doméstica), Sudeste Asiático (domesticado)
    "Himantopus himantopus", # Pernilongo-de-costas-pretas, Europa, Ásia e África
    "Lagenorhynchus acutus", # Golfinho-roaz-de-laterais-brancas, Atlântico Norte
    "Lama glama", # Lhama (doméstica), Andes (domesticada)
    "Lanius senator", # Picanço-barreteiro, Europa, África e Ásia
    "Larus fuscus", # Gaivota-d'asa-escura, Europa e África
    "Lepus europaeus", # Lebre-comum, Europa e Ásia
    "Mabuya mabouya", # Escinco-das-ilhas, Ilhas do Caribe (nativa, mas frequentemente citada como introduzida em outros contextos ou com distribuição ambígua devido à similaridade com outras espécies)
    "Milvus migrans", # Milhafre-preto, Europa, África, Ásia e Austrália
    "Plegadis falcinellus", # Íbis-preta, Europa, África, Ásia e Austrália
    "Platalea leucorodia", # Colhereiro-comum, Europa, África e Ásia
    "Psittacula eupatria", # Periquito-alexandrino, Sul da Ásia
    "Spheniscus humboldti", # Pinguim-de-Humboldt, Costa do Pacífico da América do Sul (apesar de ser da América do Sul, a espécie foi listada como "não nativa das Américas" na sua lista inicial, e não é nativa de todas as Américas, apenas da região costeira do Peru e Chile)
    "Spheniscus magellanicus" # Pinguim-de-Magalhães, América do Sul (similar ao pinguim-de-Humboldt, nativo apenas da região austral da América do Sul)
  )

#Check introduced species
introduced_df <- df_gbif %>% filter(species %in% introduced)

#Keep only natives
df_gbif_nat <- df_gbif %>% filter(!(species %in% introduced))
df5 <- df4 %>% filter(!(species %in% introduced))

#Merge info on class, order and family
df_gbif_nat <- df_gbif_nat %>% left_join(., df5[,c("species", "class", "order", "family")])

#Save temporary file
fwrite(df_gbif_nat, "Data/n_records_in_gbif.gz", compress = "gzip")

#Save table
write.csv(df5, "Data/Species_list_checked_gbif.csv", row.names = FALSE)

table(df_gbif_nat$class)
df_gbif_nat %>% group_by(class) %>% summarise(n = sum(occ_with_coordinates, na.rm = TRUE))
df_gbif_nat %>% group_by(order) %>% summarise(n = sum(occ_with_coordinates, na.rm = T)) %>% View()
