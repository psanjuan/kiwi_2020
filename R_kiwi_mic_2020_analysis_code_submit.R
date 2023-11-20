############################################################################
#   Author:       Priscilla A San Juan
#   Topic:        Microbiome analysis of kiwi gut microbiome - edited script
#   Manuscript:   "Kiwi development through a microbial lens: How time, 
#                 local environment, and disease history influence the 
#                 Brown Kiwi (Apteryx mantelli) gut microbiome"

############################################################################

#       Table of contents   
# -------------------------------------------------------------------------
# 1  -  Load packages and color vectors 
# 2  -  Load files
# 3  -  Decontam
# 4  -  Subset
# 5  -  Alpha diversity
# 6  -  Microbial Abundance 
# 7  -  Ordination
# 8  -  PermANOVA
# 9  -  Beta diversity
# 10 -  clamtest
# 11 -  Temporal betadiversity 

############################################################################

# 1 LOAD PACKAGES  
# -------------------------------------------------------------------------
req_pkg <- c("readr","microbiome","dplyr","phyloseq","vegan","ggplot2","utils", 
             "microbiomeutilities", "DESeq2","ape","forcats","DescTools","scales",
             "tidyr","ResourceSelection", "gridExtra","ggbeeswarm","mvabund",
             "RColorBrewer","randomcoloR","MASS", "viridis","ggpubr","decontam",
             "ggbeeswarm","fantaxtic","plotly", "microbial") 

# Load all required packages and show version
for(i in req_pkg){
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}
library(SRS)
library(performance)
library(ade4)        # aux
library(agricolae)   # aux
library(lmerTest)    # aux

# Color palettes ----------------------------------------------------------
n<-22
z<-34
qual_col_pals=brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector=unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
palette_qual <- distinctColorPalette(n)
pie(rep(1, n), col=palette_qual)
palette_qual_all <- distinctColorPalette(z)
pie(rep(1, z), col=palette_qual_all)
colors_scheme_life_stage<-c("Hatch room"="#a1dab4","Brooder room"="#41b6c4","Runs"="#225ea8")
colors_scheme_what<-c("soil"="#fdb462","kiwi poo"="#b3de69")
color_scheme_cohort<-c("young"="#e29580", "old"="#9f5a47")
color_scheme_cocc <- c("positive"="#d95f02","negative"="#1b9e77")
theme_set(theme_light())

###########################################################################
# 2 LOAD FILES 
# -------------------------------------------------------------------------
bacteria_otu_matrix<-as.matrix((read_csv("data/bac_otu.csv", col_names=T)[,-1])) 
bacteria_tax_matrix<-as.matrix((read_csv("data/bac_tax.csv", col_names=T)[,-1])) 
# Add in rownames
rownames(bacteria_otu_matrix)<-paste0("OTU", 1:nrow(bacteria_otu_matrix))
rownames(bacteria_tax_matrix)<-paste0("OTU", 1:nrow(bacteria_tax_matrix))
# Reading the sample data file
sampledata_kiwi <- data.frame(read_csv("data/sample_datasheet_kiwi_2020.csv"))
# Naming the the first rows using the first column
row.names(sampledata_kiwi) <- sampledata_kiwi$samplename
# Add column to sample data - dates from hatch
sampledata_kiwi$daysfromhatch <- sampledata_kiwi$julian..collection.-sampledata_kiwi$julian..hatch.
# Add column to sample data - cohort category
sampledata_kiwi <- sampledata_kiwi %>%
  mutate(cohort = case_when(daysfromhatch < 670 ~ "young", 
                            daysfromhatch > 670 ~ "old"))
# Add column to sample data - cohort category
sampledata_kiwi <- sampledata_kiwi %>%
  mutate(agebin = case_when(age..days. <= 50 ~ "young", 
                            age..days. > 50 ~ "old"))
# Count how many samples in each category
table(sampledata_kiwi$cohort)
# Remove blank
sampledata_kiwi <- droplevels(sampledata_kiwi[!sampledata_kiwi$what. == 'blank',])
# Add column to sample data - control versus sample
sampledata_kiwi$sample_or_control<-forcats::fct_collapse(
  sampledata_kiwi$what., 
  Control=c("negative control"),
  Sample=c("kiwi poo","possum poo","soil","positive control"),
  group_other = T)

sampledata_kiwi$history.of.positive.results <- sampledata_kiwi$history.of.positive.results %>% replace_na('None')

# add column to sample_data
sampledata_kiwi$disease_history <- forcats::fct_collapse(
  sampledata_kiwi$history.of.positive.results, 
  History=c("diarrhea", "coccidiosis","coccidiosis, worms","operation, coccidiosis, worms",  
            "dehydrated","coccidiosis, worms, eyes", "emphysema, coccidiosis",      
            "coccidiosis, hernia", "emphysema, coccidiosis, worms", "dehydrated, failed weight gain",
            "coccidiosis, fungal"),
  None=c("None"),
  group_other = T)

# Read file into phyloseq
OTU_kiwi_bac<-otu_table(bacteria_otu_matrix, taxa_are_rows=T)
TAX_kiwi_bac<-tax_table(bacteria_tax_matrix)
SAM_kiwi_bac<-sample_data(sampledata_kiwi)
# Combining OTU, TAX, SAM data
bacteria_kiwi_2020<-phyloseq(OTU_kiwi_bac, TAX_kiwi_bac, SAM_kiwi_bac) # 9816 otus, 1137 samples

###########################################################################
# 3 DECONTAM: FILTER CONTAMS USING NEG CONTROL  
# -------------------------------------------------------------------------
# Inspect library sizes
df_bac <- as.data.frame(sample_data(bacteria_kiwi_2020))
df_bac$LibrarySize <- sample_sums(bacteria_kiwi_2020)
df_bac <- df_bac[order(df_bac$LibrarySize),]
df_bac$Index <- seq(nrow(df_bac))
tem_bac <- ggplot(data=df_bac, aes(x=Index, y=LibrarySize, color=sample_or_control))+geom_point(); tem_bac
# Use prevalence method to filter out suspected contaminants
sample_data(bacteria_kiwi_2020)$is.neg <- sample_data(bacteria_kiwi_2020)$sample_or_control == c("Control")
contamdf.prev.bac <- isContaminant(bacteria_kiwi_2020, method="prevalence", neg="is.neg", threshold = 0.4)
table(contamdf.prev.bac$contaminant)
head(which(contamdf.prev.bac$contaminant))
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa.bac <- transform_sample_counts(bacteria_kiwi_2020, function(abund) 1*(abund>0))
ps.pa.neg.bac <- prune_samples(sample_data(ps.pa.bac)$sample_or_control == "Control", ps.pa.bac)
ps.pa.pos.bac <- prune_samples(sample_data(ps.pa.bac)$sample_or_control == "Sample", ps.pa.bac)
# Make data.frame of prevalence in positive and negative samples
df.pa.bac <- data.frame(
  pa.pos=taxa_sums(ps.pa.pos.bac),pa.neg=taxa_sums(ps.pa.neg.bac),contaminant=contamdf.prev.bac$contaminant)
ggplotly(ggplot(data=df.pa.bac, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
           xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)"))
# Remove contams from phyloseq obj
bac.noncontam <- prune_taxa(!contamdf.prev.bac$contaminant, bacteria_kiwi_2020)
bac.noncontam # Use this phyloseq object for the following steps, 8738 otus, 1137 samples
# normalize first then subset
bac.noncontam.norm <- normalize(bac.noncontam, method="TMM") # 8738 otus, 1137 samples

###########################################################################
# Trimming data 
# -------------------------------------------------------------------------

# Remove chloroplast (because plant DNA)
class <- as.vector(data.frame(tax_table(bac.noncontam))$class)
class <- (!(class%in%c("Chloroplast")))
class[is.na(class)]=FALSE
Bac.data.prune_kiwi=prune_taxa(class, bac.noncontam)
# Remove chloroplast (because plant DNA) TMM
class <- as.vector(data.frame(tax_table(bac.noncontam.norm))$class)
class <- (!(class%in%c("Chloroplast")))
class[is.na(class)]=FALSE
Bac.data.prune_kiwi.TMM=prune_taxa(class, bac.noncontam.norm) # 8738 otus, 1137 samples
# Remove chloroplast (because plant DNA)
phyla <- as.vector(data.frame(tax_table(Bac.data.prune_kiwi))$phylum)
phyla <- (!(phyla%in%c("Cyanobacteria/Chloroplast")))
phyla[is.na(phyla)]=FALSE
Bac.data.prune_kiwi=prune_taxa(phyla, Bac.data.prune_kiwi) # 8616 otus, 1137 samples
# Remove chloroplast (because plant DNA) TMM
phyla <- as.vector(data.frame(tax_table(Bac.data.prune_kiwi.TMM))$phylum)
phyla <- (!(phyla%in%c("Cyanobacteria/Chloroplast")))
phyla[is.na(phyla)]=FALSE
Bac.data.prune_kiwi.TMM=prune_taxa(phyla, Bac.data.prune_kiwi.TMM) # 8616 otus, 1137 samples
# Remove controls
control <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi))$sample_or_control)
control <- (!(control%in%c("Control")))
control[is.na(control)]=FALSE
Bac.data.prune_kiwi=prune_samples(control, Bac.data.prune_kiwi) # 8616 otus, 1114 samples 
# Remove controls TMM
control <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi.TMM))$sample_or_control)
control <- (!(control%in%c("Control")))
control[is.na(control)]=FALSE
Bac.data.prune_kiwi.TMM=prune_samples(control, Bac.data.prune_kiwi.TMM) # 8616 otus, 1114 samples 
Bac.data.prune_kiwi.TMM.keepwild=Bac.data.prune_kiwi.TMM # 8616 otus, 1114 samples
# Remove wild samples
captive <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi))$wild.cap)
captive <- (!(captive%in%c("wild")))
captive[is.na(captive)] = FALSE
Bac.data.prune_kiwi = prune_samples(captive, Bac.data.prune_kiwi) # 8616 otus, 962 samples
# Remove wild samples TMM
captive <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi.TMM))$wild.cap)
captive <- (!(captive%in%c("wild")))
captive[is.na(captive)] = FALSE
Bac.data.prune_kiwi.TMM = prune_samples(captive, Bac.data.prune_kiwi.TMM) # 8616 otus, 962 samples
# Remove samples with less than 100 reads
Bac.data.prune_kiwi=prune_samples(sample_sums(Bac.data.prune_kiwi)>=100, Bac.data.prune_kiwi) # 8616 otus, 946 samples
Bac.data.prune_kiwi.TMM=prune_samples(sample_sums(Bac.data.prune_kiwi.TMM)>=100, Bac.data.prune_kiwi.TMM) # 8616 otus, 962 samples
# Remove OTUs with zero counts
Bac.data.prune_kiwi=prune_taxa(taxa_sums(Bac.data.prune_kiwi)>0, Bac.data.prune_kiwi) # 7547 otus, 946 samples
Bac.data.prune_kiwi.TMM=prune_taxa(taxa_sums(Bac.data.prune_kiwi.TMM)>0, Bac.data.prune_kiwi.TMM) # 8616 otus, 962 samples

###########################################################################
# 4 SUBSET DATA
# -------------------------------------------------------------------------
kiwi <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi))$what.)
kiwi <- kiwi%in%c("kiwi poo","soil")
kiwi[is.na(kiwi)] = FALSE
Bac.data.prune_kiwi = prune_samples(kiwi, Bac.data.prune_kiwi) # 7547 otus, 832 samples
# TMM
kiwi <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi.TMM))$what.)
kiwi <- kiwi%in%c("kiwi poo","soil")
kiwi[is.na(kiwi)] = FALSE
Bac.data.prune_kiwi.TMM = prune_samples(kiwi, Bac.data.prune_kiwi.TMM) # 8616 otus, 847 samples
# wild and cap TMM
kiwi_wild_cap <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi.TMM.keepwild))$wild.cap)
kiwi_wild_cap <- kiwi_wild_cap%in%c("cap","wild")
kiwi_wild_cap[is.na(kiwi_wild_cap)] = FALSE
Bac.data.prune_kiwi.TMM.keepwild = prune_samples(kiwi_wild_cap, Bac.data.prune_kiwi.TMM.keepwild) # 8616 otus, 1004 samples
# Subset by poo 
poo <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi))$what.)
poo <- poo%in%c("kiwi poo")
poo[is.na(poo)]=FALSE
poodata.bac <- prune_samples(poo, Bac.data.prune_kiwi) # 7547 otus, 752 samples
# Subset by poo TMM
poo <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi.TMM))$what.)
poo <- poo%in%c("kiwi poo")
poo[is.na(poo)]=FALSE
poodata.bac.TMM <- prune_samples(poo, Bac.data.prune_kiwi.TMM) # 8616 otus, 767 samples
# Subset by soil
soil <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi))$what.)
soil <- soil%in%c("soil")
soil[is.na(soil)]=FALSE
soildata.bac <- prune_samples(soil, Bac.data.prune_kiwi) # 7547 otus, 80 samples
# Subset by soil
soil <- as.vector(data.frame(sample_data(Bac.data.prune_kiwi.TMM))$what.)
soil <- soil%in%c("soil")
soil[is.na(soil)]=FALSE
soildata.bac.TMM <- prune_samples(soil, Bac.data.prune_kiwi.TMM) # 8616 otus, 80 samples
# Subset by younger birds 
youngerdata.bac <- poodata.bac %>%
  subset_samples(daysfromhatch<670) %>%
  subset_samples(name=!NA)
youngerdata.bac=prune_taxa(taxa_sums(youngerdata.bac)>0, youngerdata.bac) # 3815 otus, 381 samples
# Subset by younger birds TMM
youngerdata.bac.TMM <- poodata.bac.TMM %>%
  subset_samples(daysfromhatch<670) %>%
  subset_samples(name=!NA)
youngerdata.bac.TMM=prune_taxa(taxa_sums(youngerdata.bac.TMM)>0, youngerdata.bac.TMM) # 8616 otus, 389 samples
# Subset out Frosty (the oldest kiwi - outlier)
no_Frostdata.bac <- poodata.bac %>%
  subset_samples(daysfromhatch<7000) %>%
  subset_samples(name=!NA)
no_Frostdata.bac=prune_taxa(taxa_sums(no_Frostdata.bac)>0, no_Frostdata.bac) # 5510 otus, 748 samples
# Subset out Frosty (the oldest kiwi - outlier)
no_Frostdata.bac.TMM <- poodata.bac.TMM %>%
  subset_samples(daysfromhatch<7000) %>%
  subset_samples(name=!NA)
no_Frostdata.bac.TMM=prune_taxa(taxa_sums(no_Frostdata.bac.TMM)>0, no_Frostdata.bac.TMM) # 8616 otus, 763 samples
# remove NA
to_remove <- c(NA)
no_Frostdata.bac <- prune_samples(!(sample_data(no_Frostdata.bac)$name %in% to_remove), no_Frostdata.bac) # 5510 otus, 746 samples
# remove NA TMM
no_Frostdata.bac.TMM <- prune_samples(!(sample_data(no_Frostdata.bac.TMM)$name %in% to_remove), no_Frostdata.bac.TMM) # 8616 otus, 761 samples
# Subset by older birds minus Frosty 
olderdata.bac <- no_Frostdata.bac %>%
  subset_samples(daysfromhatch>670) %>%
  subset_samples(name=!NA)
olderdata.bac=prune_taxa(taxa_sums(olderdata.bac)>0, olderdata.bac) # 4139 otus, 367 samples
# Subset by older birds minus Frosty 
olderdata.bac.TMM <- no_Frostdata.bac.TMM %>%
  subset_samples(daysfromhatch>670) %>%
  subset_samples(name=!NA)
olderdata.bac.TMM=prune_taxa(taxa_sums(olderdata.bac.TMM)>0, olderdata.bac.TMM) # 8616 otus, 374 samples
# Subset no_Frostdata.bac with coccidiosis data
coccid <- as.vector(data.frame(sample_data(poodata.bac))$coccidiosis.status)
coccid <- coccid%in%c("negative","positive")
coccid[is.na(coccid)]=FALSE
poodata.bac.coccid <- prune_samples(coccid, poodata.bac) # 7547 otus, 325 samples (all life stages)
# Subset no_Frostdata.bac with coccidiosis data TMM
coccid <- as.vector(data.frame(sample_data(poodata.bac.TMM))$coccidiosis.status)
coccid <- coccid%in%c("negative","positive")
coccid[is.na(coccid)]=FALSE
poodata.bac.coccid.TMM <- prune_samples(coccid, poodata.bac.TMM) # 8616 otus, 336 samples (all life stages)
# Subset poodata.fun.coccid by life stage (only Runs) 
coccid <- as.vector(data.frame(sample_data(poodata.bac.coccid))$location)
coccid <- coccid%in%c("Runs")
coccid[is.na(coccid)]=FALSE
poodata.bac.coccid.runsonly <- prune_samples(coccid, poodata.bac.coccid) # 7547 otus, 62 samples
# Subset poodata.fun.coccid by life stage (only Runs) TMM
coccid <- as.vector(data.frame(sample_data(poodata.bac.coccid.TMM))$location)
coccid <- coccid%in%c("Runs")
coccid[is.na(coccid)]=FALSE
poodata.bac.coccid.runsonly.TMM <- prune_samples(coccid, poodata.bac.coccid.TMM) # 8616 otus, 65 samples
# Subset poodata.fun.coccid by life stage (only Runs) 
current <- as.vector(data.frame(sample_data(poodata.bac))$location)
current <- current%in%c("Runs")
current[is.na(current)]=FALSE
poodata.bac.runsonly.current <- prune_samples(current, poodata.bac) # 7547 otus, 488 samples
# Subset poodata.fun.coccid by life stage (only Runs) TMM
current <- as.vector(data.frame(sample_data(poodata.bac.TMM))$location)
current <- current%in%c("Runs")
current[is.na(current)]=FALSE
poodata.bac.runsonly.current.TMM <- prune_samples(current, poodata.bac.TMM) # 8616 otus, 495 samples
# Samples (post-decontam, post-clean, includes poo and soil)
Bac.data.prune_kiwi
Bac.data.prune_kiwi.TMM
Bac.data.prune_kiwi.TMM.keepwild
# Samples (subset of Bac.data.prune_kiwi, just poo)
poodata.bac
poodata.bac.TMM
# Samples (subset of Bac.data.prune_kiwi, just soil)
soildata.bac
soildata.bac.TMM
# Samples (subet of poodata.bac, just poo minus the resident kiwi)
no_Frostdata.bac
no_Frostdata.bac.TMM
# Samples (subet of poodata.bac, just younger cohort)
youngerdata.bac
youngerdata.bac.TMM
# Samples (subet of poodata.bac, just older cohort)
olderdata.bac
olderdata.bac.TMM
# collapse replicates that are TMM normalized
bac_poo_soil_collapse <- collapse_replicates(Bac.data.prune_kiwi.TMM,method = "sample", replicate_fields = c("name", "daysfromhatch")) # 661 samples
bac_poo_collapse <- collapse_replicates(no_Frostdata.bac.TMM,method = "sample", replicate_fields = c("name", "daysfromhatch")) # 624 samples
bac_young_collapse <- collapse_replicates(youngerdata.bac.TMM,method = "sample", replicate_fields = c("name", "daysfromhatch")) # 337 samples
bac_old_collapse <- collapse_replicates(olderdata.bac.TMM,method = "sample", replicate_fields = c("name", "daysfromhatch")) # 289 samples

# The most abundant OTUs within all samples - bacteria
common.taxa.bac <- sort(taxa_sums(no_Frostdata.bac),T)
common.taxa.bac <- tax_table(no_Frostdata.bac)[names(common.taxa.bac[1:50])]; common.taxa.bac

###########################################################################
# 5 Alpha diversity
# -------------------------------------------------------------------------

# what (i.e. sample type) -------------------------------------------------
alpha.bac.shan.what <- plot_richness(
  Bac.data.prune_kiwi, x="what.",measures=c("Shannon"), color="what.") + 
  geom_boxplot() + 
  scale_colour_manual(values = colors_scheme_what) + 
  geom_quasirandom(size = 3.0)
alpha.bac.shan.what$layers <- alpha.bac.shan.what$layers[-1]; alpha.bac.shan.what
# ANOVA and Tukey test | anova better for categorical so use lm
aov.what.bac <- aov(value~what., data=alpha.bac.shan.what$data)
summary(aov.what.bac)
HSD.what.bac <- HSD.test(aov.what.bac, "what.", group=T);HSD.what.bac

# soil location -----------------------------------------------------------
alpha.soil <- plot_richness(
  soildata.bac, x="location", measures=c("Shannon"), color="location") + 
  geom_boxplot() + 
  geom_quasirandom(size = 3.0) + 
  scale_colour_manual(values=colors_scheme_life_stage)
alpha.soil$layers <- alpha.soil$layers[-1];alpha.soil 
# ANOVA and Tukey test | anova better for categorical so use lm
aov.soil.bac <- aov(value~location, data=alpha.soil$data)
summary(aov.soil.bac)
HSD.soil.bac <- HSD.test(aov.soil.bac, "location", group=T);HSD.soil.bac


# cohort ------------------------------------------------------------------
alpha.bac.shan.cohort <- plot_richness(
  no_Frostdata.bac, x="cohort", measures=c("Shannon"), color="cohort") + 
  geom_boxplot() +
  scale_colour_manual(values = color_scheme_cohort) + 
  geom_quasirandom(size = 3.0)
alpha.bac.shan.cohort$data$cohort <- factor(alpha.bac.shan.cohort$data$cohort,levels=c("young", "old"))
alpha.bac.shan.cohort$layers <- alpha.bac.shan.cohort$layers[-1];alpha.bac.shan.cohort 
# ANOVA and Tukey test | anova better for categorical so use lm
aov.cohort.bac <- aov(value~cohort, data=alpha.bac.shan.cohort$data)
summary(aov.cohort.bac)
HSD.cohort.bac <- HSD.test(aov.cohort.bac, "cohort", group=T);HSD.cohort.bac

# kiwi location aka life stage --------------------------------------------
alpha.bac.shan.all <- plot_richness(
  no_Frostdata.bac, x="location", measures=c("Shannon"), color="location") + 
  geom_boxplot() + 
  scale_colour_manual(values = colors_scheme_life_stage) + 
  geom_quasirandom(size = 3.0)
alpha.bac.shan.all$data$location <- factor(alpha.bac.shan.all$data$location,levels=c("Hatch room", "Brooder room", "Runs"))
alpha.bac.shan.all$layers <- alpha.bac.shan.all$layers[-1];alpha.bac.shan.all
# ANOVA and Tukey test | anova better for categorical so use lm
aov.kiwi.bac.all <- aov(value~location, data=alpha.bac.shan.all$data)
summary(aov.kiwi.bac.all)
HSD.kiwi.bac.all <- HSD.test(aov.kiwi.bac.all, "location", group=T);HSD.kiwi.bac.all


# kiwi origin -------------------------------------------------------------
alpha.bac.shan.origin <- plot_richness(
  no_Frostdata.bac, x="origin", measures=c("Shannon"), color="origin") + 
  geom_boxplot() + 
  geom_quasirandom(size = 3.0)
alpha.bac.shan.origin$layers <- alpha.bac.shan.origin$layers[-1];alpha.bac.shan.origin
# ANOVA and Tukey test | anova better for categorical so use lm
aov.kiwi.bac.ori <- aov(value~origin, data=alpha.bac.shan.origin$data)
summary(aov.kiwi.bac.ori)
HSD.kiwi.bac.ori <- HSD.test(aov.kiwi.bac.ori, "origin", group=T);HSD.kiwi.bac.ori


# coccidiosis status ------------------------------------------------------
alpha.bac.shan.coccid <- plot_richness(
  poodata.bac.coccid.runsonly, x="coccidiosis.status", measures=c("Shannon"), color="coccidiosis.status") +
  geom_boxplot() + 
  geom_quasirandom(size = 3.0) 
alpha.bac.shan.coccid$layers <- alpha.bac.shan.coccid$layers[-1]; alpha.bac.shan.coccid
# ANOVA and Tukey test | anova better for categorical so use lm
aov.kiwi.bac.cocc <- aov(value~coccidiosis.status, data=alpha.bac.shan.coccid$data)
summary(aov.kiwi.bac.cocc)
HSD.kiwi.bac.cocc <- HSD.test(aov.kiwi.bac.cocc, "coccidiosis.status", group=T);HSD.kiwi.bac.cocc


# days from hatch younger cohort ------------------------------------------
alpha.bac.shan.daysfromhatch.younger <- plot_richness(
  youngerdata.bac, x="daysfromhatch", measures=c("Shannon"), color="daysfromhatch")
ggscatter(alpha.bac.shan.daysfromhatch.younger$data, x="daysfromhatch", y="value", 
          add = "reg.line",  # Add loess
          conf.int = TRUE,   # Add confidence interval
          cor.coef = T,      # Add correlation coefficient. see ?stat_cor
          size = 3,          # Size of dots
          color="#a09f9e",
          show.legend.text = NA) + 
  theme_light() + 
  theme(legend.position = "none")
ggscatter(alpha.bac.shan.daysfromhatch.younger$data, x="daysfromhatch", y="value",
          add = "reg.line",   # Add loess
          add.params = list(color = "name", fill = "name"), # Customize reg. line
          conf.int = TRUE,    # Add confidence interval
          cor.coef = T,       # Add correlation coefficient. see ?stat_cor
          size = 1,           # Size of dots
          color = "#a09f9e",
          show.legend.text = NA) + 
  facet_wrap(name~.) + 
  theme_light() +
  theme(legend.position = "none") +
  ylim(0,5)


# alpha diversity by days from hatch older --------------------------------
alpha.bac.shan.daysfromhatch.older <- plot_richness(
  olderdata.bac, x="daysfromhatch", measures=c("Shannon"), color="daysfromhatch") 
ggscatter(alpha.bac.shan.daysfromhatch.older$data, x="daysfromhatch", y="value", 
          add = "reg.line",  # Add loess
          conf.int = TRUE,   # Add confidence interval
          cor.coef = T,      # Add correlation coefficient. see ?stat_cor
          size = 3,          # Size of dots
          color = "#a09f9e",
          show.legend.text = NA) + 
  theme(legend.position = "none")
ggscatter(alpha.bac.shan.daysfromhatch.older$data, x="daysfromhatch", y="value", color = "cohort",
          add = "reg.line",  # Add loess
          palette = c("#9f5a47","#e29580"),
          add.params = list(color = "cohort", fill = "cohort"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          size = 1, # Size of dots
          show.legend.text = NA) + 
  facet_wrap(name~.) + 
  theme_light() +
  theme(legend.position = "none") +
  ylim(0,6)


# alpha diversity by collection date --------------------------------------
alpha.bac.shan.colldat <- plot_richness(
  poodata.bac, x="julian..collection.", measures=c("Shannon"), color="julian..collection.")
ggscatter(alpha.bac.shan.colldat$data, x="julian..collection.", y="value", 
          add = "reg.line",  # Add loess
          conf.int = TRUE,   # Add confidence interval
          cor.coef = T,      # Add correlation coefficient. see ?stat_cor
          size = 3,          # Size of dots
          color = "#a09f9e",
          show.legend.text = F) + 
  theme_light() +
  theme(legend.position = "none")
ggscatter(alpha.bac.shan.colldat$data, x="julian..collection.", y="value", color = "cohort",
          palette = c("#9f5a47","#e29580"),
          add = "reg.line",  # Add loess
          add.params = list(color = "cohort", fill = "cohort"), # Customize reg. line
          conf.int = TRUE,   # Add confidence interval
          cor.coef = T,      # Add correlation coefficient. see ?stat_cor
          size = 1,          # Size of dots
          show.legend.text = F) + 
  facet_wrap(name~.) + 
  ylim(0,6) + 
  theme_light() +
  theme(legend.position = "none") + 
  theme(strip.background =element_rect(fill="black")) +
  theme(strip.text = element_text(colour = 'white')) +
  xlab("Collection date") +
  ylab("Shannon diversity value")


# alpha diversity by age --------------------------------------------------
alpha.bac.shan.age <- plot_richness(
  no_Frostdata.bac, x="age..days.", measures=c("Shannon")) 
alpha.bac.shan.age.all <- ggscatter(alpha.bac.shan.age$data, x="age..days.", y="value", color="grey",
                                    add = "reg.line",  # Add loess
                                    conf.int = TRUE,   # Add confidence interval
                                    add.params = list(color = "black", fill = "black"), # Customize reg. line
                                    #cor.coef = T, # Add correlation coefficient. see ?stat_cor
                                    size = 3,          # Size of dots
                                    alpha=0.7,
                                    #color = "black",
                                    show.legend.text = NA) + 
  theme_light() +
  theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  xlab("Age (days)") + 
  ylab("Shannon diversity value")
ggscatter(alpha.bac.shan.age$data, x="age..days.", y="value", color = "#a09f9e",
          add = "reg.line",  # Add loess
          palette = c("#9f5a47","#e29580"),
          add.params = list(color = "cohort", fill = "cohort"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          size = 1, # Size of dots
          show.legend.text = F) + 
  facet_wrap(name~.) + 
  theme(legend.position = "none") + 
  ylim(0,6) + 
  theme_light() + 
  theme(strip.background =element_rect(fill="black")) +
  theme(strip.text = element_text(colour = 'white'))

# linear mod with age  ----------------------------------------------------
hist(alpha.bac.shan.age.all$data$value)
mod_age_bac <- lm(alpha.bac.shan.age.all$data$value~alpha.bac.shan.age.all$data$age..days.)
summary(mod_age_bac)
r2(mod_age_bac)
alpha.bac.shan.age.all
# glmm with age and name as random effect
library(lme4)
mod_age_ran_bac <- lmer(value ~ age..days. + (1|name), data=alpha.bac.shan.age.all$data)
fm2 <- lme(value ~ age..days., data=alpha.bac.shan.age.all$data, random = ~ 1|name)
summary(fm2)
library(nlme)
newdat <- expand.grid(name=unique(alpha.bac.shan.age.all$data$name),
                      age..days.=c(min(alpha.bac.shan.age.all$data$age..days.),
                                   max(alpha.bac.shan.age.all$data$age..days.)))
ggplot(alpha.bac.shan.age.all$data, aes(x=age..days.,y=value, colour=name)) +
  geom_point(size=3) +
  geom_line(aes(y=predict(fm2), group=name, size="name")) +
  geom_line(data=newdat, aes(y=predict(fm2, level=0, newdata=newdat), size="All")) +
  scale_size_manual(name="Predictions", values=c("name"=0.5, "All"=3)) +
  theme_bw(base_size=22) 
print(p)

summary(mod_age_ran_bac)
r2(mod_age_ran_bac)
# extract the estimates of the fixed effects
fixef(mod_age_ran_bac)
# extract the estimates of the random effects
ranef(mod_age_ran_bac)
r2(mod_age_ran_bac, tolerance=1e-10)
check_singularity(mod_age_ran_bac)
plot(mod_age_ran_bac)
qqnorm(resid(mod_age_ran_bac))

mod_age_bac_check <- lmer(value ~ 1 + (1|name), data=alpha.bac.shan.age$data)
summary(mod_age_bac_check)
r2(mod_age_bac_check)
anova(mod_age_bac_check, mod_age_ran_bac)
AIC(mod_age_bac_check, mod_age_ran_bac)


glm_mod_age_ran_bac <- glmr(value ~ age..days. + (1|name), data=alpha.bac.shan.age.all$data)
summary(glm_mod_age_ran_bac)
check_singularity(glm_mod_age_ran_bac)
r2(glm_mod_age_ran_bac)
plot(mod_age_ran_bac)
qqnorm(resid(mod_age_ran_bac))


###########################################################################
# 6 Bacterial Abundance 
# -------------------------------------------------------------------------

# Area plot
bac_poo_collapse_no_norm <- collapse_replicates(
  no_Frostdata.bac, method = "sample", replicate_fields = c("name", "daysfromhatch"))

# colors for taxa
bac_fam_col <- c(
  "Bradyrhizobiaceae"="#A6CEE3", "Comamonadaceae"="#1F78B4", "Enterobacteriaceae"="#B2DF8A", 
  "Lachnospiraceae"="#33A02C", "Lactobacillaceae"="#FB9A99", "Moraxellaceae"="#E31A1C",
  "Other"="#FDBF6F", "Pseudomonadaceae"="#FF7F00", "Ruminococcaceae"="#CAB2D6", 
  "Sphingomonadaceae"="#6A3D9A","Unknown"="#FFFF99","#B15928")
# chose specific participants to plot data
pts <- c(
  "Monty","Kaitiaki","Python","Sonic","Sojourn","Palindrome","Manawa","Pukukino","Ludo",
  "Palomita","Gizmo","Whakaroau", "Moses","KJ","Gummy","Kerrigan","Ata","Adieu", "Knuckles",
  "Loki","Frenchie","Sammy","Pippen","Kauri", "Helios","Etrick","Ki Tua","Pakiki","Lonestar",
  "Grawp", "Valentin","Leon","Maui","Pihara")
pts_young <- c(
  "Monty","Kaitiaki","Python","Sonic","Sojourn","Palindrome", "Manawa","Pukukino","Gizmo",
  "Whakaroau","Gummy","Ata","Knuckles","Loki","Frenchie","Sammy","Grawp", "Valentin","Leon",
  "Maui","Pihara")
pts_old <- c(
  "Ludo","Palomita","Moses","KJ","Kerrigan","Adieu", "Pippen","Kauri","Helios","Etrick",
  "Ki Tua","Pakiki","Lonestar")

bac_poo_collapse.namesub <- subset_samples(bac_poo_collapse, name %in% pts)
bac_poo_collapse.ready <- microbiome::transform(bac_poo_collapse.namesub, "compositional") 

bac_poo_collapse.namesub.no.norm <- subset_samples(bac_poo_collapse_no_norm, name %in% pts)
bac_poo_collapse.ready.no.norm <- microbiome::transform(bac_poo_collapse.namesub.no.norm, "compositional") 

# all cohorts by time (collection date, age, or days from hatch)
area_plot_bac_all_age <- plot_area(
  bac_poo_collapse.ready.no.norm, 
  xvar="age..days.",
  level = "phylum", facet.by = "name",
  abund.thres = 0.05, prev.thres=0.1, 
  fill.colors=brewer.pal(12,"Paired"), 
  ncol=4, nrow=9)
area_plot_bac_all_age + 
  xlab("Age (days)") +
  ylab("Relative Abundance of Bacterial Phyla") + 
  scale_x_continuous(limits = c(1,132), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult=c(0,0))) +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(strip.text.x = element_text(size = 12))

area_plot_bac_all_age_genus <- plot_area(ps.rel, xvar="daysfromhatch", 
                                         level = "genus",
                                         facet.by = "name",
                                         abund.thres = 0.1,
                                         prev.thres=0.1,
                                         fill.colors=brewer.pal(12,"Paired"),
                                         ncol=6,
                                         nrow=6)
area_plot_bac_all_age_genus + 
  xlab("Days from Hatch") +
  ylab("Relative Abundance of Fungal Family") + 
  scale_x_continuous(limits = c(677,764), expand = c(0, 0)) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(mult=c(0,0))) +
  theme(strip.background =element_rect(fill="#9f5a47"))+
  theme(strip.text = element_text(colour = 'white')) + 
  theme(strip.text.x = element_text(size = 12))

###########################################################################################
# 7 ORDINATION
# -------------------------------------------------------------------------

plotbeta(Bac.data.prune_kiwi.TMM.keepwild, group = "wild.cap", ellipse=T, method="PCoA")
set.seed(723)
# use plot_beta function from microbial package, only for categorical
# beta diversity, default is bray
# soil vs poo
bac.pcoa.what.TMM <- plotbeta(bac_poo_soil_collapse, group="what.", ellipse=T, method="PCoA")
bac.pcoa.what.TMM + scale_colour_manual(values = colors_scheme_what)
bac.pcoa.soil.TMM <- plotbeta(soildata.bac.TMM, group="location", ellipse=T, method="PCoA")
bac.pcoa.soil.TMM + scale_colour_manual(values = colors_scheme_life_stage)
# just poo, young vs old
bac.pcoa.cohort.TMM <- plotbeta(bac_poo_collapse, group="cohort", ellipse=T, method="PCoA")
bac.pcoa.cohort.TMM + scale_colour_manual(values = color_scheme_cohort)
plot_ly(x=bac.pcoa.cohort.TMM$data$Axis.1, y=bac.pcoa.cohort.TMM$data$Axis.2, z=bac.pcoa.cohort.TMM$data$Axis.3, 
        type="scatter3d", mode="markers", color=bac.pcoa.cohort.TMM$data$cohort)
# just poo, life stage
bac.pcoa.lifestage.TMM <- plotbeta(bac_poo_collapse, group="location", ellipse=T, method="PCoA")
bac.pcoa.lifestage.TMM + scale_colour_manual(values = colors_scheme_life_stage)
# ordination
older.ordi.bac <- ordinate(bac_old_collapse, method="PCoA", distance="bray", k=2); beep()
young.ordi.bac <- ordinate(bac_young_collapse, method="PCoA", distance="bray", k=2); beep()
all.poo.ordi.bac <- ordinate(bac_poo_collapse, method="PCoA", distance="bray", k=2); beep()
soil.ordi.bac <- ordinate(soildata.bac.TMM, method="PCoA", distance="bray", k=2); beep()
coccid.ordi.bac <- ordinate(poodata.bac.coccid.runsonly.TMM, method="PCoA", distance="bray", k=2); beep()

sample_data(bac_poo_collapse)$disease_history <- forcats::fct_collapse(
  sample_data(bac_poo_collapse)$history.of.positive.results, 
  History=c("diarrhea", "coccidiosis", "coccidiosis, worms","operation, coccidiosis, worms",  "dehydrated",                    
            "emphysema, coccidiosis", "coccidiosis, worms, eyes","coccidiosis, hernia",  
            "emphysema, coccidiosis, worms",  "dehydrated, failed weight gain", "coccidiosis, fungal"), 
  None=c(NA),group_other = T, stringsAsFactors = FALSE)
levels(sample_data(bac_poo_collapse)$disease_history) = c("History", "None")
sample_data(bac_poo_collapse) <- sample_data(bac_poo_collapse) %>% 
  replace_na(list(current='None'))

# NMDS all samples colored by INSERT VARIABLE  
PCoA.time.bac <- plot_ordination(
  (bac_poo_collapse),all.poo.ordi.bac,color="age..days.",shape="disease_history") +
  geom_point(size=4, alpha=0.7, na.rm = T) +
  scale_color_viridis(option = "D") +
  scale_shape_manual(values=c(15,19,18))+
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black")); PCoA.time.bac 
# NMDS all samples colored by INSERT VARIABLE  
PCoA.coccidiosis.bac <- plot_ordination(
  (poodata.bac.coccid.runsonly.TMM),coccid.ordi.bac,color="current.positive.results.",shape="coccidiosis.status") +
  geom_point(size=4, alpha=0.7, na.rm = T) +
  scale_shape_manual(values=c(15,19,18))+
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black")); PCoA.coccidiosis.bac 
# remove young birds

# NMDS all samples colored by INSERT VARIABLE  
PCoA.coccidiosis.bac <- plot_ordination((poodata.bac.coccid.runsonly.TMM),
                                        coccid.ordi.bac,color="coccidiosis.status",
                                        shape="coccidiosis.status") +
  geom_point(size=4, alpha=0.7, na.rm = T) +
  scale_shape_manual(values=c(15,19,18))+
  stat_ellipse(aes(col=coccidiosis.status), size = 1) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()); PCoA.coccidiosis.bac 
# remove young birds

# remove outliers
# AGO529.0773, AGO529.0742
poodata.bac.coccid.runsonly.TMM
outlier_samples <- c("AG0529.0773", "AG0529.0742")
poodata.bac.coccid.runsonly.TMM.remove.out <- 
  prune_samples(!(sample_data(poodata.bac.coccid.runsonly.TMM)$samplename %in% outlier_samples), 
                poodata.bac.coccid.runsonly.TMM)
coccid.ordi.bac.no.out <- ordinate(poodata.bac.coccid.runsonly.TMM.remove.out, method="PCoA", distance="bray", k=2); beep()

PCoA.coccidiosis.bac.no.out <- plot_ordination((poodata.bac.coccid.runsonly.TMM.remove.out),
                                               coccid.ordi.bac.no.out,color="coccidiosis.status",
                                               shape="coccidiosis.status") +
  geom_point(size=4, alpha=0.7, na.rm = T) +
  scale_shape_manual(values=c(15,19,18))+
  stat_ellipse(aes(col=coccidiosis.status), size = 1) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()); PCoA.coccidiosis.bac.no.out 

###########################################################################
# 8 PermANOVA
# -------------------------------------------------------------------------

# dealing with NAs
no_NAs <- as.vector(data.frame(sample_data(no_Frostdata.bac))$weight.of.food.consumed..grams.)
no_NAs <- (!(no_NAs%in%c(NA)))
no_NAs[is.na(no_NAs)] = FALSE
poodata.bac.food.consumed = prune_samples(no_NAs, no_Frostdata.bac)
# variables to look at sample type, location, cohort, coccidiosis, collection date, food weight
# Pick relative abundances (compositional) and sample metadata 
otu <- abundances(soildata.bac.TMM)
meta <- meta(soildata.bac.TMM)
permanova <- adonis(t(otu) ~ location, data = meta, permutations=999, method = "bray")
# P-value
print(as.data.frame(permanova$aov.tab)["what.", "Pr(>F)"])
temp_data_age=subset_samples(Bac.data.prune_kiwi.TMM, what.=="kiwi poo")
temp_data_age=subset_samples(temp_data_age,!is.na(cohort))
otu_age <- abundances(temp_data_age)
meta_age <- meta(temp_data_age)
permanova_age <- adonis(t(otu_age) ~ age..days., 
                        data = meta_age, 
                        permutations=999, method = "bray")
permanova_age
try_mod <- lm(t(otu_age) ~ age..days., 
              data = meta_age)
summary(try_mod)
plot(try_mod)

###########################################################################################
# 9 BETA DIVERSITY  
# -------------------------------------------------------------------------

# remove NA
to_remove <- c(NA, "Tiggy Winkle")
z <- prune_samples(!(sample_data(bac_poo_collapse)$name %in% to_remove), bac_poo_collapse)
y <- subset_samples(Bac.data.prune_kiwi.TMM, !is.na(cohort))
# OTU table after normalization
otu.beta <- otu_table(y)

# Convert to data frame and transpose and calculate using betadiver 
# betadiver is for presence/absence may need to try vegdist or dist 
data.beta.pa <- betadiver(as.data.frame(t(otu.beta)), "hk")         # presence/absence
data.beta <- vegdist(as.data.frame(t(otu.beta)), method = "bray")   # abundance

# Calculate beta diversity using distances from above and matching to sample data
data.dispersion <- with(sample_data(y), betadisper(data.beta, age..days., type="centroid"))
data.dispersion.cohort <- with(sample_data(y), betadisper(data.beta, cohort, type="centroid"))
boxplot(data.dispersion.cohort)

# Permutation test
perm <- permutest(data.dispersion.cohort, pairwise=T)

# max age of kiwi
max_age_of_individuals<-read_csv("data/max_age_of_individuals.csv")
max_age_of_individuals<- max_age_of_individuals %>%
  arrange(max_age)
max_age_of_individuals$name
# Putting data in long form so ggplot can read
distance.name <- data.frame(distance=data.dispersion.cohort$distances, name=sample_data(y)$name, cohort=sample_data(y)$cohort, age=sample_data(y)$age..days.)
data.dispersion.cohort
# Plotting in ggplot
kiwi.beta <- ggplot(distance.name, aes(x=cohort, y=distance, col=cohort)) + 
  scale_color_manual(values = color_scheme_cohort) + 
  xlab("") + ylab("Distance to Centroid") + 
  geom_boxplot() +
  # geom_smooth(method="lm") + 
  geom_quasirandom(size=3) + 
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none");kiwi.beta 
kiwi.beta$data$name <- factor(kiwi.beta$data$name,levels=as.vector(max_age_of_individuals$name))
kiwi.beta$data$cohort <- factor(kiwi.beta$data$cohort,levels=as.vector(c("young","old")))
# anova to see if statistically significantly different 
temp.aov <- aov(distance~name, data=distance.name)
summary(temp.aov)
temp.HSD <- HSD.test(temp.aov, "cohort", group = T)

###########################################################################
# 10 clamtest
# -------------------------------------------------------------------------

# Bacteria
# Standard simper code
community.kiwi.bac.t2 <- as.matrix(t(otu_table(Bac.data.prune_kiwi)))
environment_data.kiwi.bac2 <- as.matrix(sample_data(Bac.data.prune_kiwi))
write.csv(environment_data.kiwi.bac2, "environment_data_2020_2.csv")
environment_data.kiwi.bac2 <- data.frame(read_csv("environment_data_2020_2.csv"))
row.names(environment_data.kiwi.bac2) <- environment_data.kiwi.bac2$submit_form_code
write.csv(community.kiwi.bac.t2, "community.kiwi.bac2.csv")
community.kiwi.bac2 <- as.matrix(read_csv("community.kiwi.bac2.csv", col_names=T)[,-1])
row.names(community.kiwi.bac2) <-  environment_data.kiwi.bac2$submit_form_code
community.kiwi.bac2 <- community.kiwi.bac2[,-1]
# Remove samples with no sequences in it
community.kiwi.bac.no.zero.OTU2 <- community.kiwi.bac2@.Data[, colSums(community.kiwi.bac2@.Data != 0) > 0]
# clamtest
clam_kiwi2 <- clamtest(community.kiwi.bac.no.zero.OTU2, environment_data.kiwi.bac2$what.)
summary(clam_kiwi2)
saveRDS(clam_kiwi2, "clam_analysis_kiwi_2020_2.RDS")
readRDS("clam_analysis_kiwi_2020_2.RDS")
#c("soil"="#fdb462","kiwi poo"="#b3de69")
base_plot_kiwi <- plot(clam_kiwi2, "Soil OTU Abundances", "Kiwi OTU Abundances", "Bacterial Species Classification", 
                       pch = c(19,19,19,19), 
                       col.points = c("tan3","#b3de69","#fdb462","grey"),cex = 1, las = 1,
                       col.lines = c("#bf812d","#b3de69","grey"), lty = c(2,2,2), lwd = 3, position = NULL)
legend("bottomright", inset=c(-0.135,0.00),legend = c("Rare", "Kiwi specialist", "Soil specialist", "Generalist"),
       col = c("grey","#b3de69","#bf812d","tan3"), pch = c(19,19,19,19), cex = 0.8)

clam_kiwi2_result_df <- as.data.frame(clam_kiwi2)
colnames(clam_kiwi2_result_df) <- c("OTUs", "Kiwi", "Soil", "Classes")
gg_plot_kiwi <- ggplot(
  clam_kiwi2_result_df, aes(Soil, Kiwi, col=Classes)) + 
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=c("Specialist_soil"="#fdb462",
                              "Specialist_kiwi poo"="#b3de69",
                              "Generalist"="tan3",
                              "Too_rare"="grey")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("Soil OTU Abundances") + 
  ylab("Kiwi OTU Abundances") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

###########################################################################
# Calculate Bray-Curtis ---------------------------------------------------

dist.matrix       <- phyloseq::distance(no_Frostdata.bac,"bray")
bray_curtis_vegan <- vegdist(community.kiwi.bac, method="bray", binary=F, diag=F, upper=F, na.rm=F)

# Do I need to calculate bray-curtis values per individual?
# Not sure how this is being calculated
bray_values <- data.frame(bray=with(sample_data(youngerdata.bac), as.numeric(as.matrix(as.dist(vegdist(t(otu_table(youngerdata.bac)), method="bray"), upper = TRUE)))), 
                          name=sample_data(youngerdata.bac)$name,
                          daysfromhatch=sample_data(youngerdata.bac)$daysfromhatch)

bray_values <- data.frame(name=meta(youngerdata.bac)$name, daysfromhatch=meta(youngerdata.bac)$daysfromhatch, 
                          bray=as.numeric(phyloseq::distance(youngerdata.bac, method="bray")))
ggplot(bray_values,aes(daysfromhatch, bray, fill=name, color=name))+geom_point()+facet_wrap(name~.)


###########################################################################################
# 11 TEMPORAL BETA DIVERSITY
# -------------------------------------------------------------------------

library(MicrobeDS)
library(microbiome)
library(dplyr)
library(vegan)

# Pick the metadata for this subject and sort the
# samples by time
# Pick the data and modify variable names
pseq <- no_Frostdata.bac
# change kiwi name
s <- "Kauri" # Selected subject
# Let us pick a subset, change kiwi name
pseq <- subset_samples(no_Frostdata.bac, name == "Kauri") 
# Rename variables
sample_data(pseq)$subject <- sample_data(pseq)$name
sample_data(pseq)$sample <- sample_data(pseq)$samplename
# Order the entries by time
df <- meta(pseq) %>% arrange(daysfromhatch)

# Calculate the beta diversity between each time point and
# the baseline (first) time point
beta <- c() # Baseline similarity
s0 <- subset(df, daysfromhatch == 716)$sample # change the daysfromhatch to the min value of this specific bird
# Let us transform to relative abundance for Bray-Curtis calculations
a <-microbiome::abundances(microbiome::transform(pseq, "compositional")) 
for (tp in df$daysfromhatch[-1]) {
  # Pick the samples for this subject
  # If the same time point has more than one sample,
  # pick one at random
  st <- sample(subset(df, daysfromhatch == tp)$sample, 1)
  # Beta diversity between the current time point and baseline
  b <- vegdist(rbind(a[, s0], a[, st]), method = "bray")
  # Add to the list
  beta <- rbind(beta, c(tp, b))
}
colnames(beta) <- c("daysfromhatch", "beta")
beta <- as.data.frame(beta)

# kiwi individual beta diversity
p <- ggplot(beta, aes(x = daysfromhatch, y = beta)) +
  geom_point() +
  geom_line() +
  geom_smooth() +
  labs(x = "Time (Days)", y = "Beta diversity (Bray-Curtis)")
print(p)

