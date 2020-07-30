#======================================================================================================================
# Title		: Systematic literature review of Mortality in VL clinical trials
# Data version	: 24-Oct-2019
# Script Date	: 15-Jan-2020
# Data source	: Supplemental file 2 
#======================================================================================================================
#rm(list=ls())
#options(digits=10)

library(tidyverse)
library(readxl)
library(tableone)
library(meta)
library(binom)

setwd("C:/Users/pdahl/Dropbox/_VL AE Working Folder/AE MS/Final AE MS and codes/Data")

####################################################################################
# Section I: Basic description of the studies used included in the meta-analysis
####################################################################################
#--------------------------------------------
# Read data with details per treatment arm
#--------------------------------------------
dat <- read_excel("Supplemental file 3_analyses data.xlsx", 
		sheet = "meta_data")
	
#---------------------------------
# Number of articles by region
#---------------------------------
CreateTableOne(
		vars 		= c("Study.Region"), 
		factorVars 	= c("Study.Region"), 
		#strata 	= "Region", 
		data 		= dat[!duplicated(dat$Tag),], 
		test		= FALSE
		)
#-----------------------------
# Breakdown by Randomisation
#-----------------------------
CreateTableOne(
		vars 		= c("Patient.Allocation"), 
		factorVars 	= c("Patient.Allocation"), 
		data 		= dat[!duplicated(dat$Tag),], 
		test		= FALSE
		)
#---------------------------------------------
# Usage of HIV Status as exclusion criteria 
#---------------------------------------------
CreateTableOne(
		vars 		= c("HIV"), 
		factorVars 	= c("HIV"), 
		data 		= dat[!duplicated(dat$Tag),], 
		test		= FALSE
		)
dat %>% 
	dplyr::group_by(HIV) %>%
		dplyr::summarise(
		n_studies = length(unique(Tag))
	)

#------------------------------------------------------------
# Table 1: Description of drug arms and patients treated
#------------------------------------------------------------

# read the excel tab with relevent information 
arms <- read_excel("Supplemental file 3_analyses data.xlsx", 
		sheet = "VL_SAEs_per_arm_04_12_2019")

# Merge the arms details with study meta-data
arms <- dplyr::inner_join(
				arms, 
				subset(dat,select=c("Tag","Follow.up.duration","age_range","HIV","Allegery_history","Pregnant","Blinding")),
		by="Tag"
	)
#----------------------
# Recode the drug names
#----------------------
arms$drug_group <- ifelse(arms$drug_group=="AMBd","AMBd",
				ifelse(arms$drug_group=="Amphotericin b fat/lipid/colloid/cholesterol","AmBb-lipid",
				ifelse(arms$drug_group=="L-AmB","L-AmB (multiple)",
				ifelse(arms$drug_group=="L-AmB (Single dose)","L-AmB (Single dose)",
				ifelse(arms$drug_group=="L-AmB (Single) + Miltefosine","L-AmB (Single)combination regimen",
				ifelse(arms$drug_group=="L-AmB (Single) + PA","L-AmB (Single)combination regimen",
				ifelse(arms$drug_group=="L-AmB (Single) + Paromomycin","L-AmB (Single)combination regimen",
				ifelse(arms$drug_group=="L-AmB + Miltefosine","L-AmB (multiple) combination regimen",
				ifelse(arms$drug_group=="L-AmB + PA","L-AmB (multiple) combination regimen",
				ifelse(arms$drug_group=="L-AmB + Paromomycin","L-AmB (multiple) combination regimen",
				ifelse(arms$drug_group=="Miltefosine","Miltefosine",
				ifelse(arms$drug_group=="Paromomycin","Paromomycin",
				ifelse(arms$drug_group=="Pentamidine","Pentamidine",
				ifelse(arms$drug_group=="Miltefosine + Paromomycin","Miltefosine + Paromomycin",
				ifelse(arms$drug_group=="Paramomycin + Miltefosine","Miltefosine + Paromomycin",
				ifelse(arms$drug_group=="PA","PA",
				ifelse(arms$drug_group=="PA + Allopurinol","PA combination regimen",
				ifelse(arms$drug_group=="Aminosidine/paromomycin + SSG","PA combination regimen",
				ifelse(arms$drug_group=="Aminosidine sulphate or Paromomycin + SSG","PA combination regimen",
				ifelse(arms$drug_group=="PA + Ketoconazole","PA combination regimen",
				ifelse(arms$drug_group=="PA + Levamisole","PA combination regimen",
				ifelse(arms$drug_group=="PA + Paromomycin","PA combination regimen",
				ifelse(arms$drug_group=="PA + Pentamidine","PA combination regimen",
				ifelse(arms$drug_group=="PA + Verapamil","PA combination regimen",
				ifelse(arms$drug_group=="Paromomycin + PA","PA combination regimen", 
				ifelse(arms$drug_group=="Sitamaquine","Sitamaquine" ,
				"Other"
				))))))))))))))))))))))))))

#----------------------------------------------
# number of study arms a drug was tested in
#----------------------------------------------
CreateTableOne(
			vars 		= c("drug_group"), 
			factorVars 	= c("drug_group"), 
			#strata	="Study.Region",
			data 		= arms , 
			test		= FALSE
			)
--------------------------------------------------
# Number of patients treated per arm by drug group
#--------------------------------------------------
arms %>% 
		dplyr::group_by(drug_group) %>%
		dplyr::summarise(
			n_arms = length(Tag),
			n_patients = sum(n_treated, na.rm=TRUE)
		)
#------------------------------------------------
# Number of patients treated per arm by region
#------------------------------------------------
arms %>% 
	dplyr::group_by(drug_group,Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE)
	)

# Number of patients treated per region
arms %>% 
			dplyr::group_by(Study.Region) %>%
			dplyr::summarise(
				n_patients = sum(n_treated, na.rm=TRUE),
				n_prop = (sum(n_treated, na.rm=TRUE)/sum(arms$n_treated) )*100
			)
#----------------------------------
# Adverse Event Monitoring System
#---------------------------------
CreateTableOne(
		vars 		= c("Was.there.an.Adverse.Event.Monitoring.System."), 
		factorVars 	= c("Was.there.an.Adverse.Event.Monitoring.System."), 
		data 		= dat[!duplicated(dat$Tag),], 
		test		= FALSE
		)
####################################################################################
# Section II: Hypersensitivity to study drugs and discontinuation from trial
####################################################################################

#--------------------------------------------
# Read data on dose testing information
#--------------------------------------------
dose_testing_information <- read_excel("Supplemental file 3_analyses data.xlsx", 
						sheet = "dose_testing"
						)
#----------------------
# Recode the drug names
#----------------------
dose_testing_information$drug_group <- ifelse(dose_testing_information$drug_group=="AMBd","AMBd",
				ifelse(dose_testing_information$drug_group=="Amphotericin b fat/lipid/colloid/cholesterol","AmBb-lipid",
				ifelse(dose_testing_information$drug_group=="L-AmB","L-AmB (multiple)",
				ifelse(dose_testing_information$drug_group=="L-AmB (Single dose)","L-AmB (Single dose)",
				ifelse(dose_testing_information$drug_group=="L-AmB (Single) + Miltefosine","L-AmB (Single)combination regimen",
				ifelse(dose_testing_information$drug_group=="L-AmB (Single) + PA","L-AmB (Single)combination regimen",
				ifelse(dose_testing_information$drug_group=="L-AmB (Single) + Paromomycin","L-AmB (Single)combination regimen",
				ifelse(dose_testing_information$drug_group=="L-AmB + Miltefosine","L-AmB (multiple) combination regimen",
				ifelse(dose_testing_information$drug_group=="L-AmB + PA","L-AmB (multiple) combination regimen",
				ifelse(dose_testing_information$drug_group=="L-AmB + Paromomycin","L-AmB (multiple) combination regimen",
				ifelse(dose_testing_information$drug_group=="Miltefosine","Miltefosine",
				ifelse(dose_testing_information$drug_group=="Paromomycin","Paromomycin",
				ifelse(dose_testing_information$drug_group=="Pentamidine","Pentamidine",
				ifelse(dose_testing_information$drug_group=="Miltefosine + Paromomycin","Miltefosine + Paromomycin",
				ifelse(dose_testing_information$drug_group=="Paramomycin + Miltefosine","Miltefosine + Paromomycin",
				ifelse(dose_testing_information$drug_group=="PA","PA",
				ifelse(dose_testing_information$drug_group=="PA + Allopurinol","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="Aminosidine/paromomycin + SSG","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="Aminosidine sulphate or Paromomycin + SSG","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="PA + Ketoconazole","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="PA + Levamisole","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="PA + Paromomycin","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="PA + Pentamidine","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="PA + Verapamil","PA combination regimen",
				ifelse(dose_testing_information$drug_group=="Paromomycin + PA","PA combination regimen", 
				ifelse(dose_testing_information$drug_group=="Sitamaquine","Sitamaquine" ,
				"Other"
				))))))))))))))))))))))))))

hypersensitivity <- subset(dat, select=c("Tag","Allegery_history"))

# Merge hypersensitivity as exclusion criteria with dose testing data
dose_testing_information <- dplyr::inner_join(
						dose_testing_information,
						hypersensitivity, 
					by="Tag"
					)
#-------------------------------------------------------------
# Hypersenstiivity to study drug as exclusion criteria
#-------------------------------------------------------------
table(dose_testing_information$dose_testing)
table(dose_testing_information$Allegery_history)
table(dose_testing_information$drug_group,dose_testing_information$Allegery_history)

#-------------------------------------------------------------
# For dose testing; keep only the Amphotericin B formulations
#-------------------------------------------------------------
amph_b_dt <- dose_testing_information[which(dose_testing_information$drug_group %in% 	
					c("AmBb-lipid","AMBd","L-AmB (multiple)","L-AmB (multiple) combination regimen","L-AmB (Single dose)","L-AmB (Single)combination regimen")),]

# Keep the arms where dose testing was conducted
amph_b_dt  <- amph_b_dt   %>% 
		filter(dose_testing=="Yes")

amph_b_dt <- droplevels(amph_b_dt )

# summarise analysis dataset
amph_b_dt  %>% 
		dplyr::summarise(
			n_patients 	= sum(n_treated, na.rm=TRUE),
			n_arms 	= length(n_treated)
		)

#-----------------------------------------------------------------------------------------------------
# Estimation of rate of discontinuation of treatment after dose testing
#-----------------------------------------------------------------------------------------------------
#  Assumption made:
# 	When dose testing was carried out, but no information was reported regarding discontinuation
#	then, it is asssumed no discountination acutally occured
#-----------------------------------------------------------------------------------------------------
amph_b_dt$n_discontinuation[is.na(amph_b_dt$n_discontinuation)]<- 0
sum(amph_b_dt$n_discontinuation)
sum(amph_b_dt$n_treated)

# Estimate of treatment discontinuation from meta-analysis
metaprop( 
	n_discontinuation, 
	n_treated, 
	Study,
	sm="PLOGIT",
	data=amph_b_dt
	)

#####################################################################
# Section III: Description of the SAEs reported at study arm levels
#####################################################################

saes <- read_excel("Supplemental file 3_analyses data.xlsx", 
		sheet = "VL_SAEs_list_04_12_2019")

# Number of SAEs by time
CreateTableOne(
		vars 		= c("Time_group"), 
		factorVars 	= c("Time_group"), 
		data 		= saes, 
		test		= FALSE
		)

# Number of SAEs by grouping of disorders
CreateTableOne(
		vars 		= c("SAEs_group"), 
		factorVars 	= c("SAEs_group"), 
		data 		= saes, 
		test		= FALSE
		)
	
# Adjudicated relationship of the SAEs to the study drug
CreateTableOne(
		vars 		= c("Relationship"), 
		factorVars 	= c("Relationship"), 
		data 		= saes, 
		test		= FALSE
		)

# Death
CreateTableOne(
		vars 		= c("Outcome"), 
		factorVars 	= c("Outcome"), 
		data 		= saes, 
		#strata	="Time_group",
		test		= FALSE
		)

####################################################################
# Section IV: Description of deaths reported at study arm levels
####################################################################
death <- read_excel("Supplemental file 3_analyses data.xlsx", 
		sheet = "VL_SAEs_per_arm_04_12_2019")

# Time of death
death %>% 
			dplyr::summarise(
				total_deaths 		= sum(number_of_deaths, na.rm=TRUE),
				total_deaths_0_30 	= sum(n_deaths_month, na.rm=TRUE),
				total_deaths_30_plus	= sum(after_30_days, na.rm=TRUE),
				total_deaths_unclear 	= sum(unclear_time, na.rm=TRUE)
			)

# Proportion of deaths for each time period
death %>% 
			dplyr::summarise(
				total_deaths 		= sum(number_of_deaths, na.rm=TRUE)/sum(death$number_of_deaths, na.rm=TRUE),
				total_deaths_0_30 	= sum(n_deaths_month, na.rm=TRUE)/sum(death$number_of_deaths, na.rm=TRUE),
				total_deaths_30_plus	= sum(after_30_days, na.rm=TRUE)/sum(death$number_of_deaths, na.rm=TRUE),
				total_deaths_unclear 	= sum(unclear_time, na.rm=TRUE)/sum(death$number_of_deaths, na.rm=TRUE)
			)
#--------------------------------------------------------
# Adjudicated relationship of deaths to the study drug
#--------------------------------------------------------
CreateTableOne(
		vars 		= c("Relationship"), 
		factorVars 	= c("Relationship"), 
		data 		= saes[which(saes$Outcome=="Death"),], 
		test		= FALSE
		)

# End script (Not run)
