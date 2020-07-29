############################################################################################################
# Title		:	Systematic literature review of Mortality in VL clinical trials
# Data version	:	24-Oct-2019
# Script Date	: 	15-Jan-2020
# Data source	: 	Supplemental file 2 
# Task		: 	Estimate of rates within 30 days of treatment initation stratfied by drug and reigon
############################################################################################################
#rm(list=ls())
library(tidyverse)
library(readxl)

# Meta analysis packages
library(meta)
library(metafor )
library(metasens) # Copas selection model for adjusting for small study effects

setwd("C:/Users/pdahl/Dropbox/_VL AE Working Folder/AE MS/Final AE MS and codes/Data")

#--------------------------------------------
# Read data with details per treatment arm
#--------------------------------------------
dat <- read_excel("Supplemental file 3_analyses data.xlsx", 
		sheet = "VL_SAEs_per_arm_04_12_2019")
	
dat$drug_group_final <- ifelse(dat$drug_group=="AMBd","AMBd",
				ifelse(dat$drug_group=="Amphotericin b fat/lipid/colloid/cholesterol","AmBb-lipid",
				ifelse(dat$drug_group=="L-AmB","L-AmB (multiple)",
				ifelse(dat$drug_group=="L-AmB (Single dose)","L-AmB (Single dose)",
				ifelse(dat$drug_group=="L-AmB (Single) + Miltefosine","L-AmB (Single)combination regimen",
				ifelse(dat$drug_group=="L-AmB (Single) + PA","L-AmB (Single)combination regimen",
				ifelse(dat$drug_group=="L-AmB (Single) + Paromomycin","L-AmB (Single)combination regimen",
				ifelse(dat$drug_group=="L-AmB + Miltefosine","L-AmB (multiple) combination regimen",
				ifelse(dat$drug_group=="L-AmB + PA","L-AmB (multiple) combination regimen",
				ifelse(dat$drug_group=="L-AmB + Paromomycin","L-AmB (multiple) combination regimen",
				ifelse(dat$drug_group=="Miltefosine","Miltefosine",
				ifelse(dat$drug_group=="Paromomycin","Paromomycin",
				ifelse(dat$drug_group=="Pentamidine","Pentamidine",
				ifelse(dat$drug_group=="Miltefosine + Paromomycin","Miltefosine + Paromomycin",
				ifelse(dat$drug_group=="Paramomycin + Miltefosine","Miltefosine + Paromomycin",
				ifelse(dat$drug_group=="PA","PA",
				ifelse(dat$drug_group=="PA + Allopurinol","PA combination regimen",
				ifelse(dat$drug_group=="Aminosidine/paromomycin + SSG","PA combination regimen",
				ifelse(dat$drug_group=="Aminosidine sulphate or Paromomycin + SSG","PA combination regimen",
				ifelse(dat$drug_group=="PA + Ketoconazole","PA combination regimen",
				ifelse(dat$drug_group=="PA + Levamisole","PA combination regimen",
				ifelse(dat$drug_group=="PA + Paromomycin","PA combination regimen",
				ifelse(dat$drug_group=="PA + Pentamidine","PA combination regimen",
				ifelse(dat$drug_group=="PA + Verapamil","PA combination regimen",
				ifelse(dat$drug_group=="Paromomycin + PA","PA combination regimen", 
				ifelse(dat$drug_group=="Sitamaquine","Sitamaquine" ,
				"Other"
				))))))))))))))))))))))))))

#--------------------------------------------
# Read meta-data with details per treatment arm
#--------------------------------------------
meta_data <- read_excel("Supplemental file 3_analyses data.xlsx", 
		sheet = "meta_data")

#------------------------------------------------
# Merge the arms details with study meta-data
#------------------------------------------------
dat <- dplyr::inner_join(
			dat, 
			subset(meta_data,select=c("Tag","Follow.up.duration","age_range","HIV","Allegery_history","Pregnant","Blinding")),
			by="Tag"
		)

dat$Follow.up.duration <- as.numeric(as.character(dat$Follow.up.duration))

# Study had a max follow-up of 14 months on one of the patients
dat$Follow.up.duration[dat$Tag==114]<- 420 
dat$Follow.up.duration[dat$Follow.up.duration <0]<- NA

#==================================================================================
# Estimating incidence rate of mortality for overall dataset within first month
#==================================================================================
dat1 <- dat[which(!is.na(dat$n_deaths_month)),]

# For studies with missing follow-up, it was assumed that 30 days is the min follow-up
# For studies with <30 days, total person-time calculated using true follow-up time
# For studies with >=30 days follow-up, total person-time cut at 30 days

dat1$follow_up <- ifelse(is.na(dat1$Follow.up.duration), 30, dat1$Follow.up.duration )

dat1$person_time <- ifelse(
			dat1$follow_up <30,dat1$n_treated*dat1$follow_up,
			dat1$n_treated*30
		)

# Using suggestion from the metafor vignette to enchance convergence
# http://127.0.0.1:30486/library/metafor/html/rma.uni.html

#====================================================================================
# ESTIMATION OF OVERALL INCIDENCE RATE OF DEATH FOR EACH DRUG REGIMEN & BY REGION 
#====================================================================================

#========================================
# AMB Deoxycholate : Regional breakdown
#========================================
AMBd <- dat1[which(dat1$drug_group_final=="AMBd"),]

AMBd %>% 
	dplyr::group_by(drug_group,Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

# Calculation of rates for each of the geographical locations seaprately 
region  <- as.data.frame(table(AMBd$Study.Region))
region  <- region[which(region$Var1!="South America"),]		# no events
region  <- region[which(region$Var1!="Eastern Africa"),]		# only 1 study arm
region  <- droplevels(region )

table(AMBd$Study.Region, AMBd$n_deaths_month)

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=AMBd[which(AMBd$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# Egger's test
		e_t <- metabias(results, method.bias="linreg", k.min=2) 

		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	

region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
					region_output[,2:ncol(region_output)], 
					function(x) as.numeric(as.character(x))
					)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#==========================================================
# AMB fat/cholestoral/lipid/emulsion : Regional breakdown
#==========================================================
AMB_fat <- dat1[which(dat1$drug_group_final=="AmBb-lipid"),]

AMB_fat %>% 
	dplyr::group_by(drug_group,Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

# Calculation of rates for each of the geographical locations seaprately 
region  <- as.data.frame(table(AMB_fat$Study.Region))
region  <- region[which(region$Var1!="Mediterranean"),]	# no events
region  <- region[which(region$Var1!="South America"),]	# no events
region  <- droplevels(region )

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=AMB_fat[which(AMB_fat$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# Egger's test
		e_t <- metabias(results, method.bias="linreg", k.min=2) 

		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
						region_output[,2:ncol(region_output)], 
						function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#======================================
# LAMB single dose: Regional breakdown
#======================================
LAmB_single <- dat1[which(dat1$drug_group_final=="L-AmB (Single dose)"),]

LAmB_single %>% 
	dplyr::group_by(drug_group,Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

# Calculation of rates for each of the geographical locations seaprately 
region  <- as.data.frame(table(LAmB_single$Study.Region))
region  <- region[which(region$Var1!="Eastern Africa"),]	# only 1 arm
region  <- droplevels(region )

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=LAmB_single[which(LAmB_single$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
					region_output[,2:ncol(region_output)], 
					function(x) as.numeric(as.character(x))
					)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
					select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
				)
print(region_output)

#==============================================================
# LAMB single dose in a combination therapy: Regional breakdown
#==============================================================
LAmB_single_comb <- dat1[which(dat1$drug_group_final=="L-AmB (Single)combination regimen"),]

LAmB_single_comb %>% 
	dplyr::group_by(drug_group,Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

# Calculation of rates for each of the geographical locations seaprately 
region  <- as.data.frame(table(LAmB_single_comb$Study.Region))
region  <- region[which(region$Var1!="South America"),] 		# only 1 arm
region  <- droplevels(region )

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=LAmB_single_comb[which(LAmB_single_comb$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms	=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
					region_output[,2:ncol(region_output)], 
					function(x) as.numeric(as.character(x))
					)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#==============================================================
# LAMB multiple dose : Regional breakdown
#==============================================================
LAmB_multiple <- dat1[which(dat1$drug_group_final=="L-AmB (multiple)"),]

LAmB_multiple %>% 
	dplyr::group_by(drug_group,Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

# Calculation of rates for each of the geographical locations seaprately 
region  <- as.data.frame(table(LAmB_multiple$Study.Region))
region  <- region[which(region$Var1!="Not stated"),]		# 1 arm
region  <- droplevels(region )
region  

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=LAmB_multiple[which(LAmB_multiple$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
						region_output[,2:ncol(region_output)], 
						function(x) as.numeric(as.character(x))
					)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#==================================================================
# LAMB multiple dose in a combination therapy : Regional breakdown
#==================================================================
LAmB_multiple_comb <- dat1[which(dat1$drug_group_final=="L-AmB (multiple) combination regimen"),]

LAmB_multiple_comb %>% 
	dplyr::group_by(Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

metarate(n_deaths_month, person_time, Study, 
				data=LAmB_multiple_comb[which(LAmB_multiple_comb$Study.Region=="Eastern Africa"),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
# Cannot fit

#==============================================================
# Miltefosine  : Regional breakdown
#==============================================================
Miltefosine <- dat1[which(dat1$drug_group_final=="Miltefosine"),]

Miltefosine %>% 
	dplyr::group_by(Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

region  <- as.data.frame(table(Miltefosine$Study.Region))

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=Miltefosine [which(Miltefosine$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test

		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
					region_output[,2:ncol(region_output)], 
					function(x) as.numeric(as.character(x))
					)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#==============================================================
# Pentavalent Antimony  : Regional breakdown
#==============================================================
PA <- dat1[which(dat1$drug_group_final=="PA"),]

PA %>% 
	dplyr::group_by(Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

region  <- as.data.frame(table(PA$Study.Region))
region  <- region[which(region$Var1!="Not stated"),]		# 1 arm
region  
table(PA$Study.Region, PA$n_deaths_month)

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=PA[which(PA$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
					region_output[,2:ncol(region_output)], 
					function(x) as.numeric(as.character(x))
					)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#=============================================================================
# Pentavalent Antimony in a non-LAmB combination therapy : Regional breakdown
#=============================================================================
PA_combination <- dat1[which(dat1$drug_group_final=="PA combination regimen"),]

PA_combination %>% 
	dplyr::group_by(Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

region  <- as.data.frame(table(PA_combination$Study.Region))

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=PA_combination[which(PA_combination$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
					region_output[,2:ncol(region_output)], 
					function(x) as.numeric(as.character(x))
					)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#=============================================================================
# Paromomycin : Regional breakdown
#=============================================================================
Paromomycin <- dat1[which(dat1$drug_group_final=="Paromomycin"),]
region  <- as.data.frame(table(Paromomycin$Study.Region))

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=Paromomycin[which(Paromomycin$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
						region_output[,2:ncol(region_output)], 
						function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#=============================================================================
# Miltefosine + Paromomycin: Regional breakdown
#=============================================================================
Miltefosine_Paromomycin	<- dat1[which(dat1$drug_group_final=="Miltefosine + Paromomycin"),]

Miltefosine_Paromomycin %>% 
	dplyr::group_by(Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)

# no events

#=============================================================================
# Pentamidine : Regional breakdown
#=============================================================================
Pentamidine<- dat1[which(dat1$drug_group_final=="Pentamidine"),]
region  <- as.data.frame(table(Pentamidine$Study.Region))

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=Pentamidine[which(Pentamidine$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
						region_output[,2:ncol(region_output)], 
						function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#=============================================================================
# Sitamaquine : Regional breakdown
#=============================================================================
Sitamaquine<- dat1[which(dat1$drug_group_final=="Sitamaquine"),]

Sitamaquine%>% 
	dplyr::group_by(Study.Region) %>%
	dplyr::summarise(
		n_arms = length(Tag),
		n_patients = sum(n_treated, na.rm=TRUE),
		n_deaths = sum(n_deaths_month, na.rm=TRUE)
	)


# Only 1 event from ISC
table(Sitamaquine$Study.Region, Sitamaquine$n_deaths_month)

# Calculation of rates for each of the geographical locations seaprately 
region  <- as.data.frame(table(Sitamaquine$Study.Region))
region  <- region[which(region$Var1=="India Subcontinent"),]
region  <- droplevels(region)

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=Sitamaquine[which(Sitamaquine$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers		= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			region_output <- cbind(est,region_output)
		}	
region_output <- as.data.frame(t(region_output))

region_output[,2:ncol(region_output)] <- lapply(
							region_output[,2:ncol(region_output)], 
							function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
region_output$n_p_t		<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)
#=============================================================================
# Other drugs : Regional breakdown
#=============================================================================
other<- dat1[which(dat1$drug_group_final=="Other"),]

# Only 1 event from ISC
table(other$Study.Region, other$n_deaths_month)
# No events

# End script (Not run)
