########################################################################################################################
# Title		:	Systematic literature review of Mortality in VL clinical trials
# Data version	:	24-Oct-2019
# Script Date	: 	15-Jan-2020
# Data source	: 	Supplemental file 2 
# Task		: 	Estimatating rates at the end of the study follow-up
#########################################################################################################################
rm(list=ls())
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
# Estimating incidence rate of mortality for overall dataset within study follow-up
#==================================================================================
dat1 <- dat[which(!is.na(dat$number_of_deaths)),]

# For studies with missing follow-up, it was assumed that 6-moths was the follow-up
# For studies with known follow-up , total person-time calculated using true follow-up time
# For studies with missing follow-up, total person-time cut at 6-months 

summary(dat1$Follow.up.duration)

dat1$follow_up <- ifelse(is.na(dat1$Follow.up.duration), 180, dat1$Follow.up.duration )

dat1$person_time <- dat1$n_treated*dat1$follow_up

# Using suggestion from the metafor vignette to enchance convergence
# http://127.0.0.1:30486/library/metafor/html/rma.uni.html

#==================================================
# Overall incidence rate of mortality
#==================================================
(overall_rate <- metarate(number_of_deaths, person_time, Study, 
			data=dat1 ,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			control=list(maxiter=1000),
			#control=list(stepadj=0.5, maxiter=1000),
			prediction =TRUE
		)
	)

(est <-	cbind(
			n_arms		=	length(overall_rate$data$Study),
			n_treated	=	sum(overall_rate$data$n_treated),
			n_events	=	sum(overall_rate$data$number_of_deaths__WORST),
			n_PT		=	sum(overall_rate$data$person_time),
			re		=	exp(overall_rate$TE.random)*1000 , 
			re_l95		=	exp(overall_rate$lower.random)*1000, 
			re_u95		=	exp(overall_rate$upper.random)*1000
		)
	)
est <- as.data.frame(est)

# Tidying up the results and export as a table
est$n_p_t	<- paste(est$n_arms,est$n_treated,est$n_events, sep="/")
est$random	<- paste(round(est$re,4),"[", round(est$re_l95,4),"-", round(est$re_u95,4), "]", sep="")
est <- subset(est, 
			select=c("n_p_t","random")
		)
print(est)

#===================================================================
# Calculation of rates for each of the study drugs seaprately 
#===================================================================
drug <- as.data.frame(table(dat1$drug_group_final))
drug <- drug[which(drug$Var1!="Miltefosine + Paromomycin"),]
drug <- drug[which(drug$Var1!="Other"),]
drug <- droplevels(drug )

table(dat1$drug_group_final,dat1$number_of_deaths)

# Null matrix to store results
drug_output	<- NULL

	for (i in 1: nrow(drug)){

		results <- metarate(number_of_deaths, person_time, Study, 
				data=dat1[which(dat1$drug_group_final==drug$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		#str(results)

		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))
		#str(k)

		est <-	rbind(
				drug		= 	results$data$drug_group_final[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$number_of_deaths__WORST),
				n_PT		=	sum(results$data$person_time),
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000
				)
			drug_output <- cbind(est,drug_output)
		}
drug_output <- as.data.frame(t(drug_output))

drug_output[,2:ncol(drug_output)] <- lapply(
							drug_output[,2:ncol(drug_output)], 
							function(x) as.numeric(as.character(x))
					)

# Tidying up the results and export as a table
drug_output$n_p_t	<- paste(drug_output$n_arms,drug_output$n_treated,drug_output$n_events, sep="/")
drug_output$random	<- paste(round(drug_output$re,4),"[", round(drug_output$re_l95,4),"-", round(drug_output$re_u95,4), "]", sep="")

(drug_output <- subset(drug_output, 
			select=c("drug","n_p_t","random")
			)
		)
# End
