#######################################################################################
# Title		:	Systematic literature review of Mortality in VL clinical trials
# Data version	:	24-Oct-2019
# Script Date	: 	15-Jan-2020
# Data source	: 	Supplemental file 2 
# Task		: 	Sensitivity analysis for mortality rate by risk of bias status
########################################################################################
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

#=========================================================
# Read Risk of Bias assessment data for randomised studies
#=========================================================
rob_data <- read_excel("Supplemental file 3_analyses data.xlsx", 
		sheet = "risk_of_bias_randomised")
rob_data$Study <- NULL 

dat1	<- dplyr::inner_join(
			dat1, rob_data,
			by="Tag"
		)

#===================================================================================
# Number of arms, treated, and events by domain and the risk of bias classification
#===================================================================================

# random_sequence_generation
a <- dat1	%>% 
		dplyr::group_by(random_sequence_generation) %>%
			dplyr::summarise(
				n_arms = length(Tag),
				n_patients = sum(n_treated, na.rm=TRUE),
				n_events= sum(n_deaths_month, na.rm=TRUE)
			)

a$domain <- "random_sequence"
a$risk_of_bias <- a$random_sequence_generation

# allocation_concealment
b <- dat1	%>% 
		dplyr::group_by(allocation_concealment) %>%
			dplyr::summarise(
				n_arms = length(Tag),
				n_patients = sum(n_treated, na.rm=TRUE),
				n_events= sum(n_deaths_month, na.rm=TRUE)
			)
b$domain <- "allocation"
b$risk_of_bias <- b$allocation_concealment

# blinding_of_participants_and_personnel
c <- dat1	%>% 
		dplyr::group_by(blinding_of_participants_and_personnel) %>%
			dplyr::summarise(
				n_arms = length(Tag),
				n_patients = sum(n_treated, na.rm=TRUE),
				n_events= sum(n_deaths_month, na.rm=TRUE)
			)
c$domain <- "blinding"
c$risk_of_bias <- c$blinding_of_participants_and_personnel

# blinding_of_outcome_assessment
d <- dat1	%>% 
		dplyr::group_by(blinding_of_outcome_assessment) %>%
			dplyr::summarise(
				n_arms = length(Tag),
				n_patients = sum(n_treated, na.rm=TRUE),
				n_events= sum(n_deaths_month, na.rm=TRUE)
			)
d$domain <- "outcome_blinding"
d$risk_of_bias <- d$blinding_of_outcome_assessment

# incomplete_outcome_data_addressed
e <- dat1	%>% 
		dplyr::group_by(incomplete_outcome_data_addressed) %>%
			dplyr::summarise(
				n_arms = length(Tag),
				n_patients = sum(n_treated, na.rm=TRUE),
				n_events= sum(n_deaths_month, na.rm=TRUE)
			)
e$domain <- "incomplete_outcome"
e$risk_of_bias <- e$incomplete_outcome_data_addressed

# selective_reporting
f <- dat1	%>% 
		dplyr::group_by(selective_reporting) %>%
			dplyr::summarise(
				n_arms = length(Tag),
				n_patients = sum(n_treated, na.rm=TRUE),
				n_events= sum(n_deaths_month, na.rm=TRUE)
			)
f$domain <- "selective"
f$risk_of_bias <- f$selective_reporting

# aems
g <- dat1	%>% 
		dplyr::group_by(aems_in_place) %>%
			dplyr::summarise(
				n_arms = length(Tag),
				n_patients = sum(n_treated, na.rm=TRUE),
				n_events= sum(n_deaths_month, na.rm=TRUE)
			)
g$domain <- "aems"
g$risk_of_bias <- g$aems_in_place

# Combine the number of arms, number of events and patients info
rob_groups <- rbind(a[,2:6], b[,2:6],c[,2:6], d[,2:6], e[,2:6], f[,2:6], g[,2:6])
rob_groups$n_p_t	<- paste(rob_groups$n_arms,rob_groups$n_patients,rob_groups$n_events, sep="/")

rob_groups <- subset(rob_groups,select=c("domain","risk_of_bias","n_p_t"))
print(rob_groups) 		

#===================================================================
# Calculation of overall rate in the randomised studies included
#===================================================================
(overall_rate <- metarate(n_deaths_month, person_time,Study, 
			data=dat1,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			prediction=TRUE,
			control=list(stepadj=0.5, maxiter=1000)
		)
	)

overall_est <-	rbind(
				n_arms		=	length(overall_rate$data$Study),
				n_treated	=	sum(overall_rate$data$n_treated),
				n_events	=	sum(overall_rate$data$n_deaths_month),
				n_PT		=	sum(overall_rate$data$person_time),
				re		=	exp(overall_rate$TE.random)*1000 , 
				re_l95		=	exp(overall_rate$lower.random)*1000, 
				re_u95		=	exp(overall_rate$upper.random)*1000,
				i2		=	round(overall_rate$I2*100,3)
				)
overall_est <- as.data.frame(t(overall_est))

# Tidying up the results and export as a table
overall_est$n_p_t	<- paste(overall_est$n_arms,overall_est$n_treated,overall_est$n_events, sep="/")
overall_est$random<- paste(round(overall_est$re,4),"[", round(overall_est$re_l95,4),"-", round(overall_est$re_u95,4), "]", sep="")

overall_est$domain 		<- "Overall"
overall_est$risk_of_bias	<- "Overall"

overall_est <- subset(overall_est, 
			select=c("domain","risk_of_bias","n_p_t","random","i2")
			)
print(overall_est)
#=========================================================================
# Test of moderator effects: different risk of bias domains
#=========================================================================
metareg(overall_rate,random_sequence_generation)
metareg(overall_rate,allocation_concealment)
metareg(overall_rate,blinding_of_participants_and_personnel)
metareg(overall_rate,blinding_of_outcome_assessment)
metareg(overall_rate,incomplete_outcome_data_addressed)
metareg(overall_rate,selective_reporting)
metareg(overall_rate,aems_in_place)

# Meta-regression returned convergence warnings
# Calculating the rates stratfied by each of the domains

#===================================================================
# Calculation of rates by random_sequence_generation
#===================================================================
bias <- as.data.frame(table(dat1$random_sequence_generation))
bias <- bias [which(bias$Var1!="High"),] 

sequence_output<- NULL

for (i in 1: nrow(bias)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$random_sequence_generation==bias$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				sequence	= 	results$data$random_sequence_generation[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$person_time),
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3)
				)
			sequence_output <- cbind(est,sequence_output)
		}
sequence_output <- as.data.frame(t(sequence_output))
sequence_output[,2:ncol(sequence_output)] <- lapply(
							sequence_output[,2:ncol(sequence_output)], 
							function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
sequence_output$n_p_t	<- paste(sequence_output$n_arms,sequence_output$n_treated,sequence_output$n_events, sep="/")
sequence_output$random	<- paste(round(sequence_output$re,4),"[", round(sequence_output$re_l95,4),"-", round(sequence_output$re_u95,4), "]", sep="")

sequence_output$domain 		<- "random_sequence"
sequence_output$risk_of_bias	<- sequence_output$sequence

sequence_output <- subset(sequence_output, 
					select=c("domain","risk_of_bias","n_p_t","random","i2")
				)
print(sequence_output)

#===================================================================
# Calculation of rates by allocation_concealment
#===================================================================
bias <- as.data.frame(table(dat1$allocation_concealment))

allocation_output<- NULL

for (i in 1: nrow(bias)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$allocation_concealment==bias$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				allocation	= 	results$data$allocation_concealment[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$person_time),
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3)
				)
			allocation_output <- cbind(est,allocation_output)
		}

allocation_output <- as.data.frame(t(allocation_output))
allocation_output[,2:ncol(allocation_output)] <- lapply(
							allocation_output[,2:ncol(allocation_output)], 
							function(x) as.numeric(as.character(x))
						)

# Tidying up the results and export as a table
allocation_output$n_p_t		<- paste(allocation_output$n_arms,allocation_output$n_treated,allocation_output$n_events, sep="/")
allocation_output$random	<- paste(round(allocation_output$re,4),"[", round(allocation_output$re_l95,4),"-", round(allocation_output$re_u95,4), "]", sep="")

allocation_output$domain 	<- "allocation"
allocation_output$risk_of_bias	<- allocation_output$allocation

allocation_output <- subset(allocation_output, 
					select=c("domain","risk_of_bias","n_p_t","random","i2")
				)
print(allocation_output)

#===================================================================
# Calculation of rates by blinding_of_participants_and_personnel
#===================================================================
bias<- as.data.frame(table(dat1$blinding_of_participants_and_personnel))
bias <- bias [which(bias$Var1!="Low"),] # no events

blinding_output<- NULL

for (i in 1: nrow(bias)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$blinding_of_participants_and_personnel==bias$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				blinding	= 	results$data$blinding_of_participants_and_personnel[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$person_time),
				re		=	exp(results$TE.random)*1000 , 
				re_l95		=	exp(results$lower.random)*1000, 
				re_u95		=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3)
				)
			blinding_output <- cbind(est,blinding_output)
		}

blinding_output <- as.data.frame(t(blinding_output))
blinding_output[,2:ncol(blinding_output)] <- lapply(
							blinding_output[,2:ncol(blinding_output)], 
							function(x) as.numeric(as.character(x))
						)

# Tidying up the results and export as a table
blinding_output$n_p_t	<- paste(blinding_output$n_arms,blinding_output$n_treated,blinding_output$n_events, sep="/")
blinding_output$random	<- paste(round(blinding_output$re,4),"[", round(blinding_output$re_l95,4),"-", round(blinding_output$re_u95,4), "]", sep="")

blinding_output$domain 		<- "blinding"
blinding_output$risk_of_bias	<- blinding_output$blinding

blinding_output <- subset(blinding_output, 
					select=c("domain","risk_of_bias","n_p_t","random","i2")
				)
print(blinding_output)

#===================================================================
# Calculation of rates by blinding_of_outcome assessment
#===================================================================
bias<- as.data.frame(table(dat1$blinding_of_outcome_assessment))
bias <- bias [which(bias$Var1!="High"),] # no events
bias <- droplevels(bias)

outcome_blinding_output<- NULL

for (i in 1: nrow(bias)){
		results <- metarate(n_deaths_month, person_time, Study, 

				data=dat1[which(dat1$blinding_of_outcome_assessment==bias$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				outcome_blinding	= 	results$data$blinding_of_outcome_assessment[1],
				n_arms			=	length(results$data$Study),
				n_treated		=	sum(results$data$n_treated),
				n_events		=	sum(results$data$n_deaths_month),
				n_PT			=	sum(results$data$person_time),
				re			=	exp(results$TE.random)*1000 , 
				re_l95			=	exp(results$lower.random)*1000, 
				re_u95			=	exp(results$upper.random)*1000,
				i2			=	round(results$I2*100,3)
				)
			outcome_blinding_output <- cbind(est,outcome_blinding_output)
		}

outcome_blinding_output <- as.data.frame(t(outcome_blinding_output))
outcome_blinding_output[,2:ncol(outcome_blinding_output)] <- lapply(
							outcome_blinding_output[,2:ncol(outcome_blinding_output)], 
							function(x) as.numeric(as.character(x))
						)

# Tidying up the results and export as a table
outcome_blinding_output$n_p_t	<- paste(outcome_blinding_output$n_arms,outcome_blinding_output$n_treated,outcome_blinding_output$n_events, sep="/")
outcome_blinding_output$random	<- paste(round(outcome_blinding_output$re,4),"[", round(outcome_blinding_output$re_l95,4),"-", round(outcome_blinding_output$re_u95,4), "]", sep="")

outcome_blinding_output$domain 		<- "outcome_blinding"
outcome_blinding_output$risk_of_bias	<- outcome_blinding_output$outcome_blinding

outcome_blinding_output <- subset(outcome_blinding_output, 
					select=c("domain","risk_of_bias","n_p_t","random","i2")
				)
print(outcome_blinding_output)

#===================================================================
# Calculation of rates by incomplete_outcome_data_addressed
#===================================================================
bias<- as.data.frame(table(dat1$incomplete_outcome_data_addressed))
bias <- bias [which(bias$Var1!="Unclear"),] # no events
table(dat1$incomplete_outcome_data_addressed, dat1$n_deaths_month)

incomplete_output<- NULL

for (i in 1: nrow(bias)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$incomplete_outcome_data_addressed==bias$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				incomplete_outcome	= 	results$data$incomplete_outcome_data_addressed[1],
				n_arms			=	length(results$data$Study),
				n_treated		=	sum(results$data$n_treated),
				n_events		=	sum(results$data$n_deaths_month),
				n_PT			=	sum(results$data$person_time),
				re			=	exp(results$TE.random)*1000 , 
				re_l95			=	exp(results$lower.random)*1000, 
				re_u95			=	exp(results$upper.random)*1000,
				i2			=	round(results$I2*100,3)
				)
			incomplete_output <- cbind(est,incomplete_output)
		}

incomplete_output <- as.data.frame(t(incomplete_output))
incomplete_output[,2:ncol(incomplete_output)] <- lapply(
							incomplete_output[,2:ncol(incomplete_output)], 
							function(x) as.numeric(as.character(x))
						)

# Tidying up the results and export as a table
incomplete_output$n_p_t	<- paste(incomplete_output$n_arms,incomplete_output$n_treated,incomplete_output$n_events, sep="/")
incomplete_output$random<- paste(round(incomplete_output$re,4),"[", round(incomplete_output$re_l95,4),"-", round(incomplete_output$re_u95,4), "]", sep="")

incomplete_output$domain 		<- "incomplete_outcome"
incomplete_output$risk_of_bias	<- incomplete_output$incomplete_outcome

incomplete_output <- subset(incomplete_output, 
					select=c("domain","risk_of_bias","n_p_t","random","i2")
				)
print(incomplete_output)

#===================================================================
# Calculation of rates by selective_reporting
#===================================================================
bias<- as.data.frame(table(dat1$selective_reporting))
bias <- bias [which(bias$Var1!="High"),] # no events
table(dat1$selective_reporting, dat1$n_deaths_month)

selective_output<- NULL

for (i in 1: nrow(bias)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$selective_reporting==bias$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				selective		= 	results$data$selective_reporting[1],
				n_arms			=	length(results$data$Study),
				n_treated		=	sum(results$data$n_treated),
				n_events		=	sum(results$data$n_deaths_month),
				n_PT			=	sum(results$data$person_time),
				re			=	exp(results$TE.random)*1000 , 
				re_l95			=	exp(results$lower.random)*1000, 
				re_u95			=	exp(results$upper.random)*1000,
				i2			=	round(results$I2*100,3)
				)
			selective_output <- cbind(est,selective_output)
		}

selective_output <- as.data.frame(t(selective_output))
selective_output[,2:ncol(selective_output)] <- lapply(
							selective_output[,2:ncol(selective_output)], 
							function(x) as.numeric(as.character(x))
						)

# Tidying up the results and export as a table
selective_output$n_p_t	<- paste(selective_output$n_arms,selective_output$n_treated,selective_output$n_events, sep="/")
selective_output$random<- paste(round(selective_output$re,4),"[", round(selective_output$re_l95,4),"-", round(selective_output$re_u95,4), "]", sep="")

selective_output$domain 		<- "selective"
selective_output$risk_of_bias	<- selective_output$selective

selective_output <- subset(selective_output, 
					select=c("domain","risk_of_bias","n_p_t","random","i2")
				)
print(selective_output)

#===================================================================
# Calculation of rates by adverse events monitoring system status
#===================================================================
bias<- as.data.frame(table(dat1$aems_in_place))
bias <- bias [which(bias$Var1!="Unclear"),] # no events
table(dat1$aems_in_place, dat1$n_deaths_month)

aems_output<- NULL

for (i in 1: nrow(bias)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$aems_in_place==bias$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				aems			= 	results$data$aems_in_place[1],
				n_arms			=	length(results$data$Study),
				n_treated		=	sum(results$data$n_treated),
				n_events		=	sum(results$data$n_deaths_month),
				n_PT			=	sum(results$data$person_time),
				re			=	exp(results$TE.random)*1000 , 
				re_l95			=	exp(results$lower.random)*1000, 
				re_u95			=	exp(results$upper.random)*1000,
				i2			=	round(results$I2*100,3)
				)
			aems_output <- cbind(est,aems_output)
		}

aems_output <- as.data.frame(t(aems_output))
aems_output[,2:ncol(aems_output)] <- lapply(
							aems_output[,2:ncol(aems_output)], 
							function(x) as.numeric(as.character(x))
						)

# Tidying up the results and export as a table
aems_output$n_p_t	<- paste(aems_output$n_arms,aems_output$n_treated,aems_output$n_events, sep="/")
aems_output$random<- paste(round(aems_output$re,4),"[", round(aems_output$re_l95,4),"-", round(aems_output$re_u95,4), "]", sep="")

aems_output$domain 		<- "aems"
aems_output$risk_of_bias	<- aems_output$aems

aems_output <- subset(aems_output, 
					select=c("domain","risk_of_bias","n_p_t","random","i2")
				)
print(aems_output)

#===================================================================
# Complie the results for each domain in a single table format
#===================================================================
OUTPUT <- rbind(overall_est,aems_output,allocation_output,blinding_output,incomplete_output,
			outcome_blinding_output,selective_output,sequence_output
		)
print(OUTPUT)

# End script (Not run)
