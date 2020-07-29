######################################################################################################################### Title		:	Systematic literature review of Mortality in VL clinical trials
# Data version	:	24-Oct-2019
# Script Date	: 	15-Jan-2020
# Data source	: 	Supplemental file 2 
# Task		: 	Sensitvity analysis for worst case outcomes within 30 days of randomisation
#########################################################################################################################rm(list=ls())
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

#==================================================
# Overall incidence rate of mortality
#==================================================
(overall_rate <- metarate(number_of_deaths__WORST, person_time, Study, 
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
# Egger's test for assessing publication bias
e_t <- metabias(overall_rate, method.bias="linreg") 	

(est <-	cbind(
			n_arms		=	length(overall_rate$data$Study),
			n_treated	=	sum(overall_rate$data$n_treated),
			n_events	=	sum(overall_rate$data$number_of_deaths__WORST),
			n_PT		=	sum(overall_rate$data$person_time),
			fe		=	exp(overall_rate$TE.fixed)*1000,
			fe_l95		=	exp(overall_rate$lower.fixed)*1000, 
			fe_u95		=	exp(overall_rate$upper.fixed)*1000 , 
			re		=	exp(overall_rate$TE.random)*1000 , 
			re_l95		=	exp(overall_rate$lower.random)*1000, 
			re_u95		=	exp(overall_rate$upper.random)*1000,
			i2		=	round(overall_rate$I2*100,3),
			predict_l95	=	exp(overall_rate$lower.predict)*1000, 
			predict_u95	=	exp(overall_rate$upper.predict)*1000 , 
			eggers		= 	e_t$p.value,
			copas_adj	=	exp(summary.copas(copas(overall_rate))$adjust$TE)*1000, 
			copas_l95	=	exp(summary.copas(copas(overall_rate))$adjust$lower)*1000, 
			copas_u95	=	exp(summary.copas(copas(overall_rate))$adjust$upper)*1000
		)
	)
est <- as.data.frame(est)

# Tidying up the results and export as a table
est$n_p_t	<- paste(est$n_arms,est$n_treated,est$n_events, sep="/")
est$fixed	<- paste(round(est$fe,4),"[", round(est$fe_l95,4),"-", round(est$fe_u95,4), "]", sep="")
est$random	<- paste(round(est$re,4),"[", round(est$re_l95,4),"-", round(est$re_u95,4), "]", sep="")
est$pred	<- paste(round(est$predict_l95,4),"-", round(est$predict_u95,4), sep="")
est$copas	<- paste(round(est$copas_adj,4),"[", round(est$copas_l95,4),"-", round(est$copas_u95,4), "]", sep="")

est <- subset(est, 
		select=c("n_p_t","fixed","random","i2","pred","copas")
		)
print(est)

#===================================================================
# Calculation of rates for each of the study drugs seaprately 
#===================================================================
drug <- as.data.frame(table(dat1$drug_group_final))
drug <- drug[which(drug$Var1!="Miltefosine + Paromomycin"),]
drug <- drug[which(drug$Var1!="Other"),]
drug <- droplevels(drug )

# Null matrix to store results
drug_output	<- NULL

	for (i in 1: nrow(drug)){

		results <- metarate(number_of_deaths__WORST, person_time, Study, 
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
			drug_output <- cbind(est,drug_output)
		}
drug_output <- as.data.frame(t(drug_output))

drug_output[,2:ncol(drug_output)] <- lapply(
					drug_output[,2:ncol(drug_output)], 
					function(x) as.numeric(as.character(x))
				)

# Tidying up the results and export as a table
drug_output$n_p_t	<- paste(drug_output$n_arms,drug_output$n_treated,drug_output$n_events, sep="/")
drug_output$fixed	<- paste(round(drug_output$fe,4),"[", round(drug_output$fe_l95,4),"-", round(drug_output$fe_u95,4), "]", sep="")
drug_output$random	<- paste(round(drug_output$re,4),"[", round(drug_output$re_l95,4),"-", round(drug_output$re_u95,4), "]", sep="")
drug_output$pred	<- paste(round(drug_output$predict_l95,4),"-", round(drug_output$predict_u95,4), sep="")
drug_output$copas	<- paste(round(drug_output$copas_adj,4),"[", round(drug_output$copas_l95,4),"-", round(drug_output$copas_u95,4), "]", sep="")

(drug_output <- subset(drug_output, 
			select=c("drug","n_p_t","fixed","random","i2","pred","copas")
				)
			)
#========================================================================
# Calculation of rates for each of the geographical locations seaprately 
#========================================================================
region  <- as.data.frame(table(dat1$Study.Region))
region   <- region[which(region$Var1!="Multi-Regional"),]
region  <- region[which(region$Var1!="Not stated"),]
region  <- droplevels(region )

# Null matrix to store results
region_output	<- NULL

for (i in 1: nrow(region)){
		results <- metarate(number_of_deaths__WORST, person_time, Study, 
				data=dat1[which(dat1$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		#str(results)

		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		# publication bias
		e_t <- metabias(results, method.bias="linreg", k.min=2) # Egger's test
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				region_name	= 	results$data$Study.Region[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$number_of_deaths__WORST),
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
region_output$n_p_t	<- paste(region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed	<- paste(round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random	<- paste(round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred	<- paste(round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas	<- paste(round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
			select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
			)
print(region_output)

#===================================================================
# Calculation of rates by HIV status as study exclusion criteria
#===================================================================
hiv_status  <- as.data.frame(table(dat1$HIV))

# Null matrix to store results
hiv_output	<- NULL
for (i in 1: nrow(hiv_status)){
		results <- metarate(number_of_deaths__WORST, person_time, Study, 
				data=dat1[which(dat1$HIV==hiv_status$Var1[i]),] ,
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

		est <-	rbind(
				hiv_status	= 	results$data$HIV[1],
				n_arms		=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$number_of_deaths__WORST),
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
			hiv_output <- cbind(est,hiv_output)
		}	
hiv_output <- as.data.frame(t(hiv_output))

hiv_output[,2:ncol(hiv_output)] <- lapply(
					hiv_output[,2:ncol(hiv_output)], 
					function(x) as.numeric(as.character(x))
				)
# Tidying up the results and export as a table
hiv_output$n_p_t		<- paste(hiv_output$n_arms,hiv_output$n_treated,hiv_output$n_events, sep="/")
hiv_output$fixed		<- paste(round(hiv_output$fe,4),"[", round(hiv_output$fe_l95,4),"-", round(hiv_output$fe_u95,4), "]", sep="")
hiv_output$random		<- paste(round(hiv_output$re,4),"[", round(hiv_output$re_l95,4),"-", round(hiv_output$re_u95,4), "]", sep="")
hiv_output$pred		<- paste(round(hiv_output$predict_l95,4),"-", round(hiv_output$predict_u95,4), sep="")
hiv_output$copas		<- paste(round(hiv_output$copas_adj,4),"[", round(hiv_output$copas_l95,4),"-", round(hiv_output$copas_u95,4), "]", sep="")

hiv_output <- subset(hiv_output, 
			select=c("hiv_status","n_p_t","fixed","random","i2","pred","copas")
			)
print(hiv_output)

# End script (Not run)
