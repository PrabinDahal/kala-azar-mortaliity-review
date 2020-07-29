#######################################################################################
# Title		:	Systematic literature review of Mortality in VL clinical trials
# Data version	:	24-Oct-2019
# Data source	: 	Supplemental file 2 
# Task		: 	Estimating mortality rate within 30 days of treatment initiation
########################################################################################
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

# randomisation status 
dat1$randomisation_status <- ifelse(dat1$randomisation=="Randomised","Randomised","Otherwise")

#==================================================
# Overall incidence rate of mortality
#==================================================
(overall_rate <- metarate(n_deaths_month, person_time, Study, 
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
			n_arms	=	length(overall_rate$data$ Study),
			n_treated	=	sum(overall_rate$data$n_treated),
			n_events	=	sum(overall_rate$data$n_deaths_month),
			n_PT		=	sum(overall_rate$data$ person_time),
			fe		=	exp(overall_rate$TE.fixed)*1000,
			fe_l95	=	exp(overall_rate$lower.fixed)*1000, 
			fe_u95	=	exp(overall_rate$upper.fixed)*1000 , 
			re		=	exp(overall_rate$TE.random)*1000 , 
			re_l95	=	exp(overall_rate$lower.random)*1000, 
			re_u95	=	exp(overall_rate$upper.random)*1000,
			i2		=	round(overall_rate$I2*100,3),
			predict_l95	=	exp(overall_rate$lower.predict)*1000, 
			predict_u95	=	exp(overall_rate$upper.predict)*1000 , 
			eggers	= 	e_t$p.value,
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
# Meta regression analysis to assess sub-group heterogenity 
#===================================================================
metareg(overall_rate,HIV)		
metareg(overall_rate,age_range)
metareg(overall_rate,randomisation_status)
metareg(overall_rate,Study.Region) 	# returns convergence warning

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

		results <- metarate(n_deaths_month, person_time, Study, 
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
		e_t <- metabias(results, method.bias="linreg", k.min=2) 
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))
		#str(k)

		est <-	rbind(
				drug		= 	results$data$drug_group_final[1],
				n_arms	=	length(results $data$ Study),
				n_treated	=	sum(results $data$n_treated),
				n_events	=	sum(results $data$n_deaths_month),
				n_PT		=	sum(results $data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				#se_se	=	results$seTE.fixed,
				fe_l95	=	exp(results$lower.fixed)*1000, 
				fe_u95	=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				#re_se	=	results$seTE.random, 
				re_l95	=	exp(results$lower.random)*1000, 
				re_u95	=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers	= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			drug_output <- cbind(est,drug_output)
		}
#[1] 4
#[1] 6
#[1] 8
#[1] 11
#[1] 12

drug_output <- as.data.frame(t(drug_output))

drug_output[,2:ncol(drug_output)] <- lapply(
							drug_output[,2:ncol(drug_output)], 
							function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
drug_output$n_p_t		<- paste(	drug_output$n_arms,drug_output$n_treated,drug_output$n_events, sep="/")
drug_output$fixed		<- paste(	round(drug_output$fe,4),"[", round(drug_output$fe_l95,4),"-", round(drug_output$fe_u95,4), "]", sep="")
drug_output$random	<- paste(	round(drug_output$re,4),"[", round(drug_output$re_l95,4),"-", round(drug_output$re_u95,4), "]", sep="")
drug_output$pred		<- paste(	round(drug_output$predict_l95,4),"-", round(drug_output$predict_u95,4), sep="")
drug_output$copas		<- paste(	round(drug_output$copas_adj,4),"[", round(drug_output$copas_l95,4),"-", round(drug_output$copas_u95,4), "]", sep="")

drug_output <- subset(drug_output, 
					select=c("drug","n_p_t","fixed","random","i2","pred","copas")
				)
print(drug_output)

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
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$Study.Region==region$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		#str(results)

		# Egger's test
		e_t <- metabias(results, method.bias="linreg", k.min=2) 
	
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
				fe_l95	=	exp(results$lower.fixed)*1000, 
				fe_u95	=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95	=	exp(results$lower.random)*1000, 
				re_u95	=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers	= 	e_t$p.value,
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
region_output$n_p_t		<- paste(	region_output$n_arms,region_output$n_treated,region_output$n_events, sep="/")
region_output$fixed		<- paste(	round(region_output$fe,4),"[", round(region_output$fe_l95,4),"-", round(region_output$fe_u95,4), "]", sep="")
region_output$random		<- paste(	round(region_output$re,4),"[", round(region_output$re_l95,4),"-", round(region_output$re_u95,4), "]", sep="")
region_output$pred		<- paste(	round(region_output$predict_l95,4),"-", round(region_output$predict_u95,4), sep="")
region_output$copas		<- paste(	round(region_output$copas_adj,4),"[", round(region_output$copas_l95,4),"-", round(region_output$copas_u95,4), "]", sep="")

region_output <- subset(region_output, 
					select=c("region_name","n_p_t","fixed","random","i2","pred","copas")
				)
print(region_output)

#===================================================================
# Calculation of rates by HIV status as study exclusion criteria
#===================================================================
metareg(overall_rate,HIV)

# Null matrix to store results
hiv_output	<- NULL
for (i in 1: nrow(hiv_status)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$HIV==hiv_status$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
			#str(results)

		# Egger's test
		e_t <- metabias(results, method.bias="linreg", k.min=2) 
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				hiv_status	= 	results$data$HIV[1],
				n_arms	=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95	=	exp(results$lower.fixed)*1000, 
				fe_u95	=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95	=	exp(results$lower.random)*1000, 
				re_u95	=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers	= 	e_t$p.value,
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

#============================================================================================
# Comparision of the overall rate of death in 30 days in east africa vs Indian Sub-Continent
#============================================================================================
east_africa <- dat1[which(dat1$Study.Region=="Eastern Africa"),]
east_africa  <-droplevels(east_africa)
indian_subcontinent <- dat1[which(dat1$Study.Region=="India Subcontinent"),]
indian_subcontinent <-droplevels(indian_subcontinent)

(rate_east_africa  <- metarate(n_deaths_month, person_time, Study, 
			data=east_africa,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			prediction=TRUE,
			control=list(stepadj=0.5, maxiter=1000)
		)
	)

(rate_ISC  <- metarate(n_deaths_month, person_time, Study, 
			data=indian_subcontinent,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			control=list(stepadj=0.5, maxiter=1000),
			prediction=TRUE
		)
	)

# Rate difference between ISC and East Africa
rate_diff <- ( exp(rate_east_africa$TE.random) - exp(rate_ISC$TE.random))

# Delta method for estimating difference and express per 1000 person days
var_rate1 <-  ((exp(rate_east_africa$TE.random)^2)*(rate_east_africa$seTE.random ^2))
var_rate2 <-  ((exp(rate_ISC$TE.random)^2)*(rate_ISC$seTE.random ^2))
se_rate_diff <- sqrt ( var_rate1  + var_rate2 )

rate_diff_lower  <- (rate_diff - 1.96*se_rate_diff )
rate_diff_upper  <- (rate_diff +1.96*se_rate_diff)

# Rate difference and 95% CI
cbind(rate_diff*1000,rate_diff_lower*1000,rate_diff_upper*1000)

#============================================================================
# Comparision of the overall rate of death between HIV vs non-HIV study arms
#============================================================================
HIV <- dat1[which(dat1$HIV=="Included"),]
HIV_ex <- dat1[which(dat1$HIV=="Excluded"),]

(rate_HIV <- metarate(n_deaths_month, person_time, Study, 
			data=HIV ,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			prediction=TRUE,
			control=list(stepadj=0.5, maxiter=1000)
		)
	)

(rate_HIV_ex <- metarate(n_deaths_month, person_time, Study, 
			data=HIV_ex,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			control=list(stepadj=0.5, maxiter=1000),
			prediction=TRUE
			)
		)

# Rate difference between HIV and non-HIV study arms
rate_diff_hiv <- (exp(rate_HIV$TE.random) - exp(rate_HIV_ex$TE.random))

# Delta method for estimating difference and express per 1000 person days
var_rate1_hiv <-  ((exp(rate_HIV$TE.random)^2)*(rate_HIV$seTE.random ^2))
var_rate2_hiv <-  ((exp(rate_HIV_ex$TE.random)^2)*(rate_HIV_ex$seTE.random ^2))
se_rate_diff_hiv <- sqrt ( var_rate1_hiv  + var_rate2_hiv )

rate_diff_hiv_lower  <- (rate_diff_hiv - 1.96*se_rate_diff_hiv )
rate_diff_hiv_upper  <- (rate_diff_hiv + 1.96*se_rate_diff_hiv)

# Rate difference and 95% CI
cbind(rate_diff_hiv*1000,rate_diff_hiv_lower*1000,rate_diff_hiv_upper*1000)

#=====================================================================================
# Comparision between ISC v Eastern Africa in study arms which excluded HIV patients
#=====================================================================================
ISC_HIV_ex <- dat1[which(dat1$HIV=="Excluded" & dat1$Study.Region=="India Subcontinent"),]
EA_HIV_ex <- dat1[which(dat1$HIV=="Excluded" & dat1$Study.Region=="Eastern Africa"),]

(rate_HIV_ex_ISC <- metarate(n_deaths_month, person_time, Study, 
			data=ISC_HIV_ex,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			control=list(stepadj=0.5, maxiter=1000),
			prediction=TRUE
			)
		)

(rate_HIV_ex_EA<- metarate(n_deaths_month, person_time, Study, 
			data=EA_HIV_ex,
			method="GLMM", 
			sm="IRLN",
			irscale=1000,
			irunit="person-days",
			control=list(stepadj=0.5, maxiter=1000),
			prediction=TRUE
			)
		)

table( dat1$Study.Region,dat1$HIV)

#========================================================================
# Calculation of rates for RCts vs non-RCts 
#========================================================================
dat1$rand <- ifelse(dat1$randomisation=="Randomised","Randomised",
			ifelse(dat1$randomisation=="Not Specified","Not Specified","Otherwise"))
random <- as.data.frame(table(dat1$rand))

# Null matrix to store results
rand_output	<- NULL

for (i in 1: nrow(random )){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$rand==random$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		#str(results)

		# Egger's test
		e_t <- metabias(results, method.bias="linreg", k.min=2) 
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				randomisation= 	results$data$rand[1],
				n_arms	=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95	=	exp(results$lower.fixed)*1000, 
				fe_u95	=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95	=	exp(results$lower.random)*1000, 
				re_u95	=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers	= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			rand_output <- cbind(est,rand_output)
		}	
rand_output <- as.data.frame(t(rand_output))

rand_output[,2:ncol(rand_output)] <- lapply(
							rand_output[,2:ncol(rand_output)], 
							function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
rand_output$n_p_t		<- paste(rand_output$n_arms,rand_output$n_treated,rand_output$n_events, sep="/")
rand_output$fixed		<- paste(round(rand_output$fe,4),"[", round(rand_output$fe_l95,4),"-", round(rand_output$fe_u95,4), "]", sep="")
rand_output$random	<- paste(round(rand_output$re,4),"[", round(rand_output$re_l95,4),"-", round(rand_output$re_u95,4), "]", sep="")
rand_output$pred		<- paste(round(rand_output$predict_l95,4),"-", round(rand_output$predict_u95,4), sep="")
rand_output$copas		<- paste(round(rand_output$copas_adj,4),"[", round(rand_output$copas_l95,4),"-", round(rand_output$copas_u95,4), "]", sep="")

rand_output <- subset(rand_output, 
					select=c("randomisation","n_p_t","fixed","random","i2","pred","copas")
				)
print(rand_output)

#=================================================
# Calculation of rates by different age categories
#=================================================
age	<- as.data.frame(table(dat1$age_range))

# Null matrix to store results
age_output	<- NULL

for (i in 1: nrow(age)){
		results <- metarate(n_deaths_month, person_time, Study, 
				data=dat1[which(dat1$age_range==age$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)
		#str(results)

		# Egger's test
		e_t <- metabias(results, method.bias="linreg", k.min=2) 
	
		# Copas selection model
		# Sometimes returns warnings when estimating the adjusted solution
		tryCatch({k <- 	summary.copas(copas(results))},warning=function(w) print(i))

		est <-	rbind(
				age_range	= 	results$data$age_range[1],
				n_arms	=	length(results$data$Study),
				n_treated	=	sum(results$data$n_treated),
				n_events	=	sum(results$data$n_deaths_month),
				n_PT		=	sum(results$data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				fe_l95	=	exp(results$lower.fixed)*1000, 
				fe_u95	=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				re_l95	=	exp(results$lower.random)*1000, 
				re_u95	=	exp(results$upper.random)*1000,
				i2		=	round(results$I2*100,3),
				predict_l95	=	exp(results$lower.predict)*1000, 
				predict_u95	=	exp(results$upper.predict)*1000 , 
				eggers	= 	e_t$p.value,
				copas_adj	=	exp(summary.copas(copas(results))$adjust$TE)*1000, 
				copas_l95	=	exp(summary.copas(copas(results))$adjust$lower)*1000, 
				copas_u95	=	exp(summary.copas(copas(results))$adjust$upper)*1000
				)
			age_output <- cbind(est,age_output)
		}	
age_output <- as.data.frame(t(age_output))

age_output[,2:ncol(age_output)] <- lapply(
							age_output[,2:ncol(age_output)], 
							function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
age_output$n_p_t		<- paste(age_output$n_arms,age_output$n_treated,age_output$n_events, sep="/")
age_output$fixed		<- paste(round(age_output$fe,4),"[", round(age_output$fe_l95,4),"-", round(age_output$fe_u95,4), "]", sep="")
age_output$random		<- paste(round(age_output$re,4),"[", round(age_output$re_l95,4),"-", round(age_output$re_u95,4), "]", sep="")
age_output$pred		<- paste(round(age_output$predict_l95,4),"-", round(age_output$predict_u95,4), sep="")
age_output$copas		<- paste(round(age_output$copas_adj,4),"[", round(age_output$copas_l95,4),"-", round(age_output$copas_u95,4), "]", sep="")

age_output <- subset(age_output, 
					select=c("age_range","n_p_t","fixed","random","i2","pred","copas")
				)
print(age_output)

#======================================================================================
# Calculation of rates for each of the study drugs seaprately from randomised studies
#======================================================================================
randomised_studies <- dat1[which(dat1$randomisation=="Randomised"),]

drug <- as.data.frame(table(randomised_studies$drug_group_final))
drug <- drug[which(drug$Var1!="L-AmB (Single dose)"),]
drug <- drug[which(drug$Var1!="Miltefosine + Paromomycin"),]
drug <- drug[which(drug$Var1!="Other"),]
drug <- droplevels(drug )

# Null matrix to store results
drug_output	<- NULL

for (i in 1: nrow(drug)){

		results <- metarate(n_deaths_month, person_time, Study, 
				data=randomised_studies[which(randomised_studies$drug_group_final==drug$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				drug		= 	results$data$drug_group_final[1],
				n_arms	=	length(results $data$ Study),
				n_treated	=	sum(results $data$n_treated),
				n_events	=	sum(results $data$n_deaths_month),
				n_PT		=	sum(results $data$ person_time),
				re		=	exp(results$TE.random)*1000 , 
				re_l95	=	exp(results$lower.random)*1000, 
				re_u95	=	exp(results$upper.random)*1000
				)
			drug_output <- cbind(est,drug_output)
		}
drug_output <- as.data.frame(t(drug_output))

drug_output[,2:ncol(drug_output)] <- lapply(
							drug_output[,2:ncol(drug_output)], 
							function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
drug_output$n_p_t		<- paste(	drug_output$n_arms,drug_output$n_treated,drug_output$n_events, sep="/")
drug_output$random	<- paste(	round(drug_output$re,4),"[", round(drug_output$re_l95,4),"-", round(drug_output$re_u95,4), "]", sep="")

drug_output <- subset(drug_output, 
					select=c("drug","n_p_t","random")
				)
print(drug_output)

#======================================================================================
# Calculation of rates for each of the study drugs seaprately from non-randomised studies
#======================================================================================
n_randomised_studies <- dat1[which(dat1$randomisation!="Randomised"),]

drug <- as.data.frame(table(n_randomised_studies$drug_group_final))
drug <- drug[which(drug$Var1!="AmBb-lipid"),]
drug <- drug[which(drug$Var1!="Paromomycin"),]
drug <- drug[which(drug$Var1!="Sitamaquine"),]
drug <- drug[which(drug$Var1!="Miltefosine + Paromomycin"),]
drug <- drug[which(drug$Var1!="Other"),]
drug <- droplevels(drug )

# Null matrix to store results
drug_output	<- NULL

for (i in 1: nrow(drug)){

		results <- metarate(n_deaths_month, person_time, Study, 
				data=n_randomised_studies[which(n_randomised_studies$drug_group_final==drug$Var1[i]),] ,
				method="GLMM", 
				sm="IRLN",
				irscale=1000,
				irunit="person-days",
				prediction=TRUE,
				control=list(stepadj=0.5, maxiter=1000)
				)

		est <-	rbind(
				drug		= 	results$data$drug_group_final[1],
				n_arms	=	length(results $data$ Study),
				n_treated	=	sum(results $data$n_treated),
				n_events	=	sum(results $data$n_deaths_month),
				n_PT		=	sum(results $data$ person_time),
				re		=	exp(results$TE.random)*1000 , 
				re_l95	=	exp(results$lower.random)*1000, 
				re_u95	=	exp(results$upper.random)*1000
				)
			drug_output <- cbind(est,drug_output)
		}
drug_output <- as.data.frame(t(drug_output))

drug_output[,2:ncol(drug_output)] <- lapply(
							drug_output[,2:ncol(drug_output)], 
							function(x) as.numeric(as.character(x))
						)
# Tidying up the results and export as a table
drug_output$n_p_t		<- paste(	drug_output$n_arms,drug_output$n_treated,drug_output$n_events, sep="/")
drug_output$random	<- paste(	round(drug_output$re,4),"[", round(drug_output$re_l95,4),"-", round(drug_output$re_u95,4), "]", sep="")

drug_output <- subset(drug_output, 
					select=c("drug","n_p_t","random")
				)
print(drug_output)

# End script (Not run)