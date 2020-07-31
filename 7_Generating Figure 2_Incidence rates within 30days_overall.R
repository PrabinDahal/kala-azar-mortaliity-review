#######################################################################################
# Title		:	Systematic literature review of Mortality in VL clinical trials
# Data version	:	24-Oct-2019
# Script Date	: 	15-Jan-2020
# Data source	: 	Supplemental file 2 
# Task		: 	Generate Figure on the manuscript
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
				n_arms		=	length(results $data$ Study),
				n_treated	=	sum(results $data$n_treated),
				n_events	=	sum(results $data$n_deaths_month),
				n_PT		=	sum(results $data$ person_time),
				fe		=	exp(results$TE.fixed)*1000,
				#se_se		=	results$seTE.fixed,
				fe_l95		=	exp(results$lower.fixed)*1000, 
				fe_u95		=	exp(results$upper.fixed)*1000 , 
				re		=	exp(results$TE.random)*1000 , 
				#re_se		=	results$seTE.random, 
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

#================================================================
# Tidying up drug names for plot
#================================================================
drug_output$drug<- drug_output$drug

drug_output$drug<- recode(drug_output$drug, "PA"="Pentavalent antimony")
drug_output$drug<- recode(drug_output$drug, "AMBd"="Amphotericin B deoxycholate")
drug_output$drug<- recode(drug_output$drug, "AmBb-lipid"= "Amphotericin B fat/lipid/colloid/cholesterol")
drug_output$drug<- recode(drug_output$drug, "L-AmB (multiple)"="L-AmB (multiple dose)")
drug_output$drug<- recode(drug_output$drug, "L-AmB (Single dose)"="L-AmB (single dose)")
drug_output$drug<- recode(drug_output$drug, "L-AmB (Single)combination regimen"="L-AmB (single dose) combination regimen")
drug_output$drug<- recode(drug_output$drug, "L-AmB (multiple) combination regimen"="L-AmB (multiple dose) combination regimen")
drug_output$drug<- recode(drug_output$drug, "PA combination regimen"="Pentavalent antimony combination regimen")

#================================================================
# Fold-difference between Copas adjusted estiamte and from the random effects analysis
#================================================================
drug_output$fold_difference <- drug_output$copas_adj/drug_output$re
drug_output[,c(1,9,15,16,19)]

#================================================================
# Plot the estimates for each of the drugs as a dot plot
#================================================================
library(ggalt)   

a<- ggplot(drug_output, aes(x=reorder(drug, re), y=re)) + 
  	geom_pointrange(aes(ymin=re_l95, ymax=re_u95),size=1,stroke =1,col="#00AFBB")+
	ylab("Incidence rate of mortality (per 1000 person-days)") +
	ylim(0,2)+
	xlab("")+
	ggtitle("Random effects meta-analysis")+
	coord_flip()

b<- ggplot(drug_output, aes(x=reorder(drug, copas_adj), y=copas_adj)) + 
  	geom_pointrange(aes(ymin=copas_l95, ymax=copas_u95),size=1,stroke =1,col="#00AFBB")+
	ylab("Incidence rate of mortality (per 1000 person-days)") +
	ylim(0,5)+
	xlab("")+
	ggtitle("Bias adjusted estimates")+
	coord_flip()

#=====================================================================
# Use cowplot libray to create panel & export as high resolution graph
#=====================================================================
library(cowplot)
setwd("C:/Users/pdahl/Dropbox/_VL AE Working Folder/AE MS/Final AE MS and codes/Results")

tiff(file="Figure_2.tiff", 
            	width=36, 
		height=14, 
		units="cm", 
            	pointsize="13", 
		compression = "lzw+p", 
            	bg="white",
		res=600, 
		antialias = "none"
	)
plot_grid(a,b)
dev.off() # End export

# End script (Not run)
