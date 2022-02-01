##install.packages("devtools")
devtools::install_github('HonghyokKim/CGPSspatialmatch')
library(dplyr) ## CGPSsm needs %>% function in dplyr
library(mgcv) ## CGPSsm needs mgcv package for Generalized Additive Model
## DO NOT LOAD library(xgboost) because xgboost package has its own "slice" function and it colides with slice function in dplyr
library(wCorr) ## CGPSsm needs wCorr package to compute correlation
library(gnm) ## CGPSsm needs gnm package for conditional regression for disease model
library(CGPSspatialmatch) ## Load CGPSsm package

data(StrokePPR_CT) ##Load sample data. This is dataset sampled from the dataset used in the study. As such, this sample code will not give you the same results of the paper. 



### Dataset for fitting Propensity Score Model (restrict <=10km from refineries to improve unmeasured spatial confounding adjustment)
PSdataset<-subset(StrokePPR_CT,exposure_10km==1)


## See functions in CGPSspatialmatch
ls("package:CGPSspatialmatch")

########One-to-One Nearest Neighbor Matching with Replacement
########One-to-One Nearest Neighbor Matching with Replacement
### CGPSsm GAM ###


### ?estimateATT will give you more details
CGPSsm_GAM_Nearest<-estimateATT(dataset=StrokePPR_CT,  ###original dataset###
                                PSdataset=PSdataset,   ###dataset for PS model###
                                bexp="exposure_5km",   ###binary exposure status variable###
                                cexp="logMetric_5km",  ###continuous exposure variable for CGPS, log-transformed###
                                fmethod.replace=TRUE,  ###one-to-n distance matching: replacement? TRUE=replacement (USE THIS); FALSE=without replacement, which has not been theorized in the paper###
                                distbuf=0.1,           ###spatial distance threshold for one-to-n distance matching
                                exp.included=FALSE,    ###USE FALSE; TRUE=include exposed locations in one-to-n distance matching (DO NOT USE THIS)
                                long="Longitude",lat="Latitude", ###Indicate coordinate variables 
                                PS.method="mgcv.GAM",  ###Select PS model fitting method
                                
                                ####Input PS model formula like y=a+b+c+d... if PS.method=mgcv.GAM
                                PS.formula="exposure_5km~Pop_65_84_18over_pt +
                                                         Hispanic_pt+NH_White_pt+EduLHigh_18over_pt+Insured_19over_pt+Poverty_18over_pt+Pop_45_64_18over_pt+Pop_85over_18over_pt+CSMOKING+s(Longitude,Latitude)",
                                
                                CGPS.method="mgcv.GAM", ###Select CGPS model fitting method
                                ###Input CGPS model formula like y=a+b+c+d...if CGPS.method=mgcv.GAM
                                CGPS.formula="logMetric_5km~s(Longitude,Latitude)",
                                
                                
                                smethod="nearest", ###GPS matching method
                                caliper_bw=NULL, ###if nearest neighbor matching, use NULL; if nearestcaliper matching, indicate caliper width
                                smethod.replace=TRUE, ###GPS matching method with replacement (TRUE)? or without replacement (FALSE)?
                                
                                ###Input Disease model formula like y=a+b+c...
                                formulaDisease="STROKE~Metric_5km",
                                ###Family for disease model
                                family="gaussian",
                                bs.N=1000, ###number of bootstrapped samples
                                bs.replace=TRUE, ### USE TRUE for bootstrapping (with replacement) Do not use FALSE
                                ###List measured potential confounders
                                varilist=c("Female_18over_pt","Pop_20_24_18over_pt","Pop_25_44_18over_pt","Pop_45_64_18over_pt","Pop_65_84_18over_pt","Pop_85over_18over_pt",
                                           "Hispanic_pt","NH_White_pt","NH_Black_pt","EduLHigh_18over_pt","Poverty_18over_pt","Insured_19over_pt","MHI","CSMOKING"),
                                modelinfo=TRUE ### FALSE=not showing details in disease model objects; TRUE=showing them
)


#### coefficient estimate with 95% CI
CGPSsm_GAM_Nearest$summary
IQR<-IQR(subset(StrokePPR_CT,exposure_5km==1)$Metric_5km)
round(CGPSsm_GAM_Nearest$summary*IQR,3)

#### Information of how many exposed units are matched
CGPSsm_GAM_Nearest$match.info
#### Standardized mean difference
CGPSsm_GAM_Nearest$smd.org ### Original dataset
CGPSsm_GAM_Nearest$smd.matched ### Matched dataset
#### Correlation
CGPSsm_GAM_Nearest$correxp.org  ### Original dataset; no correlation exists in matched dataset

#### Disease model objects (exist only if modelinfo=TRUE)
summary(CGPSsm_GAM_Nearest$modelfit[[1]]) ###first bootstrapped sample
summary(CGPSsm_GAM_Nearest$modelfit[[2]]) ###second bootstrapped sample


#### Distribution from bootstrapped samples
par(mfrow=c(2,1))
hist(CGPSsm_GAM_Nearest$coef_dist,breaks=50,xlab="Estimate",main="Coefficient estimates by boostrapping")
abline(v=mean(CGPSsm_GAM_Nearest$coef_dist),col="red")
abline(v=median(CGPSsm_GAM_Nearest$coef_dist),col="blue")
hist(CGPSsm_GAM_Nearest$se_dist,breaks=50,xlab="Estimate",main="Standard error estimates by boostrapping")
abline(v=mean(CGPSsm_GAM_Nearest$se_dist),col="red")
abline(v=median(CGPSsm_GAM_Nearest$se_dist),col="blue")



########One-to-One Nearest Neighbor Caliper Matching with Replacement
########One-to-One Nearest Neighbor Caliper Matching with Replacement
### CGPSsm GAM ###

### ?estimateATT will give you more details
CGPSsm_GAM_NearestCaliper1<-estimateATT(dataset=StrokePPR_CT,  ###original dataset###
                                        PSdataset=PSdataset,   ###dataset for PS model###
                                        bexp="exposure_5km",   ###binary exposure status variable###
                                        cexp="logMetric_5km",  ###continuous exposure variable for CGPS, log-transformed###
                                        fmethod.replace=TRUE,  ###one-to-n distance matching: replacement? TRUE=replacement (USE THIS); FALSE=without replacement, which has not been theorized in the paper###
                                        distbuf=0.1,           ###spatial distance threshold for one-to-n distance matching
                                        exp.included=FALSE,    ###USE FALSE; TRUE=include exposed locations in one-to-n distance matching (DO NOT USE THIS)
                                        long="Longitude",lat="Latitude", ###Indicate coordinate variables 
                                        PS.method="mgcv.GAM",  ###Select PS model fitting method
                                        
                                        ####Input PS model formula like y=a+b+c+d... if PS.method=mgcv.GAM
                                        PS.formula="exposure_5km~Pop_65_84_18over_pt +
                                                         Hispanic_pt+NH_White_pt+EduLHigh_18over_pt+Insured_19over_pt+Poverty_18over_pt+Pop_45_64_18over_pt+Pop_85over_18over_pt+CSMOKING+s(Longitude,Latitude)",
                                        
                                        CGPS.method="mgcv.GAM", ###Select CGPS model fitting method
                                        ###Input CGPS model formula like y=a+b+c+d...if CGPS.method=mgcv.GAM
                                        CGPS.formula="logMetric_5km~s(Longitude,Latitude)",
                                        
                                        
                                        smethod="nearestcaliper", ###GPS matching method
                                        caliper_bw=0.4, ###if nearest neighbor matching, use NULL; if nearestcaliper matching, indicate caliper width
                                        smethod.replace=TRUE, ###GPS matching method with replacement (TRUE)? or without replacement (FALSE)?
                                        
                                        ###Input Disease model formula like y=a+b+c...
                                        formulaDisease="STROKE~Metric_5km",
                                        ###Family for disease model
                                        family="gaussian",
                                        bs.N=1000, ###number of bootstrapped samples
                                        bs.replace=TRUE, ### USE TRUE for bootstrapping (with replacement) Do not use FALSE
                                        ###List measured potential confounders
                                        varilist=c("Female_18over_pt","Pop_20_24_18over_pt","Pop_25_44_18over_pt","Pop_45_64_18over_pt","Pop_65_84_18over_pt","Pop_85over_18over_pt",
                                                   "Hispanic_pt","NH_White_pt","NH_Black_pt","EduLHigh_18over_pt","Poverty_18over_pt","Insured_19over_pt","MHI","CSMOKING"),
                                        modelinfo=TRUE ### FALSE=not showing details in disease model objects; TRUE=showing them
)


#### coefficient estimate with 95% CI
CGPSsm_GAM_NearestCaliper1$summary
IQR<-IQR(subset(StrokePPR_CT,exposure_5km==1)$Metric_5km)
round(CGPSsm_GAM_NearestCaliper1$summary*IQR,3)
#### Information of how many exposed units are matched
CGPSsm_GAM_NearestCaliper1$match.info
#### Standardized mean difference
CGPSsm_GAM_NearestCaliper1$smd.org ### Original dataset
CGPSsm_GAM_NearestCaliper1$smd.matched ### Matched dataset
#### Correlation
CGPSsm_GAM_NearestCaliper1$correxp.org  ### Original dataset; no correlation exists in matched dataset

#### Disease model objects (exist only if modelinfo=TRUE)
summary(CGPSsm_GAM_NearestCaliper1$modelfit[[1]]) ###first bootstrapped sample
summary(CGPSsm_GAM_NearestCaliper1$modelfit[[2]]) ###second bootstrapped sample


#### Distribution from bootstrapped samples
par(mfrow=c(2,1))
hist(CGPSsm_GAM_NearestCaliper1$coef_dist,breaks=50,xlab="Estimate",main="Coefficient estimates by boostrapping")
abline(v=mean(CGPSsm_GAM_NearestCaliper1$coef_dist),col="red")
abline(v=median(CGPSsm_GAM_NearestCaliper1$coef_dist),col="blue")
hist(CGPSsm_GAM_NearestCaliper1$se_dist,breaks=50,xlab="Estimate",main="Standard error estimates by boostrapping")
abline(v=mean(CGPSsm_GAM_NearestCaliper1$se_dist),col="red")
abline(v=median(CGPSsm_GAM_NearestCaliper1$se_dist),col="blue")





###Create graphs for covariates balance

par(mfrow=c(1,1))
plot(x=c(CGPSsm_GAM_Nearest$smd.org,NA,NA),y=seq(16),bty="n",yaxt="n",ylab="",xlab="Standardized mean difference",xlim=c(-1,1),pch=1,col="black")
axis(2,at=seq(14),las = 2,
     
     c("% Female", "% 20-24 yrs","% 25-44 yrs","% 45-64 yrs","% 65-84 yrs","% 85+ yrs",
       "% Hispanic","% Non-Hispanic White","% Non-Hispanic Black",
       "% Low education","% Poverty", "% Insured","Median Household Income",
       "Smoking"),)
lines(x=c(-0.25,-0.25),y=c(0,14),col="grey",lty=3)
lines(x=c(0.25,0.25),y=c(0,14),col="grey",lty=3)
lines(x=c(-0.1,-0.1),y=c(0,14),col="grey",lty=2)
lines(x=c(0.1,0.1),y=c(0,14),col="grey",lty=2)

points(x=CGPSsm_GAM_Nearest$smd.matched,y=seq(14),col=rgb(0,0,1,0.7,maxColorValue = 1),pch=2)
points(x=CGPSsm_GAM_NearestCaliper1$smd.matched,y=seq(14),col=rgb(1,0,0,0.7,maxColorValue = 1),pch=2)
legend("topright",c("Original (Unmatched)","Nearest neighbor","Nearest caliper")
       ,col=c("black","blue","red"),pch=c(1,2),bty="n")
