JR_EEG_Attn <- function(x) {

#Load Packages
library(XML)
library(methods)
library(plyr) 
library(dplyr)
library(fuzzyjoin)
library(zoo)
library(tidyverse)
library(R.matlab)
library(eegkit)
library(reshape2)
library(purrr)
library(lmtest) 
library(naniar)


cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix(NA, mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}
###########################################

parent.folder <- "~/Desktop/Old EEG/JR_Proc/"
#locate epoched EEG folder
psd.folder <- "~/Desktop/Old EEG/JR_Proc/psd_NoRJ"
#locate epoched mff folder with markers
mff.folder <- "~/Desktop/Old EEG/JR_Proc/mff_JR"
#locate epoched EEG folder
ecg.folder <- "~/Desktop/Old EEG/JR_Proc/ecg_JR"
#locate epoched ICA info folder
ica.folder <- "~/Desktop/Old EEG/JR_Proc/ica"
output.folder <- "~/Desktop/Old EEG/JR_Proc/output_RR_2022/"


###DEFINE EEG Fz Bins
low_theta=4
high_theta=6

low_alpha=6
high_alpha=9

#Define frontal electrode #s
FE_R=c('2','3','5','57','59','60')
FE_L=c('9','10', '11','12','13','14')
FE_C=c('6')
FE_total=c('2','3','5','6', '9','10', '11','12','13','14','57','59','60')
Other=c( "18",	"20",	"24",	"28",	"30",	"34",	"35",	"39",	"42",	"44",	"50",	"52",	"58")

### FIND MATCHING EEG IBI FILES ##########
l1=substr(list.files (mff.folder), start = 7, stop = 9)
l2=substr(list.files (ecg.folder), start = 7, stop = 9)
l3=substr(list.files (psd.folder), start = 7, stop = 9)
eeg_ids=intersect(l1,l2)
eeg_ids=intersect(eeg_ids,l3)
eeg_ids=eeg_ids[-1]

############### IBI #####################################

eeg.file=data.frame()   
time.series=data.frame()
topo_attn_sub_theta=data.frame()
topo.file=data.frame()

for(i in eeg_ids){
  setwd(ecg.folder)
  id=substr(i, 7, 9)
  print(i)

  bb.ibitime.file <- list.files(pattern = i)

  bb.ibi.file=bb.ibitime.file[grepl (bb.ibitime.file,pattern = "ibi", ignore.case = T)]
  ibi1 <- map_df(bb.ibi.file,read.table)
  ibi1 = dplyr::rename(ibi1, IBI=V1)


  bb.time.file= bb.ibitime.file[grepl (bb.ibitime.file,pattern = "time",ignore.case = T)]
  time1 <- map_df(bb.time.file,read.table)

  time1=head(time1,-1)
  names(time1)[names(time1)=="V1"]<- "time_real"
  
  if(nrow(ibi1)==nrow(time1)){ 
    print("equal ibi and time files")
    merged=cbind(ibi1,time1)}
  
  if(nrow(ibi1)!=nrow(time1)){ 
    print("unequal ibi and time files")
    merged=cbindPad(ibi1,time1) }

  data=merged

  #Epoch data
  data$time = round(data$time_real,digits = 0)
  data=aggregate(.~time, data=data, mean)
  
  
# Parsing for looks (Ylok or Nlok)
  setwd(mff.folder)
  setwd(list.files(pattern=i))
  
  pattern=intersect(list.files(pattern = "Event", ignore.case = T),
                    list.files(pattern = "User", ignore.case = T))
  xml=xmlToDataFrame(list.files(pattern=pattern))
  
  
  evt_file=xml[c("relativeBeginTime","code")]
  evt_file=evt_file[c(-1:-2),]
  evt_file$relativeBeginTime=as.numeric(evt_file$relativeBeginTime)/1000000
  
  evt_file$time=as.numeric(evt_file$relativeBeginTime)
  
  new=difference_left_join(data,evt_file,by="time",max_dist=.5)
  
  #Create Attention Phases
  
  data2=new
  
  #Add dummy rows of at least 5 sec of early data in case first look comes before that
  data2=data2 %>% 
    add_row(.before=1) %>%
    add_row(.before=1) %>%
    add_row(.before=1) %>%
    add_row(.before=1) %>%
  replace_with_na(replace = list(code = c("Rich","SEIZ","net")))
           
  
  #create forward moving medians of IBI - the good one! 
  #median of 5 preceeding IBIs
  data2$rollmed_preceed=rollapply(data2$IBI,width=5,median,fill=NA,na.rm = TRUE,align="right")
  #create backward moving medians of IBI
  #median of 5 subsequent IBIs
  data2$rollmed_subseq=rollapply(data2$IBI,width=5,median,fill=NA,align="left")
  
  #create variable of looking periods
  data2=data2 %>% mutate (code2 = na.locf(code,na.rm = FALSE))
  
  #Ok - now delete dummy rows
  data2 = data2[-(1:4),]
  
  data2 = data2 %>% select (-c(relativeBeginTime, time.y))
  
  #create variable for pre-stim value which is the median of 5 preceeding IBIs
  #following a Ylook
  data2$prestim=ifelse(data2$code=="Ylok",lag(data2$rollmed_preceed),NA)
  #data2$prestim_alt=ifelse(data2$code=="Ylok" & data2$time.x > 5 ,data2$rollmed_preceed,NA)
  
  data2$poststim=ifelse(data2$code=="Ylok",data2$rollmed_subseq,NA)
  
  #prestim rolling
  data2=data2 %>% mutate (prestim = na.locf(prestim,na.rm = FALSE),
                          poststim = na.locf(poststim,na.rm = FALSE))

  
  #variable that counts every instance during a Ylok period that the corresponding 
  #IBI is shorter than the pre-stim median IBI:for stop phases
  data2$c1=ifelse(data2$code2=="Ylok" & data2$prestim>data2$IBI, 1,0) 
  #rolling sum variable based on previous variable to get the 5 consecutive IBIs
  #shorter than pre-stim IBI:for stop phases
  data2$c2=rollapply(data2$c1,width=5,sum,fill=NA,align="right")
  
  #variable that counts every instance during a Ylok period that the corresponding 
  #IBI is shorter than the pre-stim median IBI:for start phases
  data2$c3=ifelse(data2$code2=="Ylok" & data2$prestim<data2$IBI, 1,0) 
  #rolling sum variable based on previous variable to get the 5 consecutive IBIs
  #shorter than pre-stim IBI:for start phases
  data2$c4=rollapply(data2$c3,width=5,sum,fill=NA,align="right")
  
  #variable marking attention phase based on Ylok and an IBI longer than
  #the median of 5 pre-stim IBIs and marking stop phases based on 5 consecutive IBIs
  #lower than the pre-stim (see c2), or a Nlok
  #data2$attn=ifelse(data2$code2=="Ylok" & data2$IBI>data2$prestim, "attn",
  #                  ifelse(data2$c2>=5 | data2$code=="Nlok","stop",NA)) 
  
  data2$attn=ifelse(data2$code2=="Ylok" & data2$c4>=5, "attn",
                    ifelse(data2$c2>=5 | data2$code=="Nlok","stop",NA)) 
  
  data2[1,"attn"] = "stop"
  
  prestim_ibi_avg=mean(data2$prestim, na.rm=T)
  
  
  ###Based on total median Nlok highest & lowest IBI
  data2[1,"code2"] = "Nlok"
  
  data2$med_low=ifelse(data2$code2=="Ylok" & median(subset(data2, code2=="Nlok")$IBI)
                       >data2$IBI, 1,0) 
  #rolling sum variable based on previous variable to get the 5 consecutive IBIs
  #shorter than pre-stim IBI:for stop phases
  data2$med_cnt=rollapply(data2$med_low,width=5,sum,fill=NA,align="right")
  
  #variable that counts every instance during a Ylok period that the corresponding 
  #IBI is shorter than median  IBI during "Ylok
  data2$med1=ifelse(data2$code2=="Ylok" & median(subset(data2, code2=="Nlok")$IBI)
                   < data2$IBI, 1,0) 
  #rolling sum variable based on previous variable to get the 5 consecutive IBIs
  #shorter than pre-stim IBI:for start phases
  data2$med2=rollapply(data2$med1,width=5,sum,fill=NA,align="right")
  
  data2$attn_med=ifelse(data2$code2=="Ylok" & data2$med2>=5, "attn",
                    ifelse(data2$med_cnt>=5 | data2$code=="Nlok","stop",NA)) 
  
  ibi_ylk_med=median(subset(data2, code2=="Nlok")$IBI)
  
  
 #change value previous to first attn to stop
  data2[1,"attn_med"] = "stop"
  
  #create rolling attention code
  attn_out=data2 %>% mutate (attn = na.locf(attn,na.rm = FALSE),
                          attn_med = na.locf(attn_med,na.rm = FALSE),
                          code2 = na.locf(code2,na.rm = FALSE),
                          attn_lks_only = ifelse(code2=="Nlok",NA, attn)) %>%
    dplyr::rename(data2,time=time.x)%>% 
    select(IBI,time,code,code2,rollmed_preceed,prestim,attn,attn_med,attn_lks_only,time_real) 
  
  ecg_points = nrow(data2)
  ecg_secs=last(data2$time)
  prestim_ibi_avg=mean(data2$prestim, na.rm=T)
  
  ####### EEG ###########################################
  setwd(psd.folder)
  
  #Read matlab file
  eeg= readMat(list.files(pattern = i ))

  #Define Frequencies
  f=eeg$f[[1]][[1]]
  theta_freqs=which(f> low_theta & f<high_theta)
  alpha_freqs=which(f> low_alpha & f<high_alpha)
  
  #Load Data
  eeg_wfp=eeg$eeg.wfp[[1]][[1]]
  eeg_wfp=eeg_wfp[,,]
  spectrum=eeg_wfp[2,,15]
  ##DEFINE MATRICES
  eeg_wfp_reorder = aperm(eeg_wfp, c(2,1,3))
  #Theta and alpha power at each 1 sec segment at each electrode
  theta_matrix=log10(colSums(eeg_wfp_reorder[theta_freqs,c(2,3,5,6,9,10,11,12,13,14,57,59,60),],dims=1))
  theta_nlog=colSums(eeg_wfp_reorder[theta_freqs,c(2,3,5,6,9,10,11,12,13,14,57,59,60),],dims=1)
  alpha_matrix=log10(colSums(eeg_wfp_reorder[alpha_freqs,c(2,3,5,6,9,10,11,12,13,14,57,59,60),],dims=1))
  alpha_nlog=colSums(eeg_wfp_reorder[alpha_freqs,c(2,3,5,6,9,10,11,12,13,14,57,59,60),],dims=1)
  
  
  #Exclude Outliers 
  outlier_mat_theta=ifelse(theta_nlog>median(theta_nlog,na.rm = T)+4*sd(theta_nlog,na.rm = T) |
                                                            theta_nlog<median(theta_nlog,na.rm = T)-4*sd(theta_nlog,na.rm = T),
                                            1,0)
theta_matrix=ifelse(outlier_mat_theta==1, NA,theta_matrix)
p_seg_reject_theta=prop.table(table(outlier_mat_theta))[2]
n_seg_reject_theta=table(outlier_mat_theta)[2]  

outlier_mat_alpha=ifelse(alpha_nlog>median(alpha_nlog,na.rm = T)+4*sd(alpha_nlog,na.rm = T) |
                     alpha_nlog<median(alpha_nlog,na.rm = T)-4*sd(alpha_nlog,na.rm = T),
                   1,0)
p_seg_reject_alpha=prop.table(table(outlier_mat_alpha))[2]
n_seg_reject_alpha=table(outlier_mat_alpha)[2]  
alpha_matrix=ifelse(outlier_mat_theta==1, NA,alpha_matrix)

  
  #3-50Hz Total Power in Each 1sec Segment At Each Electrode 
  spectrum_matrix=colSums(eeg_wfp_reorder[4:52,c(2,3,5,6,9,10,11,12,13,14,57,59,60),],dims=1)
  
  theta_matrix_rel=theta_nlog/spectrum_matrix
  alpha_matrix_rel=alpha_nlog/spectrum_matrix
  
  #Reshape and label electrode value
  theta_matrix_t=as.data.frame(t(theta_matrix))
  names(theta_matrix_t) = paste0(FE_total)
  alpha_matrix_t=as.data.frame(t(alpha_matrix))
  names(alpha_matrix_t) = paste0(FE_total)
  

  theta_matrix_rel_t=as.data.frame(t(theta_matrix_rel))
  names(theta_matrix_rel_t) = paste0(FE_total)
  theta_matrix_rel_t$frontal_abs=rowMeans(theta_matrix_rel_t[c(FE_total)])
  names(theta_matrix_rel_t) = paste0("theta_rel_",names(theta_matrix_rel_t))
  
  
  alpha_matrix_rel_t=as.data.frame(t(alpha_matrix_rel))
  names(alpha_matrix_rel_t) = paste0(FE_total)
  alpha_matrix_rel_t$frontal_abs=rowMeans(alpha_matrix_rel_t[c(FE_total)])
  names(alpha_matrix_rel_t) = paste0("alpha_rel_",names(alpha_matrix_rel_t))
  
  theta_matrix_t$frontal_abs_R=rowMeans(theta_matrix_t[c(FE_R)])
  theta_matrix_t$frontal_abs_L=rowMeans(theta_matrix_t[c(FE_L)])
  theta_matrix_t$frontal_abs=rowMeans(theta_matrix_t[c(FE_total)])
  names(theta_matrix_t) = paste0("theta_",names(theta_matrix_t))
  
  alpha_matrix_t$frontal_abs_R=rowMeans(alpha_matrix_t[c(FE_R)])
  alpha_matrix_t$frontal_abs_L=rowMeans(alpha_matrix_t[c(FE_L)])
  alpha_matrix_t$frontal_abs=rowMeans(alpha_matrix_t[c(FE_total)])
  names(alpha_matrix_t) = paste0("alpha_",names(alpha_matrix_t))
  
  
  eeg_electrode=cbind(theta_matrix_t,alpha_matrix_t,theta_matrix_rel_t,alpha_matrix_rel_t)
  
  eeg_electrode$time=1:nrow(eeg_electrode)
  
  eeg_codes=attn_out
  
  eeg_electrode=join(eeg_codes,eeg_electrode,by="time")
  
  
  # #FOR TOPOPLOTS
   theta_matrix_all=log10(colSums(eeg_wfp_reorder[theta_freqs,c(18,20,24,28,30,34,35,39,42,44,50,52,58),],dims=1))
   alpha_matrix_all=log10(colSums(eeg_wfp_reorder[alpha_freqs,c(18,20,24,28,30,34,35,39,42,44,50,52,58),],dims=1))

   theta_matrix_all_t=as.data.frame(t(theta_matrix_all))
   names(theta_matrix_all_t)=paste0(Other)
   alpha_matrix_all_t=as.data.frame(t(alpha_matrix_all))
   names(alpha_matrix_all_t)=paste0(Other)


   eeg_electrode=cbind(theta_matrix_t,alpha_matrix_t,theta_matrix_rel_t,alpha_matrix_rel_t,theta_matrix_all_t,alpha_matrix_all_t)

  eeg_electrode$time=1:nrow(eeg_electrode)

  eeg_codes=attn_out

  eeg_electrode=join(eeg_codes,eeg_electrode,by="time")


  #Topos
  eeg_avgs_alpha_attn_topo= eeg_electrode %>%
    group_by(attn) %>%
    summarize(FP2 = mean(theta_5, na.rm = TRUE),
              Fz = mean(theta_6, na.rm = TRUE),
              FP1 = mean(theta_10, na.rm = TRUE),
              F3 = mean(theta_12, na.rm =TRUE),
              F7 = mean(`18`, na.rm =TRUE),
              C3 = mean(`20`, na.rm=TRUE),
              T5 = mean(`24`, na.rm = TRUE),
              P3 = mean(`28`, na.rm =TRUE),
              T3 = mean(`30`, na.rm =TRUE),
              PZ = mean(`34`, na.rm=TRUE),
              O1 = mean(`35`, na.rm = TRUE),
              O2 = mean(`39`, na.rm =TRUE),
              P4 = mean(`42`, na.rm =TRUE),
              T4 = mean(`44`, na.rm=TRUE),
              C4 = mean(`50`, na.rm = TRUE),
              T6 = mean(`52`, na.rm =TRUE),
              F8 = mean(`58`, na.rm =TRUE),
              F4= mean(theta_60, na.rm=TRUE)) %>%
    pivot_wider(names_from = attn, 
                values_from = c(FP2:F4) ) %>%
    rename_all(function(x) paste0(x,"_alpha"))
  
  eeg_avgs_theta_attn_topo= eeg_electrode %>%
    group_by(attn) %>%
    summarize(FP2 = mean(theta_5, na.rm = TRUE),
              Fz = mean(theta_6, na.rm = TRUE),
              FP1 = mean(theta_10, na.rm = TRUE),
              F3 = mean(theta_12, na.rm =TRUE),
              F7 = mean(`18`, na.rm =TRUE),
              C3=mean(`20`, na.rm=TRUE),
              T5 = mean(`24`, na.rm = TRUE),
              P3 = mean(`28`, na.rm =TRUE),
              T3 = mean(`30`, na.rm =TRUE),
              PZ=mean(`34`, na.rm=TRUE),
              O1 = mean(`35`, na.rm = TRUE),
              O2 = mean(`39`, na.rm =TRUE),
              P4 = mean(`42`, na.rm =TRUE),
              T4=mean(`44`, na.rm=TRUE),
              C4 = mean(`50`, na.rm = TRUE),
              T6 = mean(`52`, na.rm =TRUE),
              F8 = mean(`58`, na.rm =TRUE),
              F4=mean(theta_60, na.rm=TRUE)) %>%
    pivot_wider(names_from = attn, 
                values_from = c(FP2:F4) ) %>%
    rename_all(function(x) paste0(x,"_theta"))

  eeg_avgs_theta_JR_topo= eeg_electrode %>%
    group_by(attn_lks_only) %>%
    summarize(FP2 = mean(theta_5, na.rm = TRUE),
              Fz = mean(theta_6, na.rm = TRUE),
              FP1 = mean(theta_10, na.rm = TRUE),
              F3 = mean(theta_12, na.rm =TRUE),
              F7 = mean(`18`, na.rm =TRUE),
              C3=mean(`20`, na.rm=TRUE),
              T5 = mean(`24`, na.rm = TRUE),
              P3 = mean(`28`, na.rm =TRUE),
              T3 = mean(`30`, na.rm =TRUE),
              PZ=mean(`34`, na.rm=TRUE),
              O1 = mean(`35`, na.rm = TRUE),
              O2 = mean(`39`, na.rm =TRUE),
              P4 = mean(`42`, na.rm =TRUE),
              T4=mean(`44`, na.rm=TRUE),
              C4 = mean(`50`, na.rm = TRUE),
              T6 = mean(`52`, na.rm =TRUE),
              F8 = mean(`58`, na.rm =TRUE),
              F4=mean(theta_60, na.rm=TRUE)) %>%
    pivot_wider(names_from = attn_lks_only, 
                values_from = c(FP2:F4) ) %>%
    rename_all(function(x) paste0(x,"_theta"))

  eeg_avgs_alpha_JR_topo= eeg_electrode %>%
    group_by(attn_lks_only) %>%
    summarize(FP2 = mean(alpha_5, na.rm = TRUE),
              Fz = mean(alpha_6, na.rm = TRUE),
              FP1 = mean(alpha_10, na.rm = TRUE),
              F3 = mean(alpha_12, na.rm =TRUE),
              F7 = mean(`18`, na.rm =TRUE),
              C3=mean(`20`, na.rm=TRUE),
              T5 = mean(`24`, na.rm = TRUE),
              P3 = mean(`28`, na.rm =TRUE),
              T3 = mean(`30`, na.rm =TRUE),
              PZ=mean(`34`, na.rm=TRUE),
              O1 = mean(`35`, na.rm = TRUE),
              O2 = mean(`39`, na.rm =TRUE),
              P4 = mean(`42`, na.rm =TRUE),
              T4=mean(`44`, na.rm=TRUE),
              C4 = mean(`50`, na.rm = TRUE),
              T6 = mean(`52`, na.rm =TRUE),
              F8 = mean(`58`, na.rm =TRUE),
              F4=mean(alpha_60, na.rm=TRUE)) %>%
    pivot_wider(names_from = attn_lks_only, 
              values_from = c(FP2:F4) ) %>%
    rename_all(function(x) paste0(x,"_alpha"))
  

 
  
  ######### MERGE EEG & IBI FOR TIMESERIES #####################
  # Create Avgs Based on Phases
  
  eeg_avgs_theta= eeg_electrode %>% 
    #dplyr::filter(dummy_outlier_theta==0) %>%
    group_by(attn) %>%
    summarize(ibi_attn_mean = mean(IBI, na.rm = TRUE),
              ibi_attn_sd = sd(IBI, na.rm = TRUE),
              abs_power_theta = mean(theta_frontal_abs, na.rm = TRUE),
              abs_power_theta_R = mean(theta_frontal_abs_R, na.rm =TRUE),
              abs_power_theta_L = mean(theta_frontal_abs_L, na.rm =TRUE),
              rel_power_theta=mean(theta_rel_frontal_abs, na.rm=TRUE)) %>%
    pivot_wider(names_from = attn, 
                values_from = c("abs_power_theta", "abs_power_theta_R","abs_power_theta_L", "ibi_attn_mean",
                                "ibi_attn_sd","rel_power_theta")) 
 
  eeg_avgs_alpha= eeg_electrode %>% 
    group_by(attn) %>%
    summarize(
              abs_power_alpha = mean(alpha_frontal_abs, na.rm = TRUE),
              abs_power_alpha_R = mean(alpha_frontal_abs_R, na.rm =TRUE),
              abs_power_alpha_L = mean(alpha_frontal_abs_L, na.rm =TRUE),
              rel_power_alpha=mean(alpha_rel_frontal_abs,na.rm=TRUE)) %>%
    pivot_wider(names_from = attn, 
                values_from = c("abs_power_alpha", "abs_power_alpha_R","abs_power_alpha_L","rel_power_alpha")) 
  
  # Create Avgs Based on Phases
  
  eeg_avgs_theta_med= eeg_electrode %>% 
    group_by(attn_med) %>%
    summarize(ibi_attn_mean_med = mean(IBI, na.rm = TRUE),
              ibi_attn_sd_med = sd(IBI, na.rm = TRUE),
              abs_power_theta_med = mean(theta_frontal_abs, na.rm = TRUE),
              abs_power_theta_R_med = mean(theta_frontal_abs_R, na.rm =TRUE),
              abs_power_theta_L_med = mean(theta_frontal_abs_L, na.rm =TRUE),
              rel_power_theta_med=mean(theta_rel_frontal_abs, na.rm=TRUE)) %>%
    pivot_wider(names_from = attn_med, 
                values_from = c("abs_power_theta_med", "abs_power_theta_R_med","abs_power_theta_L_med", "ibi_attn_mean_med",
                                "ibi_attn_sd_med","rel_power_theta_med")) 
  
  eeg_avgs_alpha_med= eeg_electrode %>% 
    group_by(attn_med) %>%
    summarize(
      abs_power_alpha_med = mean(alpha_frontal_abs, na.rm = TRUE),
      abs_power_alpha_R_med = mean(alpha_frontal_abs_R, na.rm =TRUE),
      abs_power_alpha_L_med = mean(alpha_frontal_abs_L, na.rm =TRUE),
      rel_power_alpha_med=mean(alpha_rel_frontal_abs,na.rm=TRUE))%>%
    pivot_wider(names_from = attn_med, 
                values_from = c("abs_power_alpha_med", "abs_power_alpha_R_med","abs_power_alpha_L_med","rel_power_alpha_med")) 
  
  
  #Create Averages Based on Looks
  eeg_looks_theta= eeg_electrode %>% 
    group_by(code2) %>%
    summarize(ibi_attn_looks = mean(IBI, na.rm = TRUE),
              abs_power_theta_looks = mean(theta_frontal_abs, na.rm = TRUE),
              abs_power_theta_looks_R = mean(theta_frontal_abs_R, na.rm =TRUE),
              abs_power_theta_looks_L = mean(theta_frontal_abs_L, na.rm =TRUE),
              rel_power_theta_looks=mean(theta_rel_frontal_abs,na.rm=TRUE)) %>%
    pivot_wider(names_from = code2, 
                values_from = c("abs_power_theta_looks", "abs_power_theta_looks_R","abs_power_theta_looks_L", "ibi_attn_looks",
                                "rel_power_theta_looks")) 
  
  #Create Averages Based on Looks
  eeg_looks_alpha= eeg_electrode %>% 
    group_by(code2) %>%
    summarize(
              abs_power_alpha_looks = mean(alpha_frontal_abs, na.rm = TRUE),
              abs_power_alpha_looks_R = mean(alpha_frontal_abs_R, na.rm =TRUE),
              abs_power_alpha_looks_L = mean(alpha_frontal_abs_L, na.rm =TRUE),
              rel_power_alpha_looks=mean(alpha_rel_frontal_abs,na.rm=TRUE)) %>%
    pivot_wider(names_from = code2, 
                values_from = c("abs_power_alpha_looks", "abs_power_alpha_looks_R","abs_power_alpha_looks_L",
                                "rel_power_alpha_looks"))
  
  eeg_attn_lks_only_theta = eeg_electrode %>% 
    group_by(attn_lks_only) %>%
    summarize(ibi_attn_lks_only_mean = mean(IBI, na.rm = TRUE),
              ibi_attn_lks_only_sd = sd(IBI, na.rm = TRUE),
              abs_power_theta_lks_only = mean(theta_frontal_abs, na.rm = TRUE),
              abs_power_theta_lks_only_R = mean(theta_frontal_abs_R, na.rm =TRUE),
              abs_power_theta_lks_only_L = mean(theta_frontal_abs_L, na.rm =TRUE),
              rel_power_theta_lks_only = mean(theta_rel_frontal_abs, na.rm=TRUE)) %>%
    pivot_wider(names_from = attn_lks_only, 
                values_from = c("abs_power_theta_lks_only", "abs_power_theta_lks_only_R","abs_power_theta_lks_only_L", 
                                "ibi_attn_lks_only_mean",
                                "ibi_attn_lks_only_sd","rel_power_theta_lks_only")) 
  
  eeg_attn_lks_only_alpha= eeg_electrode %>% 
    group_by(attn_lks_only) %>%
    summarize(
      abs_power_alpha_lks_only = mean(alpha_frontal_abs, na.rm = TRUE),
      abs_power_alpha_lks_only_R = mean(alpha_frontal_abs_R, na.rm =TRUE),
      abs_power_alpha_lks_only_L = mean(alpha_frontal_abs_L, na.rm =TRUE),
      rel_power_alpha_lks_only=mean(alpha_rel_frontal_abs,na.rm=TRUE)) %>%
    pivot_wider(names_from = attn_lks_only, 
                values_from = c("abs_power_alpha_lks_only", "abs_power_alpha_lks_only_R","abs_power_alpha_lks_only_L",
                                "rel_power_alpha_lks_only")) 
  
  
  eeg_electrode=eeg_electrode%>%mutate(attn_lks_only_Nlks=ifelse((is.na(attn_lks_only) & code2=="Nlok") ,"Nlok",attn_lks_only))
  
  
   data3=eeg_electrode %>% mutate(
    D = lag(attn),
    Dlooks= lag(code2),
    Dmed=lag(attn_med),
    
    Dlks_only=lag(attn_lks_only_Nlks),
    
    switch_dummy=ifelse(attn==D,0,1),
    switch_lk_dummy=ifelse(code2==Dlooks,0,1),
    switch_attnmed_dummy=ifelse(attn_med==Dmed,0,1),
    
    switch_attn_lks_only_dummy=ifelse(attn_lks_only_Nlks==Dlks_only,0,1), 
    
    attn_switch=ifelse(attn=="attn" & D=="stop","stop",
                       ifelse(attn =="stop" & D=="attn","attn",NA)),
    ylk_switch=ifelse(code2=="Ylok" & Dlooks=="Nlok","Nlok",
                      ifelse(code2 =="Nlok" & Dlooks=="Ylok","Ylok",NA)),
    attn_med_switch=ifelse(attn_med=="attn" & Dmed=="stop","stop",
                       ifelse(attn_med =="stop" & Dmed=="attn","attn",NA)),
    
    attn_lks_only_switch=ifelse(attn_lks_only_Nlks=="attn" & Dlks_only!="attn",attn_lks_only_Nlks,
                           ifelse(attn_lks_only_Nlks =="stop" & Dlks_only!="stop",attn_lks_only_Nlks,NA)),
    
    switch_attn_time=ifelse(attn_switch=="attn"|attn_switch=="stop",time,NA),
    switch_lk_time=ifelse(ylk_switch=="Ylok"|ylk_switch=="Nlok",time,NA),
    switch_attnmed_time=ifelse(attn_med_switch=="attn"|attn_med_switch=="stop",time,NA),
    
    switch_attn_lksonly_time=ifelse(attn_lks_only_switch=="attn"|attn_lks_only_switch=="stop",time,NA))

  
#force start and end of data
data3[nrow(data3), "attn_switch"] = data3[nrow(data3), "attn"]
data3[nrow(data3), "switch_attn_time"] = tail(data3$time,1)
data3[1, "switch_attn_time"] = head(data3$time,1)

data3[nrow(data3), "attn_med_switch"] = data3[nrow(data3), "attn"]
data3[nrow(data3), "switch_attnmed_time"] = tail(data3$time,1)
data3[1, "switch_attnmed_time"] = head(data3$time,1)
data3[nrow(data3), "attn_med_switch"] = data3[nrow(data3), "attn"]

data3[nrow(data3), "ylk_switch"] = data3[nrow(data3), "code2"]
data3[nrow(data3), "switch_lk_time"] = tail(data3$time,1)
data3[1, "switch_lk_time"] = head(data3$time,1)

data3[nrow(data3), "attn_lks_only_switch"] = data3[nrow(data3), "attn_lks_only_Nlks"]
data3[1, "switch_attn_lksonly_time"] = head(data3$time,1)
data3[nrow(data3), "switch_attn_lksonly_time"] = tail(data3$time,1)


  #create bins variable which calculates time of start when a switch from a nlok to a look occurs,
  #groups by looking and non looking to calculate summary stats
  descrip_attn=data3 %>% select(switch_attn_time,attn_switch) %>%
  drop_na(attn_switch)  %>% mutate(bins = switch_attn_time - lag(switch_attn_time, 1,default = 0)) %>%
    drop_na() %>%
    group_by(attn_switch,.drop=T) %>%
    summarize(count_phases=n(),
              avg_sec = mean(bins,na.rm=T),
              min = min(bins,na.rm=T),
              max = max (bins,na.rm=T),
              median= median(bins,na.rm=T),
              sum=sum(bins,na.rm=T),
              stdev=sd(bins,na.rm=T)) %>%
    pivot_wider(names_from = "attn_switch", values_from=c(count_phases:stdev))
  
  descrip_attn_med=data3 %>% select(switch_attnmed_time,attn_med_switch) %>%
    drop_na(attn_med_switch)  %>% mutate(bins = switch_attnmed_time - lag(switch_attnmed_time, 1,default = 0)) %>%
    drop_na() %>%
    group_by(attn_med_switch,.drop=T) %>%
    summarize(count_phases_med=n(),
              avg_sec_med = mean(switch_attnmed_time,na.rm=T),
              min_med = min(switch_attnmed_time,na.rm=T),
              max_med = max (switch_attnmed_time,na.rm=T),
              median_med= median(switch_attnmed_time,na.rm=T),
              sum_med=sum(switch_attnmed_time,na.rm=T),
              stdev_med=sd(switch_attnmed_time,na.rm=T)) %>%
    pivot_wider(names_from = "attn_med_switch", values_from=c(count_phases_med:stdev_med))
  
  
  #create bins variable which calculates time of start when a switch from a nlok to a look occurs,
  #groups by looking and non looking to calculate summary stats
  descrip_looks=data3 %>% 
    select(ylk_switch,switch_lk_time) %>%
    drop_na(ylk_switch) %>% mutate(bins = switch_lk_time - lag(switch_lk_time, 1,default = 0)) %>%
    drop_na() %>%
    group_by(ylk_switch,.drop=T) %>%
    summarize(count=n(),
              avg = mean(bins,na.rm=T),
              min = min(bins,na.rm=T),
              max = max (bins,na.rm=T),
              median= median(bins,na.rm=T),
              sum=sum(bins,na.rm=T),
              stdev=sd(bins,na.rm=T)) %>%
    pivot_wider(names_from = "ylk_switch", values_from=c(count:stdev))
  
  
  #
  
  descrip_lks_only=data3 %>% 
    select(attn_lks_only_switch,switch_attn_lksonly_time) %>%
    drop_na(attn_lks_only_switch) %>% mutate(bins = switch_attn_lksonly_time - lag(switch_attn_lksonly_time, 1,default = 0)) %>%
    drop_na() %>%
    group_by(attn_lks_only_switch,.drop=T) %>%
    summarize(count_phases_JR=n(),
              avg_JR = mean(bins,na.rm=T),
              min_JR = min(bins,na.rm=T),
              max_JR = max (bins,na.rm=T),
              median_JR= median(bins,na.rm=T),
              sum_JR=sum(bins,na.rm=T),
              stdev_JR=sd(bins,na.rm=T)) %>%
    pivot_wider(names_from = "attn_lks_only_switch", values_from=c(count_phases_JR:stdev_JR))
  
  
  #Timeseries 
  data3 = data3 %>% 
    mutate(code_num=dplyr::recode(code2, "Nlok"=0,"Ylok"=1))
  
  
  sub_file=data3 %>% select (IBI, time,time_real,code2,theta_frontal_abs,alpha_frontal_abs)
  
  #Massive timeseries file
  long_data=cbind(data.frame(ID=i),data3)
  time.series = rbind.fill(time.series,long_data)
  
  sub_file=data3 %>% select (IBI, time,time_real,code2,theta_frontal_abs,alpha_frontal_abs)
  
  ## output files
  
  outlier_theta_count=sum(eeg_electrode$dummy_outlier_theta,na.rm = T)
  outlier_alpha_count=sum(eeg_electrode$dummy_outlier_alpha, na.rm = T)
  
  eeg_secs=NROW(eeg_electrode)
  
  
  #REPORT ICA
  setwd(ica.folder)
  ica=readMat(list.files(pattern = i ))$file.proc.info
  IC_REJ=c((ica[[15]][3])[[1]][[1]][[1]])
  IC_Mean_Art_Prob_Kept=c((ica[[15]][5])[[1]][[1]][[1]])
  IC_Med_Art_Prob_Kept=c((ica[[15]][6])[[1]][[1]][[1]])
 
  eeg.mean=cbind(data.frame(ID=i),eeg_avgs_theta,eeg_avgs_alpha, eeg_avgs_theta_med,eeg_avgs_alpha_med,
                 eeg_looks_theta,eeg_looks_alpha,eeg_attn_lks_only_theta,eeg_attn_lks_only_alpha,
                 descrip_attn,descrip_attn_med,descrip_looks,descrip_lks_only,
                 outlier_theta_count,outlier_alpha_count,
                 #attn_count, stop_count, attn_med_count, stop_med_count,
                 #ylook_count,nlook_count, 
                 ecg_points,ecg_secs, prestim_ibi_avg, eeg_secs,ibi_ylk_med,p_seg_reject_theta,n_seg_reject_theta,
                 p_seg_reject_alpha,n_seg_reject_alpha,
                 IC_REJ,IC_Mean_Art_Prob_Kept,IC_Med_Art_Prob_Kept)


  eeg.file=rbind.fill(eeg.file, eeg.mean)
  
  elecs_attn=eeg_avgs_theta_JR_topo %>% select (matches("attn"))
  elects_stop=eeg_avgs_theta_JR_topo %>% select (matches("stop"))
  
  attn_m_stop=elecs_attn-elects_stop
  elecs_attn_m_stop=(cbind(data.frame(ID=i),attn_m_stop))
  
  topo_attn_sub_theta=rbind.fill(topo_attn_sub_theta,elecs_attn_m_stop)
  
  all.elects=cbind(data.frame(ID=i),eeg_avgs_alpha_JR_topo,eeg_avgs_theta_JR_topo)
  all.elects=all.elects %>% select(matches("ID|attn|stop"))
  
  topo.file=rbind.fill(topo.file,all.elects)  
    
  
  #Output 

  setwd(paste0(output.folder,"/JR_fulldata_new_RR"))
  write.csv(topo_attn_sub_theta,"topo_attn_sub_theta.csv")
  write.csv( topo.file, "topo.file.csv")
}

}
