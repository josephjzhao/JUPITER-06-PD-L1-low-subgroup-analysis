JUPITER-06 PD-L1 TPS \<1% subgroup analysis
================
Joseph J Zhao, MBBS<sup>1,2</sup>; Filippo Pietrantonio, MD<sup>3</sup>;
Raghav Sundar, MBBS, PhD<sup>1,2</sup>

<sup>1</sup> Yong Loo Lin School of Medicine, National University of
Singapore, Singapore

<sup>2</sup> Department of Haematology-Oncology, National University
Cancer Institute, Singapore

<sup>3</sup> Medical Oncology Department, Fondazione IRCCS Istituto
Nazionale dei Tumori, Milan, Italy

### Letter in Reply to Gao *et al* - Supporting documents

Regarding concerns raised pertaining to our original article: Yap, D. W.
T. *et al*. Effectiveness of Immune Checkpoint Inhibitors in Patients
With Advanced Esophageal Squamous Cell Carcinoma: A Meta-analysis
Including Low PD-L1 Subgroups. JAMA Oncology,
<doi:10.1001/jamaoncol.2022.5816> (2022).

This Github repository contains supporting workings, data, and figures
for our Letter in Reply.

### Load packages

``` r
package.name=c(# data reconstruction
  "IPDfromKM", 
  # data manipulation
  "tidyverse", "readr", "readxl", "dplyr", "tidyr", "lubridate", "tibble", "plyr", "devtools", "stringr", "stringi",
  # survival analysis
  "survminer", "survival", "rstpm2", "survRM2", "prodlim",
  # plotting
  "ggplot2"
)

for (package.name in package.name){
  tryCatch(
  {
  if (!require(package.name, character.only = TRUE)){ 
    install.packages(package.name, character.only = TRUE) 
    library(package.name, character.only = TRUE)} else {library(package.name, character.only = TRUE)}
  }, 
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
```

### Set working directory

``` r
wd="C:/Users/jzhao/OneDrive/Research_Cloud/NUH_NCIS/pdl1_escc_kms/"
setwd(wd)
```

### Load reconstructed time-to-event data of JUPITER-06 PD-L1 TPS \<1% subgroups reported by Wu *et al* (figures 1B &D)

Reconstruction of survival data was conducted separately with the
*IPDfromKM* package. Reconstructed data is provided in the
“reconstructed_data” folder.

Reference: Wu, H. X. *et al*. Clinical Benefit of First-Line Programmed
Death-1 Antibody Plus Chemotherapy in Low Programmed Cell Death Ligand
1-Expressing Esophageal Squamous Cell Carcinoma: A Post Hoc Analysis of
JUPITER-06 and Meta-Analysis. J Clin Oncol, Jco2201490,
<doi:10.1200/jco.22.01490> (2022).

``` r
files <- list.files(path = paste0(wd, "data/recon/"), pattern = paste(".csv",sep=""), full.names = T)
df <- sapply(files, read_csv, simplify=F) %>% 
          bind_rows(.id = "id") %>% 
          dplyr::select(time, status, arm, outcome)

head(df)
```

    ## # A tibble: 6 x 4
    ##    time status arm             outcome
    ##   <dbl>  <dbl> <chr>           <chr>  
    ## 1 0.340      1 Placebo plus TP OS     
    ## 2 1.04       0 Placebo plus TP OS     
    ## 3 1.62       1 Placebo plus TP OS     
    ## 4 1.66       1 Placebo plus TP OS     
    ## 5 1.86       0 Placebo plus TP OS     
    ## 6 2.82       1 Placebo plus TP OS

### Figure 1 Compare reconstructed KM plots to original plots

``` r
for (i.outcome in unique(df$outcome)){
  
  cox=coxph(Surv(time, status)~arm, data=subset(df, df$outcome==i.outcome))
  sum.cox=summary(cox)
  
  km=survfit(Surv(time, status)~arm, data=subset(df, df$outcome==i.outcome))
  
  plot_km = ggsurvplot(km,
                      data = subset(df, df$outcome==i.outcome), 
                      size=0.8,
                      risk.table = TRUE,    
                      censor.shape="|",
                      censor.size = 1.2,
                      palette = c("black","steelblue4"),
                      xlim = c(0,max(subset(df, df$outcome==i.outcome)$time)),
                      xlab = "Time, months",
                      break.time.by = 3,   
                      ggtheme = theme_bw(),
                      risk.table.y.text = F,
                      risk.table.fontsize=2.7,
                      legend.labs = (cox$xlevels$arm),
                      title=paste0("Reconstructed JUPITER-06 PD-L1 TPS<1%: ", i.outcome)
                      )
  
  print(plot_km)
}
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

### Calculate follow-up time with reverse KM method with the *prodlim* package

``` r
quantile(prodlim(Hist(time=time, event=status)~arm, data=subset(df, df$outcome=="OS"), reverse=TRUE))
```

    ## Quantiles of the potential follow up time distribution based on the Kaplan-Meier method 
    ## applied to the censored times reversing the roles of event status and censored.
    ## 
    ## Table of quantiles and corresponding confidence limits:
    ## 
    ##  arm=Placebo plus TP 
    ## 
    ##      q quantile lower upper
    ## 1 0.00       NA    NA    NA
    ## 2 0.25     12.8  10.4  18.0
    ## 3 0.50      8.1   7.1   9.7
    ## 4 0.75      5.3   4.6   6.6
    ## 5 1.00      1.0   1.0   3.2
    ## Median time (IQR):8.06 (5.28;12.80)
    ## 
    ## Table of quantiles and corresponding confidence limits:
    ## 
    ##  arm=Toripalimab plus TP 
    ## 
    ##      q quantile lower upper
    ## 1 0.00       NA    NA    NA
    ## 2 0.25     13.8  10.3  19.2
    ## 3 0.50      8.3   7.6   9.8
    ## 4 0.75      6.4   5.3   7.2
    ## 5 1.00      3.1   3.1   3.7
    ## Median time (IQR):8.27 (6.42;13.76)

### Figure 2 Evaluate for proportional hazards assumption with the Grambsch-Therneau test and plotted Schoenfeld residuals

``` r
for (i.outcome in unique(df$outcome)){print(ggcoxzph(cox.zph(coxph(Surv(time, status)~arm, data=subset(df, df$outcome==i.outcome))), title=i.outcome))}
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

### Figure 3 Hazard difference over time

``` r
for (i.outcome in unique(df$outcome)){
  
  df_temp=subset(df, df$outcome==i.outcome)
  df_temp$arm=ifelse(df_temp$arm=="Placebo plus TP", 0, 1)
  
  # fit model for stpm2 
  fit=stpm2(Surv(time, status) ~ arm, data=df_temp, df=3)
  pred=predict(fit, newdata=data.frame(arm=0), type="hdiff", var="arm", exposed=function(data) transform(data, arm=1), grid=TRUE, full=TRUE, se.fit=TRUE)
  
  plot_timevar=ggplot(pred, aes(x=time, y=Estimate, ymin=lower, ymax=upper)) +
            geom_ribbon(alpha=0.4, fill="steelblue4") + 
            geom_line(size=1, color="darkgoldenrod4") +
            theme(plot.title= element_text(face="bold", size=14), axis.text=element_text(color="black"))+
            theme_bw()+
            coord_cartesian(xlim=c(0, 0.9*max(df_temp$time)))+
            geom_hline(yintercept=0, linetype="dashed")+
            labs(color="black", x = "", y = "Hazard difference\nfavours Toripalimab plus TP | favours Placebo plus TP", 
                 title = paste0("JUPITER-06 PD-L1 TPS<1% subgroup\nHazard-difference over time from reconstructed ", i.outcome, " KM plots"), face="bold") +
            geom_vline(xintercept = pred$time[pred$Estimate==min(pred$Estimate)], linetype="dashed", size=0.8) +
            annotate("text", x=pred$time[pred$Estimate==min(pred$Estimate)]+0.5, y=min(pred$Estimate)-0.005, hjust=0, label=paste0("Maximal hazard difference: ", format(round(pred$time[pred$Estimate==min(pred$Estimate)], 3), nsmall=3), " months"))
  
  # historam for censorship counts
  df_temp_censor=subset(df, df$status==0 & df$outcome==i.outcome)
  
  plot_censor=ggplot(aes(x=time, fill=arm), data=df_temp_censor)+
    geom_histogram( position="identity", width=0.2, alpha=0.3)+
    theme_bw()+
    coord_cartesian(xlim=c(0, 0.9*max(df_temp$time)))+
    labs(y="censorship count", x="Time, months", fill="")+
    scale_fill_manual(values=c("black", "steelblue4")) +
    geom_vline(xintercept = median(df_temp_censor$time), linetype="dashed", size=0.8) +
            annotate("text", x=median(df_temp_censor$time)+0.3, y=4, hjust=0, label=paste0("Median: ", format(round(median(df_temp_censor$time), 3), nsmall=3), " months"))
  
  ggarrange(plot_timevar, plot_censor, ncol = 1, nrow = 2,
            heights = c(2, 0.5),
            common.legend = TRUE, legend = "bottom", align='v') %>% print
    
}
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

Abbreviations: PD-L1, programmed death ligand 1; OS, overall survival;
PFS, progression-free survival; RMST, restricted mean survival time;
TPS, Tumor Proportion Score; TP, paclitaxel plus cisplatin.
