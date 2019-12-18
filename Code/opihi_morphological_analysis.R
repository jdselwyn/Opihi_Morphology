#### Computer set up ####
rm(list=ls())

#### Libraries ####
library(tidyverse)
library(magrittr)
library(nlme)
library(mgcv)
library(lmtest)
library(cowplot)
library(purrrlyr)
library(brms)
library(tidybayes)
library(future)

plan(multiprocess)
#
#### Functions ####
SurfArea <- function(A,B,H) {
  #Function caluclates the Lateral surface area of an elliptical cone
  require(cubature)
  integrand<-function(t,a,b,h) {
    #https://rechneronline.de/pi/elliptic-cone.php
    sqrt(a^2*b^2 + h^2*(A^2*sin(t)^2 + b^2 * cos(t)^2))
  }
  
  pmap_dbl(.l=list(A,B,H),
           .f=function(x,y,z) adaptIntegrate(integrand, 
                                             lowerLimit=c(0),
                                             upperLimit=c(2*pi), 
                                             a=x,b=y,h=z)$integral/2)
}

triangulate <- function(IJData) {
  #JS - I haven't looked into how this works or how to improve.
  
  #define triangulateY3grY2
  triangulateY3grY2 <- function(IJData){
    IJData["X4"]<- (IJData["X2"]) #creating pt 4 using value X2 and Y3
    IJData["Y4"]<- (IJData["Y3"])
    IJData["L13"]<- sqrt(((IJData["X3"]-IJData["X1"])^2)+((IJData["Y3"]-IJData["Y1"])^2)) #length of L13 using pythag, posterior length, spire to aperture
    IJData["Slope_L13"]<- ((IJData["Y1"]-IJData["Y3"])/(IJData["X1"]-IJData["X3"])) #find slope using rise over run
    IJData["b_L13"]<- (IJData["Y3"]-(IJData["Slope_L13"]*IJData["X3"])) #find y-int to plug into point-slope formula
    IJData["X5"]<- IJData["X2"] #define x of point 5
    IJData["Y5"]<- (IJData["X5"]*IJData["Slope_L13"])+IJData["b_L13"] #calulate y5 using y=mx + b 
    IJData["L35"]<- sqrt(((IJData["X5"]-IJData["X3"])^2)+((IJData["Y5"]-IJData["Y3"])^2)) #length of L35 using pythag
    IJData["L23"]<- sqrt(((IJData["X2"]-IJData["X3"])^2)+((IJData["Y2"]-IJData["Y3"])^2)) #length of L23 using pythag
    IJData["L34"]<- abs(IJData["X4"]-IJData["X3"]) #length of L34 using diffs bw x of 3 and 4 b/c they have the same y values
    IJData["A435"]<- acos(IJData["L34"]/IJData["L35"]) #angle of point 3 on triangle 435 
    IJData["A432"]<- acos(IJData["L34"]/IJData["L23"]) #angle of point 3 on triangle 432 
    IJData["A231"]<- IJData["A435"]- IJData["A432"]#angle of point 3 on triangle 231
    IJData["L36"]<- IJData["L13"]*(cos(IJData["A231"])) #biologically relevant thing!!! (posterior apex position)
    IJData["L26"]<- IJData["L23"]-IJData["L36"] #anterior apex position
    
    #ceb Oct 1 2018 Additional measures
    IJData["L12"]<- sqrt(((IJData["X2"]-IJData["X1"])^2)+((IJData["Y2"]-IJData["Y1"])^2)) #anterior length, spire to aperture
    IJData["L16"]<- sqrt(IJData["L12"]^2 - IJData["L26"]^2) #shell height
    IJData["A213"]<- (180/pi) * (acos((IJData["L12"]^2 + IJData["L13"]^2 - IJData["L23"]^2) / (2*IJData["L12"]*IJData["L13"])))
    
    return(IJData[,c("ShellID","L13","L26","L36","L12","L16","A213")])
  } #end of triangulateY3grY2
  
  #define triangulateY3lsY2
  triangulateY3lsY2 <- function(IJData){
    names(IJData) <- gsub("2$", "2orig", names(IJData)) #rename columns so that point 2 is now called point 3 and vice versa
    names(IJData) <- gsub("3$", "2", names(IJData))
    names(IJData) <- gsub("2orig$", "3", names(IJData))
    IJData["X4"]<- (IJData["X2"]) #creating pt 4 using value X2 and Y3
    IJData["Y4"]<- (IJData["Y3"])
    IJData["L13"]<- sqrt(((IJData["X3"]-IJData["X1"])^2)+((IJData["Y3"]-IJData["Y1"])^2)) #length of L13 using pythag
    IJData["L12"]<- sqrt(((IJData["X2"]-IJData["X1"])^2)+((IJData["Y2"]-IJData["Y1"])^2)) #length of L13 using pythag
    IJData["Slope_L13"]<- ((IJData["Y3"]-IJData["Y1"])/(IJData["X3"]-IJData["X1"])) #find slope using rise over run #DIFFERENT
    IJData["b_L13"]<- (IJData["Y3"]-(IJData["Slope_L13"]*IJData["X3"])) #find y-int to plug into point-slope formula
    IJData["X5"]<- IJData["X2"] #define x of point 5
    IJData["Y5"]<- (IJData["X5"]*IJData["Slope_L13"])+IJData["b_L13"] #calulate y5 using y=mx + b 
    IJData["L35"]<- sqrt(((IJData["X5"]-IJData["X3"])^2)+((IJData["Y5"]-IJData["Y3"])^2)) #length of L35 using pythag
    IJData["L23"]<- sqrt(((IJData["X2"]-IJData["X3"])^2)+((IJData["Y2"]-IJData["Y3"])^2)) #length of L23 using pythag
    IJData["L34"]<- abs(IJData["X4"]-IJData["X3"]) #length of L34 using diffs bw x of 3 and 4 b/c they have the same y values
    IJData["A435"]<- acos(IJData["L34"]/IJData["L35"]) #angle of point 3 on triangle 435 
    IJData["A432"]<- acos(IJData["L34"]/IJData["L23"]) #angle of point 3 on triangle 432 
    IJData["A231"]<- IJData["A435"]- IJData["A432"]#angle of point 3 on triangle 231
    IJData["L36"]<- IJData["L13"]*(cos(IJData["A231"])) #biologically relevant thing!!! (apex position)
    IJData["L26"]<- IJData["L23"]-IJData["L36"] 
    
    #ceb Oct 1 2018 additional measures
    #IJData["L12"]<- sqrt(((IJData["X2"]-IJData["X1"])^2)+((IJData["Y2"]-IJData["Y1"])^2))
    IJData["L16"]<- sqrt(IJData["L12"]^2 - IJData["L26"]^2) #shell height
    IJData["A312"]<- (180/pi) * (acos((IJData["L12"]^2 + IJData["L13"]^2 - IJData["L23"]^2) / (2*IJData["L12"]*IJData["L13"])))
    
    names(IJData) <- gsub("2", "2orig", names(IJData))
    names(IJData) <- gsub("3", "2", names(IJData))
    names(IJData) <- gsub("2orig", "3", names(IJData))
    return(IJData[,c("ShellID","L13","L26","L36","L12","L16","A213")])
  } #end of triangulateY3lsY2
  
  # define triangulateY3eqY2  # PLACEHOLDER FOR NOW, NEED TO FIGURE OUT REAL FUNCTION FOR THIS CASE
  triangulateY3eqY2 <- function (IJData) {
    IJData["L13"]<- sqrt(((IJData["X3"]-IJData["X1"])^2)+((IJData["Y3"]-IJData["Y1"])^2)) #length of L13 using pythag
    IJData["Slope_L13"]<- ((IJData["Y1"]-IJData["Y3"])/(IJData["X1"]-IJData["X3"])) #find slope using rise over run
    IJData["b_L13"]<- (IJData["Y3"]-(IJData["Slope_L13"]*IJData["X3"])) #find y-int to plug into point-slope formula
    IJData["X5"]<- IJData["X2"] #define x of point 5
    IJData["Y5"]<- (IJData["X5"]*IJData["Slope_L13"])+IJData["b_L13"] #calulate y5 using y=mx + b 
    IJData["L35"]<- sqrt(((IJData["X5"]-IJData["X3"])^2)+((IJData["Y5"]-IJData["Y3"])^2)) #length of L35 using pythag
    IJData["L23"]<- sqrt(((IJData["X2"]-IJData["X3"])^2)+((IJData["Y2"]-IJData["Y3"])^2)) #length of L23 using pythag
    IJData["A235"]<- acos(IJData["L23"]/IJData["L35"]) #angle of point 3 on triangle 435 
    IJData["L36"]<- IJData["L13"]*(cos(IJData["A235"])) #biologically relevant thing!!! (apex position)
    IJData["L26"]<- IJData["L23"]-IJData["L36"] 
    
    #ceb Oct 1 2018 Calculate the apex angle, lengthwise
    IJData["L12"]<- sqrt(((IJData["X2"]-IJData["X1"])^2)+((IJData["Y2"]-IJData["Y1"])^2))
    IJData["L16"]<- sqrt(IJData["L12"]^2 - IJData["L26"]^2) #shell height
    IJData["A213"]<- (180/pi) * (acos((IJData["L12"]^2 + IJData["L13"]^2 - IJData["L23"]^2) / (2*IJData["L12"]*IJData["L13"])))
    
    return(IJData[,c("ShellID","L13","L26","L36","L12","L16","A213")])
  } #end of triangulateY3eqY2
  
  if(IJData["Y3"]>IJData["Y2"]){
    IJDataSub <- triangulateY3grY2(IJData)
  }else if(IJData["Y3"]<IJData["Y2"]){
    IJDataSub <- triangulateY3lsY2(IJData)
  }else {
    IJDataSub <- triangulateY3eqY2(IJData)
  } # end of if statements
  
  # IJData <- cbind(IJData[,1:7],IJDataSub[,2:7])
  # print("woohoo, it worked!")
  # return(IJData)
  IJDataSub[-1]
}

cld_display<-function(y){
  #Modified algorithm from here:
  #https://stackoverflow.com/questions/27293770/compact-letter-display-from-logical-matrix
  require(igraph)
  
  output<-y %$%
    group %>% 
    as_tibble(rownames='comparison') %>%
    separate(comparison,sep='-',into=c('A','B')) %>%
    dplyr::select(1,2) %>%
    gather %>% 
    dplyr::select(value) %>%
    distinct()
  
  cliq<-y %$%
    group %>%
    as_tibble(rownames='comparison') %>%
    separate(comparison,sep='-',into=c('A','B')) %>%
    mutate(sig=`p adj`<0.05) %>%
    dplyr::select(1,2,'sig') %>%
    filter(!sig) %>%
    dplyr::select(-sig) %>%
    
    graph.data.frame(directed = F,vertices = output$value) %>% 
    maximal.cliques
  
  
  lab.txt <- vector(mode="list", nrow(output)) # empty list
  lab <- letters[seq(cliq)] # clique labels
  for(i in seq(cliq)){ # loop to concatenate clique labels
    for(j in cliq[[i]]){
      lab.txt[[j]] <- str_c(lab.txt[[j]], lab[i])
    }
  }
  output %>%
    mutate(grouping=unlist(lab.txt))
}

ggdiag_plots<-function(model,independent,dependent){
  
  continuous_vars<-function(r,f){
    #r<-model_data$.resid;f<-model_data$.fitted
    
    continue_data<-tibble(r=r,f=f)
    
    if(length(unique(f))<=4){sub.text='Too few points for spline'}
    if(length(unique(f))>4){
      if (length(unique(f)) < 10) {k <- round(length(unique(f))/2)} else {k <- 10}
      
      rf.model<-gam(r~s(f, k=k), data=data.frame(r=as.double(r),f=as.double(f)), se.fit=T)
      spline_text <- paste("Linearity p-value =", signif(summary(rf.model)$s.table[1,4],4),'(spline)')
    }
    
    bp_text <- str_c('Homoscedasticity p-value =',signif(bptest(lm(r~f))$p.value,4),'(Breush-Pagan test)',sep=' ')
    
    continue_data %>%
      nest %>%
      mutate(plot=map(data,~.x %>%
                        ggplot(aes(x=f,y=r)) +
                        geom_point() +
                        geom_hline(yintercept = 0,colour='blue') +
                        geom_smooth(method='loess',colour='red',fill='red') +
                        xlab(str_c('Fitted',spline_text,bp_text,sep='\n')) +
                        ylab('Residuals'))) %>%
      dplyr::select(plot)
  }
  
  catagorical_vars<-function(r,f){
    cat_data<-tibble(r=r,f=f) %>%
      group_by(f) %>%
      mutate(robust=abs(r-median(r))) %>%
      ungroup %>%
      mutate(plot_resid=resid(aov(r~f))) %>%
      nest %>%
      mutate(BFL = map(data,~aov(robust~f,data=.x))) %>%
      mutate(plot=map2(data,BFL,~.x %>%
                         ggplot(aes(x=f,y=plot_resid)) +
                         geom_boxplot() +
                         xlab(str_c('Homoscedasticity p-value =',signif(anova(.y)[1,5],4),'(BFL test)',sep=' ')) +
                         ylab('Residuals'))) %>%
      dplyr::select(plot)
    cat_data
  }
  
  qqPlot<-function(r){
    n<-length(r)
    y.mx<-matrix(0,5000,n)
    
    for(i in 1:5000) {y.mx[i,]<-sort(rnorm(n,median(r),mad(r)))}
    
    tibble(x=qnorm(ppoints(n)),y=sort(r),
           meds=apply(y.mx,2,median),
           lwr99=apply(y.mx,2, function(x) {quantile(x,0.005)}),
           lwr95=apply(y.mx,2, function(x) {quantile(x,0.025)}),
           upr95=apply(y.mx,2, function(x) {quantile(x,0.975)}),
           upr99=apply(y.mx,2, function(x) {quantile(x,0.995)})) %>%
      nest %>%
      mutate(SWtest=map(data,~shapiro.test(.x$y))) %>%
      mutate(plot=map2(data,SWtest,~.x %>%
                         ggplot(aes(x=x,y=y)) +
                         geom_ribbon(aes(ymin=lwr99,ymax=upr99),alpha=0.2) +
                         geom_ribbon(aes(ymin=lwr95,ymax=upr95),alpha=0.4) +
                         geom_smooth(aes(y=meds),lty=2,method='lm') +
                         geom_point() +
                         xlab(str_c("Theoretical Quantiles",
                                    str_c('Normality p-value =',signif(.y$p.value,4),'(SW test)',sep=' '),
                                    str_c('Number of observations:',nrow(.x),sep=' '),sep='\n')) +
                         ylab("Sample Quantiles"))) %>%
      dplyr::select(plot)
  }
  
  cookPlot<-function(cook){
    cook %>%
      as.tibble %>%
      mutate(index=row_number()) %>%
      nest %>%
      mutate(plot=map(data,~.x %>%
                        ggplot(aes(x=index,ymax=value)) +
                        geom_linerange(aes(ymin=0)) +
                        xlab(str_c('Observation',str_c('Most influencial point',signif(max(.$value),4),"(Cook's distance)",sep=' '),sep='\n')) +
                        ylab("Cook's Distance"))) %>%
      dplyr::select(plot)
  }
  
  # Goal to include logic to pick out the independent and dependent variables
  
  if(any(class(model)=='gnls')) {
    model_data<-tibble(.fitted=fitted(model),.resid=resid(model, type="normalized"))
  } else {
    model_data<-augment(model) %>%
      as_tibble()
  }
  
  
  #Create continuous plots
  continuous_plots<-model_data %>% 
    dplyr::select(one_of(c(independent,'.fitted','.resid'))) %>%
    summarize_if(.predicate = function(x) {is.numeric(x) | is.integer(x)},
                 .funs=function(y) continuous_vars(f=y,r=.$.resid)) %>%
    dplyr::select(-matches('.resid')) %>%
    map(function(x) x$plot[[1]])
  
  
  #Create QQ plot
  qq_plot<-model_data %>% 
    dplyr::select(one_of(c('.resid'))) %>%
    summarize_all(qqPlot) %>%
    map(function(x) x$plot[[1]])
  
  #Create cook's distance plot - not yet working
  cook_plot<-NULL
  if(".cooksd" %in% names(model_data)) {
    cook_plot<-model_data %>%
      dplyr::select(one_of(c('.cooksd'))) %>%
      summarise_all(cookPlot) %>%
      map(function(x) x$plot[[1]])
  }
  
  #If factor create BFL plot
  factor_plots<-model_data %>%
    dplyr::select(one_of(c(independent,'.fitted',dependent))) %>%
    summarize_if(.predicate = function(x) {is.factor(x) | is.character(x)},
                 .funs=function(y) catagorical_vars(f=y,r=.[[dependent]])) %>%
    map(function(x) x$plot[[1]])
  
  output<-list(QQ=qq_plot,Continuous=continuous_plots,Catagorical=factor_plots,Cooks=cook_plot) %>%
    unlist(recursive = F)
  
  removal<-(output %>%
              map_lgl(is.null) %>%
              which(!.) %>% 
              names())
  
  if(length(removal)>0) {
    output<-unlist(output,recursive = F) %>% 
      list_modify(!!removal:=NULL)
  }
  
  output<-plot_grid(plotlist=output,labels=names(output))
}

cld_display_bayes<-function(x){
  #Modified algorithm from here:
  #https://stackoverflow.com/questions/27293770/compact-letter-display-from-logical-matrix
  suppressMessages(library(igraph))
  
  output <- x %>%
    dplyr::select(1,2) %>%
    gather %>%
    dplyr::select(value) %>%
    distinct 
  
  cliq<-x %>%
    ungroup %>%
    mutate(sig = !(.[,4] < 0 & .[,5] > 0)) %>%
    dplyr::select(1,2,sig) %>%
    filter(!sig) %>%
    dplyr::select(-sig) %>%
    
    graph.data.frame(directed = FALSE, vertices = output$value) %>%
    maximal.cliques()
  
  lab.txt <- vector(mode = 'list', nrow(output))
  lab <- letters[seq(cliq)]
  for(i in seq(cliq)){
    for(j in cliq[[i]]){
      lab.txt[[j]] <- str_c(lab.txt[[j]], lab[i])
    }
  }
  
  output %>%
    mutate(grouping = unlist(lab.txt)) %>%
    rename(!!word(colnames(x)[1], 1, sep = '\\.') := value)
  
}

#### Read in data ####
columns_for_distance<-c("Width_normalized","Height_normalized","ApertureArea","Posterior_Apex_Position_normalized","Curve_normalized","L13_normalized","Volume","SurfaceArea","Circumference")

island_names<-c("GP","LPP",'MMM','NIH','KAU','OAH','MAUI','BI')

caliper_measures <- read_csv(file = "../Data/9Nov18_WS.csv", na = c("N/A", "n/a"))

ImageJ_Data<-read_csv('../Data/ImageJData.csv')

#### Sampling Sites ####
Sampling_Sites <- tibble(Island=c("GP","LPP","MMM", "NIH", 
                                  "KAU", 'KAU', "OAH",'OAH', "MAUI",'BI'),
                         SamplingSite = c('GP','LPP',"MMM",'Nihoa',
                                          'Miloloii','PMRF','Kewalo','Magic_Island',"MAUI",'Hilo'),
                         
                         longitude=c(-167.99915277777777,-166.261, -164.70361944444443, -161.92190555555555, 
                                     -159.7243861111111, -159.78777222222223, -157.8630472222222,-157.8462777777778, -156.0308777777778, -155.0523777777778),
                         
                         latitude=c(24.998308333333334, 23.76872222222222, 23.577730555555554, 23.058647222222223, 
                                    22.147141666666666, 22.03068888888889, 21.291522222222223,21.281875, 20.674702777777778, 19.733288888888886),
                         Location=c(0,222.8,382.8,668.4,
                                    911.3, 911.3, 1028.3, 1028.3, 1083.6, 1130.4)) %>%
  mutate(island_group = case_when(Island == 'BI' | Island == 'MAUI' | Island == 'OAH' | Island == 'KAU' ~ 'MHI',
                                  TRUE ~ 'NWHI')) %>%
  mutate(latitude.c = latitude - mean(latitude),
         Location.c = if_else(island_group == 'MHI', Location - 911.3, Location))

#### Impute missing length/width ####
caliper_measures <- caliper_measures %>%
  #Convert to mm
  mutate(Length=Length/10,
         Width=Width/10,
         Height=HeightLengthwise/10) %>%
  
  #Impute missing lengths and widths from the other one
  nest(data = everything()) %>%
  mutate(length_imputation = map(data, ~ .x %>%
                                   dplyr::select(Length,Width) %>%
                                   na.omit() %>%
                                   gnls(Length ~ a*Width^b, 
                                        start = c(a=1,b=1),
                                        weights = varExp(form=~fitted(.)),
                                        data=.)),
         length_imputation_gnls = map(data, ~ .x %>%
                                        dplyr::select(Length,Width) %>%
                                        na.omit() %>%
                                        gnls(Length ~ a*Width^b, 
                                             start = c(a=1,b=1),
                                             data=.)),
         length_imputation_lm = map(data, ~ .x %>%
                                      dplyr::select(Length,Width) %>%
                                      na.omit() %>%
                                      lm(Length ~ Width, data=.)),
         length_imputation_null = map(data, ~ .x %>%
                                        dplyr::select(Length,Width) %>%
                                        na.omit() %>%
                                        gnls(Length ~ a, 
                                             start = c(a=1),
                                             data=.)),
         
         width_imputation = map(data, ~ .x %>%
                                  dplyr::select(Length,Width) %>%
                                  na.omit() %>%
                                  gnls(Width ~ a*Length^b, 
                                       start = c(a=1,b=1),
                                       weights = varExp(form=~fitted(.)),
                                       data=.,
                                       na.action = na.exclude)),
         width_imputation_gnls = map(data, ~ .x %>%
                                       dplyr::select(Length,Width) %>%
                                       na.omit() %>%
                                       gnls(Width ~ a*Length^b, 
                                            start = c(a=1,b=1),
                                            data=.)),
         width_imputation_lm = map(data, ~ .x %>%
                                     dplyr::select(Length,Width) %>%
                                     na.omit() %>%
                                     lm(Width ~ Length, data=.)),
         width_imputation_null = map(data, ~ .x %>%
                                       dplyr::select(Length,Width) %>%
                                       na.omit() %>%
                                       gnls(Width ~ a, 
                                            start = c(a=1),
                                            data=.))) %>%
  
  mutate(length_summaries = pmap(list(full_model = length_imputation, 
                                      basic_nls = length_imputation_gnls,
                                      basic_linear = length_imputation_lm,
                                      null = length_imputation_null),
                                 function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                         as_tibble(rownames='Model') %>%
                                                                                         mutate(RSE = c(null$sigma,summary(basic_linear)$sigma,basic_nls$sigma,full_model$sigma)),
                                                                                       
                                                                                       anova(null,full_model))),
         
         width_summaries = pmap(list(full_model = width_imputation, 
                                     basic_nls = width_imputation_gnls,
                                     basic_linear = width_imputation_lm,
                                     null = width_imputation_null),
                                function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                        as_tibble(rownames='Model') %>%
                                                                                        mutate(RSE = c(null$sigma,summary(basic_linear)$sigma,basic_nls$sigma,full_model$sigma)),
                                                                                      anova(null,full_model))) ) %>% 
  
  
  mutate(data = pmap(.l=list(data,length_imputation,width_imputation), 
                     .f=function(x,y,z) x %>%
                       mutate(Imputed=if_else(is.na(Length) & !is.na(Width),'Length',
                                              if_else(is.na(Width) & !is.na(Length),'Width','No')),
                              Length=ifelse(is.na(Length),predict(y,x,na.action = na.pass),Length),
                              Width=ifelse(is.na(Width),predict(z,x,na.action = na.pass),Width)))) %>%
  
  #Remove immature shells - those less than 29 mm and those with missing lengths
  #mutate(data = map(data, ~filter(.x,Length > 29))) %>%
  mutate(data = map(data, ~filter(.x,!is.na(Length)))) %>%
  
  #Make plots to check imputation
  mutate(length_imputation_plot = map2(data,length_imputation, ~.x %>%
                                         ggplot(aes(x=Width,y=Length)) +
                                         geom_point(aes(colour=Imputed)) +
                                         geom_smooth(method='lm',colour='blue',se=F) + 
                                         geom_line(data=(.x %$%
                                                           seq(min(Width),max(Width),length.out = 1000) %>%
                                                           tibble(Width=.) %>%
                                                           mutate(Length=predict(.y,.))),colour='red') +
                                         scale_colour_manual(values=c('red','gray50','blue'))),
         width_imputation_plot = map2(data,width_imputation, ~.x %>%
                                        ggplot(aes(x=Length,y=Width)) +
                                        geom_point(aes(colour=Imputed)) +
                                        geom_smooth(method='lm',colour='blue') + 
                                        geom_line(data=(.x %$%
                                                          seq(min(Length),max(Length),length.out = 1000) %>%
                                                          tibble(Length=.) %>%
                                                          mutate(Width=predict(.y,.))),colour='red') +
                                        scale_colour_manual(values=c('red','gray50','blue')))) %>%
  
  #Check imputation assumptuons
  mutate(length_imputation_assumptions = map(length_imputation,ggdiag_plots,independent = c('Width'), dependent = c('Length')),
         width_imputation_assumptions = map(length_imputation,ggdiag_plots,independent = c('Length'), dependent = c('Width')))

#### Process ImageJ Data ####
ImageJ_Data<-ImageJ_Data %>%
  mutate(Curve=ifelse(Curve==0,NA,Curve)) %>%
  dplyr::select(ShellID,Point,X,Y,Curve) %>%
  gather(Axis,position,-ShellID,-Point,-Curve) %>%
  mutate(Point=if_else(Point==1,'X1',
                       if_else(Point==2,'X2',
                               if_else(Point==3,'X3',
                                       if_else(Point==4,'ScalePt1X',
                                               if_else(Point==5,'ScalePt2X','PosteriorLengthX'))))),
         Point=if_else(Axis=='X',Point,str_replace(Point,'X','Y'))) %>%
  dplyr::select(-Axis) %>%
  spread(Point,-ShellID) %>% 
  group_by(ShellID) %>%
  summarise_all(mean,na.rm=T) %>%
  
  #Get distances from points
  by_row(triangulate,.collate = 'cols') %>%
  rename(L13=L131,L26=L261,L36=L361,L12=L121,L16=L161,A213=A2131) %>% #by_row adds a 1 to the columns it makes. Not sure why or how to prevent yet so I remove that here
  
  #Scale distances to mm
  mutate(ScalingXY = sqrt((ScalePt2X - ScalePt1X)^2 + (ScalePt2Y - ScalePt1Y)^2),
         Posterior_Apex_Position = (L36*5)/ScalingXY,
         L26 = L26*5/ScalingXY,
         L16 = L16*5/ScalingXY,
         Curve = Curve*5/ScalingXY,
         L12 = L12*5/ScalingXY,
         L13 = L13*5/ScalingXY,
         Length_Picture = L26+Posterior_Apex_Position) 

#### Normalize by Length ####
opihi_data<-caliper_measures %>%
  dplyr::select(data) %>%
  unnest %>%
  inner_join(ImageJ_Data) %>% 
  
  #Calibrate Image lengths
  mutate(Calibration_coef = Length/Length_Picture, #Is this the proper way to calibrate for differences between the picture and caliper measures
         Posterior_Apex_Position = Posterior_Apex_Position * Calibration_coef,
         L26 = L26 * Calibration_coef,
         L16 = L16 * Calibration_coef,
         L12 = L12 * Calibration_coef,
         L13 = L13 * Calibration_coef,
         Curve = Curve * Calibration_coef) %>%
  nest %>%
  
  #Normalize measurements to account for allometry
  mutate(width_normalization = map(data, ~ .x %>%
                                     dplyr::select(Length,Width) %>%
                                     na.omit() %>%
                                     gnls(Width ~ a*Length^b, data=., start = c(a=1,b=1), weights = varExp(form=~fitted(.)))),
         width_normalization_gnls = map(data, ~ .x %>%
                                          dplyr::select(Length,Width) %>%
                                          na.omit() %>%
                                          gnls(Width ~ a*Length^b, data=., start = c(a=1,b=1))),
         width_normalization_lm = map(data, ~ .x %>%
                                        dplyr::select(Length,Width) %>%
                                        na.omit() %>%
                                        lm(Width ~ Length, data=.)),
         width_normalization_null = map(data, ~ .x %>%
                                          dplyr::select(Length,Width) %>%
                                          na.omit() %>%
                                          gnls(Width ~ a, 
                                               start = c(a=1),
                                               data=.)),
         width_normalization_summary = pmap(list(full_model = width_normalization, 
                                                 basic_nls = width_normalization_gnls,
                                                 basic_linear = width_normalization_lm,
                                                 null = width_normalization_null),
                                            function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                                    as_tibble(rownames='Model') %>%
                                                                                                    mutate(RSE = c(null$sigma,
                                                                                                                   summary(basic_linear)$sigma,
                                                                                                                   basic_nls$sigma,
                                                                                                                   full_model$sigma)),
                                                                                                  anova(null,full_model))),
         
         
         
         
         height_normalization = map(data, ~ .x %>%
                                      dplyr::select(Length,Height) %>%
                                      na.omit() %>%
                                      mutate(Length=log(Length),Height=log(Height)) %>%
                                      gnls(Height ~ a*Length^b, data=., start = c(a=1,b=1))),
         height_normalization_gnls = map(data, ~ .x %>%
                                           dplyr::select(Length,Height) %>%
                                           na.omit() %>%
                                           mutate(Length=log(Length),Height=log(Height)) %>%
                                           gnls(Height ~ a*Length^b, data=., start = c(a=1,b=1))),
         height_normalization_lm = map(data, ~ .x %>%
                                         dplyr::select(Length,Height) %>%
                                         na.omit() %>%
                                         mutate(Length=log(Length),Height=log(Height)) %>%
                                         lm(Height ~ Length, data=.)),
         height_normalization_null = map(data, ~ .x %>%
                                           dplyr::select(Length,Height) %>%
                                           na.omit() %>%
                                           mutate(Length=log(Length),Height=log(Height)) %>%
                                           gnls(Height ~ a, 
                                                start = c(a=1),
                                                data=.)),
         height_normalization_summary = pmap(list(full_model = height_normalization, 
                                                  basic_nls = height_normalization_gnls,
                                                  basic_linear = height_normalization_lm,
                                                  null = height_normalization_null),
                                             function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                                     as_tibble(rownames='Model') %>%
                                                                                                     mutate(RSE = c(null$sigma,
                                                                                                                    summary(basic_linear)$sigma,
                                                                                                                    basic_nls$sigma,
                                                                                                                    full_model$sigma)),
                                                                                                   anova(null,full_model))),
         
         posterioir_apex_position_normalization = map(data, ~ .x %>%
                                                        dplyr::select(Length,Posterior_Apex_Position) %>%
                                                        na.omit() %>%
                                                        mutate(Length=log(Length),Posterior_Apex_Position=log(Posterior_Apex_Position)) %>%
                                                        gnls(Posterior_Apex_Position ~ a*Length^b, data=., start = c(a=1,b=1))),
         posterioir_apex_position_normalization_gnls = map(data, ~ .x %>%
                                                             dplyr::select(Length,Posterior_Apex_Position) %>%
                                                             na.omit() %>%
                                                             mutate(Length=log(Length),Posterior_Apex_Position=log(Posterior_Apex_Position)) %>%
                                                             gnls(Posterior_Apex_Position ~ a*Length^b, data=., start = c(a=1,b=1))),
         posterioir_apex_position_normalization_lm = map(data, ~ .x %>%
                                                           dplyr::select(Length,Posterior_Apex_Position) %>%
                                                           na.omit() %>%
                                                           mutate(Length=log(Length),Posterior_Apex_Position=log(Posterior_Apex_Position)) %>%
                                                           lm(Posterior_Apex_Position ~ Length, data=.)),
         posterioir_apex_position_normalization_null = map(data, ~ .x %>%
                                                             dplyr::select(Length,Posterior_Apex_Position) %>%
                                                             na.omit() %>%
                                                             mutate(Length=log(Length),Posterior_Apex_Position=log(Posterior_Apex_Position)) %>%
                                                             gnls(Posterior_Apex_Position ~ a, 
                                                                  start = c(a=1),
                                                                  data=.)),
         posterioir_apex_position_normalization_summary = pmap(list(full_model = posterioir_apex_position_normalization, 
                                                                    basic_nls = posterioir_apex_position_normalization_gnls,
                                                                    basic_linear = posterioir_apex_position_normalization_lm,
                                                                    null = posterioir_apex_position_normalization_null),
                                                               function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                                                       as_tibble(rownames='Model') %>%
                                                                                                                       mutate(RSE = c(null$sigma,
                                                                                                                                      summary(basic_linear)$sigma,
                                                                                                                                      basic_nls$sigma,
                                                                                                                                      full_model$sigma)),
                                                                                                                     anova(null,full_model))),
         
         curve_normalization = map(data, ~ .x %>%
                                     dplyr::select(Length,Curve) %>%
                                     na.omit() %>%
                                     mutate(Length=log(Length),Curve=log(Curve)) %>%
                                     gnls(Curve ~ a*Length^b, data=., start = c(a=1,b=1), weights = varExp(form=~fitted(.)))),
         curve_normalization_gnls = map(data, ~ .x %>%
                                          dplyr::select(Length,Curve) %>%
                                          na.omit() %>%
                                          mutate(Length=log(Length),Curve=log(Curve)) %>%
                                          gnls(Curve ~ a*Length^b, data=., start = c(a=1,b=1))),
         curve_normalization_lm = map(data, ~ .x %>%
                                        dplyr::select(Length,Curve) %>%
                                        na.omit() %>%
                                        mutate(Length=log(Length),Curve=log(Curve)) %>%
                                        lm(Curve ~ Length, data=.)),
         curve_normalization_null = map(data, ~ .x %>%
                                          dplyr::select(Length,Curve) %>%
                                          na.omit() %>%
                                          mutate(Length=log(Length),Curve=log(Curve)) %>%
                                          gnls(Curve ~ a, 
                                               start = c(a=1),
                                               data=.)),
         curve_normalization_summary = pmap(list(full_model = curve_normalization, 
                                                 basic_nls = curve_normalization_gnls,
                                                 basic_linear = curve_normalization_lm,
                                                 null = curve_normalization_null),
                                            function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                                    as_tibble(rownames='Model') %>%
                                                                                                    mutate(RSE = c(null$sigma,
                                                                                                                   summary(basic_linear)$sigma,
                                                                                                                   basic_nls$sigma,
                                                                                                                   full_model$sigma)),
                                                                                                  anova(null,full_model))),
         
         L13_normalization = map(data, ~ .x %>%
                                   dplyr::select(Length,L13) %>%
                                   na.omit() %>% 
                                   mutate(Length=log(Length),L13=log(L13)) %>%
                                   gnls(L13 ~ a*Length^b, data=., start = c(a=1,b=1), weights = varExp(form=~fitted(.)))),
         L13_normalization_gnls = map(data, ~ .x %>%
                                        dplyr::select(Length,L13) %>%
                                        na.omit() %>%
                                        mutate(Length=log(Length),L13=log(L13)) %>%
                                        gnls(L13 ~ a*Length^b, data=., start = c(a=1,b=1))),
         L13_normalization_lm = map(data, ~ .x %>%
                                      dplyr::select(Length,L13) %>%
                                      na.omit() %>%
                                      mutate(Length=log(Length),L13=log(L13)) %>%
                                      lm(L13 ~ Length, data=.)),
         L13_normalization_null = map(data, ~ .x %>%
                                        dplyr::select(Length,L13) %>%
                                        na.omit() %>%
                                        mutate(Length=log(Length),L13=log(L13)) %>%
                                        gnls(L13 ~ a, 
                                             start = c(a=1),
                                             data=.)),
         L13_normalization_summary = pmap(list(full_model = L13_normalization, 
                                               basic_nls = L13_normalization_gnls,
                                               basic_linear = L13_normalization_lm,
                                               null = L13_normalization_null),
                                          function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                                  as_tibble(rownames='Model') %>%
                                                                                                  mutate(RSE = c(null$sigma,
                                                                                                                 summary(basic_linear)$sigma,
                                                                                                                 basic_nls$sigma,
                                                                                                                 full_model$sigma)),
                                                                                                anova(null,full_model))),
         
         L16_normalization = map(data, ~ .x %>%
                                   dplyr::select(Length,L16) %>%
                                   na.omit() %>%
                                   mutate(Length=log(Length),L16=log(L16)) %>%
                                   gnls(L16 ~ a*Length^b, data=., start = c(a=1,b=1))),
         L16_normalization_gnls = map(data, ~ .x %>%
                                        dplyr::select(Length,L16) %>%
                                        na.omit() %>%
                                        mutate(Length=log(Length),L16=log(L16)) %>%
                                        gnls(L16 ~ a*Length^b, data=., start = c(a=1,b=1))),
         L16_normalization_lm = map(data, ~ .x %>%
                                      dplyr::select(Length,L16) %>%
                                      na.omit() %>%
                                      mutate(Length=log(Length),L16=log(L16)) %>%
                                      lm(L16 ~ Length, data=.)),
         L16_normalization_null = map(data, ~ .x %>%
                                        dplyr::select(Length,L16) %>%
                                        na.omit() %>%
                                        mutate(Length=log(Length),L16=log(L16)) %>%
                                        gnls(L16 ~ a, 
                                             start = c(a=1),
                                             data=.)),
         L16_normalization_summary = pmap(list(full_model = L16_normalization, 
                                               basic_nls = L16_normalization_gnls,
                                               basic_linear = L16_normalization_lm,
                                               null = L16_normalization_null),
                                          function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                                  as_tibble(rownames='Model') %>%
                                                                                                  mutate(RSE = c(null$sigma,
                                                                                                                 summary(basic_linear)$sigma,
                                                                                                                 basic_nls$sigma,
                                                                                                                 full_model$sigma)),
                                                                                                anova(null,full_model))),
         
         L26_normalization = map(data, ~ .x %>%
                                   dplyr::select(Length,L26) %>%
                                   na.omit() %>%
                                   mutate(Length=log(Length),L26=log(L26)) %>%
                                   gnls(L26 ~ a*Length^b, data=., start = c(a=1,b=1))),
         L26_normalization_gnls = map(data, ~ .x %>%
                                        dplyr::select(Length,L26) %>%
                                        na.omit() %>%
                                        mutate(Length=log(Length),L26=log(L26)) %>%
                                        gnls(L26 ~ a*Length^b, data=., start = c(a=1,b=1))),
         L26_normalization_lm = map(data, ~ .x %>%
                                      dplyr::select(Length,L26) %>%
                                      na.omit() %>%
                                      mutate(Length=log(Length),L26=log(L26)) %>%
                                      lm(L26 ~ Length, data=.)),
         L26_normalization_null = map(data, ~ .x %>%
                                        dplyr::select(Length,L26) %>%
                                        na.omit() %>%
                                        mutate(Length=log(Length),L26=log(L26)) %>%
                                        gnls(L26 ~ a, 
                                             start = c(a=1),
                                             data=.)),
         L26_normalization_summary = pmap(list(full_model = L26_normalization, 
                                               basic_nls = L26_normalization_gnls,
                                               basic_linear = L26_normalization_lm,
                                               null = L26_normalization_null),
                                          function(full_model,basic_nls,basic_linear,null) list(AIC(null,basic_linear,basic_nls,full_model) %>%
                                                                                                  as_tibble(rownames='Model') %>%
                                                                                                  mutate(RSE = c(null$sigma,
                                                                                                                 summary(basic_linear)$sigma,
                                                                                                                 basic_nls$sigma,
                                                                                                                 full_model$sigma)),
                                                                                                anova(null,full_model)))) %>%
  
  
  
  
  mutate(data = pmap(.l=list(data,width_normalization,height_normalization,
                             posterioir_apex_position_normalization,curve_normalization,
                             L13_normalization,
                             L16_normalization,
                             L26_normalization),
                     .f=function(x,y,z,A,B,C,E,f) x %>% #D
                       mutate(Width_normalized = Width * (mean(Length)/Length)^coef(y)[2],
                              Height_normalized = exp(log(Height) * (log(mean(Length))/log(Length))^coef(z)[2]),
                              Posterior_Apex_Position_normalized = exp(log(Posterior_Apex_Position) * (log(mean(Length))/log(Length))^coef(A)[2]),
                              Curve_normalized = exp(log(Curve) * (log(mean(Length))/log(Length))^coef(B)[2]),
                              L13_normalized = exp(log(L13) * (log(mean(Length))/log(Length))^coef(C)[2]),
                              L16_normalized = exp(log(L16) * (log(mean(Length))/log(Length))^coef(E)[2]),
                              L26_normalized = exp(log(L26) * (log(mean(Length))/log(Length))^coef(f)[2])))) %>%
  
  #Make graphs from normalization process
  mutate(width_normalization_plot = map2(data,width_normalization, ~.x %>%
                                           dplyr::select(Island,ShellID,Length,Width,Width_normalized) %>%
                                           gather(Type,Width,-ShellID,-Length,-Island) %>%
                                           ggplot(aes(x=Length,y=Width)) +
                                           geom_point(aes(shape=Type)) +
                                           geom_line(data=(.x %$%
                                                             seq(min(Length),max(Length),length.out = 1000) %>%
                                                             tibble(Length=.) %>%
                                                             mutate(Width=twNlme::varPredictNlmeGnls(.y,.)[,'fit']))) +
                                           geom_ribbon(data=(.x %$%
                                                               seq(min(Length),max(Length),length.out = 1000) %>%
                                                               tibble(Length=.) %>%
                                                               mutate(Width=twNlme::varPredictNlmeGnls(.y,.)[,'fit'],
                                                                      sdPop=twNlme::varPredictNlmeGnls(.y,.)[,'sdPop'])),aes(ymin=Width-sdPop,ymax=Width+sdPop),alpha=0.4) +
                                           scale_shape_manual(name='Type',labels=c('Non-normalized','Normalized'),values=c(1,19)) +
                                           geom_vline(aes(xintercept=mean(Length)),lty=2)),
         height_normalization_plot = map2(data,height_normalization, ~.x %>%
                                            dplyr::select(Island,ShellID,Length,Height,Height_normalized) %>%
                                            gather(Type,Height,-ShellID,-Length,-Island) %>%
                                            ggplot(aes(x=Length,y=Height)) +
                                            geom_point(aes(shape=Type)) +
                                            geom_line(data=(.x %$%
                                                              seq(min(Length),max(Length),length.out = 1000) %>%
                                                              tibble(Length=.) %>%
                                                              mutate(Length=log(Length)) %>%
                                                              mutate(Height=twNlme::varPredictNlmeGnls(.y,.)[,'fit']) %>%
                                                              mutate(Length=exp(Length),
                                                                     Height=exp(Height)))) +
                                            geom_ribbon(data=(.x %$%
                                                                seq(min(Length),max(Length),length.out = 1000) %>%
                                                                tibble(Length=.) %>%
                                                                mutate(Length=log(Length)) %>%
                                                                mutate(Height=twNlme::varPredictNlmeGnls(.y,.)[,'fit'],
                                                                       sdPop=twNlme::varPredictNlmeGnls(.y,.)[,'sdPop']) %>%
                                                                mutate(Length=exp(Length),
                                                                       Height=exp(Height),
                                                                       sdPop=exp(sdPop))),aes(ymin=Height-sdPop,ymax=Height+sdPop),alpha=0.4) +
                                            scale_shape_manual(name='Type',labels=c('Non-normalized','Normalized'),values=c(1,19)) +
                                            geom_vline(aes(xintercept=mean(Length)),lty=2)),
         Posterior_Apex_Position_normalization_plot = map2(data,posterioir_apex_position_normalization, ~.x %>%
                                                             dplyr::select(ShellID,Length,Posterior_Apex_Position,Posterior_Apex_Position_normalized,Island) %>%
                                                             gather(Type,Posterior_Apex_Position,-ShellID,-Length,-Island) %>%
                                                             ggplot(aes(x=Length,y=Posterior_Apex_Position)) +
                                                             geom_point(aes(shape=Type)) +
                                                             geom_line(data=(.x %$%
                                                                               seq(min(Length),max(Length),length.out = 1000) %>%
                                                                               tibble(Length=.) %>%
                                                                               mutate(Length=log(Length)) %>%
                                                                               mutate(Posterior_Apex_Position=twNlme::varPredictNlmeGnls(.y,.)[,'fit']) %>%
                                                                               mutate(Length=exp(Length),
                                                                                      Posterior_Apex_Position=exp(Posterior_Apex_Position)))) +
                                                             geom_ribbon(data=(.x %$%
                                                                                 seq(min(Length),max(Length),length.out = 1000) %>%
                                                                                 tibble(Length=.) %>%
                                                                                 mutate(Length=log(Length)) %>%
                                                                                 mutate(Posterior_Apex_Position=twNlme::varPredictNlmeGnls(.y,.)[,'fit'],
                                                                                        sdPop=twNlme::varPredictNlmeGnls(.y,.)[,'sdPop']) %>%
                                                                                 mutate(Length=exp(Length),
                                                                                        Posterior_Apex_Position=exp(Posterior_Apex_Position),
                                                                                        sdPop=exp(sdPop))),aes(ymin=Posterior_Apex_Position-sdPop,ymax=Posterior_Apex_Position+sdPop),alpha=0.4) +
                                                             scale_shape_manual(name='Type',labels=c('Non-normalized','Normalized'),values=c(1,19)) +
                                                             geom_vline(aes(xintercept=mean(Length)),lty=2)),
         Curve_normalization_plot = map2(data,curve_normalization, ~.x %>%
                                           dplyr::select(ShellID,Length,Curve,Curve_normalized,Island) %>%
                                           gather(Type,Curve,-ShellID,-Length,-Island) %>%
                                           ggplot(aes(x=Length,y=Curve)) +
                                           geom_point(aes(shape=Type)) +
                                           geom_line(data=(.x %$%
                                                             seq(min(Length),max(Length),length.out = 1000) %>%
                                                             tibble(Length=.) %>%
                                                             mutate(Length=log(Length)) %>%
                                                             mutate(Curve=twNlme::varPredictNlmeGnls(.y,.)[,'fit']) %>%
                                                             mutate(Length=exp(Length),
                                                                    Curve=exp(Curve)))) +
                                           geom_ribbon(data=(.x %$%
                                                               seq(min(Length),max(Length),length.out = 1000) %>%
                                                               tibble(Length=.) %>%
                                                               mutate(Length=log(Length)) %>%
                                                               mutate(Curve=twNlme::varPredictNlmeGnls(.y,.)[,'fit'],
                                                                      sdPop=twNlme::varPredictNlmeGnls(.y,.)[,'sdPop']) %>%
                                                               mutate(Length=exp(Length),
                                                                      Curve=exp(Curve),
                                                                      sdPop=exp(sdPop))),aes(ymin=Curve-sdPop,ymax=Curve+sdPop),alpha=0.4) +
                                           scale_shape_manual(name='Type',labels=c('Non-normalized','Normalized'),values=c(1,19)) +
                                           geom_vline(aes(xintercept=mean(Length)),lty=2)),
         L13_normalization_plot = map2(data,L13_normalization, ~.x %>%
                                         dplyr::select(ShellID,Length,L13,L13_normalized,Island) %>%
                                         gather(Type,L13,-ShellID,-Length,-Island) %>%
                                         ggplot(aes(x=Length,y=L13)) +
                                         geom_point(aes(shape=Type)) +
                                         geom_line(data=(.x %$%
                                                           seq(min(Length),max(Length),length.out = 1000) %>%
                                                           tibble(Length=.) %>%
                                                           mutate(Length=log(Length)) %>%
                                                           mutate(L13=twNlme::varPredictNlmeGnls(.y,.)[,'fit']) %>%
                                                           mutate(Length=exp(Length),
                                                                  L13=exp(L13)))) +
                                         geom_ribbon(data=(.x %$%
                                                             seq(min(Length),max(Length),length.out = 1000) %>%
                                                             tibble(Length=.) %>%
                                                             mutate(Length=log(Length)) %>%
                                                             mutate(L13=twNlme::varPredictNlmeGnls(.y,.)[,'fit'],
                                                                    sdPop=twNlme::varPredictNlmeGnls(.y,.)[,'sdPop']) %>%
                                                             mutate(Length=exp(Length),
                                                                    L13=exp(L13),
                                                                    sdPop=exp(sdPop))),aes(ymin=L13-sdPop,ymax=L13+sdPop),alpha=0.4) +
                                         scale_shape_manual(name='Type',labels=c('Non-normalized','Normalized'),values=c(1,19)) +
                                         geom_vline(aes(xintercept=mean(Length)),lty=2)),
         L16_normalization_plot = map2(data,L16_normalization, ~.x %>%
                                         dplyr::select(ShellID,Length,L16,L16_normalized,Island) %>%
                                         gather(Type,L16,-ShellID,-Length,-Island) %>%
                                         ggplot(aes(x=Length,y=L16)) +
                                         geom_point(aes(shape=Type)) +
                                         geom_line(data=(.x %$%
                                                           seq(min(Length),max(Length),length.out = 1000) %>%
                                                           tibble(Length=.) %>%
                                                           mutate(Length=log(Length)) %>%
                                                           mutate(L16=twNlme::varPredictNlmeGnls(.y,.)[,'fit']) %>%
                                                           mutate(Length=exp(Length),
                                                                  L16=exp(L16)))) +
                                         geom_ribbon(data=(.x %$%
                                                             seq(min(Length),max(Length),length.out = 1000) %>%
                                                             tibble(Length=.) %>%
                                                             mutate(Length=log(Length)) %>%
                                                             mutate(L16=twNlme::varPredictNlmeGnls(.y,.)[,'fit'],
                                                                    sdPop=twNlme::varPredictNlmeGnls(.y,.)[,'sdPop']) %>%
                                                             mutate(Length=exp(Length),
                                                                    L16=exp(L16),
                                                                    sdPop=exp(sdPop))),aes(ymin=L16-sdPop,ymax=L16+sdPop),alpha=0.4) +
                                         scale_shape_manual(name='Type',labels=c('Non-normalized','Normalized'),values=c(1,19)) +
                                         geom_vline(aes(xintercept=mean(Length)),lty=2)),
         L26_normalization_plot = map2(data,L26_normalization, ~.x %>%
                                         dplyr::select(ShellID,Length,L26,L26_normalized,Island) %>%
                                         gather(Type,L26,-ShellID,-Length,-Island) %>%
                                         ggplot(aes(x=Length,y=L26)) +
                                         geom_point(aes(shape=Type)) +
                                         geom_line(data=(.x %$%
                                                           seq(min(Length),max(Length),length.out = 1000) %>%
                                                           tibble(Length=.) %>%
                                                           mutate(Length=log(Length)) %>%
                                                           mutate(L26=twNlme::varPredictNlmeGnls(.y,.)[,'fit']) %>%
                                                           mutate(Length=exp(Length),
                                                                  L26=exp(L26)))) +
                                         geom_ribbon(data=(.x %$%
                                                             seq(min(Length),max(Length),length.out = 1000) %>%
                                                             tibble(Length=.) %>%
                                                             mutate(Length=log(Length)) %>%
                                                             mutate(L26=twNlme::varPredictNlmeGnls(.y,.)[,'fit'],
                                                                    sdPop=twNlme::varPredictNlmeGnls(.y,.)[,'sdPop']) %>%
                                                             mutate(Length=exp(Length),
                                                                    L26=exp(L26),
                                                                    sdPop=exp(sdPop))),aes(ymin=L26-sdPop,ymax=L26+sdPop),alpha=0.4) +
                                         scale_shape_manual(name='Type',labels=c('Non-normalized','Normalized'),values=c(1,19)) +
                                         geom_vline(aes(xintercept=mean(Length)),lty=2))) %>%
  
  #Check normalization assumptuons
  mutate(width_normalization_assumptions = map(width_normalization,ggdiag_plots,dependent = c('Width'), independent = c('Length')),
         height_normalization_assumptions = map(height_normalization,ggdiag_plots,dependent = c('Length'), independent = c('Length')),
         Posterior_Apex_Position_normalization_assumptions = map(posterioir_apex_position_normalization,ggdiag_plots,dependent = c('Posterior_Apex_Position'), independent = c('Length')),
         Curve_normalization_assumptions = map(curve_normalization,ggdiag_plots,dependent = c('Curve'), independent = c('Length')),
         L13_normalization_assumptions = map(L13_normalization,ggdiag_plots,dependent = c('L13'), independent = c('Length')),
         L16_normalization_assumptions = map(L16_normalization,ggdiag_plots,dependent = c('L16'), independent = c('Length')),
         L26_normalization_assumptions = map(L26_normalization,ggdiag_plots,dependent = c('L26'), independent = c('Length'))) %>%
  
  #Make new indicies
  mutate(data = map(data, ~.x %>%
                      mutate(Circ.hvalue=((mean(Length)/2-Width_normalized/2)^2) / (mean(Length)/2 + Width_normalized/2)^2, 
                             Circumference = pi*(mean(Length)/2 + Width_normalized/2)^2 * (1+ 3*Circ.hvalue / 10 + sqrt(4-3*Circ.hvalue)), 
                             ApertureArea = pi * mean(Length)/2 * Width_normalized/2,
                             Volume = (pi/3) * (mean(Length)/2 * Width_normalized/2 * Height_normalized),
                             SurfaceArea = SurfArea(mean(Length),Width_normalized,Height_normalized))))

#### Bayesian Analysis of Shell Morphology ####
opihi_shells<-opihi_data %>%
  dplyr::select(data) %>%
  unnest(cols = c(data)) %>%
  dplyr::select(ShellID, Island, SamplingSite, Length, SurfaceArea, Curve_normalized, 
                L13_normalized, Height_normalized, Width_normalized,
                Imputed, starts_with('RibTipColor')) %>%
  filter(Length > 29) %>%
  mutate(C13 = Curve_normalized/L13_normalized,
         Height_Index = Height_normalized/mean(Length),
         SamplingSite = if_else(is.na(SamplingSite) | SamplingSite=='', Island, SamplingSite),
         SurfaceArea = SurfaceArea/100) %>%
  mutate(RibTipColor= (RibTipColor_AnteriorSinistrall + RibTipColor_PosteriorSinistral + 
                         RibTipColor_PosteriorDextral + RibTipColor_AnteriorDextral)/4,
         Morphotype=if_else(RibTipColor <= 50, 'White','Black')) %>%
  filter(!is.na(Morphotype)) %>%
  ungroup %>%
  dplyr::select(-RibTipColor_AnteriorSinistrall:-RibTipColor_AnteriorDextral) %>%
  inner_join(Sampling_Sites)


#### Model SA ####
opihi_shells %>%
  ggplot(aes(x=SurfaceArea)) +
  geom_histogram(bins = 30)

the_model_SA_null<- opihi_shells %>%
  filter(!is.na(SurfaceArea)) %>%
  brm(data = ., family = gaussian,
      bf(SurfaceArea ~ 1
      ),
      iter = 2000, warmup = 1000, chains = 4,
      sample_prior = TRUE,
      control = list(adapt_delta = 0.8, max_treedepth = 10),
      future = TRUE, seed = 12, refresh = -1,
      file = '../Results/SurfaceArea_null') 
summary(the_model_SA_null)
plot(the_model_SA_null) 
pp_check(the_model_SA_null, nsamples = 100)
the_model_SA_null <- add_criterion(the_model_SA_null, "loo")
the_model_SA_null <- add_criterion(the_model_SA_null, "waic")

the_model_SA_basic<- opihi_shells %>%
  filter(!is.na(SurfaceArea)) %>%
  brm(data = ., family = gaussian,
      bf(SurfaceArea ~ latitude.c*island_group + 
           (1 | Island) + (1 | Island:SamplingSite)
      ),
      iter = 2000, warmup = 1000, chains = 4,
      sample_prior = TRUE,
      control = list(adapt_delta = 0.99, stepsize = 1, max_treedepth = 10),
      future = TRUE, seed = 12, refresh = -1,
      file = '../Results/SurfaceArea_basic') 
summary(the_model_SA_basic)
plot(the_model_SA_basic) 
pp_check(the_model_SA_basic, nsamples = 100)
the_model_SA_basic <- add_criterion(the_model_SA_basic, "loo")
the_model_SA_basic <- add_criterion(the_model_SA_basic, "waic")

the_model_SA_basic2<- opihi_shells %>%
  filter(!is.na(SurfaceArea)) %>%
  brm(data = ., family = gaussian,
      bf(SurfaceArea ~ latitude.c*island_group + 
           (1 |ID1| Island) + (1 |ID2| Island:SamplingSite),
         sigma ~ (1 |ID1| Island) + (1 |ID2| Island:SamplingSite)
      ),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.99, stepsize = 0.1, max_treedepth = 10),
      future = TRUE, seed = 12, refresh = -1,
      sample_prior = TRUE,
      file = '../Results/SurfaceArea_basic2')
summary(the_model_SA_basic2)
plot(the_model_SA_basic2) 
pp_check(the_model_SA_basic2, nsamples = 100)
the_model_SA_basic2 <- add_criterion(the_model_SA_basic2, "loo")
the_model_SA_basic2 <- add_criterion(the_model_SA_basic2, "waic")

the_model_SA_delux <- opihi_shells %>%
  filter(!is.na(SurfaceArea)) %>%
  brm(data = ., family = gaussian,
      bf(SurfaceArea ~ latitude.c * island_group + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite),
         sigma ~ island_group + Location.c : island_group + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite)
      ),
      prior = c(prior(gamma(1, 0.5), class = sd),
                prior(gamma(1, 0.5), class = sd, dpar = sigma),
                prior(lkj(2), class = cor)),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.99, stepsize = 0.1, max_treedepth = 20),
      future = TRUE, seed = 12, refresh = -1,
      inits = 0,
      sample_prior = TRUE,
      file = '../Results/SurfaceArea_delux')
summary(the_model_SA_delux) 
plot(the_model_SA_delux) 
pp_check(the_model_SA_delux, nsamples = 100)
the_model_SA_delux <- add_criterion(the_model_SA_delux, "loo")
the_model_SA_delux <- add_criterion(the_model_SA_delux, "waic")

## Model Comparison ##
tibble(model = 1:4, 
       loo_ic = c(the_model_SA_null$loo$estimates[3,1], the_model_SA_basic$loo$estimates[3,1], 
                  the_model_SA_basic2$loo$estimates[3,1], the_model_SA_delux$loo$estimates[3,1]),
       loo_weights = model_weights(the_model_SA_null, the_model_SA_basic, 
                                   the_model_SA_basic2, the_model_SA_delux, 
                                   criterion = 'loo') %>%
         round(3),
       R2_lwr = c(bayes_R2(the_model_SA_null)[3],bayes_R2(the_model_SA_basic)[3],
              bayes_R2(the_model_SA_basic2)[3],bayes_R2(the_model_SA_delux)[3]) %>% 
         round(3),
       R2_upr = c(bayes_R2(the_model_SA_null)[4],bayes_R2(the_model_SA_basic)[4],
                  bayes_R2(the_model_SA_basic2)[4],bayes_R2(the_model_SA_delux)[4]) %>% 
         round(3))

#### Model Curve/L13 ####
opihi_shells %>%
  ggplot(aes(x=C13)) +
  geom_histogram(bins = 50)

the_model_Curve_null <- opihi_shells %>%
  filter(!is.na(C13)) %>%
  brm(data = ., family = gaussian,
      bf(C13 ~ 1
      ),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.9, max_treedepth = 10),
      sample_prior = TRUE,
      future = TRUE, seed = 12, refresh = -1,
      file = '../Results/C13_null') 
summary(the_model_Curve_null)
plot(the_model_Curve_null) 
pp_check(the_model_Curve_null, nsamples = 100)
the_model_Curve_null <- add_criterion(the_model_Curve_null, "loo")
the_model_Curve_null <- add_criterion(the_model_Curve_null, "waic")

the_model_Curve_basic <- opihi_shells %>%
  filter(!is.na(C13)) %>%
  brm(data = ., family = gaussian,
      bf(C13 ~ latitude.c * island_group + 
           (1 | Island) + (1 | Island:SamplingSite)
      ),
      prior = c(prior(gamma(1, 0.5), class = sd)),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20),
      sample_prior = TRUE,
      future = TRUE, seed = 12, refresh = -1,
      file = '../Results/C13_basic') 
summary(the_model_Curve_basic)
plot(the_model_Curve_basic) 
pp_check(the_model_Curve_basic, nsamples = 100)
the_model_Curve_basic <- add_criterion(the_model_Curve_basic, "loo")
the_model_Curve_basic <- add_criterion(the_model_Curve_basic, "waic")

the_model_Curve_basic2 <- opihi_shells %>%
  filter(!is.na(C13)) %>%
  brm(data = ., family = gaussian,
      bf(C13 ~ latitude.c * island_group + 
           (1 |ID1| Island) + (1 |ID2| Island:SamplingSite),
         sigma ~ (1 |ID1| Island) + (1 |ID2| Island:SamplingSite)
      ),
      prior = c(prior(gamma(1, 0.5), class = sd),
                prior(gamma(1, 0.5), class = sd, dpar = sigma),
                prior(lkj(2), class = cor)),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 20),
      future = TRUE, seed = 12, refresh = -1,
      sample_prior = TRUE,
      file = '../Results/C13_basic2')
summary(the_model_Curve_basic2)
plot(the_model_Curve_basic2) 
pp_check(the_model_Curve_basic2, nsamples = 100)
the_model_Curve_basic2 <- add_criterion(the_model_Curve_basic2, "loo")
the_model_Curve_basic2 <- add_criterion(the_model_Curve_basic2, "waic")

the_model_C13_delux<- opihi_shells %>%
  filter(!is.na(C13)) %>%
  brm(data = ., family = gaussian,
      bf(C13 ~ latitude.c * island_group + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite),
         sigma ~ island_group + Location.c : island_group + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite)
      ),
      prior = c(prior(gamma(1, 0.5), class = sd),
                prior(gamma(1, 0.5), class = sd, dpar = sigma),
                prior(lkj(2), class = cor)),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 10),
      future = TRUE, seed = 12, refresh = -1,
      inits = 0,
      sample_prior = TRUE,
      file = '../Results/C13_delux')
summary(the_model_C13_delux)
plot(the_model_C13_delux) 
pp_check(the_model_C13_delux, nsamples = 100)
the_model_C13_delux <- add_criterion(the_model_C13_delux, "loo", reloo = TRUE)
the_model_C13_delux <- add_criterion(the_model_C13_delux, "waic")

## Model Comparison ##
tibble(model = 1:4, 
       loo_ic = c(the_model_Curve_null$loo$estimates[3,1], the_model_Curve_basic$loo$estimates[3,1], 
                  the_model_Curve_basic2$loo$estimates[3,1], the_model_C13_delux$loo$estimates[3,1]),
       loo_weights = model_weights(the_model_Curve_null, the_model_Curve_basic, 
                                   the_model_Curve_basic2, the_model_C13_delux, 
                                   criterion = 'loo') %>%
         round(3),
       R2_lwr = c(bayes_R2(the_model_Curve_null)[3],bayes_R2(the_model_Curve_basic)[3],
                  bayes_R2(the_model_Curve_basic2)[3],bayes_R2(the_model_C13_delux)[3]) %>% 
         round(3),
       R2_upr = c(bayes_R2(the_model_Curve_null)[4],bayes_R2(the_model_Curve_basic)[4],
                  bayes_R2(the_model_Curve_basic2)[4],bayes_R2(the_model_C13_delux)[4]) %>% 
         round(3))

#### Model Height Index ####
opihi_shells %>%
  ggplot(aes(x=Height_Index)) +
  geom_histogram(bins = 50)

the_model_Height_null <- opihi_shells %>%
  filter(!is.na(Height_Index)) %>%
  brm(data = ., family = gaussian,
      bf(Height_Index ~ 1
      ),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.8, max_treedepth = 10),
      sample_prior = TRUE,
      future = TRUE, seed = 12, refresh = -1,
      file = '../Results/Height_null') 
summary(the_model_Height_null)
plot(the_model_Height_null) 
pp_check(the_model_Height_null, nsamples = 100)
the_model_Height_null <- add_criterion(the_model_Height_null, "loo")
the_model_Height_null <- add_criterion(the_model_Height_null, "waic")

the_model_Height_basic <- opihi_shells %>%
  filter(!is.na(Height_Index)) %>%
  brm(data = ., family = gaussian,
      bf(Height_Index ~ latitude.c * island_group + 
           (1 | Island) + (1 | Island:SamplingSite)
      ),
      prior = c(prior(gamma(1, 0.5), class = sd)),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 20),
      sample_prior = TRUE,
      future = TRUE, seed = 12, refresh = -1,
      file = '../Results/Height_basic') 
summary(the_model_Height_basic)
plot(the_model_Height_basic) 
pp_check(the_model_Height_basic, nsamples = 100)
the_model_Height_basic <- add_criterion(the_model_Height_basic, "loo")
the_model_Height_basic <- add_criterion(the_model_Height_basic, "waic")

the_model_Height_basic2<- opihi_shells %>%
  filter(!is.na(Height_Index)) %>%
  brm(data = ., family = gaussian,
      bf(Height_Index ~ latitude.c * island_group + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite),
         sigma ~  + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite)
      ),
      prior = c(prior(gamma(1, 0.5), class = sd),
                prior(gamma(1, 0.5), class = sd, dpar = sigma),
                prior(lkj(2), class = cor)),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 20),
      sample_prior = TRUE,
      future = TRUE, seed = 12, refresh = -1,
      file = '../Results/Height_basic2')
summary(the_model_Height_basic2)
plot(the_model_Height_basic2) 
pp_check(the_model_Height_basic2, nsamples = 100)
the_model_Height_basic2 <- add_criterion(the_model_Height_basic2, "loo")
the_model_Height_basic2 <- add_criterion(the_model_Height_basic2, "waic")

the_model_Height_delux<- opihi_shells %>%
  filter(!is.na(Height_Index)) %>%
  brm(data = ., family = gaussian,
      bf(Height_Index ~ latitude.c * island_group + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite),
         sigma ~ island_group + Location.c : island_group + (1 |ID1| Island) + (1 |ID2| Island:SamplingSite)
      ),
      prior = c(prior(gamma(1, 0.5), class = sd),
                prior(gamma(1, 0.5), class = sd, dpar = sigma),
                prior(lkj(2), class = cor)),
      iter = 2000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 20),
      future = TRUE, seed = 12, refresh = -1,
      inits = 0,
      sample_prior = TRUE,
      file = '../Results/Height_delux')
summary(the_model_Height_delux)
plot(the_model_Height_delux) 
pp_check(the_model_Height_delux, nsamples = 100)
the_model_Height_delux <- add_criterion(the_model_Height_delux, "loo")
the_model_Height_delux <- add_criterion(the_model_Height_delux, "waic")


## Model Comparison ##
tibble(model = 1:4, 
       loo_ic = c(the_model_Height_null$loo$estimates[3,1], the_model_Height_basic$loo$estimates[3,1], 
                  the_model_Height_basic2$loo$estimates[3,1], the_model_Height_delux$loo$estimates[3,1]),
       loo_weights = model_weights(the_model_Height_null, the_model_Height_basic, 
                                   the_model_Height_basic2, the_model_Height_delux, 
                                   criterion = 'loo') %>%
         round(3),
       R2_lwr = c(bayes_R2(the_model_Height_null)[3],bayes_R2(the_model_Height_basic)[3],
                  bayes_R2(the_model_Height_basic2)[3],bayes_R2(the_model_Height_delux)[3]) %>% 
         round(3),
       R2_upr = c(bayes_R2(the_model_Height_null)[4],bayes_R2(the_model_Height_basic)[4],
                  bayes_R2(the_model_Height_basic2)[4],bayes_R2(the_model_Height_delux)[4]) %>% 
         round(3))

#### Fit best models with many longer chains to improve convergence of credible intervals ####
plan(multiprocess, workers = 20)
the_model_SA_mega <- update(the_model_SA_delux, iter = 11000, warmup = 1000, chains = 20,
                            future = TRUE, seed = 12, refresh = -1,
                            inits = 0,
                            sample_prior = TRUE, file = '../Results/SurfaceArea_mega')
summary(the_model_SA_mega)
rm(list = c('the_model_SA_null', 'the_model_SA_basic', 'the_model_SA_basic2', 'the_model_SA_delux'))

the_model_Curve_mega <- update(the_model_C13_delux, iter = 11000, warmup = 1000, chains = 20,
                               future = TRUE, seed = 12, refresh = -1,
                               inits = 0,
                               sample_prior = TRUE, file = '../Results/C13_mega')
summary(the_model_Curve_mega)
rm(list = c('the_model_Curve_null', 'the_model_Curve_basic', 'the_model_Curve_basic2', 'the_model_C13_delux'))

the_model_Height_mega <- update(the_model_Height_delux, iter = 11000, warmup = 1000, chains = 20,
                                future = TRUE, seed = 12, refresh = -1,
                                inits = 0,
                                sample_prior = TRUE, file = '../Results/Height_mega')
summary(the_model_Height_mega)
rm(list = c('the_model_Height_null', 'the_model_Height_basic', 'the_model_Height_basic2', 'the_model_Height_delux'))

#### Hypotheses Testing ################################################################################
## SA means#####
the_model_SA_mega$fit@sim$samples %>%
  map(~as_tibble(.x) %>%
        slice(-1:-1000)) %>%
  bind_rows() %>% 
  mutate(NWHI_fit_Slope = b_latitude.c + `b_latitude.c:island_groupNWHI`, 
         MHI_fit_Slope = b_latitude.c,
         NWHI_sigma_slope = `b_sigma_island_groupNWHI:Location.c`, 
         MHI_sigma_Slope = `b_sigma_island_groupMHI:Location.c`) %>%
  dplyr::select(NWHI_fit_Slope, MHI_fit_Slope, NWHI_sigma_slope, MHI_sigma_Slope) %>%
  gather() %>%
  separate(key, into = c('Island_Group','Parameter','Slope')) %>%
  dplyr::select(-Slope) %>%
  group_by(Island_Group, Parameter) %>%
  summarise(median_slope = -1 * median(value),
            lwr95 = -1*hdci(value, 0.95)[1],
            upr95 = -1*hdci(value, 0.95)[2],
            
            pct.less.zero = sum(value < 0)/n(),
            pct.grtr.zero = sum(value > 0)/n())

hypothesis(the_model_SA_mega, 
           c("b_latitude.c + b_latitude.c:island_groupNWHI < 0",
             "b_latitude.c < 0",
             "b_latitude.c + b_latitude.c:island_groupNWHI < b_latitude.c",
             "b_Intercept + b_island_groupNWHI > b_Intercept",
             "r_Island:SamplingSite[KAU_Miloloii,Intercept] = r_Island:SamplingSite[KAU_PMRF,Intercept]",
             "r_Island:SamplingSite[OAH_Kewalo,Intercept] = r_Island:SamplingSite[OAH_Magic_Island,Intercept]",
             "r_Island:SamplingSite__sigma[KAU_Miloloii,Intercept] = r_Island:SamplingSite__sigma[KAU_PMRF,Intercept]",
             "r_Island:SamplingSite__sigma[OAH_Kewalo,Intercept] = r_Island:SamplingSite__sigma[OAH_Magic_Island,Intercept]",
             "b_sigma_island_groupNWHI:Location.c = 0",
             "b_sigma_island_groupMHI:Location.c = 0",
             "b_sigma_island_groupNWHI:Location.c > b_sigma_island_groupMHI:Location.c",
             "b_sigma_island_groupNWHI > b_sigma_Intercept"),
           class = "",
           alpha = 0.05)

## Curve##############################################################################
the_model_Curve_mega$fit@sim$samples %>%
  map(~as_tibble(.x) %>%
        slice(-1:-1000)) %>%
  bind_rows() %>% 
  mutate(NWHI_fit_Slope = b_latitude.c + `b_latitude.c:island_groupNWHI`, 
         MHI_fit_Slope = b_latitude.c,
         NWHI_sigma_slope = `b_sigma_island_groupNWHI:Location.c`, 
         MHI_sigma_Slope = `b_sigma_island_groupMHI:Location.c`) %>%
  dplyr::select(NWHI_fit_Slope, MHI_fit_Slope, NWHI_sigma_slope, MHI_sigma_Slope) %>%
  gather() %>%
  separate(key, into = c('Island_Group','Parameter','Slope')) %>%
  dplyr::select(-Slope) %>%
  group_by(Island_Group, Parameter) %>%
  summarise(median_slope = -1 * median(value),
            lwr95 = -1*hdci(value, 0.95)[1],
            upr95 = -1*hdci(value, 0.95)[2],
            
            pct.less.zero = sum(value < 0)/n(),
            pct.grtr.zero = sum(value > 0)/n())

hypothesis(the_model_Curve_mega, 
           c("b_latitude.c + b_latitude.c:island_groupNWHI < 0",
             "b_latitude.c < 0",
             "b_latitude.c + b_latitude.c:island_groupNWHI < b_latitude.c",
             "b_Intercept + b_island_groupNWHI > b_Intercept",
             "r_Island:SamplingSite[KAU_Miloloii,Intercept] = r_Island:SamplingSite[KAU_PMRF,Intercept]",
             "r_Island:SamplingSite[OAH_Kewalo,Intercept] = r_Island:SamplingSite[OAH_Magic_Island,Intercept]",
             "r_Island:SamplingSite__sigma[KAU_Miloloii,Intercept] = r_Island:SamplingSite__sigma[KAU_PMRF,Intercept]",
             "r_Island:SamplingSite__sigma[OAH_Kewalo,Intercept] = r_Island:SamplingSite__sigma[OAH_Magic_Island,Intercept]",
             "b_sigma_island_groupNWHI:Location.c = 0",
             "b_sigma_island_groupMHI:Location.c = 0",
             "b_sigma_island_groupNWHI:Location.c > b_sigma_island_groupMHI:Location.c",
             "b_sigma_island_groupNWHI > b_sigma_Intercept"),
           class = "",
           alpha = 0.05)

## Height##########################################################
the_model_Height_mega$fit@sim$samples %>%
  map(~as_tibble(.x) %>%
        slice(-1:-1000)) %>%
  bind_rows() %>% 
  mutate(NWHI_fit_Slope = b_latitude.c + `b_latitude.c:island_groupNWHI`, 
         MHI_fit_Slope = b_latitude.c,
         NWHI_sigma_slope = `b_sigma_island_groupNWHI:Location.c`, 
         MHI_sigma_Slope = `b_sigma_island_groupMHI:Location.c`) %>%
  dplyr::select(NWHI_fit_Slope, MHI_fit_Slope, NWHI_sigma_slope, MHI_sigma_Slope) %>%
  gather() %>%
  separate(key, into = c('Island_Group','Parameter','Slope')) %>%
  dplyr::select(-Slope) %>%
  group_by(Island_Group, Parameter) %>%
  summarise(median_slope = -1 * median(value),
            lwr95 = -1*hdci(value, 0.95)[1],
            upr95 = -1*hdci(value, 0.95)[2],
            
            pct.less.zero = sum(value < 0)/n(),
            pct.grtr.zero = sum(value > 0)/n())

hypothesis(the_model_Height_mega, 
           c("b_latitude.c + b_latitude.c:island_groupNWHI < 0",
             "b_latitude.c < 0",
             "b_latitude.c + b_latitude.c:island_groupNWHI < b_latitude.c",
             "b_Intercept + b_island_groupNWHI > b_Intercept",
             "r_Island:SamplingSite[KAU_Miloloii,Intercept] = r_Island:SamplingSite[KAU_PMRF,Intercept]",
             "r_Island:SamplingSite[OAH_Kewalo,Intercept] = r_Island:SamplingSite[OAH_Magic_Island,Intercept]",
             "r_Island:SamplingSite__sigma[KAU_Miloloii,Intercept] = r_Island:SamplingSite__sigma[KAU_PMRF,Intercept]",
             "r_Island:SamplingSite__sigma[OAH_Kewalo,Intercept] = r_Island:SamplingSite__sigma[OAH_Magic_Island,Intercept]",
             "b_sigma_island_groupNWHI:Location.c = 0",
             "b_sigma_island_groupMHI:Location.c = 0",
             "b_sigma_island_groupNWHI:Location.c > b_sigma_island_groupMHI:Location.c",
             "b_sigma_island_groupNWHI > b_sigma_Intercept"),
           class = "",
           alpha = 0.05)

#### Plots ####
#### SA ####
SurfaceArea_island_grouping<-opihi_shells %>%
  group_by(Island) %>%
  summarise() %>%
  expand(Island, Island) %>%
  filter(Island < Island1) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               #add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma', re_formula = ~ (1 |ID1| Island)) %>%
               #mutate(SamplingSite = 'new') %>%
               add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(Island, .draw, fit, sigma),
             by = c('Island' = 'Island')) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               #add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma', re_formula = ~ (1 |ID1| Island)) %>%
               #mutate(SamplingSite = 'new') %>%
               add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(Island, .draw, fit, sigma) ,
             by = c('Island1' = 'Island', '.draw' = '.draw')) %>%
  mutate(difference_fit = fit.x - fit.y,
         difference_sigma = sigma.x - sigma.y) %>%
  dplyr::select(Island, Island1, difference_fit, difference_sigma) %>%
  gather(dpar, value, -Island, -Island1) %>%
  separate(dpar, into = c('remove','dpar')) %>%
  dplyr::select(-remove) %>%
  group_by(dpar, Island, Island1) %>%
  summarise(mean_diff = mean(value),
            lwr95_diff = hdci(value)[1],
            upr95_diff = hdci(value)[2]) %>%
  ungroup %>%
  group_by(dpar) %>%
  nest %>%
  mutate(groupings = map(data, cld_display_bayes)) %>%
  dplyr::select(-data) %>%
  unnest %>%
  spread(dpar, grouping) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group) %>%
               summarise(Location = mean(Location),
                         latitude=mean(latitude))) %>%
  dplyr::rename(fit_island_group = fit, sigma_island_group = sigma)

SurfaceArea_SamplingSite_grouping<-opihi_shells %>%
  group_by(SamplingSite) %>%
  summarise() %>%
  expand(SamplingSite, SamplingSite) %>%
  filter(SamplingSite < SamplingSite1) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               #add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma', re_formula = ~ (1 |ID1| Island)) %>%
               #mutate(SamplingSite = 'new') %>%
               add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(SamplingSite, .draw, fit, sigma),
             by = c('SamplingSite' = 'SamplingSite')) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               #add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma', re_formula = ~ (1 |ID1| Island)) %>%
               #mutate(SamplingSite = 'new') %>%
               add_fitted_draws(the_model_SA_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(SamplingSite, .draw, fit, sigma) ,
             by = c('SamplingSite1' = 'SamplingSite', '.draw' = '.draw')) %>%
  mutate(difference_fit = fit.x - fit.y,
         difference_sigma = sigma.x - sigma.y) %>%
  dplyr::select(SamplingSite, SamplingSite1, difference_fit, difference_sigma) %>%
  gather(dpar, value, -SamplingSite, -SamplingSite1) %>%
  separate(dpar, into = c('remove','dpar')) %>%
  dplyr::select(-remove) %>%
  group_by(dpar, SamplingSite, SamplingSite1) %>%
  summarise(mean_diff = mean(value),
            lwr95_diff = hdci(value)[1],
            upr95_diff = hdci(value)[2]) %>%
  ungroup %>%
  group_by(dpar) %>%
  nest %>%
  mutate(groupings = map(data, cld_display_bayes)) %>%
  dplyr::select(-data) %>%
  unnest %>%
  spread(dpar, grouping) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite) %>%
               summarise(Location = mean(Location),
                         latitude=mean(latitude))) %>%
  dplyr::rename(fit_island_group = fit, sigma_island_group = sigma)

SA_pt_1<-opihi_shells %>%
  group_by(island_group, latitude.c, latitude, Location, Location.c) %>%
  summarise() %>%
  group_by(island_group) %>%
  nest %>%
  mutate(addition_lat = map_dbl(data, ~.x %>% mutate(difference = latitude-latitude.c) %$% unique(difference)),
         addition_loc = map_dbl(data, ~.x %>% mutate(difference = Location-Location.c) %$% unique(difference)),
         data = map(data, ~.x %$%
                      tibble(latitude.c = seq(min(latitude.c)-0.1, max(latitude.c)+0.1, length.out = 100),
                             Location.c = seq(min(Location.c)-30, max(Location.c)+30,length.out = 100))
                    
                    
         )) %>%
  unnest %>%
  add_fitted_draws(the_model_SA_mega, re_formula = ~(1|island_group), dpar = 'sigma') %>%
  mutate(latitude = latitude.c + addition_lat,
         Location = Location.c + addition_loc) %>%
  group_by(island_group, latitude, Location) %>%
  dplyr::select(island_group, latitude, Location, .value, sigma) %>%
  point_interval(.width = c(0.68, 0.95), .interval = hdci)


SA_pt_2<-opihi_shells %>%
  mutate(Island = factor(Island, levels = c("GP", "LPP", "MMM", "NIH", "KAU", "OAH", "MAUI", "BI"))) %>%
  group_by(Island, island_group, SamplingSite, latitude.c, latitude, Location, Location.c) %>%
  summarise() %>%
  add_fitted_draws(the_model_SA_mega, dpar = 'sigma') %>%
  group_by(island_group, Island, SamplingSite, latitude, Location) %>%
  dplyr::select(island_group, Island, SamplingSite, latitude, Location, .value, sigma) %>%
  point_interval(.width = c(0.68, 0.95), .interval = hdci)


SurfaceArea_island_grouping <- SurfaceArea_island_grouping  %>%
  inner_join(SA_pt_2 %>%
               group_by(Island) %>%
               summarise(.value = mean(.value), sigma = mean(sigma)))

SurfaceArea_SamplingSite_grouping <- SurfaceArea_SamplingSite_grouping %>%
  inner_join(SA_pt_2 %>%
               group_by(.value, sigma, Island, SamplingSite) %>%
               summarise)

SurfaceArea_SamplingSite_grouping$fit_island_group_ceb <- c("de","ef","fg","b","g","efg","","c","a","d") #CEB hack: put labels in order, a=tallest
SurfaceArea_SamplingSite_grouping$sigma_island_group_ceb <- c("ab","bcd","de","bcd","e","bcde","de","a","bc","cd")

fit_text_x <- SurfaceArea_SamplingSite_grouping$latitude - 0.2 ; fit_text_x[10] <- fit_text_x[10]+0.25 ; fit_text_x[6] <- fit_text_x[6] - 0.05 #; text_x[6] <- text_x[6]-27
fit_text_y <- SurfaceArea_SamplingSite_grouping$.value + 0.2 ; fit_text_y[1] <- fit_text_y[1]+.2; fit_text_y[10] <- fit_text_y[10]+.3; fit_text_y[6] <- fit_text_y[6]-0.25; fit_text_y[5] <- fit_text_y[5]-0.1
SA_mean <- SA_pt_1 %>% 
  ggplot(aes(x = latitude, y = .value, group = interaction(island_group,.width), 
             ymin = .value.lower, ymax = .value.upper)) +
  geom_line(size=1) +
  geom_ribbon(alpha = 0.2, fill = 'black') +
  geom_point(data = SA_pt_2, size = 4) +
  geom_linerange(data = SA_pt_2 %>% filter(.width == 0.95), size = 0.5, linetype = 1) +
  geom_linerange(data = SA_pt_2 %>% filter(.width == 0.68), size = 2) +
  # geom_text(data = SurfaceArea_island_grouping, 
  #           aes(x = latitude, label = fit_island_group, y = .value + 2), 
  #           size= 7, show.legend = FALSE, inherit.aes = FALSE) +
  geom_text(data = SurfaceArea_SamplingSite_grouping, 
            aes(x = fit_text_x, label = fit_island_group_ceb, y = fit_text_y), 
            size= 6, show.legend = FALSE, inherit.aes = FALSE) +
  scale_x_reverse(breaks = c(25,24,23,22,21,20)) +
  xlab("Latitude ~ Degrees") +
  ylab("Surface Area ~ cm^2") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 23))
#ggsave(SA_mean, 'SA_mean.png',width = 7, height = 7)

sigma_text_x <- SurfaceArea_SamplingSite_grouping$Location - 40; sigma_text_x[2] <- sigma_text_x[2]+95; sigma_text_x[4] <- sigma_text_x[4]-20; sigma_text_x[6] <- sigma_text_x[6]+40
sigma_text_y <- SurfaceArea_SamplingSite_grouping$sigma - 0.04; sigma_text_y[2] <- sigma_text_y[2]+.1; sigma_text_y[4] <- sigma_text_y[4]+.1 ; sigma_text_y[6] <- sigma_text_y[6]+.4 ; sigma_text_y[1] <- sigma_text_y[1]+.1
SA_sigma <- SA_pt_1 %>% 
  ggplot(aes(x = Location, y = sigma, group = interaction(island_group,.width), 
             ymin = sigma.lower, ymax = sigma.upper)) +
  geom_line(size=1) +
  geom_ribbon(alpha = 0.2, fill = 'black') +
  geom_point(data = SA_pt_2, size = 4) +
  geom_linerange(data = SA_pt_2 %>% filter(.width == 0.95), size = 0.5, linetype = 1) +
  geom_linerange(data = SA_pt_2 %>% filter(.width == 0.68), size = 2) +
  # geom_text(data = SurfaceArea_island_grouping, 
  #           aes(x = Location, label = sigma_island_group, y = sigma + 0.5), 
  #           size= 7, show.legend = FALSE, inherit.aes = FALSE) +
  geom_text(data = SurfaceArea_SamplingSite_grouping, 
            #aes(x = Location+30, label = sigma_island_group_ceb, y = sigma + 0.5), 
            aes(x = sigma_text_x, label = sigma_island_group_ceb, y = sigma_text_y), 
            size= 6, show.legend = FALSE, inherit.aes = FALSE) +
  xlab("One Dimensional Location ~ km") +
  ylab("Variation in SurfaceArea ~ cm^2") +
  scale_x_continuous(breaks = c(0,300,600,900,1200)) +
  scale_y_continuous(breaks = c(0,1,2,3)) +
  theme_classic() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 23))
#ggsave(SA_sigma, 'SA_sigma.png',width = 7, height = 7)


#### Curve #######################################################################################
C13_island_grouping<-opihi_shells %>%
  group_by(Island) %>%
  summarise() %>%
  expand(Island, Island) %>%
  filter(Island < Island1) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Curve_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(Island, .draw, fit, sigma) ,
             by = c('Island' = 'Island')) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Curve_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(Island, .draw, fit, sigma) ,
             by = c('Island1' = 'Island', '.draw' = '.draw')) %>%
  mutate(difference_fit = fit.x - fit.y,
         difference_sigma = sigma.x - sigma.y) %>%
  dplyr::select(Island, Island1, difference_fit, difference_sigma) %>%
  gather(dpar, value, -Island, -Island1) %>%
  separate(dpar, into = c('remove','dpar')) %>%
  dplyr::select(-remove) %>%
  group_by(dpar, Island, Island1) %>%
  summarise(mean_diff = mean(value),
            lwr95_diff = hdci(value)[1],
            upr95_diff = hdci(value)[2]) %>%
  ungroup %>%
  group_by(dpar) %>%
  nest %>%
  mutate(groupings = map(data, cld_display_bayes)) %>%
  dplyr::select(-data) %>%
  unnest %>%
  spread(dpar, grouping) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group) %>%
               summarise(Location = mean(Location),
                         latitude=mean(latitude))) %>%
  dplyr::rename(fit_island_group = fit, sigma_island_group = sigma)

C13_SamplingSite_grouping<-opihi_shells %>%
  group_by(SamplingSite) %>%
  summarise() %>%
  expand(SamplingSite, SamplingSite) %>%
  filter(SamplingSite < SamplingSite1) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Curve_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(SamplingSite, .draw, fit, sigma),
             by = c('SamplingSite' = 'SamplingSite')) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Curve_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(SamplingSite, .draw, fit, sigma) ,
             by = c('SamplingSite1' = 'SamplingSite', '.draw' = '.draw')) %>%
  mutate(difference_fit = fit.x - fit.y,
         difference_sigma = sigma.x - sigma.y) %>%
  dplyr::select(SamplingSite, SamplingSite1, difference_fit, difference_sigma) %>%
  gather(dpar, value, -SamplingSite, -SamplingSite1) %>%
  separate(dpar, into = c('remove','dpar')) %>%
  dplyr::select(-remove) %>%
  group_by(dpar, SamplingSite, SamplingSite1) %>%
  summarise(mean_diff = mean(value),
            lwr95_diff = hdci(value)[1],
            upr95_diff = hdci(value)[2]) %>%
  ungroup %>%
  group_by(dpar) %>%
  nest %>%
  mutate(groupings = map(data, cld_display_bayes)) %>%
  dplyr::select(-data) %>%
  unnest %>%
  spread(dpar, grouping) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite) %>%
               summarise(Location = mean(Location),
                         latitude=mean(latitude))) %>%
  dplyr::rename(fit_island_group = fit, sigma_island_group = sigma)


C13_pt_1<-opihi_shells %>%
  group_by(island_group, latitude.c, latitude, Location, Location.c) %>%
  summarise() %>%
  group_by(island_group) %>%
  nest %>%
  mutate(addition_lat = map_dbl(data, ~.x %>% mutate(difference = latitude-latitude.c) %$% unique(difference)),
         addition_loc = map_dbl(data, ~.x %>% mutate(difference = Location-Location.c) %$% unique(difference)),
         data = map(data, ~.x %$%
                      tibble(latitude.c = seq(min(latitude.c)-0.1, max(latitude.c)+0.1, length.out = 100),
                             Location.c = seq(min(Location.c)-30, max(Location.c)+30,length.out = 100))
                    
                    
         )) %>%
  unnest %>%
  add_fitted_draws(the_model_Curve_mega, re_formula = NA, dpar = 'sigma') %>%
  mutate(latitude = latitude.c + addition_lat,
         Location = Location.c + addition_loc) %>%
  group_by(island_group, latitude, Location) %>%
  dplyr::select(island_group, latitude, Location, .value, sigma) %>%
  point_interval(.width = c(0.68, 0.95), .interval = hdci)


C13_pt_2<-opihi_shells %>%
  mutate(Island = factor(Island, levels = c("GP", "LPP", "MMM", "NIH", "KAU", "OAH", "MAUI", "BI"))) %>%
  group_by(Island, island_group, SamplingSite, latitude.c, latitude, Location, Location.c) %>%
  summarise() %>%
  add_fitted_draws(the_model_Curve_mega, dpar = 'sigma') %>%
  group_by(island_group, Island, SamplingSite, latitude, Location) %>%
  dplyr::select(island_group, Island, SamplingSite, latitude, Location, .value, sigma) %>%
  point_interval(.width = c(0.68, 0.95), .interval = hdci)


C13_island_grouping <- C13_island_grouping  %>%
  inner_join(C13_pt_2 %>%
               group_by(Island) %>%
               summarise(.value = mean(.value), sigma = mean(sigma)))

C13_SamplingSite_grouping <- C13_SamplingSite_grouping %>%
  inner_join(C13_pt_2 %>%
               group_by(.value, sigma, Island, SamplingSite) %>%
               summarise)
#C13_island_grouping$fit_island_group_ceb <- c("d", "c", "c", "b", "d", "b", "a" ,"d") #CEB hack: put labels in order, a=tallest
C13_island_grouping$sigma_island_group_ceb <- c("de",  "c",   "c",   "cd",  "e",   "b",   "a",   "cde")
C13_mean <- C13_pt_1 %>% 
  ggplot(aes(x = latitude, y = .value, group = interaction(island_group,.width), 
             ymin = .value.lower, ymax = .value.upper)) +
  geom_line(size=1) +
  geom_ribbon(alpha = 0.2, fill = 'black') +
  geom_point(data = C13_pt_2, size = 4) +
  geom_linerange(data = C13_pt_2 %>% filter(.width == 0.95), size = 0.5, linetype = 1) +
  geom_linerange(data = C13_pt_2 %>% filter(.width == 0.68), size = 2) +
  geom_text(data = C13_island_grouping, 
            aes(x = latitude, label = fit_island_group, y = .value + 0.01), 
            size= 6, show.legend = FALSE, inherit.aes = FALSE) +
  # geom_text(data = C13_SamplingSite_grouping, 
  #           aes(x = latitude, label = fit_island_group, y = .value - 0.02), 
  #           size= 7, show.legend = FALSE, inherit.aes = FALSE) +
  scale_x_reverse() +
  scale_y_continuous(breaks = c(1.00,1.02,1.04,1.06,1.08)) +
  xlab("Latitude ~ Degrees") +
  ylab("Doming Index ~ C13/L13") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 23))
#ggsave('C13_mean.png',width = 7, height = 7)

C13_sigma <- C13_pt_1 %>% 
  ggplot(aes(x = Location, y = sigma, group = interaction(island_group,.width), 
             ymin = sigma.lower, ymax = sigma.upper)) +
  geom_line(size=1) +
  geom_ribbon(alpha = 0.2, fill = 'black') +
  geom_point(data = C13_pt_2, size = 4) +
  geom_linerange(data = C13_pt_2 %>% filter(.width == 0.95), size = 0.5, linetype = 1) +
  geom_linerange(data = C13_pt_2 %>% filter(.width == 0.68), size = 2) +
  geom_text(data = C13_island_grouping, 
            aes(x = Location, label = sigma_island_group_ceb, y = sigma + 0.005), 
            size= 6, show.legend = FALSE, inherit.aes = FALSE) +
  # geom_text(data = C13_SamplingSite_grouping, 
  #           aes(x = Location, label = sigma_island_group, y = sigma - 0.01), 
  #           size= 7, show.legend = FALSE, inherit.aes = FALSE) +
  xlab("One Dimensional Location ~ km") +
  ylab("Variation in Doming Index") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 23))
#ggsave('C13_sigma.png',width = 7, height = 7)

#### Height #####################################################################################
Height_island_grouping<-opihi_shells %>%
  group_by(Island) %>%
  summarise() %>%
  expand(Island, Island) %>%
  filter(Island < Island1) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Height_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(Island, .draw, fit, sigma) ,
             by = c('Island' = 'Island')) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Height_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(Island, .draw, fit, sigma) ,
             by = c('Island1' = 'Island', '.draw' = '.draw')) %>%
  mutate(difference_fit = fit.x - fit.y,
         difference_sigma = sigma.x - sigma.y) %>%
  dplyr::select(Island, Island1, difference_fit, difference_sigma) %>%
  gather(dpar, value, -Island, -Island1) %>%
  separate(dpar, into = c('remove','dpar')) %>%
  dplyr::select(-remove) %>%
  group_by(dpar, Island, Island1) %>%
  summarise(mean_diff = mean(value),
            lwr95_diff = hdci(value)[1],
            upr95_diff = hdci(value)[2]) %>%
  ungroup %>%
  group_by(dpar) %>%
  nest %>%
  mutate(groupings = map(data, cld_display_bayes)) %>%
  dplyr::select(-data) %>%
  unnest %>%
  spread(dpar, grouping) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group) %>%
               summarise(Location = mean(Location),
                         latitude=mean(latitude))) %>%
  dplyr::rename(fit_island_group = fit, sigma_island_group = sigma)

Height_SamplingSite_grouping<-opihi_shells %>%
  group_by(SamplingSite) %>%
  summarise() %>%
  expand(SamplingSite, SamplingSite) %>%
  filter(SamplingSite < SamplingSite1) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Height_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(SamplingSite, .draw, fit, sigma),
             by = c('SamplingSite' = 'SamplingSite')) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite, latitude.c, Location.c) %>%
               summarise() %>%
               add_fitted_draws(the_model_Height_mega, value = 'fit', dpar = 'sigma') %>%
               ungroup %>%
               dplyr::select(SamplingSite, .draw, fit, sigma) ,
             by = c('SamplingSite1' = 'SamplingSite', '.draw' = '.draw')) %>%
  mutate(difference_fit = fit.x - fit.y,
         difference_sigma = sigma.x - sigma.y) %>%
  dplyr::select(SamplingSite, SamplingSite1, difference_fit, difference_sigma) %>%
  gather(dpar, value, -SamplingSite, -SamplingSite1) %>%
  separate(dpar, into = c('remove','dpar')) %>%
  dplyr::select(-remove) %>%
  group_by(dpar, SamplingSite, SamplingSite1) %>%
  summarise(mean_diff = mean(value),
            lwr95_diff = hdci(value)[1],
            upr95_diff = hdci(value)[2]) %>%
  ungroup %>%
  group_by(dpar) %>%
  nest %>%
  mutate(groupings = map(data, cld_display_bayes)) %>%
  dplyr::select(-data) %>%
  unnest %>%
  spread(dpar, grouping) %>%
  inner_join(opihi_shells %>%
               group_by(Island, island_group, SamplingSite) %>%
               summarise(Location = mean(Location),
                         latitude=mean(latitude))) %>%
  dplyr::rename(fit_island_group = fit, sigma_island_group = sigma)


Height_pt_1<-opihi_shells %>%
  group_by(island_group, latitude.c, latitude, Location, Location.c) %>%
  summarise() %>%
  group_by(island_group) %>%
  nest %>%
  mutate(addition_lat = map_dbl(data, ~.x %>% mutate(difference = latitude-latitude.c) %$% unique(difference)),
         addition_loc = map_dbl(data, ~.x %>% mutate(difference = Location-Location.c) %$% unique(difference)),
         data = map(data, ~.x %$%
                      tibble(latitude.c = seq(min(latitude.c)-0.1, max(latitude.c)+0.1, length.out = 100),
                             Location.c = seq(min(Location.c)-30, max(Location.c)+30,length.out = 100))
                    
                    
         )) %>%
  unnest %>%
  add_fitted_draws(the_model_Height_mega, re_formula = NA, dpar = 'sigma') %>%
  mutate(latitude = latitude.c + addition_lat,
         Location = Location.c + addition_loc) %>%
  group_by(island_group, latitude, Location) %>%
  dplyr::select(island_group, latitude, Location, .value, sigma) %>%
  point_interval(.width = c(0.68, 0.95), .interval = hdci)


Height_pt_2<-opihi_shells %>%
  mutate(Island = factor(Island, levels = c("GP", "LPP", "MMM", "NIH", "KAU", "OAH", "MAUI", "BI"))) %>%
  group_by(Island, island_group, SamplingSite, latitude.c, latitude, Location, Location.c) %>%
  summarise() %>%
  add_fitted_draws(the_model_Height_mega, dpar = 'sigma') %>%
  group_by(island_group, Island, SamplingSite, latitude, Location) %>%
  dplyr::select(island_group, Island, SamplingSite, latitude, Location, .value, sigma) %>%
  point_interval(.width = c(0.68, 0.95), .interval = hdci)


Height_island_grouping <- Height_island_grouping  %>%
  inner_join(Height_pt_2 %>%
               group_by(Island) %>%
               summarise(.value = mean(.value), sigma = mean(sigma)))

Height_SamplingSite_grouping <- Height_SamplingSite_grouping %>%
  inner_join(Height_pt_2 %>%
               group_by(.value, sigma, Island, SamplingSite) %>%
               summarise)
Height_island_grouping$fit_island_group_ceb <- c("e","c","d","b","e","b","a","e") #CEB hack: put labels in order, a=tallest
Height_SamplingSite_grouping$sigma_site_group_ceb <- c("ab", "e",  "",  "bc", "de", "de", "de", "a",  "bc", "cd")

Height_mean <- Height_pt_1 %>% 
  ggplot(aes(x = latitude, y = .value, group = interaction(island_group,.width), 
             ymin = .value.lower, ymax = .value.upper)) +
  geom_line(size=1) +
  geom_ribbon(alpha = 0.2, fill = 'black') +
  geom_point(data = Height_pt_2, size = 4) +
  geom_linerange(data = Height_pt_2 %>% filter(.width == 0.95), size = .5, linetype = 1) +
  geom_linerange(data = Height_pt_2 %>% filter(.width == 0.68), size = 2) +
  geom_text(data = Height_island_grouping, 
            aes(x = latitude, label = fit_island_group_ceb, y = .value + 0.027), 
            size= 6, show.legend = FALSE, inherit.aes = FALSE) +
  # geom_text(data = Height_SamplingSite_grouping,
  #           aes(x = latitude, label = fit_island_group, y = .value - 0.05),
  #           size= 7, show.legend = FALSE, inherit.aes = FALSE) +
  scale_x_reverse() +
  scale_y_continuous(breaks = c(0.3,0.4,0.5)) +
  xlab("Latitude ~ Degrees") +
  ylab("Height Index ~ H/L") +
  theme_classic() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 23))
#ggsave('Height_mean.png',width = 7, height = 7)

text_x <- Height_SamplingSite_grouping$Location - 40   ; text_x[6] <- text_x[6]+80;     text_x[2] <- text_x[2]+70
text_y <- Height_SamplingSite_grouping$sigma + 0.0025  ; text_y[6] <- text_y[6]+0.001 ; text_y[2] <- text_y[2]-0.0045 ; text_y[5] <- text_y[5]-.005 
Height_sigma <- Height_pt_1 %>% 
  ggplot(aes(x = Location, y = sigma, group = interaction(island_group,.width), 
             ymin = sigma.lower, ymax = sigma.upper)) +
  geom_line(size=1) +
  geom_ribbon(alpha = 0.2, fill = 'black') +
  geom_point(data = Height_pt_2, size = 4) +
  geom_linerange(data = Height_pt_2 %>% filter(.width == 0.95), size = 0.5, linetype = 1) +
  geom_linerange(data = Height_pt_2 %>% filter(.width == 0.68), size = 2) +
  # geom_text(data = Height_island_grouping,
  #           aes(x = Location, label = sigma_island_group, y = sigma + 0.02),
  #           size= 7, show.legend = FALSE, inherit.aes = FALSE) +
  geom_text(data = Height_SamplingSite_grouping,
            aes(x = text_x, label = sigma_site_group_ceb, y = text_y ),
            size= 6, show.legend = FALSE, inherit.aes = FALSE) +
  xlab("One Dimensional Location ~ km") +
  ylab("Variation in Height Index") +
  scale_x_continuous(breaks=c(0,300,600,900,1200)) +
  theme_classic() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 23))
#ggsave('Height_sigma.png',width = 7, height = 7)




#### Plots for output ####
library(gridExtra)

grid.arrange(SA_mean, C13_mean, Height_mean, ncol=2)
ggsave('mean_plots.png',width = 7, height = 7)

grid.arrange(SA_sigma, C13_sigma, Height_sigma, ncol=2)
ggsave('sigma_plots.png',width = 7, height = 7)
