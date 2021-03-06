# # Persistent Effect of Temperature on GDP Identified from Lower Frequency Temperature Variability
# # Bastien-Olvera, B. A.; Granella, F. and Moore, F. C.

# ## Index
#     # 1. Setup
#     # 2. Simulation (Figure 1 and 2)
#         # 2.1. Temperature fluctuations and filters - Figure 1
#         # 2.2. Global temperature profile
#         # 2.3. Perform MAIN simulation - Figure 2
#         # 2.4. Additional simulations used for Response to Reviewers
#     # 3. Empirical analysis (Figure 3 and Table 1)
#         # 3.1 Country-level regressions
#         # 3.2. Panel Regression
#         # 3.2. Categorizing - Figure 3
#         # 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2
#         # 3.4. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2
#         # 3.5. Other datasets - Table 1, columns 3 to 6; Supp Fig 3

## 1. Setup (start)
    #set working directory
    dir <- "C:/Users/bastien/Documents/GitHub/TempEffectGDP/"
    setwd(dir)
    x<-c("sandwich","TTR", "ggpubr","mFilter", "ncdf4","TSA","tidyverse","ggplot2",
        "ggpubr","dplR","reshape2","raster","dplyr",
        "RColorBrewer","colorspace","spData","sf","countrycode", "ggstatsplot",
        "ggsignif","dlnm","lfe","directlabels","splines","timeSeries",
        "MASS","stargazer","jtools","wbstats","dplR","lmtest","sandwich","dlnm","splines")
    lapply(x, require, character.only = TRUE)
    # some functions and color pallettes
    tsyears <- function(ts) as.numeric(trunc(time(ts)))
    cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
    cols=rev(c("#76d3ae","#0dcfca","#055692"))

    #Load data
    barro <- read.csv(paste(dir,"Data/Barro_UDel.csv",sep="")) #Barro-Ursua and University of Delaware
    mad <- read.csv(paste(dir,"Data/Maddison_UDel.csv",sep="")) #Maddison Project and University of Delaware
    wb <- read.csv(paste(dir,"Data/WB_UDel.csv",sep=""))   #World Bank and University of Delaware
    bhm <- read.csv("C:/Users/bastien/Box/Long Run GDP Growth/Data/BHM.csv") #BHM dataset (World Bank and University of Delaware)
    
## 1. Setup (end)

## 2. Simulation - Figures 1 and 2 (start)

    ## 2.1. Temperature fluctuations and filters - Figure 1 (start)
        
        wbUS <- wb[wb$countrycode=="USA",]
        mt <- lm(UDel_pop_temp~year, data = wbUS)
        t <- resid(mt)    
        t <- interpNA(t, method = "linear")
        t <- removeNA(t)
        tdf <- data.frame(year = mt$model[,2], temp= t+1.8, Filter = "Unfiltered", mean = 1.8)
        t3 <- pass.filt(t, W=3, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t3+1, Filter = "2-3 years", mean = 1))
        t5 <- pass.filt(t, W=5, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t5+0.5, Filter = "2-5 years", mean = 0.5))
        t10 <- pass.filt(t, W=10, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t10+0.25, Filter = "2-10 years", mean = 0.25))
        t15 <- pass.filt(t, W=15, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t15+0.1, Filter = "2-15 years", mean = 0.1))
        
        tdf$Filter <- factor(tdf$Filter, levels = c("Unfiltered", "2-3 years", "2-5 years","2-10 years","2-15 years"))
        

        
        ggplot(data=tdf[which(tdf$Filter=="Unfiltered"),], aes(x=year, y=temp, color = Filter))+geom_line()+theme_bw()+
            theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
            geom_errorbar(aes(x=1958,ymin=1,ymax=2), color = "black")+
            geom_text(aes(x=1960,y=1.5),label="1C", color="black")+
            #geom_dl(aes(label = Filter), method = list(dl.combine("last.points")), cex = 0.9) +
            xlim(1958,2018)+ ggtitle("United States temperature fluctuation")+ labs(color= "Filtered oscillations \n (periods)")+
            geom_line(data = tdf[which(tdf$Filter=="Unfiltered"),], aes(x=year, y =mean, color = Filter))
            ggsave("UStemp.png",dpi=300)   

            ggplot(data=tdf, aes(x=year, y=temp, color = Filter))+geom_line()+theme_bw()+
            theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
            geom_errorbar(aes(x=1958,ymin=1,ymax=2), color = "black")+
            geom_text(aes(x=1960,y=1.5),label="1C", color="black")+
            #geom_dl(aes(label = Filter), method = list(dl.combine("last.points")), cex = 0.9) +
            xlim(1958,2018)+ ggtitle("United States temperature fluctuation")+ labs(color= "Filtered oscillations \n (periods)")+
            geom_line(data = tdf, aes(x=year, y =mean, color = Filter))
            ggsave("UStemp_full.png",dpi=300)
    ## 2.1. Temperature fluctuations and filters - Figure 1 (end)

    ## 2.1. Temperature fluctuations and filters no demean (start)
        
        wbUS <- wb[wb$countrycode=="USA",]
        mt <- lm(UDel_pop_temp~year, data = wbUS)
        t <- wbUS$UDel_pop_temp  
        t <- interpNA(t, method = "linear")
        t <- removeNA(t)
        tdf <- data.frame(year = mt$model[,2], temp= t+1.8, Filter = "Unfiltered", mean = 1.8)
        t3 <- pass.filt(t, W=3, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t3+1, Filter = "2-3 years", mean = 1))
        t5 <- pass.filt(t, W=5, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t5+0.5, Filter = "2-5 years", mean = 0.5))
        t10 <- pass.filt(t, W=10, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t10+0.25, Filter = "2-10 years", mean = 0.25))
        t15 <- pass.filt(t, W=15, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t15+0.1, Filter = "2-15 years", mean = 0.1))
        
        tdf$Filter <- factor(tdf$Filter, levels = c("Unfiltered", "2-3 years", "2-5 years","2-10 years","2-15 years"))
        

        
        ggplot(data=tdf[which(tdf$Filter=="Unfiltered"),], aes(x=year, y=temp, color = Filter))+geom_line()+theme_bw()+
            theme(legend.position="none",axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
            geom_errorbar(aes(x=1958,ymin=15,ymax=16), color = "black")+
            geom_text(aes(x=1960,y=15.5),label="1C", color="black")+
            #geom_dl(aes(label = Filter), method = list(dl.combine("last.points")), cex = 0.9) +
            xlim(1958,2018)+ ggtitle("United States temperature fluctuation")+ labs(color= "Filtered oscillations \n (periods)")
            #geom_line(data = tdf[which(tdf$Filter=="Unfiltered"),], aes(x=year, y =mean, color = Filter))
            ggsave("UStemp_notrend.png",dpi=300)   

            ggplot(data=tdf, aes(x=year, y=temp, color = Filter))+geom_line()+theme_bw()+
            theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
            geom_errorbar(aes(x=1958,ymin=15,ymax=16), color = "black")+
            geom_text(aes(x=1960,y=15.5),label="1C", color="black")+
            #geom_dl(aes(label = Filter), method = list(dl.combine("last.points")), cex = 0.9) +
            xlim(1958,2018)+ ggtitle("United States temperature fluctuation")+ labs(color= "Filtered oscillations \n (periods)")
            #geom_line(data = tdf, aes(x=year, y =mean, color = Filter))
            ggsave("UStemp_full_nodetrend.png",dpi=300)
    ## 2.1. Temperature fluctuations and filters no demean (end)

    ## 2.1. Temperature fluctuations and filters no demean (start)
        
        wbUS <- wb[wb$countrycode=="USA",]
        mt <- lm(UDel_pop_temp~year, data = wbUS)
        t <- wbUS$gdppc  
        t <- interpNA(t, method = "linear")
        t <- removeNA(t)
        tdf <- data.frame(year = mt$model[,2], temp= t+1.8, Filter = "Unfiltered", mean = 1.8)
        t3 <- pass.filt(t, W=3, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t3+1, Filter = "2-3 years", mean = 1))
        t5 <- pass.filt(t, W=5, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t5+0.5, Filter = "2-5 years", mean = 0.5))
        t10 <- pass.filt(t, W=10, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t10+0.25, Filter = "2-10 years", mean = 0.25))
        t15 <- pass.filt(t, W=15, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t15+0.1, Filter = "2-15 years", mean = 0.1))
        
        tdf$Filter <- factor(tdf$Filter, levels = c("Unfiltered", "2-3 years", "2-5 years","2-10 years","2-15 years"))
        

        
        ggplot(data=tdf[which(tdf$Filter=="Unfiltered"),], aes(x=year, y=temp, color = Filter))+geom_line()+theme_bw()+
            theme(legend.position="none",axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
            #geom_errorbar(aes(x=1958,ymin=15,ymax=16), color = "black")+
            #geom_text(aes(x=1960,y=15.5),label="1C", color="black")+
            #geom_dl(aes(label = Filter), method = list(dl.combine("last.points")), cex = 0.9) +
            xlim(1958,2018)+ ggtitle("United States GDP")+ labs(color= "")
            #geom_line(data = tdf[which(tdf$Filter=="Unfiltered"),], aes(x=year, y =mean, color = Filter))
            ggsave("USgdp_notrend.png",dpi=300)   

            ggplot(data=tdf, aes(x=year, y=temp, color = Filter))+geom_line()+theme_bw()+
            theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
            geom_errorbar(aes(x=1958,ymin=15,ymax=16), color = "black")+
            geom_text(aes(x=1960,y=15.5),label="1C", color="black")+
            #geom_dl(aes(label = Filter), method = list(dl.combine("last.points")), cex = 0.9) +
            xlim(1958,2018)+ ggtitle("United States temperature fluctuation")+ labs(color= "Filtered oscillations \n (periods)")
            #geom_line(data = tdf, aes(x=year, y =mean, color = Filter))
            ggsave("UStemp_full_nodetrend.png",dpi=300)
    ## 2.1. Temperature fluctuations and filters no demean (end)


    ## 2.2. Global temperature profile (start)
        
        #Read and save global temperature file. NetCDF too large to be in GitHub

        #file=nc_open(paste(dir,"Data/gmt_MCruns_ensemble_full.nc",sep=""))
        #gtemp=ncvar_get(file,"gmt")
        #gtemp=apply(gtemp,MARGIN=3,FUN=mean)
        #saveRDS(gtemp, file = "Data/gtemp.rds")
        
        # Uncomment for US temp profile (end)
            file=nc_open("C:/Users/bastien/Box/Long Run GDP Growth/Data/LMR Data/air_MCruns_ensemble_mean_LMRv2.0.nc")
            gtemp=ncvar_get(file,"air",start=c(1,1,1,1),count=c(-1,-1,1,-1))
            lon <- ncvar_get(file,"lon")
            lat <- ncvar_get(file,"lat")
            time <- ncvar_get(file,"time")
            tunits <- ncatt_get(file,"time","units")
            data(world)
            geom_iso <- world$geom[world$iso_a2=="US"]
            geom_iso <- st_shift_longitude(geom_iso)
            geom_iso <- st_cast(geom_iso, "POLYGON")
            geom_iso <-as_Spatial(geom_iso)
            
            #gtemp <- brick(gtemp, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
            #ustemp_lmr <- rep(NA,1500)
            #for (i in 1:1500){
            #    gtemp_year <- raster::subset(gtemp,i)
            #   #gtemp_year <- rotate(gtemp_year)
            #  gtemp_year <- t(gtemp_year)
            # ustemp <- extract(gtemp_year, geom_iso,method='simple')
                #ustemp_lmr[i] <- mean(ustemp[[1]])

            #}
            
            #print(file)
            #ustemp_lmr_orig <- ustemp_lmr
            #ustemp_lmr=lm(ustemp_lmr_orig~x)$residuals
            #p=periodogram(ustemp_lmr)
        
        # Uncomment for US temp profile (end)
        
        # Load gtemp
        gtemp <- readRDS(file = "Data/gtemp.rds")
        
        #take first 1500 years, prior to anthropogenic influence
        gtemp=gtemp[1:1500]

        #take out linear time trend
        x=1:1500
        gtemp=lm(gtemp~x)$residuals
        p=periodogram(gtemp)
        dd=data.frame(freq=p$freq,spec=p$spec)
        dd=dd[order(-dd$spec),]

        #look at top 50 frequencies
        dd=dd[1:50,]
        dd$period=1/dd$freq

        #create randomized time series 
        randomts=function(timeseries){
            #returns a time series with the same spectral profile as the argument, but with randomly chosen phases
            ft=fft(timeseries)
            N=length(timeseries)
            rphase=runif(N,min=0,max=360)
            newft=complex(real=abs(ft)*cos(rphase),imaginary=abs(ft)*sin(rphase))
            return(fft(newft,inverse=T)/N)
        }
    ## 2.2. Global temperature profile (end)

    ## 2.3. Perform MAIN simulation (start)


            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef=-0.05 #effect size - change in growth per degree warming()
            coef2=-0.05
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,3,5,10,15)
            nsims=1000
            sims=array(dim=c(nsims,length(periods),2))
            sims_adjusted=array(dim=c(nsims,length(periods),2))
            sims_ar=array(dim=c(nsims,length(periods),2))
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                for(j in 2:time){
                #linear
                randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                }
                dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                for(k in 1:length(periods)){
                mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                if (k==1){
                        sims[i,k,]=c(mod_g,mod_l)
                        sims_adjusted[i,k,]=c(mod_g,mod_l)
                        
                    }else{
                    x <- seq(1:length(randomtemp))
                    randomtemp <- lm(randomtemp~x)$residuals
                    tempts <- pass.filt(randomtemp, W=periods[k], type="low", method="Butterworth")
                    
                    ratio <- median(randomtemp/tempts)
                    growth <- randomgrowth_g
                    #growth <- pass.filt(growth, W=periods[k], type="low", method="Butterworth")
                    level <- randomgrowth_l
                    level <- pass.filt(level, W=periods[k], type="low", method="Butterworth")
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        level = unclass(level))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]
                    dw <- dwtest(lm(growth~temp,data=filt)) #test for autocorrelation
                    sims[i,k,]=c(mod_gfilt,mod_lfilt)
                    sims_adjusted[i,k,]=c(mod_gfilt/ratio,mod_lfilt/ratio)
                    sims_ar[i,k,]=c(dw$statistic[1],dw$p.value[1])
                }
                

                }
                if(i%%50==0) print(i)
            }
            #
          
                # Test - Figure 2
                    
                    simsfiltmean=apply(sims_ar,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_ar,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("frequencies","ImpactType","dw_estimate","dw_pvalue")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "10 years", "20 years", "50 years", "100 years"))
                    simsfiltdat_ar <- simsfiltdat

                    
                    simsfiltmean=apply(sims_adjusted,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_adjusted,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    simsfiltdat_adjusted <- simsfiltdat

                    simsfiltmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    #simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    
                    
                    a=ggplot(simsfiltdat_adjusted,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + 
                    ggtitle("") #+ ylim(-10,5)
                    a
                    
                    ggsave("Fig2.pdf",dpi=600,device="pdf")
                    
    ## 2.3. Perform MAIN simulation (end)
    ## 2.3. Perform MAIN simulation with combined growth and levels effect (start)


            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef=-0.05 #effect size - change in growth per degree warming()
            coef2=-0.05
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,3,5,10,15)
            nsims=1000
            sims=array(dim=c(nsims,length(periods),3))
            sims_adjusted=array(dim=c(nsims,length(periods),3))
            sims_ar=array(dim=c(nsims,length(periods),2))
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                randomgrowth_gl1=basegr #levels impact model
                for(j in 2:time){
                #linear
                randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_gl1=c(randomgrowth_gl1, basegr+(randomtemp[j]-randomtemp[j-1])*0.02+(randomtemp[j]*-0.07)+rnorm(1,mean=0,sd=growthsd))
                }
                dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l,gl1=randomgrowth_gl1)
                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                for(k in 1:length(periods)){
                mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                mod_l=lm(l~temp+templag,data=dataset)$coef[2]
                mod_gl1=lm(gl1~temp+templag,data=dataset)$coef[2]

                if (k==1){
                        sims[i,k,]=c(mod_g,mod_l,mod_gl1)
                        sims_adjusted[i,k,]=c(mod_g,mod_l,mod_gl1)
                        
                    }else{
                    x <- seq(1:length(randomtemp))
                    randomtemp <- lm(randomtemp~x)$residuals
                    tempts <- pass.filt(randomtemp, W=periods[k], type="low", method="Butterworth")
                    
                    ratio <- median(randomtemp/tempts)
                    growth <- randomgrowth_g
                    #growth <- pass.filt(growth, W=periods[k], type="low", method="Butterworth")
                    level <- randomgrowth_l
                    #level <- pass.filt(level, W=periods[k], type="low", method="Butterworth")
                    gl1 <- randomgrowth_gl1
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        gl1=unclass(gl1),level = unclass(level))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","gl1","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]
                    mod_gl1filt=lm(gl1~temp,data=filt)$coef[2]
                    dw <- dwtest(lm(growth~temp,data=filt)) #test for autocorrelation
                    sims[i,k,]=c(mod_gfilt,mod_lfilt,mod_gl1filt)
                    sims_adjusted[i,k,]=c(mod_gfilt/ratio,mod_lfilt/ratio,mod_gl1filt/ratio)
                    sims_ar[i,k,]=c(dw$statistic[1],dw$p.value[1])
                }
                

                }
                if(i%%50==0) print(i)
            }
            #
          
                # Test - Figure 2
                    
                    # simsfiltmean=apply(sims_ar,MARGIN=c(2,3),FUN=mean)
                    # simsfiltsd=apply(sims_ar,MARGIN=c(2,3),FUN=sd)
                    # colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    # simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    # colnames(simsfiltdat)=c("frequencies","ImpactType","dw_estimate","dw_pvalue")
                    # theme_set(theme_bw(base_size = 20))
                    
                    # simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    # simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    #     levels = c("Unfiltered", "10 years", "20 years", "50 years", "100 years"))
                    # simsfiltdat_ar <- simsfiltdat

                    
                    simsfiltmean=apply(sims_adjusted,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_adjusted,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth only \n (\u03b3 = -0.05, \u03b2 = 0) \n","Growth and Level \n (\u03b3 = -0.07, \u03b2 = 0.02)\n","Level only \n (\u03b3 = 0, \u03b2 = -0.05)\n");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    simsfiltdat_adjusted <- simsfiltdat

                    levels(simsfiltdat_adjusted$ImpactType)
                    simsfiltdat_adjusted$ImpactType <- factor(simsfiltdat_adjusted$ImpactType, 
                        levels=c("Level only \n (\u03b3 = 0, \u03b2 = -0.05)\n",
                         "Growth only \n (\u03b3 = -0.05, \u03b2 = 0) \n",
                         "Growth and Level \n (\u03b3 = -0.07, \u03b2 = 0.02)\n"))
                    
                    a=ggplot(simsfiltdat_adjusted,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType,alpha=ImpactType))+
                    geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model",alpha="Impact Model")
                    a=a+scale_color_manual(values=c("#7375a4","#d3818c","#d3818c")) + 
                    scale_alpha_manual(values=c(1,1,0.3))+
                    ggtitle("") #+ ylim(-10,5)
                    a
                    
                    ggsave("Fig2_version2.pdf",dpi=600,device="pdf")
                    ggsave("Fig2_version2.png",dpi=600)
                    
    ## 2.3. Perform MAIN simulation with combined growth and levels effect (end)

    ## 2.4. Additional simulations and figures (Used for Response to Reviewers)
        ## 2.4a. Perform simulation - WITH QUADRATIC TERM (start)
            
            # file=nc_open("C:/Users/bastien/Box/Long Run GDP Growth/Data/LMR Data/air_MCruns_ensemble_mean_LMRv2.0.nc")
            # gtemp=ncvar_get(file,"air",start=c(1,1,1,1),count=c(-1,-1,1,-1))
            # lon <- ncvar_get(file,"lon")
            # lat <- ncvar_get(file,"lat")
            # time <- ncvar_get(file,"time")
            # tunits <- ncatt_get(file,"time","units")
            # data(world)
            # geom_iso <- world$geom[world$iso_a2=="US"]
            # geom_iso <- st_shift_longitude(geom_iso)
            # geom_iso <- st_cast(geom_iso, "POLYGON")
            # geom_iso <-as_Spatial(geom_iso)
            
            # gtemp <- brick(gtemp, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
            # ustemp_lmr <- rep(NA,1500)
            # for (i in 1:1500){
            #    gtemp_year <- raster::subset(gtemp,i)
            #   #gtemp_year <- rotate(gtemp_year)
            #  gtemp_year <- t(gtemp_year)
            # ustemp <- extract(gtemp_year, geom_iso,method='simple')
            #     ustemp_lmr[i] <- mean(ustemp[[1]])

            # }
            
            # print(file)
            # ustemp_lmr_orig <- ustemp_lmr
            # x <- seq(length(ustemp_lmr_orig))
            # ustemp_lmr=lm(ustemp_lmr_orig~x)$residuals
            # p=periodogram(ustemp_lmr)
            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef=0.0127# BHM coefficient effect size - change in growth per degree warming()
            coef2=-0.0005 #BHM coefficient
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,3,5,10,15)
            nsims=500
            sims=array(dim=c(nsims,length(periods),2))
            for(i in 1:nsims){
                #randomtemp=Re(randomts(gtemp))[1:time]+25
                #randomtemp=Re(randomts(ustemp_lmr))[1:time] +25
                randomtemp = rnorm(time,mean=25,sd=0.5)
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                for(j in 2:time){
                #quadratic standard form
                randomgrowth_g=c(randomgrowth_g, basegr+coef*(randomtemp[j])+coef2*(randomtemp[j])^2+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l=c(randomgrowth_l, basegr+randomtemp[j]*coef+coef2*randomtemp[j]^2-coef*randomtemp[j-1]-coef2* randomtemp[j-1]^2+rnorm(1,mean=0,sd=growthsd))
               }
                x <- seq(length(randomtemp))
                randomtemp <- lm(randomtemp~x)$residuals
                dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                dataset <- dataset[-c(1),]
                for(k in 1:length(periods)){
                mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                if (k==1){sims[i,k,]=c(mod_g,mod_l)
                    }else{
                    tempts <- pass.filt(randomtemp, W=periods[k], type="low", method="Butterworth")
                    ratio <- median(randomtemp/tempts)
                    growth <- randomgrowth_g
                    level <- randomgrowth_l
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        level = unclass(level))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]/ratio
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]/ratio
                    sims[i,k,]=c(mod_gfilt,mod_lfilt)
                }
                

                }
                if(i%%50==0) print(i)
            }
            # Comparing quadratic
                #Unfiltered Growth (start)
                    linreg <- lm(g~temp,data=dataset)
                    nonlinreg <- lm(g~temp+I(temp^2),data=dataset)
                   
                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated2 <- gestimated
                    dataset$ci12 <- ci12
                    dataset$ci22 <- ci22
                    dataset$x <- xmat[,2]

                    model <- linreg
                    sigma <- vcov(model)
                    beta.hat <- coef(model)
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated <- gestimated
                    dataset$ci1<- ci1
                    dataset$ci2 <- ci2

                    cols=c("#d3818c","#7375a4")
                    u_g <- ggplot(data=dataset, aes(x=temp,y=g))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature")+
                    ylab("Growth") + xlim(-1.2,1.2)+
                    ggtitle("Unfiltered Growth World")
                #Unfiltered Growth (end)
                
                #Unfiltered levels (start)
                    linreg <- lm(l~temp,data=dataset)
                    nonlinreg <- lm(l~temp+I(temp^2),data=dataset)
                    summary(linreg)            
                    summary(nonlinreg)

                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated2 <- gestimated
                    dataset$ci12 <- ci12
                    dataset$ci22 <- ci22
                    dataset$x <- xmat[,2]

                    model <- linreg
                    Sigma <- vcov(model)
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==coefT)
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated <- gestimated
                    dataset$ci12 <- ci1
                    dataset$ci22 <- ci2

                    cols=c("#d3818c","#7375a4")
                    u_l <- ggplot(data=dataset, aes(x=temp,y=l))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature")+
                    ylab("Growth") + xlim(-1.2,1.2)+
                    ggtitle("Unfiltered Levels World")
                #Unfiltered levels (end)

                #15y-filtered Growth (start)
                    linreg <- lm(growth~temp,data=filt)
                    nonlinreg <- lm(growth~temp+I(temp^2),data=filt)
                    summary(linreg)            
                    summary(nonlinreg)

                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated2 <- gestimated
                    filt$ci12 <- ci12
                    filt$ci22 <- ci22
                    filt$x <- xmat[,2]

                    model <- linreg
                    sigma <- vcov(model)
                    beta.hat <- coef(model)
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated <- gestimated
                    filt$ci1 <- ci1
                    filt$ci2 <- ci2

                    cols=c("#d3818c","#7375a4")
                    g_15 <- ggplot(data=filt, aes(x=temp,y=growth))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature")+
                    ylab("Growth") + xlim(-1.2,1.2)+
                    ggtitle("15-y filtered Growth World")
                #15y-filtered Growth (end)
                
                #15y-filtered levels (start)
                    linreg <- lm(levels~temp,data=filt)
                    nonlinreg <- lm(levels~temp+I(temp^2),data=filt)
                    summary(linreg)            
                    summary(nonlinreg)

                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated2 <- gestimated
                    filt$ci12 <- ci12
                    filt$ci22 <- ci22
                    filt$x <- xmat[,2]

                    model <- linreg
                    Sigma <- vcov(model)
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==coefT)
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated <- gestimated
                    filt$ci1 <- ci1
                    filt$ci2 <- ci2

                    cols=c("#d3818c","#7375a4")
                    l_15 <- ggplot(data=filt, aes(x=temp,y=levels))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature")+
                    ylab("Growth") + xlim(-1.2,1.2)+
                    ggtitle("15-y filtered Levels World")
                #15y-filteredd levels (end)
                
                # Test 
                    simsfiltmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                     
                    a=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth\n Impact (% per Degree)",col="Impact Model")
                    a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Impacts estimation using low pass filters") #+ ylim(-10,5)
                    a
                # Comparing quadratic (end)
                
                ggarrange(ggarrange(u_l,l_15,u_g,g_15),a,heights = c(2, 0.7),nrow=2,ncol=1)
                #ggsave("Figures/FigResponse_Quadratic.png",dpi=600)
        ## 2.4a. Perform simulation - WITH QUADRATIC TERM (end)

        ## 2.4b.Perform simulation - Combination Gammas and Betas (start)


            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef_gamma=-0.05 #effect size - change in growth per degree warming()
            coef_beta=-0.05
            coef_beta2=0.05
            coef_beta3=0.07
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,3,5,10,15)
            nsims=500
            sims=array(dim=c(nsims,length(periods),2))
            sims_adjusted=array(dim=c(nsims,length(periods),3))
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                randomgrowth_l2=basegr #levels impact model
                for(j in 2:time){
                randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef_gamma)+(randomtemp[j]-randomtemp[j-1])*coef_beta+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]*coef_gamma)+(randomtemp[j]-randomtemp[j-1])*coef_beta2+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l2=c(randomgrowth_l2, basegr+(randomtemp[j]*coef_gamma)+(randomtemp[j]-randomtemp[j-1])*coef_beta3+rnorm(1,mean=0,sd=growthsd))
                
                }
                dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l,l2=randomgrowth_l2)
                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                for(k in 1:length(periods)){
                mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                mod_l=lm(l~temp+templag,data=dataset)$coef[2]
                mod_l2=lm(l2~temp+templag,data=dataset)$coef[2]

                if (k==1){
                        sims[i,k,]=c(mod_g,mod_l)
                        sims_adjusted[i,k,]=c(mod_g,mod_l,mod_l2)
                        
                    }else{
                    x <- seq(1:length(randomtemp))
                    randomtemp <- lm(randomtemp~x)$residuals
                    tempts <- pass.filt(randomtemp, W=periods[k], type="low", method="Butterworth")
                    
                    ratio <- median(randomtemp/tempts)
                    #ratio <- median(randomtemp / tempts)
                    #ratio <- median( abs(randomtemp) / abs(tempts))
                    growth <- randomgrowth_g
                    level <- randomgrowth_l
                    level2 <- randomgrowth_l2
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        level = unclass(level),level2 = unclass(level2))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","levels2","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]
                    mod_lfilt2=lm(levels2~temp,data=filt)$coef[2]
                    sims[i,k,]=c(mod_gfilt,mod_lfilt)
                    sims_adjusted[i,k,]=c(mod_gfilt/ratio,mod_lfilt/ratio,mod_lfilt2/ratio)
                }
                

                }
                if(i%%50==0) print(i)
            }
            # Comparing quadratic
          
                # Test 
                    simsfiltmean=apply(sims_adjusted,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_adjusted,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Attenuating","Not Converging","Changing Sign");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    simsfiltdat_adjusted <- simsfiltdat

                    # simsfiltmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                    # simsfiltsd=apply(sims,MARGIN=c(2,3),FUN=sd)
                    # colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    # simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    # colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    # theme_set(theme_bw(base_size = 20))
                    
                    # #simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    # simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    # simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    #     levels = c("Unfiltered", "10 years", "20 years", "50 years", "100 years"))
                    
                    
                    a=ggplot(simsfiltdat_adjusted,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    a=a+scale_color_manual(values=c("#7fc97f","#beaed4","#fdc086"),
                        labels=c(expression(paste(beta,"=-0.05; ",gamma,"=-0.05")),expression(paste(beta,"=0.05; ",gamma,"=-0.05")),expression(paste(beta,"=0.07; ",gamma,"=-0.05")))) + 
                    ggtitle("Combined effects") #+ ylim(-10,5)
                    a


                    b=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    b=b+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    b=b+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    b=b+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    b=b+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Impacts estimation using low pass filters") #+ ylim(-10,5)
                    b
                # Comparing quadratic (end)
                
                ggarrange(b,a,nrow=2,ncol=1)
        ## 2.4b. Perform simulation - Combination Gammas and Betas (end)

        ## 2.4c. Simulation comparing low-pass filter with distributed lags
            
                time = 120
                basegr=0.01 #1% per year baseline growth rate
                #basegr=0.025 #2.5% USA per year baseline growth rate
                
                start=100
                baseline=rep(basegr,time)
                coef=-0.05 #effect size - change in growth per degree warming()
                growthsd=0.005 #standard deviation of growth variability unexplained by temperature
                #growthsd=0.01 #USA standard deviation of growth variability unexplained by temperature
                periods <- c(0,5,10,15,50)
                # arglagbothlist <- list(list(knots=logknots(periods[k],2)),
                #     list(fun="poly",degree = 4),
                #     list(fun="strata",breaks = 4),
                #     list(fun="ns",df = 4),
                #     list(fun="bs",df = 4))
                # #nsims=176
                nsims=1000
                sims=array(dim=c(nsims,length(periods),2))
                sims_lp=array(dim=c(nsims,length(periods),2))
                sims_lp_se=array(dim=c(nsims,length(periods),2))
                sims_se=array(dim=c(nsims,length(periods),2))
                sims_dlnm=array(dim=c(nsims,length(periods),2))
                sims_dlnm_se=array(dim=c(nsims,length(periods),2))
                for(i in 1:nsims){
                    randomtemp=Re(randomts(gtemp))[1:time]
                    randomgrowth_g=basegr #growth impact model
                    randomgrowth_l=basegr #levels impact model
                    for(j in 6:time){
                    randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                    randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                    #randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-0.25*randomtemp[j-1]-0.25*randomtemp[j-2]-
                    #  0.25*randomtemp[j-3]-0.25*randomtemp[j-4]-0.25*randomtemp[j-5])*coef+rnorm(1,mean=0,sd=growthsd))
                    }
                    dataset=data.frame(years=1:(time-5),temp=randomtemp[6:(time)],
                    g=randomgrowth_g[2:(time-4)],l=randomgrowth_l[2:(time-4)])
                    #dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                    dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                    for(k in 1:length(periods)){
                        mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                        mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                    if (k==1){
                            m_g <- lm(g~temp+templag,data=dataset)
                            m_l <- lm(l~temp+templag,data=dataset)
                            mod_g <- m_g$coefficients[2]
                            mod_l <- m_l$coefficients[2]
                            sims[i,k,]=c(mod_g,mod_l)
                            se_l <- summary(m_l)$coefficients[2,2]
                            se_g <- summary(m_g)$coefficients[2,2]
                            sims_se[i,k,]=c(se_g,se_l)
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            sims_lp[i,k,]=c(mod_g,mod_l)

                            sims_dlnm[i,k,] <- c(mod_g,mod_l)
                            sims_dlnm_se[i,k,]=c(se_g,se_l)

                        }else{
                            tempts <- pass.filt(randomtemp[6:(time)], W=periods[k], type="low", method="Butterworth")
                            ratio <- median(randomtemp[6:(time)]/tempts)
                            growth <- randomgrowth_g[2:(time-4)]#[6:time]
                            level <- randomgrowth_l[2:(time-4)]#[6:time]
                            filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                                level = unclass(level))
                            filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                            names(filt) <- c("temp","growth","levels","templag")
                            mod_gfilt=lm(growth~temp,data=filt)
                            mod_lfilt=lm(levels~temp,data=filt)
                            sims_lp[i,k,]=c(mod_gfilt$coef[2]/ratio,mod_lfilt$coef[2]/ratio)
                            se_l <- summary(mod_lfilt)$coefficients[2,2]
                            se_g <- summary(mod_gfilt)$coefficients[2,2]
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            
                            num_lags <- c(0,periods[k])
                            arglagboth <- list(fun="ns",df = 3)
                            cbtemp <- crossbasis(dataset$temp,
                            lag=num_lags,
                            argvar=list("poly",degree=1,scale=1),
                            arglag=arglagboth)
                            #model_g <-  felm((randomgrowth_g[2:(time-4)]) ~ cbtemp)
                            model_g <-  lm(dataset$g ~ cbtemp )
                            redvar_g <- crossreduce(cbtemp,model_g,cen=0)
                            #plot(redvar_g)
                            #summary(redvar_g)
                            beta.hat_g <- redvar_g$coefficients
                            sigma_g <- redvar_g$vcov
                            #model_l <-  felm((randomgrowth_l[2:(time-4)]) ~ cbtemp)
                            model_l <-  felm(dataset$l ~ cbtemp)
                            redvar_l <- crossreduce(cbtemp,model_l,cen=0)
                            beta.hat_l <- redvar_l$coefficients
                            sigma_l <- redvar_l$vcov
                            dynlag_l <-  crosspred(cbtemp, model_l,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            dynlag_g <-  crosspred(cbtemp, model_g,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            
                            dynlag_l$allfit
                            sims_dlnm[i,k,] <- c(dynlag_g$allfit[[2]],beta.hat_l[[1]])
                            sims_dlnm_se[i,k,]=c(dynlag_g$allse[[2]],sigma_l[[1]]^0.5)

                            if (k>1){
                                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                                dataset$templag2 =c(NA,dataset$templag[1:(dim(dataset)[1]-1)])
                                dataset$templag3 =c(NA,dataset$templag2[1:(dim(dataset)[1]-1)])
                                dataset$templag4 =c(NA,dataset$templag3[1:(dim(dataset)[1]-1)])
                                dataset$templag5 =c(NA,dataset$templag4[1:(dim(dataset)[1]-1)])
                                m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5,data=dataset)
                                m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5,data=dataset)
                                mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                sims[i,k,]=c(mod_g,mod_l)
                                if (k>2) {
                                    dataset$templag6 =c(NA,dataset$templag5[1:(dim(dataset)[1]-1)])
                                    dataset$templag7 =c(NA,dataset$templag6[1:(dim(dataset)[1]-1)])
                                    dataset$templag8 =c(NA,dataset$templag7[1:(dim(dataset)[1]-1)])
                                    dataset$templag9 =c(NA,dataset$templag8[1:(dim(dataset)[1]-1)])
                                    dataset$templag10 =c(NA,dataset$templag9[1:(dim(dataset)[1]-1)])
                                    m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10,data=dataset)
                                    m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10,data=dataset)
                                    mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                    mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                    if (k>3){
                                        dataset$templag11 =c(NA,dataset$templag10[1:(dim(dataset)[1]-1)])
                                        dataset$templag12 =c(NA,dataset$templag11[1:(dim(dataset)[1]-1)])
                                        dataset$templag13 =c(NA,dataset$templag12[1:(dim(dataset)[1]-1)])
                                        dataset$templag14 =c(NA,dataset$templag13[1:(dim(dataset)[1]-1)])
                                        dataset$templag15 =c(NA,dataset$templag14[1:(dim(dataset)[1]-1)])
                                        m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15,data=dataset)
                                        m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15,data=dataset)

                                        if (k>4){
                                            dataset$templag16 =c(NA,dataset$templag15[1:(dim(dataset)[1]-1)])
                                            dataset$templag17 =c(NA,dataset$templag16[1:(dim(dataset)[1]-1)])
                                            dataset$templag18 =c(NA,dataset$templag17[1:(dim(dataset)[1]-1)])
                                            dataset$templag19 =c(NA,dataset$templag18[1:(dim(dataset)[1]-1)])
                                            dataset$templag20 =c(NA,dataset$templag19[1:(dim(dataset)[1]-1)])
                                            dataset$templag21 =c(NA,dataset$templag20[1:(dim(dataset)[1]-1)])
                                            dataset$templag22 =c(NA,dataset$templag21[1:(dim(dataset)[1]-1)])
                                            dataset$templag23 =c(NA,dataset$templag22[1:(dim(dataset)[1]-1)])
                                            dataset$templag24 =c(NA,dataset$templag23[1:(dim(dataset)[1]-1)])
                                            dataset$templag25 =c(NA,dataset$templag24[1:(dim(dataset)[1]-1)])
                                            dataset$templag26 =c(NA,dataset$templag25[1:(dim(dataset)[1]-1)])
                                            dataset$templag27 =c(NA,dataset$templag26[1:(dim(dataset)[1]-1)])
                                            dataset$templag28 =c(NA,dataset$templag27[1:(dim(dataset)[1]-1)])
                                            dataset$templag29 =c(NA,dataset$templag28[1:(dim(dataset)[1]-1)])
                                            dataset$templag30 =c(NA,dataset$templag29[1:(dim(dataset)[1]-1)])
                                            dataset$templag31 =c(NA,dataset$templag30[1:(dim(dataset)[1]-1)])
                                            dataset$templag32 =c(NA,dataset$templag31[1:(dim(dataset)[1]-1)])
                                            dataset$templag33 =c(NA,dataset$templag32[1:(dim(dataset)[1]-1)])
                                            dataset$templag34 =c(NA,dataset$templag33[1:(dim(dataset)[1]-1)])
                                            dataset$templag35 =c(NA,dataset$templag34[1:(dim(dataset)[1]-1)])
                                            dataset$templag36 =c(NA,dataset$templag35[1:(dim(dataset)[1]-1)])
                                            dataset$templag37 =c(NA,dataset$templag36[1:(dim(dataset)[1]-1)])
                                            dataset$templag38 =c(NA,dataset$templag37[1:(dim(dataset)[1]-1)])
                                            dataset$templag39 =c(NA,dataset$templag38[1:(dim(dataset)[1]-1)])
                                            dataset$templag40 =c(NA,dataset$templag39[1:(dim(dataset)[1]-1)])
                                            dataset$templag41 =c(NA,dataset$templag40[1:(dim(dataset)[1]-1)])
                                            dataset$templag42 =c(NA,dataset$templag41[1:(dim(dataset)[1]-1)])
                                            dataset$templag43 =c(NA,dataset$templag42[1:(dim(dataset)[1]-1)])
                                            dataset$templag44 =c(NA,dataset$templag43[1:(dim(dataset)[1]-1)])
                                            dataset$templag45 =c(NA,dataset$templag44[1:(dim(dataset)[1]-1)])
                                            dataset$templag46 =c(NA,dataset$templag45[1:(dim(dataset)[1]-1)])
                                            dataset$templag47 =c(NA,dataset$templag46[1:(dim(dataset)[1]-1)])
                                            dataset$templag48 =c(NA,dataset$templag47[1:(dim(dataset)[1]-1)])
                                            dataset$templag49 =c(NA,dataset$templag48[1:(dim(dataset)[1]-1)])
                                            dataset$templag50 =c(NA,dataset$templag49[1:(dim(dataset)[1]-1)])
                                            m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15+templag16+templag17+templag18+templag19+templag20+
                                            templag21+templag22+templag23+templag24+templag25+templag26+templag27+templag28+templag29+templag30+
                                            templag31+templag32+templag33+templag34+templag35+templag36+templag37+templag38+templag39+templag40+
                                            templag41+templag42+templag43+templag44+templag45+templag46+templag47+templag48+templag49+templag50,data=dataset)
                                            m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15+templag16+templag17+templag18+templag19+templag20+
                                            templag21+templag22+templag23+templag24+templag25+templag26+templag27+templag28+templag29+templag30+
                                            templag31+templag32+templag33+templag34+templag35+templag36+templag37+templag38+templag39+templag40+
                                            templag41+templag42+templag43+templag44+templag45+templag46+templag47+templag48+templag49+templag50,data=dataset)
                                        }
                                    }
                                    
                                }
                                
                                


                            }
                            mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                            summary(m_l)
                            #sum(m_g$coefficients[2:length(m_g$coefficients)])
                            mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                            sigma <- vcov(m_l)
                            a <- rep(1,length(m_l$coefficients))
                            a[1] <- 0
                            se_sum_l <- c(sqrt(t(a) %*% sigma %*% a))
                            sigma <- vcov(m_g)
                            se_sum_g <- c(sqrt(t(a) %*% sigma %*% a))
                            sims_se[i,k,]=c(se_sum_g,se_sum_l)
                            sims[i,k,]=c(mod_g,mod_l)
                    }
                    

                    }
                    if(i%%50==0) print(i)
                }


                ## Plotting Figure 2  - Simulation exercise (start)
                
                
                library(ggbreak) 
                simsfiltmean=apply(sims_dlnm,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_dlnm_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsdlnmdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsdlnmdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsdlnmdat$periodsregationPeriod <- c("No lags", "5 years", "10 years", "15 years", "50 years","No lags", "5 years", "10 years", "15 years", "50 years")
                simsdlnmdat$periodsregationPeriod <- factor(simsdlnmdat$periodsregationPeriod,
                    levels = c("No lags", "5 years", "10 years", "15 years", "50 years"))
                glimpse(simsdlnmdat)

                
                
                a=ggplot(simsdlnmdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Number of Lags",y="",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("") +ylim(-50,45)
                c_110 <- a
                 #c_110   + scale_y_break(c(-50,-10))  + scale_y_break(c(10, 40))   + theme_minimal()        





                simslagmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean)=c("Growth","Level");rownames(simslagmean)=periods
                simslagmean_se=apply(sims_se,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean_se)=c("Growth","Level");rownames(simslagmean_se)=periods
                simslagdat=cbind(melt(simslagmean),melt(simslagmean_se)[,3])
                colnames(simslagdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simslagdat$periodsregationPeriod <- c("No lags", "5 years", "10 years", "15 years", "50 years","No lags", "5 years", "10 years", "15 years", "50 years")
                simslagdat$periodsregationPeriod <- factor(simslagdat$periodsregationPeriod,
                    levels = c("No lags", "5 years", "10 years", "15 years", "50 years"))
                glimpse(simslagdat)
                
                
                a=ggplot(simslagdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Number of lags",y="",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("") +ylim(-50,45)
                b_110 <- a
                b_110



                simsfiltmean=apply(sims_lp,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_lp_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsfiltdat$periodsregationPeriod <- c("Unfiltered", "5 years", "10 years", "15 years", "50 years","Unfiltered", "5 years", "10 years", "15 years", "50 years")
                simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    levels = c("Unfiltered", "5 years", "10 years", "15 years", "50 years"))
                glimpse(simsfiltdat)
                
                
                a=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Filters",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("") +ylim(-50,45)
                a_110 <- a

               plot_110 <-  ggarrange(a_110,b_110,c_110,nrow=1,ncol=3,common.legend=TRUE,legend="bottom")
               plot_110 <-  annotate_figure(plot_110, top = text_grob("120-year timeseries", 
               color = "black", face = "bold", size = 14))
                plot_110

            
                time = 70
                start=100
                baseline=rep(basegr,time)
                coef=-0.05 #effect size - change in growth per degree warming()
                growthsd=0.005 #standard deviation of growth variability unexplained by temperature
                #growthsd=0.01 #USA standard deviation of growth variability unexplained by temperature
                periods <- c(0,5,10,15,30)
                # arglagbothlist <- list(list(knots=logknots(periods[k],2)),
                #     list(fun="poly",degree = 4),
                #     list(fun="strata",breaks = 4),
                #     list(fun="ns",df = 4),
                #     list(fun="bs",df = 4))
                # #nsims=176
                nsims=1000
                sims=array(dim=c(nsims,length(periods),2))
                sims_lp=array(dim=c(nsims,length(periods),2))
                sims_lp_se=array(dim=c(nsims,length(periods),2))
                sims_se=array(dim=c(nsims,length(periods),2))
                sims_dlnm=array(dim=c(nsims,length(periods),2))
                sims_dlnm_se=array(dim=c(nsims,length(periods),2))
                for(i in 1:nsims){
                    randomtemp=Re(randomts(gtemp))[1:time]
                    randomgrowth_g=basegr #growth impact model
                    randomgrowth_l=basegr #levels impact model
                    for(j in 6:time){
                    randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                    randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                    #randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-0.25*randomtemp[j-1]-0.25*randomtemp[j-2]-
                    #  0.25*randomtemp[j-3]-0.25*randomtemp[j-4]-0.25*randomtemp[j-5])*coef+rnorm(1,mean=0,sd=growthsd))
                    }
                    dataset=data.frame(years=1:(time-5),temp=randomtemp[6:(time)],
                    g=randomgrowth_g[2:(time-4)],l=randomgrowth_l[2:(time-4)])
                    #dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                    dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                    for(k in 1:length(periods)){
                        mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                        mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                    if (k==1){
                            m_g <- lm(g~temp+templag,data=dataset)
                            m_l <- lm(l~temp+templag,data=dataset)
                            mod_g <- m_g$coefficients[2]
                            mod_l <- m_l$coefficients[2]
                            sims[i,k,]=c(mod_g,mod_l)
                            se_l <- summary(m_l)$coefficients[2,2]
                            se_g <- summary(m_g)$coefficients[2,2]
                            sims_se[i,k,]=c(se_g,se_l)
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            sims_lp[i,k,]=c(mod_g,mod_l)

                            sims_dlnm[i,k,] <- c(mod_g,mod_l)
                            sims_dlnm_se[i,k,]=c(se_g,se_l)

                        }else{
                            tempts <- pass.filt(randomtemp[6:(time)], W=periods[k], type="low", method="Butterworth")
                            ratio <- median(randomtemp[6:(time)]/tempts)
                            growth <- randomgrowth_g[2:(time-4)]#[6:time]
                            level <- randomgrowth_l[2:(time-4)]#[6:time]
                            filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                                level = unclass(level))
                            filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                            names(filt) <- c("temp","growth","levels","templag")
                            mod_gfilt=lm(growth~temp,data=filt)
                            mod_lfilt=lm(levels~temp,data=filt)
                            sims_lp[i,k,]=c(mod_gfilt$coef[2]/ratio,mod_lfilt$coef[2]/ratio)
                            se_l <- summary(mod_lfilt)$coefficients[2,2]
                            se_g <- summary(mod_gfilt)$coefficients[2,2]
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            
                            num_lags <- c(0,periods[k])
                            arglagboth <- list(fun="ns",df = 3)
                            cbtemp <- crossbasis(dataset$temp,
                            lag=num_lags,
                            argvar=list("poly",degree=1,scale=1),
                            arglag=arglagboth)
                            #model_g <-  felm((randomgrowth_g[2:(time-4)]) ~ cbtemp)
                            model_g <-  lm(dataset$g ~ cbtemp )
                            redvar_g <- crossreduce(cbtemp,model_g,cen=0)
                            #plot(redvar_g)
                            #summary(redvar_g)
                            beta.hat_g <- redvar_g$coefficients
                            sigma_g <- redvar_g$vcov
                            #model_l <-  felm((randomgrowth_l[2:(time-4)]) ~ cbtemp)
                            model_l <-  felm(dataset$l ~ cbtemp)
                            redvar_l <- crossreduce(cbtemp,model_l,cen=0)
                            beta.hat_l <- redvar_l$coefficients
                            sigma_l <- redvar_l$vcov
                            dynlag_l <-  crosspred(cbtemp, model_l,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            dynlag_g <-  crosspred(cbtemp, model_g,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            
                            dynlag_l$allfit
                            sims_dlnm[i,k,] <- c(dynlag_g$allfit[[2]],beta.hat_l[[1]])
                            sims_dlnm_se[i,k,]=c(dynlag_g$allse[[2]],sigma_l[[1]]^0.5)

                            if (k>1){
                                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                                dataset$templag2 =c(NA,dataset$templag[1:(dim(dataset)[1]-1)])
                                dataset$templag3 =c(NA,dataset$templag2[1:(dim(dataset)[1]-1)])
                                dataset$templag4 =c(NA,dataset$templag3[1:(dim(dataset)[1]-1)])
                                dataset$templag5 =c(NA,dataset$templag4[1:(dim(dataset)[1]-1)])
                                m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5,data=dataset)
                                m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5,data=dataset)
                                mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                sims[i,k,]=c(mod_g,mod_l)
                                if (k>2) {
                                    dataset$templag6 =c(NA,dataset$templag5[1:(dim(dataset)[1]-1)])
                                    dataset$templag7 =c(NA,dataset$templag6[1:(dim(dataset)[1]-1)])
                                    dataset$templag8 =c(NA,dataset$templag7[1:(dim(dataset)[1]-1)])
                                    dataset$templag9 =c(NA,dataset$templag8[1:(dim(dataset)[1]-1)])
                                    dataset$templag10 =c(NA,dataset$templag9[1:(dim(dataset)[1]-1)])
                                    m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10,data=dataset)
                                    m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10,data=dataset)
                                    mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                    mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                    if (k>3){
                                        dataset$templag11 =c(NA,dataset$templag10[1:(dim(dataset)[1]-1)])
                                        dataset$templag12 =c(NA,dataset$templag11[1:(dim(dataset)[1]-1)])
                                        dataset$templag13 =c(NA,dataset$templag12[1:(dim(dataset)[1]-1)])
                                        dataset$templag14 =c(NA,dataset$templag13[1:(dim(dataset)[1]-1)])
                                        dataset$templag15 =c(NA,dataset$templag14[1:(dim(dataset)[1]-1)])
                                        m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15,data=dataset)
                                        m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15,data=dataset)

                                        if (k>4){
                                            dataset$templag16 =c(NA,dataset$templag15[1:(dim(dataset)[1]-1)])
                                            dataset$templag17 =c(NA,dataset$templag16[1:(dim(dataset)[1]-1)])
                                            dataset$templag18 =c(NA,dataset$templag17[1:(dim(dataset)[1]-1)])
                                            dataset$templag19 =c(NA,dataset$templag18[1:(dim(dataset)[1]-1)])
                                            dataset$templag20 =c(NA,dataset$templag19[1:(dim(dataset)[1]-1)])
                                            dataset$templag21 =c(NA,dataset$templag20[1:(dim(dataset)[1]-1)])
                                            dataset$templag22 =c(NA,dataset$templag21[1:(dim(dataset)[1]-1)])
                                            dataset$templag23 =c(NA,dataset$templag22[1:(dim(dataset)[1]-1)])
                                            dataset$templag24 =c(NA,dataset$templag23[1:(dim(dataset)[1]-1)])
                                            dataset$templag25 =c(NA,dataset$templag24[1:(dim(dataset)[1]-1)])
                                            dataset$templag26 =c(NA,dataset$templag25[1:(dim(dataset)[1]-1)])
                                            dataset$templag27 =c(NA,dataset$templag26[1:(dim(dataset)[1]-1)])
                                            dataset$templag28 =c(NA,dataset$templag27[1:(dim(dataset)[1]-1)])
                                            dataset$templag29 =c(NA,dataset$templag28[1:(dim(dataset)[1]-1)])
                                            dataset$templag30 =c(NA,dataset$templag29[1:(dim(dataset)[1]-1)])
                                            
                                            m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15+templag16+templag17+templag18+templag19+templag20+
                                            templag21+templag22+templag23+templag24+templag25+templag26+templag27+templag28+templag29+templag30,data=dataset)
                                            m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15+templag16+templag17+templag18+templag19+templag20+
                                            templag21+templag22+templag23+templag24+templag25+templag26+templag27+templag28+templag29+templag30,data=dataset)
                                        }
                                    }
                                    
                                }
                                
                                


                            }
                            mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                            summary(m_l)
                            #sum(m_g$coefficients[2:length(m_g$coefficients)])
                            mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                            sigma <- vcov(m_l)
                            a <- rep(1,length(m_l$coefficients))
                            a[1] <- 0
                            se_sum_l <- c(sqrt(t(a) %*% sigma %*% a))
                            sigma <- vcov(m_g)
                            se_sum_g <- c(sqrt(t(a) %*% sigma %*% a))
                            sims_se[i,k,]=c(se_sum_g,se_sum_l)
                            sims[i,k,]=c(mod_g,mod_l)
                    }
                    

                    }
                    if(i%%50==0) print(i)
                }


                ## Plotting Figure 2  - Simulation exercise (start)
                
                
                simsfiltmean=apply(sims_dlnm,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_dlnm_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsdlnmdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsdlnmdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsdlnmdat$periodsregationPeriod <- c("No lags", "5 years", "10 years", "15 years", "30 years","No lags", "5 years", "10 years", "15 years", "30 years")
                simsdlnmdat$periodsregationPeriod <- factor(simsdlnmdat$periodsregationPeriod,
                    levels = c("No lags", "5 years", "10 years", "15 years", "30 years"))
                glimpse(simsdlnmdat)

                
                
                a=ggplot(simsdlnmdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="",y="",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Distributed lag non-linear model")+ylim(-95,90)
                c <- a
                c            





                simslagmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean)=c("Growth","Level");rownames(simslagmean)=periods
                simslagmean_se=apply(sims_se,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean_se)=c("Growth","Level");rownames(simslagmean_se)=periods
                simslagdat=cbind(melt(simslagmean),melt(simslagmean_se)[,3])
                colnames(simslagdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simslagdat$periodsregationPeriod <- c("No lags", "5 years", "10 years", "15 years", "30 years","No lags", "5 years", "10 years", "15 years", "30 years")
                simslagdat$periodsregationPeriod <- factor(simslagdat$periodsregationPeriod,
                    levels = c("No lags", "5 years", "10 years", "15 years", "30 years"))
                glimpse(simslagdat)
                
                
                a=ggplot(simslagdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="",y="",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Lagged temperatures") +ylim(-95,90)+ scale_x_break(c(10, 80))
                b <- a
                b



                simsfiltmean=apply(sims_lp,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_lp_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsfiltdat$periodsregationPeriod <- c("Unfiltered", "5 years", "10 years", "15 years", "30 years","Unfiltered", "5 years", "10 years", "15 years", "30 years")
                simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    levels = c("Unfiltered", "5 years", "10 years", "15 years", "30 years"))
                glimpse(simsfiltdat)
                
                
                a=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Low-pass filter")+ylim(-95,90)
                a

               plot_60 <-  ggarrange(a,b,c,nrow=1,ncol=3,common.legend=TRUE,legend="none")
               plot_60 <-  annotate_figure(plot_60, top = text_grob("70-year timeseries", 
               color = "black", face = "bold", size = 14))

               ggarrange(plot_60,plot_110,nrow=2,ncol=1,common.legend=TRUE)
               ggsave("SF1_lags_dlnm_filt.png",dpi=300)

            
            
            times=c(65,100,200,500)
            for (time1 in 1:3){
                time <- times[time1]
                basegr=0.01 #1% per year baseline growth rate
                #basegr=0.025 #2.5% USA per year baseline growth rate
                
                start=100
                baseline=rep(basegr,time)
                coef=-0.05 #effect size - change in growth per degree warming()
                growthsd=0.005 #standard deviation of growth variability unexplained by temperature
                #growthsd=0.01 #USA standard deviation of growth variability unexplained by temperature
                periods <- c(0,3,5,10,15)
                arglagbothlist <- list(list(knots=logknots(periods[k],2)),
                    list(fun="poly",degree = 4),
                    list(fun="strata",breaks = 4),
                    list(fun="ns",df = 4),
                    list(fun="bs",df = 4))
                #nsims=176
                nsims=200
                sims=array(dim=c(nsims,length(periods),2))
                sims_lp=array(dim=c(nsims,length(periods),2))
                sims_lp_se=array(dim=c(nsims,length(periods),2))
                sims_se=array(dim=c(nsims,length(periods),2))
                sims_dlnm=array(dim=c(nsims,length(periods),2))
                sims_dlnm_se=array(dim=c(nsims,length(periods),2))
                for(i in 1:nsims){
                    randomtemp=Re(randomts(gtemp))[1:time]
                    randomgrowth_g=basegr #growth impact model
                    randomgrowth_l=basegr #levels impact model
                    for(j in 6:time){
                    randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                    randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                    #randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-0.25*randomtemp[j-1]-0.25*randomtemp[j-2]-
                    #  0.25*randomtemp[j-3]-0.25*randomtemp[j-4]-0.25*randomtemp[j-5])*coef+rnorm(1,mean=0,sd=growthsd))
                    }
                    dataset=data.frame(years=1:(time-5),temp=randomtemp[6:(time)],
                    g=randomgrowth_g[2:(time-4)],l=randomgrowth_l[2:(time-4)])
                    #dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                    dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                    for(k in 1:length(periods)){
                        mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                        mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                    if (k==1){
                            m_g <- lm(g~temp+templag,data=dataset)
                            m_l <- lm(l~temp+templag,data=dataset)
                            mod_g <- m_g$coefficients[2]
                            mod_l <- m_l$coefficients[2]
                            sims[i,k,]=c(mod_g,mod_l)
                            se_l <- summary(m_l)$coefficients[2,2]
                            se_g <- summary(m_g)$coefficients[2,2]
                            sims_se[i,k,]=c(se_g,se_l)
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            sims_lp[i,k,]=c(mod_g,mod_l)

                            sims_dlnm[i,k,] <- c(mod_g,mod_l)
                            sims_dlnm_se[i,k,]=c(se_g,se_l)

                        }else{
                            tempts <- pass.filt(randomtemp[6:(time)], W=periods[k], type="low", method="Butterworth")
                            ratio <- median(randomtemp[6:(time)]/tempts)
                            growth <- randomgrowth_g[2:(time-4)]#[6:time]
                            level <- randomgrowth_l[2:(time-4)]#[6:time]
                            filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                                level = unclass(level))
                            filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                            names(filt) <- c("temp","growth","levels","templag")
                            mod_gfilt=lm(growth~temp,data=filt)
                            mod_lfilt=lm(levels~temp,data=filt)
                            sims_lp[i,k,]=c(mod_gfilt$coef[2]/ratio,mod_lfilt$coef[2]/ratio)
                            se_l <- summary(mod_lfilt)$coefficients[2,2]
                            se_g <- summary(mod_gfilt)$coefficients[2,2]
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            
                            num_lags <- c(0,periods[k])
                            arglagboth <- list(fun="ns",df = 3)
                            cbtemp <- crossbasis(dataset$temp,
                            lag=num_lags,
                            argvar=list("poly",degree=1,scale=1),
                            arglag=arglagboth)
                            #model_g <-  felm((randomgrowth_g[2:(time-4)]) ~ cbtemp)
                            model_g <-  lm(dataset$g ~ cbtemp )
                            redvar_g <- crossreduce(cbtemp,model_g,cen=0)
                            #plot(redvar_g)
                            #summary(redvar_g)
                            beta.hat_g <- redvar_g$coefficients
                            sigma_g <- redvar_g$vcov
                            #model_l <-  felm((randomgrowth_l[2:(time-4)]) ~ cbtemp)
                            model_l <-  felm(dataset$l ~ cbtemp)
                            redvar_l <- crossreduce(cbtemp,model_l,cen=0)
                            beta.hat_l <- redvar_l$coefficients
                            sigma_l <- redvar_l$vcov
                            dynlag_l <-  crosspred(cbtemp, model_l,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            dynlag_g <-  crosspred(cbtemp, model_g,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            
                            sims_dlnm[i,k,] <- c(dynlag_g$allfit[[2]],beta.hat_l[[1]])
                            sims_dlnm_se[i,k,]=c(dynlag_g$allse[[2]],sigma_l[[1]]^0.5)

                            if (k>1){
                                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                                dataset$templag2 =c(NA,dataset$templag[1:(dim(dataset)[1]-1)])
                                dataset$templag3 =c(NA,dataset$templag2[1:(dim(dataset)[1]-1)])
                                m_g <- lm(g~temp+templag+templag2,data=dataset)
                                m_l <- lm(l~temp+templag+templag2,data=dataset)
                                mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                sims[i,k,]=c(mod_g,mod_l)
                                if (k>2) {
                                    dataset$templag4 =c(NA,dataset$templag3[1:(dim(dataset)[1]-1)])
                                    dataset$templag5 =c(NA,dataset$templag4[1:(dim(dataset)[1]-1)])
                                    m_g <- lm(g~temp+templag+templag2+templag3+templag4,data=dataset)
                                    m_l <- lm(l~temp+templag+templag2+templag3+templag4,data=dataset)
                                    mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                    mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                    if (k>3){
                                        dataset$templag6 =c(NA,dataset$templag5[1:(dim(dataset)[1]-1)])
                                        dataset$templag7 =c(NA,dataset$templag6[1:(dim(dataset)[1]-1)])
                                        dataset$templag8 =c(NA,dataset$templag7[1:(dim(dataset)[1]-1)])
                                        dataset$templag9 =c(NA,dataset$templag8[1:(dim(dataset)[1]-1)])
                                        dataset$templag10 =c(NA,dataset$templag9[1:(dim(dataset)[1]-1)])
                                        m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9,data=dataset)
                                        m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9,data=dataset)

                                        if (k>4){
                                            dataset$templag11 =c(NA,dataset$templag10[1:(dim(dataset)[1]-1)])
                                            dataset$templag12 =c(NA,dataset$templag11[1:(dim(dataset)[1]-1)])
                                            dataset$templag13 =c(NA,dataset$templag12[1:(dim(dataset)[1]-1)])
                                            dataset$templag14 =c(NA,dataset$templag13[1:(dim(dataset)[1]-1)])
                                            dataset$templag15 =c(NA,dataset$templag14[1:(dim(dataset)[1]-1)])
                                            m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+templag11+templag12+templag13+templag14,data=dataset)
                                            m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+templag11+templag12+templag13+templag14,data=dataset)
                                        }
                                    }
                                    
                                }
                                
                                


                            }
                            mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                            mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                            sigma <- vcov(m_l)
                            a <- rep(1,length(m_l$coefficients))
                            a[1] <- 0
                            se_sum_l <- c(sqrt(t(a) %*% sigma %*% a))
                            sigma <- vcov(m_g)
                            se_sum_g <- c(sqrt(t(a) %*% sigma %*% a))
                            sims_se[i,k,]=c(se_sum_g,se_sum_l)
                            sims[i,k,]=c(mod_g,mod_l)
                    }
                    

                    }
                    if(i%%50==0) print(i)
                }


                ## Plotting Figure 2  - Simulation exercise (start)
                
                
                simsfiltmean=apply(sims_dlnm,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_dlnm_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsdlnmdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsdlnmdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsdlnmdat$periodsregationPeriod <- c("No lags", "3 years", "5 years", "10 years", "15 years","No lags", "3 years", "5 years", "10 years", "15 years")
                simsdlnmdat$periodsregationPeriod <- factor(simsdlnmdat$periodsregationPeriod,
                    levels = c("No lags", "3 years", "5 years", "10 years", "15 years"))
                glimpse(simsdlnmdat)

                
                
                a=ggplot(simsdlnmdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Number of Lags",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Distributed lag non-linear model")+ylim(-20,15)
                c <- a
                c            





                simslagmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean)=c("Growth","Level");rownames(simslagmean)=periods
                simslagmean_se=apply(sims_se,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean_se)=c("Growth","Level");rownames(simslagmean_se)=periods
                simslagdat=cbind(melt(simslagmean),melt(simslagmean_se)[,3])
                colnames(simslagdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simslagdat$periodsregationPeriod <- c("No lags", "3 years", "5 years", "10 years", "15 years","No lags", "3 years", "5 years", "10 years", "15 years")
                simslagdat$periodsregationPeriod <- factor(simslagdat$periodsregationPeriod,
                    levels = c("No lags", "3 years", "5 years", "10 years", "15 years"))
                glimpse(simslagdat)
                
                
                a=ggplot(simslagdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Number of lags",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Lagged temperatures") +ylim(-20,15)
                b <- a



                simsfiltmean=apply(sims_lp,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_lp_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                glimpse(simsfiltdat)
                
                
                a=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Filters",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Low-pass filter") +ylim(-20,15)
                a

                if (time1==1){
                    
                    compdlnm <- cbind(rbind(cbind(simsfiltdat,method="filter"),cbind(simslagdat,method="lags"),cbind(simsdlnmdat,method="dlnm")),time)
                }else{
                    temporal_compdlnm <- cbind(rbind(cbind(simsfiltdat,method="filter"),cbind(simslagdat,method="lags"),cbind(simsdlnmdat,method="dlnm")),time)
                    compdlnm <- rbind(compdlnm,temporal_compdlnm)

                }
            }

            ggplot(compdlnm[which(compdlnm$periodsregationPeriod=="15 years"),],
            aes(x=time,y=SDEffect,color=method,shape=ImpactType))+
            geom_point()+xlab("length of simulated timeseries") + ylab("Standard deviation")
            #ggsave("Figures/Methods_Comparison_increasinglogknots.png",dpi=600)

            times=c(120,300,500)
            for (time1 in 1:3){
                time <- times[time1]
                basegr=0.01 #1% per year baseline growth rate
                #basegr=0.025 #2.5% USA per year baseline growth rate
                
                start=100
                baseline=rep(basegr,time)
                coef=-0.05 #effect size - change in growth per degree warming()
                growthsd=0.005 #standard deviation of growth variability unexplained by temperature
                #growthsd=0.01 #USA standard deviation of growth variability unexplained by temperature
                periods <- c(0,3,5,10,50)
                arglagbothlist <- list(list(knots=logknots(periods[k],2)),
                    list(fun="poly",degree = 4),
                    list(fun="strata",breaks = 4),
                    list(fun="ns",df = 4),
                    list(fun="bs",df = 4))
                #nsims=176
                nsims=200
                sims=array(dim=c(nsims,length(periods),2))
                sims_lp=array(dim=c(nsims,length(periods),2))
                sims_lp_se=array(dim=c(nsims,length(periods),2))
                sims_se=array(dim=c(nsims,length(periods),2))
                sims_dlnm=array(dim=c(nsims,length(periods),2))
                sims_dlnm_se=array(dim=c(nsims,length(periods),2))
                for(i in 1:nsims){
                    randomtemp=Re(randomts(gtemp))[1:time]
                    randomgrowth_g=basegr #growth impact model
                    randomgrowth_l=basegr #levels impact model
                    for(j in 6:time){
                    randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                    randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                    #randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-0.25*randomtemp[j-1]-0.25*randomtemp[j-2]-
                    #  0.25*randomtemp[j-3]-0.25*randomtemp[j-4]-0.25*randomtemp[j-5])*coef+rnorm(1,mean=0,sd=growthsd))
                    }
                    dataset=data.frame(years=1:(time-5),temp=randomtemp[6:(time)],
                    g=randomgrowth_g[2:(time-4)],l=randomgrowth_l[2:(time-4)])
                    #dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                    dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                    for(k in 1:length(periods)){
                        mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                        mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                    if (k==1){
                            m_g <- lm(g~temp+templag,data=dataset)
                            m_l <- lm(l~temp+templag,data=dataset)
                            mod_g <- m_g$coefficients[2]
                            mod_l <- m_l$coefficients[2]
                            sims[i,k,]=c(mod_g,mod_l)
                            se_l <- summary(m_l)$coefficients[2,2]
                            se_g <- summary(m_g)$coefficients[2,2]
                            sims_se[i,k,]=c(se_g,se_l)
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            sims_lp[i,k,]=c(mod_g,mod_l)

                            sims_dlnm[i,k,] <- c(mod_g,mod_l)
                            sims_dlnm_se[i,k,]=c(se_g,se_l)

                        }else{
                            tempts <- pass.filt(randomtemp[6:(time)], W=periods[k], type="low", method="Butterworth")
                            ratio <- median(randomtemp[6:(time)]/tempts)
                            growth <- randomgrowth_g[2:(time-4)]#[6:time]
                            level <- randomgrowth_l[2:(time-4)]#[6:time]
                            filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                                level = unclass(level))
                            filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                            names(filt) <- c("temp","growth","levels","templag")
                            mod_gfilt=lm(growth~temp,data=filt)
                            mod_lfilt=lm(levels~temp,data=filt)
                            sims_lp[i,k,]=c(mod_gfilt$coef[2]/ratio,mod_lfilt$coef[2]/ratio)
                            se_l <- summary(mod_lfilt)$coefficients[2,2]
                            se_g <- summary(mod_gfilt)$coefficients[2,2]
                            sims_lp_se[i,k,]=c(se_g,se_l)
                            
                            num_lags <- c(0,periods[k])
                            arglagboth <- list(fun="ns",df = 3)
                            cbtemp <- crossbasis(dataset$temp,
                            lag=num_lags,
                            argvar=list("poly",degree=1,scale=1),
                            arglag=arglagboth)
                            #model_g <-  felm((randomgrowth_g[2:(time-4)]) ~ cbtemp)
                            model_g <-  lm(dataset$g ~ cbtemp )
                            redvar_g <- crossreduce(cbtemp,model_g,cen=0)
                            #plot(redvar_g)
                            #summary(redvar_g)
                            beta.hat_g <- redvar_g$coefficients
                            sigma_g <- redvar_g$vcov
                            #model_l <-  felm((randomgrowth_l[2:(time-4)]) ~ cbtemp)
                            model_l <-  felm(dataset$l ~ cbtemp)
                            redvar_l <- crossreduce(cbtemp,model_l,cen=0)
                            beta.hat_l <- redvar_l$coefficients
                            sigma_l <- redvar_l$vcov
                            dynlag_l <-  crosspred(cbtemp, model_l,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            dynlag_g <-  crosspred(cbtemp, model_g,from=-1,to=1,by=1,cumul=TRUE, cen = 0)
                            
                            sims_dlnm[i,k,] <- c(dynlag_g$allfit[[2]],beta.hat_l[[1]])
                            sims_dlnm_se[i,k,]=c(dynlag_g$allse[[2]],sigma_l[[1]]^0.5)

                            if (k>1){
                                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                                dataset$templag2 =c(NA,dataset$templag[1:(dim(dataset)[1]-1)])
                                dataset$templag3 =c(NA,dataset$templag2[1:(dim(dataset)[1]-1)])
                                m_g <- lm(g~temp+templag+templag2,data=dataset)
                                m_l <- lm(l~temp+templag+templag2,data=dataset)
                                mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                sims[i,k,]=c(mod_g,mod_l)
                                if (k>2) {
                                    dataset$templag4 =c(NA,dataset$templag3[1:(dim(dataset)[1]-1)])
                                    dataset$templag5 =c(NA,dataset$templag4[1:(dim(dataset)[1]-1)])
                                    m_g <- lm(g~temp+templag+templag2+templag3+templag4,data=dataset)
                                    m_l <- lm(l~temp+templag+templag2+templag3+templag4,data=dataset)
                                    mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                                    mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                                    if (k>3){
                                        dataset$templag6 =c(NA,dataset$templag5[1:(dim(dataset)[1]-1)])
                                        dataset$templag7 =c(NA,dataset$templag6[1:(dim(dataset)[1]-1)])
                                        dataset$templag8 =c(NA,dataset$templag7[1:(dim(dataset)[1]-1)])
                                        dataset$templag9 =c(NA,dataset$templag8[1:(dim(dataset)[1]-1)])
                                        dataset$templag10 =c(NA,dataset$templag9[1:(dim(dataset)[1]-1)])
                                        m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9,data=dataset)
                                        m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9,data=dataset)

                                        if (k>4){
                                            dataset$templag11 =c(NA,dataset$templag10[1:(dim(dataset)[1]-1)])
                                            dataset$templag12 =c(NA,dataset$templag11[1:(dim(dataset)[1]-1)])
                                            dataset$templag13 =c(NA,dataset$templag12[1:(dim(dataset)[1]-1)])
                                            dataset$templag14 =c(NA,dataset$templag13[1:(dim(dataset)[1]-1)])
                                            dataset$templag15 =c(NA,dataset$templag14[1:(dim(dataset)[1]-1)])
                                            dataset$templag16 =c(NA,dataset$templag15[1:(dim(dataset)[1]-1)])
                                            dataset$templag17 =c(NA,dataset$templag16[1:(dim(dataset)[1]-1)])
                                            dataset$templag18 =c(NA,dataset$templag17[1:(dim(dataset)[1]-1)])
                                            dataset$templag19 =c(NA,dataset$templag18[1:(dim(dataset)[1]-1)])
                                            dataset$templag20 =c(NA,dataset$templag19[1:(dim(dataset)[1]-1)])
                                            dataset$templag21 =c(NA,dataset$templag20[1:(dim(dataset)[1]-1)])
                                            dataset$templag22 =c(NA,dataset$templag21[1:(dim(dataset)[1]-1)])
                                            dataset$templag23 =c(NA,dataset$templag22[1:(dim(dataset)[1]-1)])
                                            dataset$templag24 =c(NA,dataset$templag23[1:(dim(dataset)[1]-1)])
                                            dataset$templag25 =c(NA,dataset$templag24[1:(dim(dataset)[1]-1)])
                                            dataset$templag26 =c(NA,dataset$templag25[1:(dim(dataset)[1]-1)])
                                            dataset$templag27 =c(NA,dataset$templag26[1:(dim(dataset)[1]-1)])
                                            dataset$templag28 =c(NA,dataset$templag27[1:(dim(dataset)[1]-1)])
                                            dataset$templag29 =c(NA,dataset$templag28[1:(dim(dataset)[1]-1)])
                                            dataset$templag30 =c(NA,dataset$templag29[1:(dim(dataset)[1]-1)])
                                            dataset$templag31 =c(NA,dataset$templag30[1:(dim(dataset)[1]-1)])
                                            dataset$templag32 =c(NA,dataset$templag31[1:(dim(dataset)[1]-1)])
                                            dataset$templag33 =c(NA,dataset$templag32[1:(dim(dataset)[1]-1)])
                                            dataset$templag34 =c(NA,dataset$templag33[1:(dim(dataset)[1]-1)])
                                            dataset$templag35 =c(NA,dataset$templag34[1:(dim(dataset)[1]-1)])
                                            dataset$templag36 =c(NA,dataset$templag35[1:(dim(dataset)[1]-1)])
                                            dataset$templag37 =c(NA,dataset$templag36[1:(dim(dataset)[1]-1)])
                                            dataset$templag38 =c(NA,dataset$templag37[1:(dim(dataset)[1]-1)])
                                            dataset$templag39 =c(NA,dataset$templag38[1:(dim(dataset)[1]-1)])
                                            dataset$templag40 =c(NA,dataset$templag39[1:(dim(dataset)[1]-1)])
                                            dataset$templag41 =c(NA,dataset$templag40[1:(dim(dataset)[1]-1)])
                                            dataset$templag42 =c(NA,dataset$templag41[1:(dim(dataset)[1]-1)])
                                            dataset$templag43 =c(NA,dataset$templag42[1:(dim(dataset)[1]-1)])
                                            dataset$templag44 =c(NA,dataset$templag43[1:(dim(dataset)[1]-1)])
                                            dataset$templag45 =c(NA,dataset$templag44[1:(dim(dataset)[1]-1)])
                                            dataset$templag46 =c(NA,dataset$templag45[1:(dim(dataset)[1]-1)])
                                            dataset$templag47 =c(NA,dataset$templag46[1:(dim(dataset)[1]-1)])
                                            dataset$templag48 =c(NA,dataset$templag47[1:(dim(dataset)[1]-1)])
                                            dataset$templag49 =c(NA,dataset$templag48[1:(dim(dataset)[1]-1)])
                                            dataset$templag50 =c(NA,dataset$templag49[1:(dim(dataset)[1]-1)])
                                            m_g <- lm(g~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15+templag16+templag17+templag18+templag19+templag20+
                                            templag21+templag22+templag23+templag24+templag25+templag26+templag27+templag28+templag29+templag30+
                                            templag31+templag32+templag33+templag34+templag35+templag36+templag37+templag38+templag39+templag40+
                                            templag41+templag42+templag43+templag44+templag45+templag46+templag47+templag48+templag49+templag50,data=dataset)
                                            m_l <- lm(l~temp+templag+templag2+templag3+templag4+templag5+templag6+templag7+templag8+templag9+templag10+
                                            templag11+templag12+templag13+templag14+templag15+templag16+templag17+templag18+templag19+templag20+
                                            templag21+templag22+templag23+templag24+templag25+templag26+templag27+templag28+templag29+templag30+
                                            templag31+templag32+templag33+templag34+templag35+templag36+templag37+templag38+templag39+templag40+
                                            templag41+templag42+templag43+templag44+templag45+templag46+templag47+templag48+templag49+templag50,data=dataset)
                                        }
                                    }
                                    
                                }
                                
                                


                            }
                            mod_g <- sum(m_g$coefficients[2:length(m_g$coefficients)])
                            mod_l <- sum(m_l$coefficients[2:length(m_l$coefficients)])
                            sigma <- vcov(m_l)
                            a <- rep(1,length(m_l$coefficients))
                            a[1] <- 0
                            se_sum_l <- c(sqrt(t(a) %*% sigma %*% a))
                            sigma <- vcov(m_g)
                            se_sum_g <- c(sqrt(t(a) %*% sigma %*% a))
                            sims_se[i,k,]=c(se_sum_g,se_sum_l)
                            sims[i,k,]=c(mod_g,mod_l)
                    }
                    

                    }
                    if(i%%50==0) print(i)
                }


                ## Plotting Figure 2  - Simulation exercise (start)
                
                
                simsfiltmean=apply(sims_dlnm,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_dlnm_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsdlnmdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsdlnmdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsdlnmdat$periodsregationPeriod <- c("No lags", "3 years", "5 years", "10 years", "15 years","No lags", "3 years", "5 years", "10 years", "15 years")
                simsdlnmdat$periodsregationPeriod <- factor(simsdlnmdat$periodsregationPeriod,
                    levels = c("No lags", "3 years", "5 years", "10 years", "15 years"))
                glimpse(simsdlnmdat)

                
                
                a=ggplot(simsdlnmdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Number of Lags",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Distributed lag non-linear model")+ylim(-20,15)
                c <- a
                c            





                simslagmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean)=c("Growth","Level");rownames(simslagmean)=periods
                simslagmean_se=apply(sims_se,MARGIN=c(2,3),FUN=mean)
                colnames(simslagmean_se)=c("Growth","Level");rownames(simslagmean_se)=periods
                simslagdat=cbind(melt(simslagmean),melt(simslagmean_se)[,3])
                colnames(simslagdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simslagdat$periodsregationPeriod <- c("No lags", "3 years", "5 years", "10 years", "15 years","No lags", "3 years", "5 years", "10 years", "15 years")
                simslagdat$periodsregationPeriod <- factor(simslagdat$periodsregationPeriod,
                    levels = c("No lags", "3 years", "5 years", "10 years", "15 years"))
                glimpse(simslagdat)
                
                
                a=ggplot(simslagdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Number of lags",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Lagged temperatures") +ylim(-20,15)
                b <- a



                simsfiltmean=apply(sims_lp,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                simsfiltmean_se=apply(sims_lp_se,MARGIN=c(2,3),FUN=mean)
                colnames(simsfiltmean_se)=c("Growth","Level");rownames(simsfiltmean_se)=periods
                simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltmean_se)[,3])
                colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                theme_set(theme_bw(base_size = 20))
                
                simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                glimpse(simsfiltdat)
                
                
                a=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                a=a+labs(x="Filters",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Low-pass filter") +ylim(-20,15)
                a

                if (time1==1){
                    
                    compdlnm <- cbind(rbind(cbind(simsfiltdat,method="filter"),cbind(simslagdat,method="lags"),cbind(simsdlnmdat,method="dlnm")),time)
                }else{
                    temporal_compdlnm <- cbind(rbind(cbind(simsfiltdat,method="filter"),cbind(simslagdat,method="lags"),cbind(simsdlnmdat,method="dlnm")),time)
                    compdlnm <- rbind(compdlnm,temporal_compdlnm)

                }
            }

            ggplot(compdlnm[which(compdlnm$periodsregationPeriod=="15 years"),],
            aes(x=time,y=SDEffect,color=method,shape=ImpactType))+
            geom_point()+xlab("length of simulated timeseries") + ylab("Standard error")
            #ggsave("Figures/Methods_Comparison_increasinglogknots.png",dpi=600)

            #ggarrange(a,b,c, common.legend=TRUE,ncol=3,nrow=1)
            #ggsave("Figures/Filters_Lags_DLNM.png",dpi=600)
             #ggsave("Figures/Filters_Lags_DLNM_1000points_200sims.png",dpi=600)
        ## 2.4c. Simulation comparing low-pass filter with distributed lags
        

    ## 2.4. Additional simulations and figures (Used for Response to Reviewers)


## 2. Simulation - Figures 1 and 2 (end)

## 3. Empirical analysis - Figure 3 and Table 1 (start)

    ## 3.1. DLNM Country-level regressions (start)
        i <- 1
        numlags <- c(3,5,10,15)
        GDPdatasets <- c("wb","barro","mad") 
        Tempdatasets <- c("UDel")
        weights <- c("pop")
        librayr()
        iter <- 1
        for (n in 1:length(GDPdatasets)){
            for (i in 1:length(numlags)){

                    if(i ==1){
                            df_lags_all <- c(2)
                        }else if (i==2){
                            df_lags_all <- c(2,3)
                        } else if(i>2){df_lags_all <- c(2,3,4)}

                    for (dfl in 1:length(df_lags_all)){
                        df_lags <- df_lags_all[dfl]
                        arglagbothlist <- list(list(fun="poly",degree = (df_lags)),
                            list(fun="bs",df = (df_lags)))

                        arglagbothlist_names<- c("poly","b-spline")

                        for (p in 1:length(arglagbothlist)){

                            DATA <- get(GDPdatasets[n])
                            num_lags <- c(0,numlags[i])
                                tempname <- "UDel_pop_temp"
                                preciname <- "UDel_pop_preci"

                            arglagboth <- arglagbothlist[[p]]
                            if(GDPdatasets[n]=="mad"){

                            DATA <- DATA[which(!is.na(DATA$growth) & !is.na(DATA$UDel_pop_temp) & !is.na(DATA$countrycode)),]
                            #DATA <- DATA[which( !is.na(DATA$countrycode)),]
                            #DATA <- setDT(DATA)[, growth := cumsum(c(0, diff(year)) > 1), by = countrycode][, if (.N > 4) .SD, by = .(countrycode, growth)][, growth := NULL][]
                            DATA <- as.data.frame(DATA)

                            # tapply(!is.na(DATA$growth), DATA$countrycode, function(x) {
                            #    max.vec <- rle(x)$length[rle(x)$values==T]
                            #    ifelse(length(max.vec) > 0, max(max.vec), 0)
                            #    })
                            }
                            
                            cbtemp <- crossbasis(DATA[,which(names(DATA)%in%tempname)],
                                lag=num_lags,argvar=list(fun="poly",degree = 2 ,scale= 1), 
                                arglag=arglagboth, group = DATA$countrycode)
                            
                        
                            # Construct Cross-bases (end)

                                model <-  felm(I(growth*100) ~ cbtemp + UDel_pop_preci +  
                                countrycode:year + countrycode:I(year^2)|
                                countrycode + year | 0 | 0, data=DATA)


                            redvar <- crossreduce(cbtemp,model,cen=10)
                    
                            beta.hat <- redvar$coefficients
                            sigma <- redvar$vcov
                            xx <- 1:30
                            xmat <-cbind(1,2*xx)
                            dgdt <- colSums(beta.hat*t(xmat)) 
                            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                            coldata <- cbind(xx,dgdt,ci1,ci2)
                            madcol <- data.frame(coldata)
                            madcol$x<- as.numeric(as.character(madcol$x))
                            madcol$dgdt <- as.numeric(as.character(madcol$dgdt))
                            madcol$ci1 <- as.numeric(as.character(madcol$ci1))
                            madcol$ci2 <- as.numeric(as.character(madcol$ci2))
                            madcol$lagfun <- arglagbothlist_names[p]
                            madcol$dataset <- GDPdatasets[n]
                            madcol$nlags <- numlags[i]
                            madcol$df <- df_lags_all[dfl]
                            madcol$lagfunnum <- paste(arglagbothlist_names[p],df_lags_all[dfl],"df",sep=" ")
                            madcol$ID <- paste(numlags[i]," lags. ",df_lags_all[dfl],arglagbothlist_names[p],sep="")

                            if(iter == 1){
                                marginal_eff <- madcol
                                iter <- 0
                            } else {marginal_eff <- rbind(marginal_eff,madcol)}


                        }

                    }
            }
        }

        glimpse(marginal_eff[marginal_eff$lagfunnum=="2-poly" & marginal_eff$dataset=="wb",])
        
        marginaltemp <- ggplot(data=marginal_eff[which(marginal_eff$dataset=="wb" ),], 
                aes(x=x, y=dgdt, group=ID, color=nlags)) +
                facet_wrap(~nlags)+
                    labs(x="Temp", y="Change in GDP growth rate",
                    title=paste("Marginal effect World Bank"), cex=4) +
                    geom_line(aes(linetype=lagfunnum))+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,fill=nlags),alpha=0.1,colour=NA) +
                    geom_hline(yintercept=0, linetype="dashed", size =1) +
                    theme_minimal() +
                    theme(legend.position="bottom") + 
                    guides(fill=guide_legend(title="Number of lags"),color=guide_legend(title="Number of lags"),
                    linetype=guide_legend(title="Lag structure"))
                    marginaltemp
                    ggsave("Figures/DLNM_WB.png",dpi=300)

            marginaltemp <- ggplot(data=marginal_eff[which(marginal_eff$dataset=="barro" ),], 
                aes(x=x, y=dgdt, group=ID, color=nlags)) +
                facet_wrap(~nlags)+
                    labs(x="Temp", y="Change in GDP growth rate",
                    title=paste("Marginal effect Barro-Ursua"), cex=4) +
                    geom_line(aes(linetype=lagfunnum))+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,fill=nlags),alpha=0.1,colour=NA) +
                    geom_hline(yintercept=0, linetype="dashed", size =1) +
                    theme_minimal() +
                    theme(legend.position="bottom") + 
                    guides(fill=guide_legend(title="Number of lags"),color=guide_legend(title="Number of lags"),
                    linetype=guide_legend(title="Lag structure"))
                    marginaltemp

            marginaltemp <- ggplot(data=marginal_eff[which(marginal_eff$dataset=="mad" ),], 
                aes(x=x, y=dgdt, group=ID, color=nlags)) +
                facet_wrap(~nlags)+
                    labs(x="Temp", y="Change in GDP growth rate",
                    title=paste("Marginal effect Maddison"), cex=4) +
                    geom_line(aes(linetype=lagfunnum))+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,fill=nlags),alpha=0.1,colour=NA) +
                    geom_hline(yintercept=0, linetype="dashed", size =1) +
                    theme_minimal() +
                    theme(legend.position="bottom") + 
                    guides(fill=guide_legend(title="Number of lags"),color=guide_legend(title="Number of lags"),
                    linetype=guide_legend(title="Lag structure"))
                    marginaltemp

            cols <- c("#66c2a5","#fc8d62","#8da0cb") 
            marginaltemp <- ggplot(data=marginal_eff, 
                aes(x=x, y=dgdt, group=interaction(ID,dataset), color=factor(dataset))) +
                facet_wrap(~nlags)+
                    labs(x="Temp", y="Change in GDP growth rate",
                    title=paste("Global cumulative marginal effect"), cex=4) +
                    geom_line(aes(linetype=lagfunnum))+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,fill=factor(dataset)),alpha=0.05,colour=NA) +
                    geom_hline(yintercept=0, linetype="dashed", size =1) +
                    theme_minimal() +
                    theme(legend.position="bottom") + 
                    guides(fill=guide_legend(title="Number of lags"),color=guide_legend(title="Number of lags"),
                    linetype=guide_legend(title="Lag structure"))+
                    scale_color_manual(values=cols)+
                    scale_fill_manual(values=cols)
                    marginaltemp
                    #annotate("text", label = Agg, x = 20, y = max(madcol_plot$dgdt), size = 8, colour = "black")
                    ggsave(marginaltemp,"")

                file=nc_open("C:/Users/bastien/Box/Long Run GDP Growth/Data/LMR Data/air_MCruns_ensemble_mean_LMRv2.0.nc")
                gtemp=ncvar_get(file,"air",start=c(1,1,1,1),count=c(-1,-1,1,-1))
                lon <- ncvar_get(file,"lon")
                lat <- ncvar_get(file,"lat")
                time <- ncvar_get(file,"time")
                tunits <- ncatt_get(file,"time","units")
                data(world)
                isos_all <- unique(world$iso_a2)
                gtemp <- brick(gtemp, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
                firstcountry=1   
                for(j in 1:20){
                    geom_iso <- world$geom[world$iso_a2==isos_all[j]]
                    geom_iso <- st_shift_longitude(geom_iso)
                    geom_iso <- st_cast(geom_iso, "POLYGON")
                    geom_iso <-as_Spatial(geom_iso)
                
                    temp_lmr <- rep(NA,100)  
                    growth_lmr_panel <- rep(NA,100) 
                    growth_lmr_nopanel <- rep(NA,100) 
                        beta <- 0.0100# BHM coefficient effect size - change in growth per degree warming()
                        beta2<- -0.0005 #BHM coefficient
                        gamma <- 0.0027
                        gamma2_a <- -0.0005
                        gamma2_b <- 0   

                        meanT <- mean(bhm$temp[which(bhm$iso==countrycode(isos_all[j], origin="iso2c", destination="iso3c"))],na.rm=TRUE)
                        if(length(meanT)==0){next}
                        
                    for (i in 1:60){
                        gtemp_year <- raster::subset(gtemp,(i+1399))
                        #gtemp_year <- rotate(gtemp_year)
                        gtemp_year <- t(gtemp_year)
                        temp <- extract(gtemp_year, geom_iso,method='simple')
                        temp_lmr[i] <- mean(temp[[1]]) + meanT
                        if(i > 2){
                            growth_lmr_panel[i] <- 0.02 + temp_lmr[i]*gamma + gamma2_a*temp_lmr[i]^2+ temp_lmr[i]*beta + beta2*temp_lmr[i]^2 - temp_lmr[i-1]*beta - beta2*temp_lmr[i-1]^2 + +rnorm(1,mean=0,sd=0.005)
                            growth_lmr_nopanel[i] <- 0.02 + temp_lmr[i]*gamma + gamma2_b*temp_lmr[i]^2+ temp_lmr[i]*beta + beta2*temp_lmr[i]^2 - temp_lmr[i-1]*beta - beta2*temp_lmr[i-1]^2 + +rnorm(1,mean=0,sd=0.005)
                            }
                    }
                    df <- data.frame(temp=temp_lmr,growth_panel=growth_lmr_panel,growth_nopanel=growth_lmr_nopanel,country=isos_all[j])
                    if(firstcountry==1){
                        testing_panel <- df
                        firstcountry <- 0
                    }else{
                        testing_panel <- rbind(testing_panel,df)

                    }
                    print(j)


                    }
                    glimpse(testing_panel)
                    testing_panel$year <- seq(1:100)
                    arglagbothlist <- list(list(fun="poly",degree = (2)))

                    arglagbothlist_names<- c("poly")

                    DATA <- testing_panel
                            num_lags <- c(0,5)
                                tempname <- "temp"

                            arglagboth <- arglagbothlist[[1]]
                            
                            cbtemp <- crossbasis(DATA[,which(names(DATA)%in%tempname)],
                                lag=num_lags,argvar=list(fun="poly",degree = 2 ,scale= 1), 
                                arglag=arglagboth, group = DATA$country)
                            
                        
                            # Construct Cross-bases (end)
                            iter <- 1
                    for(panel in c(0,1)){
                        if(panel==0){
                            model_nop <-  felm(I(growth_nopanel*100) ~ cbtemp +  
                                country:year + country:I(year^2)|
                                country + year | 0 | 0, data=DATA)
                            
                        redvar <- crossreduce(cbtemp,model_nop,cen=10)
                        }else{
                            model_p <-  felm(I(growth_panel*100) ~ cbtemp +  
                                country:year + country:I(year^2)|
                                country + year | 0 | 0, data=DATA)
                        redvar <- crossreduce(cbtemp,model_p,cen=10)
                        }
                    
                            beta.hat <- redvar$coefficients
                            sigma <- redvar$vcov
                            xx <- 1:30
                            xmat <-cbind(1,2*xx)
                            dgdt <- colSums(beta.hat*t(xmat)) 
                            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                            coldata <- cbind(xx,dgdt,ci1,ci2)
                            madcol <- data.frame(coldata)
                            madcol$x<- as.numeric(as.character(madcol$x))
                            madcol$dgdt <- as.numeric(as.character(madcol$dgdt))
                            madcol$ci1 <- as.numeric(as.character(madcol$ci1))
                            madcol$ci2 <- as.numeric(as.character(madcol$ci2))
                            madcol$lagfun <- arglagbothlist_names[1]
                            madcol$dataset <- panel
                            madcol$nlags <- 5
                            madcol$df <- 2
                            madcol$lagfunnum <- paste(arglagbothlist_names[1],"2","df",sep=" ")
                            madcol$ID <- paste("5"," lags. ","2",arglagbothlist_names[1],sep="")

                            if(iter == 1){
                                marginal_eff <- madcol
                                iter <- 0
                            } else {marginal_eff <- rbind(marginal_eff,madcol)}
                    }

                    cols <- c("#66c2a5","#fc8d62","#8da0cb") 
            marginaltemp <- ggplot(data=marginal_eff, 
                aes(x=x, y=dgdt, group=interaction(dataset), color=factor(dataset))) +
                #facet_wrap(~nlags)+
                    labs(x="Temp", y="Change in GDP growth rate",
                    title=paste("Global cumulative marginal effect"), cex=4) +
                    geom_line(aes(linetype=lagfunnum))+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,fill=factor(dataset)),alpha=0.05,colour=NA) +
                    geom_hline(yintercept=0, linetype="dashed", size =1) +
                    theme_minimal() +
                    theme(legend.position="bottom") + 
                    guides(fill=guide_legend(title="Growth effect related \n with avg temp^2"),color=guide_legend(title="Growth effect related \n with avg temp^2"),
                    linetype=guide_legend(title="Lag structure"))+
                    scale_color_manual(values=cols)+
                    scale_fill_manual(values=cols)
                    marginaltemp
                    
                    plagged <- crosspred(cbtemp, model_p,from=0,to=29,by=1,cumul=TRUE)
                    relmax <- which(plagged$matfit[,1]==max(plagged$matfit[,1]))
                    plagged <- crosspred(cbtemp, model_p,from=0,to=29,by=1,cumul=TRUE,cen=relmax)
                    x2=as.numeric(rownames(plagged$matfit))
                    plot(x2,plagged$matfit[1:length(x2),1],type="l",lwd=2,col="#25a2a6",
                    xlab="Population-Weighted Annual Avereage Temperature (degrees C)",
                    ylab="GDP per-capita Growth Rate\n(pp Difference from Max)",main="Contemporaneous effect: Growth effect related \n with avg temp^2",bty="n")
                    abline(h=0,lwd=1.75)
                    polygon(c(x2,rev(x2)),c(plagged$matlow[1:length(x2),1],rev(plagged$mathigh[1:length(x2),1])),col=rgb(t(col2rgb("#25a2a6")),max=255,alpha=70),border=NA)
                    

                    plagged <- crosspred(cbtemp, model_nop,from=0,to=29,by=1,cumul=TRUE)
                    relmax <- which(plagged$matfit[,1]==max(plagged$matfit[,1]))
                    plagged <- crosspred(cbtemp, model_nop,from=0,to=29,by=1,cumul=TRUE,cen=relmax)
                    x2=as.numeric(rownames(plagged$matfit))
                    plot(x2,plagged$matfit[1:length(x2),1],type="l",lwd=2,col="#25a2a6",
                    xlab="Population-Weighted Annual Avereage Temperature (degrees C)",
                    ylab="GDP per-capita Growth Rate\n(pp Difference from Max)",main="Contemporaneous effect: Growth effect NOT related \n with avg temp^2",bty="n")
                    abline(h=0,lwd=1.75)
                    polygon(c(x2,rev(x2)),c(plagged$matlow[1:length(x2),1],rev(plagged$mathigh[1:length(x2),1])),col=rgb(t(col2rgb("#25a2a6")),max=255,alpha=70),border=NA)
                    
            
    ## 3.1. DLNM Country-level regressions (end)

    ## 3.1. LAGS Country-level regressions (start)
        dataset <- c("wb","barro","mad") 
        datasetweather <- c("UDel")
        periods <- c(0,3,5,10,15)
        countries <- unique(wb$countrycode)
        #periods <- c(0,10,20,25,27) Uncomment to get Supp Fig 2
        fullmods_filter=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_var=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_ar=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_p=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        panel_data <- data.frame(years = integer(),temp = double(), growth = double(), 
            preci = double(), countrycode = factor(), climdata = character(),
            econdata = character(), filter = character(), meant = double(), meanp=double(), growth_wdi=double())
        for (mm in 1:length(datasetweather)){
            tempname <- paste(datasetweather[mm],"_pop_temp",sep="")
            preciname <- paste(datasetweather[mm],"_pop_preci",sep="")
            for (jj in (1:length(dataset))){
                DATA <- get(dataset[jj])
                countries=unique(factor(DATA$countrycode))
                #fullmods=array(dim=c(length(countries),2,length(periods)))
                for(k in 1:length(periods)){
                for(i in 1:length(countries)){
                    dat=DATA[which(DATA$countrycode==countries[i]),
                        which(colnames(DATA)%in%c("countrycode","year","growth",tempname,preciname))] 
                    dat <- dat[is.finite(dat$growth),]

                        gonext <- 0
                    for (nn in (1:(dim(dat)[2]))){
                        if(sum(is.na(dat[,nn]))==dim(dat)[1]){
                            #print(paste(countries[i],i,"missing complete column"))
                            gonext <- 1 }
                        
                    }
                    if(gonext==1){
                        gonext <- 0
                        next}
                    if(sum(complete.cases(dat))<2){
                            #print(paste(countries[i],i,"not complete cases"))
                            next
                    }
                    dat <- dat %>% rename("temp"=tempname,"preci"=preciname)
                    dat <- dat[complete.cases(dat),]
                    if(max(dat$year)-min(dat$year) != (dim(dat)[1]-1)){
                        #print(paste("Missing observation.",dataset[jj]))
                        next

                    }
                    #dat <- dat[-c(1),] #first observation is missing
                    if(dim(dat)[1]<5){next}
                    mt <- lm(temp~year+I(year^2), data = dat)
                    meanT <- mean(dat$temp, na.rm = TRUE)
                    t <- resid(mt)    
                    mp <- lm(preci~year+I(year^2), data = dat)
                    p <- resid(mp)
                    meanP <- mean(dat$preci, na.rm = TRUE)
                    mg <- lm(growth~year+I(year^2), data = dat)
                    g <- resid(mg)
                    growth_wdi <- dat$growth
                    if(length(t)<(3+2*(periods[k]))){
                            #print(paste(countries[i],i,"not enough data for this filter","length=",length(t),"periods=",periods[k]))
                            next}
                    
                    tempts <- t
                    precits <- p
                    ratio <- 1
                    temp <- data.frame(year = dat$year, temp= unclass(tempts))
                    preci <- data.frame(year = dat$year, preci= unclass(precits))
                    growth <- data.frame(year = dat$year, growth = unclass(g))
                    #growth_wdi <- data.frame(year = dat$year, growth = unclass(dat$growth))
                    filterdata <- merge(temp,growth, by = "year")
                    filterdata <- merge(filterdata,preci, by = "year")
                    #filterdata <- merge(filterdata,growth_wdi, by = "year")
                    names(filterdata) <- c("years","temp","growth","preci")

                    if(k==1){
                        filterdata2 <- filterdata
                        filterdata2$templag =c(NA,filterdata2$temp[1:(dim(filterdata2)[1]-1)])
                        mod_gfilterdata=lm(growth~temp+preci+templag,data=filterdata2)
                        se <- summary(mod_gfilterdata)$coefficients[2,2]
                        coefl <- sum(summary(mod_gfilterdata)$coefficients[2,1])
                        vcov_ <-vcov(mod_gfilterdata)
                    } else if(k==2){
                        filterdata2 <- filterdata
                        filterdata2$templag =c(NA,filterdata2$temp[1:(dim(filterdata2)[1]-1)])
                        filterdata2$templag2 =c(NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-2)])
                        filterdata2$templag3 =c(NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-3)])
                        mod_gfilterdata=lm(growth~temp+templag+templag2+templag3+preci,data=filterdata2)
                        a <- c(0,1,1,1,1,0)
                        se <- c(sqrt(t(a) %*% vcov(mod_gfilterdata) %*% a))
                        coefl <- sum(summary(mod_gfilterdata)$coefficients[2:(periods[k]+1),1])
                        
                    }else if(k==3){
                        filterdata2 <- filterdata
                        filterdata2$templag =c(NA,filterdata2$temp[1:(dim(filterdata2)[1]-1)])
                        filterdata2$templag2 =c(NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-2)])
                        filterdata2$templag3 =c(NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-3)])
                        filterdata2$templag4 =c(NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-4)])
                        filterdata2$templag5 =c(NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-5)])
                        mod_gfilterdata=lm(growth~temp+templag+templag2+templag3+templag4+templag5+preci,data=filterdata2)
                        a <- c(0,1,1,1,1,1,1,0)
                        se <- c(sqrt(t(a) %*% vcov(mod_gfilterdata) %*% a))
                        coefl <- sum(summary(mod_gfilterdata)$coefficients[2:(periods[k]+1),1])
                        
                    }else if(k==4){
                        filterdata2 <- filterdata
                        filterdata2$templag =c(NA,filterdata2$temp[1:(dim(filterdata2)[1]-1)])
                        filterdata2$templag2 =c(NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-2)])
                        filterdata2$templag3 =c(NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-3)])
                        filterdata2$templag4 =c(NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-4)])
                        filterdata2$templag5 =c(NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-5)])
                        filterdata2$templag6 =c(NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-6)])
                        filterdata2$templag7 =c(NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-7)])
                        filterdata2$templag8 =c(NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-8)])
                        filterdata2$templag9 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-9)])
                        filterdata2$templag10 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-10)])
                        mod_gfilterdata=lm(growth~temp+templag+templag2+templag3+templag4+templag5+
                            templag6+templag7+templag8+templag9+templag10+preci,data=filterdata2)
                        a <- c(0,1,1,1,1,1,1,1,1,1,1,1,0)
                        se <- c(sqrt(t(a) %*% vcov(mod_gfilterdata) %*% a))
                        coefl <- sum(summary(mod_gfilterdata)$coefficients[2:(periods[k]+1),1])
                        
                    }else if(k==5){
                        filterdata2 <- filterdata
                        filterdata2$templag =c(NA,filterdata2$temp[1:(dim(filterdata2)[1]-1)])
                        filterdata2$templag2 =c(NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-2)])
                        filterdata2$templag3 =c(NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-3)])
                        filterdata2$templag4 =c(NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-4)])
                        filterdata2$templag5 =c(NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-5)])
                        filterdata2$templag6 =c(NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-6)])
                        filterdata2$templag7 =c(NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-7)])
                        filterdata2$templag8 =c(NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-8)])
                        filterdata2$templag9 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-9)])
                        filterdata2$templag10 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-10)])
                        filterdata2$templag11 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-11)])
                        filterdata2$templag12 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-12)])
                        filterdata2$templag13 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-13)])
                        filterdata2$templag14 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-14)])
                        filterdata2$templag15 =c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,filterdata2$temp[1:(dim(filterdata2)[1]-15)])
                        mod_gfilterdata=lm(growth~temp+templag+templag2+templag3+templag4+templag5+
                            templag6+templag7+templag8+templag9+templag10+
                            templag11+templag12+templag13+templag14+templag15+preci,data=filterdata2)
                        a <- c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0)
                        se <- c(sqrt(t(a) %*% vcov(mod_gfilterdata) %*% a))
                        coefl <- sum(summary(mod_gfilterdata)$coefficients[2:(periods[k]+1),1])
                        
                    }
                    
                    fullmods_filter[i,,k,jj,mm]=c(coefl,se)
                    
                    #filterdata$years <- dat$year
                    filterdata$countrycode <- rep(dat$countrycode[1],dim(filterdata)[1])
                    filterdata$climdata <- rep(datasetweather[mm],dim(filterdata)[1])
                    filterdata$econdata <- rep(dataset[jj],dim(filterdata)[1])
                    filterdata$filter <- rep(paste(periods[k],sep="-"),dim(filterdata)[1])
                    filterdata$meant <- rep(meanT,dim(filterdata)[1])
                    filterdata$meanp <- rep(meanP,dim(filterdata)[1])
                    filterdata$growth_wdi <- dat$growth 
                    panel_data <- bind_rows(panel_data,filterdata)
                    #fullmods_filter_var[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[2,1]/ratio,vcov_[2,2]/(ratio^2))
                    #fullmods_filter_ar[i,,k,jj,mm]=c(dw$statistic[1],dw$p.value[1])
                    #fullmods_filter_p[i,,k,jj,mm]=summary(mod_gfilterdata)$coefficients[3,1:2] #Not HAC s.e. If usd should be divided by the ratio of precipietation timeseries
                    }
                    
                
                }
            }
        }
        original_fullmods_filter <- fullmods_filter
        original_fullmods_filter_var <- fullmods_filter_var
        glimpse(fullmods_filter)
        
        
        
        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filterm <- melt(fullmods_filter)
        glimpse(fullmods_filterm)
        #
        names(fullmods_filterm) <- c("countrycode","variable","Lags","econdata","climdata","value")
        new_fft_l <- fullmods_filterm[order(fullmods_filterm[,2], fullmods_filterm[,3] ),] # Sort by column index [1] then [3]
        glimpse(new_fft_l)
        estimates <- new_fft_l$value[which(new_fft_l$variable=="Estimate")]        
        ses <- new_fft_l$value[which(new_fft_l$variable=="StandardError")]
        new_fft_l <- new_fft_l[c(1:length(estimates)),-c(2,6)]    
        new_fft_l$Estimate <- estimates
        new_fft_l$StandardError <- ses   
        glimpse(new_fft_l) 
        new_fft_l <- new_fft_l[order(new_fft_l[,1]),] # Sort by column index [1] then [3]     
        fullmods_filter_l <- new_fft_l
        data1 <- fullmods_filter_l[fullmods_filter_l$econdata=="wb" & fullmods_filter$climdata=="UDel",]
        glimpse(data1)
        data1[1:10,]
        data1[data1$Lags==10,]

        ggplot(data=data1)+
        geom_line(aes(x=Lags,y=Estimate,group=countrycode))
        #
        # countriesbarro <- unique(factor(barro$countrycode))
        # countries_barro <- unlist(lapply(countriesbarro , as.character))
        # countries_barro_extended <- c(countries_barro,countries_wb[(length(countries_barro)+1):length(countries_wb)*NA])
        # countriesmad <- unique(factor(mad$countrycode))
        # countries_mad <- unlist(lapply(countriesmad , as.character))
        # countries_mad_extended <- c(countries_mad,countries_wb[(length(countries_mad)+1):length(countries_wb)*NA])
        # countries_wbrep <- rep(countries_wb,2*5)
        # countries_barrorep <- rep(countries_barro_extended,2*5)
        # countries_madrep <- rep(countries_mad_extended,2*5) #2*5
        # names1 <- c(countries_wbrep,countries_barrorep,countries_madrep)
        # #names2 <- c(names1,names1)
        # names2 <- names1 #only usiong UDel and not LMR
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter$Estimate[which(is.infinite(fullmods_filter$StandardError))]=NA
        # names(fullmods_filter) <- c("countrycode","frequencies","econdata","climdata","Estimate","StandardError")
        
        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter_var)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filter_varm <- melt(fullmods_filter_var)
        #glimpse(fullmods_filter_var)
        #
        #names(fullmods_filter_varm) <- c("countrycode","variable","frequencies","econdata","climdata","value")
        #new_fft_var <- fullmods_filter_varm[order(fullmods_filter_varm[,2], fullmods_filter_varm[,3] ),] # Sort by column index [1] then [3]
        #estimates <- new_fft_var$value[which(new_fft_var$variable=="Estimate")]        
        #ses <- new_fft_var$value[which(new_fft_var$variable=="StandardError")]
        #new_fft_var <- new_fft_var[c(1:length(estimates)),-c(2,6)]    
        #new_fft_var$Estimate <- estimates
        #new_fft_var$Variance <- ses   
        #glimpse(new_fft_var) 
        #new_fft_var  <- new_fft_var[order(new_fft_var[,1]),] # Sort by column index [1] then [3]     
        #fullmods_filter_var <- new_fft_var    


        # dimnames(fullmods_filter_var)=list(countries,c("Estimate","Var"),ranges,dataset,datasetweather)
        # fullmods_filterm <- melt(fullmods_filter_var)
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter_v=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter_v$Estimate[which(is.infinite(fullmods_filter_v$Var))]=NA
        # names(fullmods_filter_v) <- c("countrycode","frequencies","econdata","climdata","Estimate","Variance")

        # dimnames(fullmods_filter_ar)=list(countries,c("dw_test","dw_pvalue"),ranges,dataset,datasetweather)
        # fullmods_filterm <- melt(fullmods_filter_ar)
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter_ar=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter_ar$dw_test[which(is.infinite(fullmods_filter_ar$Var))]=NA
        # names(fullmods_filter_ar) <- c("countrycode","frequencies","econdata","climdata","dw_test","dw_pvalue")

        
    ## 3.1. LAGS Country-level regressions LAGS (end)

            
    ## 3.3. LAGS categorizing - Plotting Figure 3   (start)

        # Categorizing statistically different estimates (start)
         
            fmod_fft <- new_fft_l
            #fmod_fft$Variance <- fullmods_filter_var$Variance
            
            #fmod_fft <- fullmods_filter
            fmod_fft <- fmod_fft[fmod_fft$econdata=="wb",]
            fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
            fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
            fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
            numcountries <- length(unique(fmod_fft$countrycode))
            #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
            fmod_fft$lowsignificant <- 0
            fmod_fft$unfilteredsignificant <- 0
            fmod_fft$signlofreq <- 0
            fmod_fft$signunfilt <- 0
                            
            uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
            `%notin%` <- Negate(`%in%`)
            fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
            
            numcountries <- length(unique(fmod_fft$countrycode))
            glimpse(fmod_fft)
            for (i in 1:numcountries){
                if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                    lastfreq <- 5
                    }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                        lastfreq <- 4
                    }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                        lastfreq <- 3
                    }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                        lastfreq <- 2
                    }else{next}
                if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                if(fmod_fft$Estimate[1+5*(i-1)] >0){fmod_fft$signunfilt[(1+5*(i-1)):(5+5*(i-1))]=1}
                m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                if (is.na(m)){next}
                if(m>0 ){
                    fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                    }
            }

                converging_stat <- 0
                converging <- 0
                not_converging <- 0
                Intensifying <- 0
                gray_area <- 0
                for (i in 1:numcountries){
                    if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                            
                            if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                            if(is.na(fmod_fft$StandardError[lastfreq+5*(i-1)])){next}
                    theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                    var <- fmod_fft$StandardError[lastfreq+5*(i-1)]^2 + fmod_fft$StandardError[1+5*(i-1)]^2
                    conf95 <-  (var^0.5)*1.65 #one-tail  95%
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                         
                        if(fmod_fft$sign[2+5*(i-1)]==0){gray_area <- c(gray_area,(1+5*(i-1)):(5+5*(i-1))) }else{
                                not_converging <- c(not_converging,(1+5*(i-1)):(5+5*(i-1)))
                        
                        if(theta+conf95 < 0 ){
                            Intensifying <- c(Intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                        }}
                    }
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                            converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                            if(theta-conf95 >= 0 ){
                            converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                        }
                    }
                    
                    
                        
                        if((lastfreq+5*(i-1)) %notin% uncertain){
                        fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                        if((1+5*(i-1)) %notin% uncertain){
                            fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                    }
                    
                    
                filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                
                fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                fmod_fft$category <- "Undefined"
                fmod_fft$category[converging] <- "converging"
                fmod_fft$category[not_converging] <- "Not Converging"
                fmod_fft$category[converging_stat] <- "converging_stat"
                fmod_fft$category[Intensifying] <- "*Intensifying"
                fmod_fft$category[gray_area] <- "gray_area"

                table(fmod_fft$category)/5
                
                fmod_fft$Category <- "Undefined"
                fmod_fft$Category[converging] <- "converging"
                fmod_fft$Category[not_converging] <- "Not Converging"
                table(fmod_fft$Category)/5
                
                fmod_fft$CoeffDifferents <- 0
                fmod_fft$CoeffDifferents[Intensifying] <- 1
                fmod_fft$CoeffDifferents[converging_stat] <- 1
                table(fmod_fft$CoeffDifferents)/5
                
                table(fmod_fft$lowsignificant)/5
                table(fmod_fft$unfilteredsignificant)/5

                
                table(fmod_fft$category[which(fmod_fft$unfilteredsignificant==1)])/5
                table(fmod_fft$category[which(fmod_fft$lowsignificant==1)])/5
                table(fmod_fft$category[which(fmod_fft$lowsignificant==0)])/5

                fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                library('rnaturalearth')
                world <- ne_countries(scale = "medium", returnclass = "sf")

                
                
                fmod_fft$iso_a3 <- fmod_fft$countrycode
                fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                

                #factor(fmod_fft$category)
                fmod_fft$category <- factor(fmod_fft$category, levels = c("*Intensifying", "Not Converging", "converging","converging_stat","gray_area"))
                df1 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                pal1 <- sequential_hcl(7, palette = "Oranges")
                pal2 <- sequential_hcl(7, palette = "Teal")
                df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                df1$lowsignificant <- factor(df1$lowsignificant)
                    
                hist_bf <- ggplot() + theme_bw() + 
                    geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                    #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                    stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                    name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                    scale_x_discrete(name="",
                        labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                        expression(paste("|",theta[f],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                df2 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                hist_bu <- ggplot() + theme_bw()+
                    geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                    #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                    scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                    name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                    scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                 
                fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(abs(fmod_fft$Estimate)>quantile(abs(fmod_fft$Estimate), na.rm=TRUE,0.99))],]
                #fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
               
               plot_fg <- ggplot(data=fmod_fft,aes(x=Lags,y=Estimate*100, group = countrycode,color=factor(category)))+
                    geom_line()+
                    scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],"dimgray",pal2[1])) +
                    theme_bw() + xlab("Number of Lags")+
                    geom_hline(yintercept=0,lty=2)+
                    #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),legend.position="none") +
                    #panel.background = element_blank()) +
                    ggtitle("")
                    plot_fg      

                table(fmod_fft$sign)/5
                table(fmod_fft$sign)/5
                table(fmod_fft$category[fmod_fft$sign==1])/5
                table(fmod_fft$category[fmod_fft$sign==0])/5
                table(fmod_fft$unfilteredsignificant[fmod_fft$sign==0])/5
                
                plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"),name="",labels=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified"))
                leg_cat <- get_legend(plot_categoriescolors)                    
                
                fig_3 <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                fig_3
                e_wb <- fmod_fft
                #ggsave("Figures/Fig3_wb_lags.png",dpi=600)
                
                # map_categories <- ggplot(data=fmod_fft_map) +
                #             geom_sf(data = fmod_fft_map_na,aes(fill=NA))+
                #             theme_minimal()+
                #             geom_sf(aes(fill = category))+
                #             scale_fill_manual(name = "",
                #                             #labels=c("Converging","Not Converging","*Intensifying"),
                #                             values=c(pal2[3],pal1[3],pal1[1]))+
                #             theme(legend.position="none")+
                #             ggtitle("Behavior of coefficients")
                #     ggarrange(map_categories,leg_cat,ncol=2,nrow=1,widths=c(5,1))
                    
                #     #ggsave("Fig3_map.png",dpi=600)

                # fmod_fft_map_low <- fmod_fft_map[which(fmod_fft_map$filters=="15 years"),]
                # missingcountries15 <- fmod_fft_map_low$countrycode[is.na(fmod_fft_map_low$Estimate)]
                # fmod_fft_map_10 <- fmod_fft_map[which(fmod_fft_map$filters=="10 years"),]
                
                # fmod_fft_map_low$Estimate[which(fmod_fft_map_low$countrycode %in% missingcountries15)] <- fmod_fft_map_10$Estimate[which(fmod_fft_map_10$countrycode %in% missingcountries15)]
                # ggplot(data = fmod_fft_map_low) +
                #                 geom_sf(data=fmod_fft_map_na,fill=NA)+theme_minimal()+
                #                 geom_sf(data=fmod_fft_map_low,aes(fill = Estimate*100))+
                #                 scale_fill_gradient2(
                #                     name = "Estimated impact \n (% per Degree)",
                #                     low = "red",
                #                     mid = "white",
                #                     high = "#00BFC4",
                #                     midpoint = 0,
                #                     space = "Lab",
                #                     na.value = "gray",
                #                     guide = "colourbar",
                #                     aesthetics = "fill")+
                #                     ggtitle("Growth effects")+
                #                 theme(legend.position="bottom")
                    #ggsave("Map_correcteddata.png",dpi=600)
                           
        # Categorizing statistically different estimates (end)
    ## 3.3. LAGS categorizing - Plotting Figure 3 (end) 
    
    ## 3.4. LAGS FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
        fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                origin = "iso3c",
                                destination = "continent")
        fmod_fft$invse <- 1/fmod_fft$StandardError
        fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
        f_unfiltpos <- fmod_fft[fmod_fft$signunfilt==1,]
        f_unfiltneg <- fmod_fft[fmod_fft$signunfilt==0,]
        
        #f_samesign <- f_samesign[f_samesign$signlofreq==1,]
        felm_1_l <- felm(Estimate ~ factor(Lags)|0|0|continent, data =f_unfiltpos,weights = f_unfiltpos$invse)
        felm_2_l <- felm(Estimate ~ factor(Lags)|0|0|continent, data =f_unfiltneg,weights = f_unfiltneg$invse)
        
        #stargazer(felm_2,felm_1, type="html", out="felm_1_2.html")
        stargazer(felm_2_l,felm_1_l,type="text")
        stargazer(felm_2_l,type="text")
        #write.csv(fmod_fft,'Growth_Estimates_Filters.csv')

        #felm_1_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft,weights = fmod_fft$invse)
        #felm_2_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft)
        #summary(felm_1_within)
        #stargazer(felm_2_within,felm_1_within,type="text")
        
        #stargazer(felm_2,felm_1,felm_2_within,felm_1_within, type="html", out="felm_1_2_within.html")
        
    ## 3.4. LAGS FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
    
    ## 3.1. Country-level regressions (start)
        dataset <- c("wb","barro","mad") 
        datasetweather <- c("UDel")
        periods <- c(0,3,5,10,15)
        countries <- unique(wb$countrycode)
        #periods <- c(0,10,20,25,27) Uncomment to get Supp Fig 2
        fullmods_filter=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_var=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_var_p=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_ar=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_p=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        panel_data <- data.frame(years = integer(),temp = double(), growth = double(), 
            preci = double(), countrycode = factor(), climdata = character(),
            econdata = character(), filter = character(), meant = double(), meanp=double(), growth_wdi=double(), tempsd=double(), numobs=double())
        for (mm in 1:length(datasetweather)){
            tempname <- paste(datasetweather[mm],"_pop_temp",sep="")
            preciname <- paste(datasetweather[mm],"_pop_preci",sep="")
            for (jj in (1:length(dataset))){
                DATA <- get(dataset[jj])
                countries=unique(factor(DATA$countrycode))
                #fullmods=array(dim=c(length(countries),2,length(periods)))
                for(k in 1:length(periods)){
                for(i in 1:length(countries)){
                    dat=DATA[which(DATA$countrycode==countries[i]),
                        which(colnames(DATA)%in%c("countrycode","year","growth",tempname,preciname))] 
                    dat <- dat[is.finite(dat$growth),]

                        gonext <- 0
                    for (nn in (1:(dim(dat)[2]))){
                        if(sum(is.na(dat[,nn]))==dim(dat)[1]){
                            #print(paste(countries[i],i,"missing complete column"))
                            gonext <- 1 }
                        
                    }
                    if(gonext==1){
                        gonext <- 0
                        next}
                    if(sum(complete.cases(dat))<2){
                            #print(paste(countries[i],i,"not complete cases"))
                            next
                    }
                    dat <- dat %>% rename("temp"=tempname,"preci"=preciname)
                    dat <- dat[complete.cases(dat),]
                    if(max(dat$year)-min(dat$year) != (dim(dat)[1]-1)){
                        #print(paste("Missing observation.",dataset[jj]))
                        next

                    }
                    #dat <- dat[-c(1),] #first observation is missing
                    if(dim(dat)[1]<5){next}
                    mt <- lm(temp~year+I(year^2), data = dat)
                    meanT <- mean(dat$temp, na.rm = TRUE)
                    t <- resid(mt)    
                    mp <- lm(preci~year+I(year^2), data = dat)
                    p <- resid(mp)
                    meanP <- mean(dat$preci, na.rm = TRUE)
                    mg <- lm(growth~year+I(year^2), data = dat)
                    g <- resid(mg)
                    growth_wdi <- dat$growth
                    if(length(t)<(2*(periods[k]))){
                            #print(paste(countries[i],i,"not enough data for this filter","length=",length(t),"periods=",periods[k]))
                            next}
                    if(k==1){
                        tempts <- t
                        precits <- p
                        gts <- g
                        ratio <- 1
                        ratio_p <- 1
                    } else{
                        tempts <- pass.filt(t, W=periods[k], type="low", method="Butterworth")
                        precits <- pass.filt(p, W=periods[k], type="low", method="Butterworth")
                        #gts <- pass.filt(g, W=periods[k], type="low", method="Butterworth")
                        gts <- g
                        x <- seq(1:length(tempts))
                        demeaned_t <- t
                        demeaned_tempts <- tempts
                        ratio <- abs(median(demeaned_t/demeaned_tempts))
                        ratio_p <- abs(median(p/precits))
                        if(median(demeaned_t/demeaned_tempts)<0){
                            #print(paste("ratio was negative for ", countries[i], "at period",toString(periods[k])))
                        }
                     
                    }
                     
                    temp <- data.frame(year = dat$year, temp= unclass(tempts))
                    preci <- data.frame(year = dat$year, preci= unclass(precits))
                    #growth <- data.frame(year = dat$year, growth = unclass(g))
                    growth <- data.frame(year = dat$year, growth = unclass(gts))
                    #growth_wdi <- data.frame(year = dat$year, growth = unclass(dat$growth))
                    filterdata <- merge(temp,growth, by = "year")
                    filterdata <- merge(filterdata,preci, by = "year")
                    #filterdata <- merge(filterdata,growth_wdi, by = "year")
                    names(filterdata) <- c("years","temp","growth","preci")

                    if(k==1){
                        filterdata2 <- filterdata
                        filterdata2$templag =c(NA,filterdata2$temp[1:(dim(filterdata2)[1]-1)])
                        mod_gfilterdata=lm(growth~temp+preci+templag,data=filterdata2)
                    } else{
                        mod_gfilterdata=lm(growth~temp+preci,data=filterdata)
                    }
                    dw <- dwtest(mod_gfilterdata) #test for autocorrelation
                    if (dw$p.value[1]<0.1){
                        vcov_ <- sandwich::vcovHAC(mod_gfilterdata)
                        hacse <- unname(sqrt(diag(vcov_)))[2] #s.e. corrected for autocorrelation
                        hacse_p<- unname(sqrt(diag(vcov_)))[3] #s.e. corrected for autocorrelation
                    } else {
                        hacse <- summary(mod_gfilterdata)$coefficients[2,2]
                        hacse_p <- summary(mod_gfilterdata)$coefficients[3,3]
                        vcov_ <-vcov(mod_gfilterdata)
                    } #We could try Newey West
                    
                    #filterdata$years <- dat$year
                    filterdata$countrycode <- rep(dat$countrycode[1],dim(filterdata)[1])
                    filterdata$climdata <- rep(datasetweather[mm],dim(filterdata)[1])
                    filterdata$econdata <- rep(dataset[jj],dim(filterdata)[1])
                    filterdata$filter <- rep(paste(periods[k],sep="-"),dim(filterdata)[1])
                    filterdata$meant <- rep(meanT,dim(filterdata)[1])
                    filterdata$meanp <- rep(meanP,dim(filterdata)[1])
                    filterdata$growth_wdi <- dat$growth 
                    filterdata$tempsd <- sd(t)
                    filterdata$numobs <- length(t)
                    panel_data <- rbind(panel_data,filterdata)
                    fullmods_filter[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[2,1]/ratio,hacse/ratio)
                    fullmods_filter_var[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[2,1]/ratio,vcov_[2,2]/(ratio^2))
                    fullmods_filter_ar[i,,k,jj,mm]=c(dw$statistic[1],dw$p.value[1])
                    fullmods_filter_p[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[3,1]/ratio_p,hacse_p/ratio_p) #Not HAC s.e. If usd should be divided by the ratio of precipietation timeseries
                    fullmods_filter_var_p[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[3,1]/ratio_p,vcov_[3,3]/(ratio_p^2))
                    }
                    
                
                }
            }
        }
        original_fullmods_filter <- fullmods_filter
        original_fullmods_filter_var <- fullmods_filter_var
        glimpse(fullmods_filter)
        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filterm <- melt(fullmods_filter)
        glimpse(fullmods_filterm)
        #
        names(fullmods_filterm) <- c("countrycode","variable","frequencies","econdata","climdata","value")
        new_fft <- fullmods_filterm[order(fullmods_filterm[,2], fullmods_filterm[,3] ),] # Sort by column index [1] then [3]
        glimpse(new_fft)
        estimates <- new_fft$value[which(new_fft$variable=="Estimate")]        
        ses <- new_fft$value[which(new_fft$variable=="StandardError")]
        new_fft <- new_fft[c(1:length(estimates)),-c(2,6)]    
        new_fft$Estimate <- estimates
        new_fft$StandardError <- ses   
        glimpse(new_fft) 
        new_fft <- new_fft[order(new_fft[,1]),] # Sort by column index [1] then [3]     
        fullmods_filter <- new_fft
        #
        # countriesbarro <- unique(factor(barro$countrycode))
        # countries_barro <- unlist(lapply(countriesbarro , as.character))
        # countries_barro_extended <- c(countries_barro,countries_wb[(length(countries_barro)+1):length(countries_wb)*NA])
        # countriesmad <- unique(factor(mad$countrycode))
        # countries_mad <- unlist(lapply(countriesmad , as.character))
        # countries_mad_extended <- c(countries_mad,countries_wb[(length(countries_mad)+1):length(countries_wb)*NA])
        # countries_wbrep <- rep(countries_wb,2*5)
        # countries_barrorep <- rep(countries_barro_extended,2*5)
        # countries_madrep <- rep(countries_mad_extended,2*5) #2*5
        # names1 <- c(countries_wbrep,countries_barrorep,countries_madrep)
        # #names2 <- c(names1,names1)
        # names2 <- names1 #only usiong UDel and not LMR
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter$Estimate[which(is.infinite(fullmods_filter$StandardError))]=NA
        # names(fullmods_filter) <- c("countrycode","frequencies","econdata","climdata","Estimate","StandardError")
        
        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter_var)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filter_varm <- melt(fullmods_filter_var)
        glimpse(fullmods_filter_var)
        #
        names(fullmods_filter_varm) <- c("countrycode","variable","frequencies","econdata","climdata","value")
        new_fft_var <- fullmods_filter_varm[order(fullmods_filter_varm[,2], fullmods_filter_varm[,3] ),] # Sort by column index [1] then [3]
        estimates <- new_fft_var$value[which(new_fft_var$variable=="Estimate")]        
        ses <- new_fft_var$value[which(new_fft_var$variable=="StandardError")]
        new_fft_var <- new_fft_var[c(1:length(estimates)),-c(2,6)]    
        new_fft_var$Estimate <- estimates
        new_fft_var$Variance <- ses   
        glimpse(new_fft_var) 
        new_fft_var  <- new_fft_var[order(new_fft_var[,1]),] # Sort by column index [1] then [3]     
        fullmods_filter_var <- new_fft_var    


        # dimnames(fullmods_filter_var)=list(countries,c("Estimate","Var"),ranges,dataset,datasetweather)
        # fullmods_filterm <- melt(fullmods_filter_var)
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter_v=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter_v$Estimate[which(is.infinite(fullmods_filter_v$Var))]=NA
        # names(fullmods_filter_v) <- c("countrycode","frequencies","econdata","climdata","Estimate","Variance")

        # dimnames(fullmods_filter_ar)=list(countries,c("dw_test","dw_pvalue"),ranges,dataset,datasetweather)
        # fullmods_filterm <- melt(fullmods_filter_ar)
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter_ar=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter_ar$dw_test[which(is.infinite(fullmods_filter_ar$Var))]=NA
        # names(fullmods_filter_ar) <- c("countrycode","frequencies","econdata","climdata","dw_test","dw_pvalue")

        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter_p)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filterm <- melt(fullmods_filter_p)
        glimpse(fullmods_filterm)
        #
        names(fullmods_filterm) <- c("countrycode","variable","frequencies","econdata","climdata","value")
        new_fft_p <- fullmods_filterm[order(fullmods_filterm[,2], fullmods_filterm[,3] ),] # Sort by column index [1] then [3]
        estimates <- new_fft_p$value[which(new_fft_p$variable=="Estimate")]        
        ses <- new_fft_p$value[which(new_fft_p$variable=="StandardError")]
        new_fft_p <- new_fft_p[c(1:length(estimates)),-c(2,6)]    
        new_fft_p$Estimate <- estimates
        new_fft_p$StandardError <- ses   
        glimpse(new_fft_p) 
        new_fft_p <- new_fft_p[order(new_fft_p[,1]),] # Sort by column index [1] then [3]     
        fullmods_filter_p <- new_fft_p

        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter_var_p)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filter_varm_p <- melt(fullmods_filter_var_p)
        glimpse(fullmods_filter_var_p)
        #
        names(fullmods_filter_varm_p) <- c("countrycode","variable","frequencies","econdata","climdata","value")
        new_fft_var_p <- fullmods_filter_varm_p[order(fullmods_filter_varm_p[,2], fullmods_filter_varm_p[,3] ),] # Sort by column index [1] then [3]
        estimates <- new_fft_var_p$value[which(new_fft_var_p$variable=="Estimate")]        
        ses <- new_fft_var_p$value[which(new_fft_var_p$variable=="StandardError")]
        new_fft_var_p <- new_fft_var_p[c(1:length(estimates)),-c(2,6)]    
        new_fft_var_p$Estimate <- estimates
        new_fft_var_p$Variance <- ses   
        glimpse(new_fft_var_p) 
        new_fft_var_p  <- new_fft_var_p[order(new_fft_var_p[,1]),] # Sort by column index [1] then [3]     
        fullmods_filter_var_p <- new_fft_var_p    
        
    ## 3.1. Country-level regressions (end)

    ## 3.3. categorizing - Plotting Figure 3  (start)

        # Categorizing statistically different estimates (start)
         
            fmod_fft <- new_fft
            fmod_fft$Variance <- fullmods_filter_var$Variance
            
            #fmod_fft <- fullmods_filter
            fmod_fft <- fmod_fft[fmod_fft$econdata=="wb",]
            fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
            fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
            fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
            numcountries <- length(unique(fmod_fft$countrycode))
            #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
            fmod_fft$lowsignificant <- 0
            fmod_fft$unfilteredsignificant <- 0
            fmod_fft$signlofreq <- 0
            fmod_fft$signunfilt <- 0
                            
            uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
            `%notin%` <- Negate(`%in%`)
            fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
            
            numcountries <- length(unique(fmod_fft$countrycode))
            glimpse(fmod_fft)
            for (i in 1:numcountries){
                if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                    lastfreq <- 5
                    }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                        lastfreq <- 4
                    }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                        lastfreq <- 3
                    }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                        lastfreq <- 2
                    }else{next}
                if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                if(fmod_fft$Estimate[1+5*(i-1)] >0){fmod_fft$signunfilt[(1+5*(i-1)):(5+5*(i-1))]=1}
                m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                if (is.na(m)){next}
                if(m>0 ){
                    fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                    }
            }

                converging_stat <- 0
                converging <- 0
                not_converging <- 0
                Intensifying <- 0
                gray_area <- 0
                for (i in 1:numcountries){
                    if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                            
                            if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                    theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                    var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                    conf95 <-  (var^0.5)*1.65 #one-tail  95%
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                         
                        if(fmod_fft$sign[2+5*(i-1)]==0){gray_area <- c(gray_area,(1+5*(i-1)):(5+5*(i-1))) }else{
                                not_converging <- c(not_converging,(1+5*(i-1)):(5+5*(i-1)))
                        
                        if(theta+conf95 < 0 ){
                            Intensifying <- c(Intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                        }}
                    }
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                            converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                            if(theta-conf95 >= 0 ){
                            converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                        }
                    }
                    
                    
                        
                        if((lastfreq+5*(i-1)) %notin% uncertain){
                        fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                        if((1+5*(i-1)) %notin% uncertain){
                            fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                    }
                    
                    
                filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                
                fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                fmod_fft$category <- "Undefined"
                fmod_fft$category[converging] <- "converging"
                fmod_fft$category[not_converging] <- "Not Converging"
                fmod_fft$category[converging_stat] <- "converging_stat"
                fmod_fft$category[Intensifying] <- "*Intensifying"
                fmod_fft$category[gray_area] <- "gray_area"

                table(fmod_fft$category)/5
                
                fmod_fft$Category <- "Undefined"
                fmod_fft$Category[converging] <- "converging"
                fmod_fft$Category[not_converging] <- "Not Converging"
                table(fmod_fft$Category)/5
                
                fmod_fft$CoeffDifferents <- 0
                fmod_fft$CoeffDifferents[Intensifying] <- 1
                fmod_fft$CoeffDifferents[converging_stat] <- 1
                table(fmod_fft$CoeffDifferents)/5
                
                table(fmod_fft$lowsignificant)/5
                table(fmod_fft$unfilteredsignificant)/5

                
                table(fmod_fft$category[which(fmod_fft$unfilteredsignificant==1)])/5
                table(fmod_fft$category[which(fmod_fft$lowsignificant==1)])/5
                table(fmod_fft$category[which(fmod_fft$lowsignificant==0)])/5
                fmod_fft_withoutliers <- 

                fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                library('rnaturalearth')
                world <- ne_countries(scale = "medium", returnclass = "sf")

                
                
                fmod_fft$iso_a3 <- fmod_fft$countrycode
                fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                

                #factor(fmod_fft$category)
                fmod_fft$category <- factor(fmod_fft$category, levels = c("*Intensifying", "Not Converging", "converging","converging_stat","gray_area"))
                df1 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                pal1 <- sequential_hcl(7, palette = "Oranges")
                pal2 <- sequential_hcl(7, palette = "Teal")
                df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                df1$lowsignificant <- factor(df1$lowsignificant)
                    
                hist_bf <- ggplot() + theme_bw() + 
                    geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                    #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                    stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                    #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"), 
                    name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                    scale_x_discrete(name="",
                        labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                        expression(paste("|",theta[f],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                df2 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                hist_bu <- ggplot() + theme_bw()+
                    geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                    #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                    scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                    #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"), 
                    name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                    scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                fmod_fft_withoutliers <- fmod_fft
                fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(abs(fmod_fft$Estimate)>quantile(abs(fmod_fft$Estimate), na.rm=TRUE,0.99))],]
                #fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
               
               plot_fg <- ggplot(data=fmod_fft,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(category)))+
                    geom_line()+
                    scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],"dimgray",pal2[1])) +
                    theme_bw() + xlab("Minimum Periodicity after Filtering")+
                    geom_hline(yintercept=0,lty=2)+
                    #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),legend.position="none") +
                    #panel.background = element_blank()) +
                    ggtitle("")
                    plot_fg      

                table(fmod_fft$sign)/5
                table(fmod_fft$sign)/5
                table(fmod_fft$category[fmod_fft$sign==1])/5
                table(fmod_fft$category[fmod_fft$sign==0])/5
                table(fmod_fft$unfilteredsignificant[fmod_fft$sign==0])/5
                
                plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"),name="",labels=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified"))
                leg_cat <- get_legend(plot_categoriescolors)                    
                
                fig_3 <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                #fig_3
                e_wb <- fmod_fft
                #ggsave("Figures/Fig3.png",dpi=600)
                
                map_categories <- ggplot(data=fmod_fft_map) +
                            geom_sf(data = fmod_fft_map_na,aes(fill=NA))+
                            theme_void()+
                            geom_sf(aes(fill = category))+
                            scale_fill_manual(name = "",
                                            #labels=c("Converging","Not Converging","*Intensifying"),
                                            values=c(pal1[1],pal2[3],"gray",pal1[3]))+
                            theme(legend.position="none")
                            #ggtitle("Behavior of coefficients")
                            #map_categories
                    
                    
                    #ggarrange(map_categories,leg_cat,ncol=2,nrow=1,widths=c(5,1))
                    
                #     #ggsave("Fig3_map.png",dpi=600)

                 fmod_fft_map_low <- fmod_fft_map[which(fmod_fft_map$filters=="15 years"),]
                 missingcountries15 <- fmod_fft_map_low$countrycode[is.na(fmod_fft_map_low$Estimate)]
                 fmod_fft_map_10 <- fmod_fft_map[which(fmod_fft_map$filters=="10 years"),]
                
                 fmod_fft_map_low$Estimate[which(fmod_fft_map_low$countrycode %in% missingcountries15)] <- fmod_fft_map_10$Estimate[which(fmod_fft_map_10$countrycode %in% missingcountries15)]
                    library("ggpattern")

                    fmod_fft_map_low$signlofreq[fmod_fft_map_low$signlofreq==0] <- "Negative"
                    fmod_fft_map_low$signlofreq[fmod_fft_map_low$signlofreq==1] <- "Positive"
                   map_cat_sign <-  ggplot(fmod_fft_map_low) +
                            geom_sf(data = fmod_fft_map_na,aes(fill=NA))+ 
                        geom_sf_pattern(data = fmod_fft_map_low, 
                        aes(fill = category, pattern = factor(signlofreq))) + # adding a mapping for the pattern
                        scale_fill_manual(name = "",
                                            #labels=c("Converging","Not Converging","*Intensifying"),
                                            values=c(pal1[1],pal2[3],"gray",pal1[3]),
                                            #values=c(pal1[1],pal2[3],pal2[1],"gray",pal1[3]),#,
                                            guide="none")+
                                            theme_void()+
                                            guides(pattern=guide_legend("Sign of Effect"))+
                                            theme(legend.position="bottom")
                                            #map_cat_sign
                        
                         #ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),map_cat_sign,ncol=1,nrow=3) 
                         #ggsave("Figures/Fig3_withMap_bothFiltered.png",dpi=400)

                         library("ggalluvial")

                        
                        fmod_fft2 <- fmod_fft_withoutliers
                        names(fmod_fft2)
                        fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                        fmod_fft2$signsign_filt <- paste(fmod_fft2$signlofreq,fmod_fft2$lowsignificant)
                        fmod_fft2$gl <- "Unclassified"
                        fmod_fft2$gl[fmod_fft2$lowsignificant ==1 ] <- "Presence of growth effects"         
                        fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                        
                        
                        alluvial_1 <- data.frame(fmod_fft2[,c(20,21,22)])
                        es_alluvial <- plyr::count(alluvial_1)
                        #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                        
                        glimpse(es_alluvial)
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 0"] <- "Negative not significant"
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 1"] <- "Negative significant"
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 0"] <- "Positive not significant"
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 1"] <- "Positive significant"

                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                        


                        glimpse(es_alluvial)
                        es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial$signsign_filt <- factor(es_alluvial$signsign_filt, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial$gl <- factor(es_alluvial$gl, levels=c("Presence of growth effects","Only level effects detected","Unclassified"))
                        

                        pal1 <- sequential_hcl(7, palette = "Purples")
                            pal2 <- sequential_hcl(7, palette = "Greens")

                        library('ggalluvial')
                        sankeyplot_pval <- ggplot(as.data.frame(es_alluvial),
                        aes(y = freq/5, axis1 = factor(signsign_unf), axis2 = factor(signsign_filt))) +
                        geom_alluvium(aes(fill = factor(gl)), width = 1/12) +
                        geom_stratum(width = 1/12, aes(fill =  gl), color = "grey")+
                        #geom_stratum() +
                        geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                        scale_x_discrete(limits = c("Unfiltered", "15-year filtered estimate"), expand = c(0.1, .1))  +
                        ylab("Number of countries")+
                        ggtitle("Evidence from point estimates")+
                        scale_fill_manual(values=c("#d3818c","#7375a4","gray"),name="Type of Impact",breaks = levels(es_alluvial$gl))+
                        theme(plot.title = element_text(hjust = 0.5)) +
                        theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
                        guides(fill = guide_legend(title.position="top", title.hjust = 0.5))
                        #ggsave("TempEffect_GvL.png",dpi=600)

                        fmod_fft3 <- fmod_fft_withoutliers
                        names(fmod_fft3)
                        #fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                        fmod_fft3$signsign_filt <- paste(fmod_fft3$signlofreq,fmod_fft3$lowsignificant)
                        fmod_fft3$gl <- "Unclassified"
                        fmod_fft3$gl[fmod_fft3$category =="*Intensifying" ] <- "Presence of growth effects detected"  
                        fmod_fft3$gl[fmod_fft3$category =="Not Converging" ] <- "Presence of growth effects detected"   
                        fmod_fft3$gl[fmod_fft3$category =="converging" ] <- "Only level effects detected"       
                        #fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                        
                        fmod_fft3$cat <- "Unclassified"
                        fmod_fft3$cat[fmod_fft3$category =="*Intensifying" ] <- "*Intensifying"  
                        fmod_fft3$cat[fmod_fft3$category =="Not Converging" ] <- "Not converging to zero"   
                        fmod_fft3$cat[fmod_fft3$category =="converging" ] <- "Converging to zero"       

                        alluvial_2 <- data.frame(fmod_fft3[,c(20,21,22)])
                        es_alluvial2 <- plyr::count(alluvial_2)
                        #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                        
                        glimpse(es_alluvial2)
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 0"] <- "Negative not significant"
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 1"] <- "Negative significant"
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 0"] <- "Positive not significant"
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 1"] <- "Positive significant"


                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                        


                        glimpse(es_alluvial2)
                        #es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial2$signsign_filt <- factor(es_alluvial2$signsign_filt, 
                            levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial2$gl <- factor(es_alluvial2$gl, 
                            levels=c("Presence of growth effects detected","Only level effects detected","Unclassified"))
                        
                        es_alluvial2$cat <- factor(es_alluvial2$cat, 
                            levels=c("*Intensifying","Not converging to zero","Converging to zero","Unclassified"))
                        factor(es_alluvial2$cat)

                        pal1 <- sequential_hcl(7, palette = "Oranges")
                        pal2 <- sequential_hcl(7, palette = "Teal")

                        library('ggalluvial')
                        glimpse(es_alluvial2)

                        sankeyplot_cat <- ggplot(as.data.frame(es_alluvial2),
                        aes(y = freq/5, axis1 = factor(signsign_filt), axis2 = factor(cat))) +
                        geom_alluvium(aes(fill = factor(cat)), width = 1/12) +
                        geom_stratum(width = 1/5, aes(fill = cat), color = "grey")+
                        #geom_stratum() +
                        geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                        scale_x_discrete(limits = c("15-year filtered estimate", "Convergence to zero category"), expand = c(0.1, .1))  +
                        ylab("")+
                        scale_fill_manual(values=c(pal1[1],pal2[3],pal1[3],"dimgray"),name="Type of Impact") +
                        ggtitle("Evidence from the magnitude change of \u03b8")+
                        theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")
                        #guides(fill = guide_legend(title.position="top", title.hjust = 0.5))

                        ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom", widths = c(2, 1) )

                        #ggsave("sankey_diag_3.pdf",device="pdf",dpi=300)

                        ggarrange(plot_fg,ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom", widths = c(2, 1) ),nrow=2,ncol=1)
                        ggsave("Figure3_estimates_sankey.png",dpi=300)
                           
        # Categorizing statistically different estimates (end)
    ## 3.3. categorizing - Plotting Figure 3 (end) 

    ## 3.3. categorizing - PRECIP  (start)

        # Categorizing statistically different estimates (start)
         
            fmod_fft <- new_fft_p
            fmod_fft$Variance <- fullmods_filter_var$Variance
            
            #fmod_fft <- fullmods_filter
            fmod_fft <- fmod_fft[fmod_fft$econdata=="wb",]
            fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
            fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
            fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
            numcountries <- length(unique(fmod_fft$countrycode))
            #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
            fmod_fft$lowsignificant <- 0
            fmod_fft$unfilteredsignificant <- 0
            fmod_fft$signlofreq <- 0
            fmod_fft$signunfilt <- 0
                            
            uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
            `%notin%` <- Negate(`%in%`)
            fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
            
            numcountries <- length(unique(fmod_fft$countrycode))
            glimpse(fmod_fft)
            for (i in 1:numcountries){
                if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                    lastfreq <- 5
                    }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                        lastfreq <- 4
                    }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                        lastfreq <- 3
                    }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                        lastfreq <- 2
                    }else{next}
                if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                if(fmod_fft$Estimate[1+5*(i-1)] >0){fmod_fft$signunfilt[(1+5*(i-1)):(5+5*(i-1))]=1}
                m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                if (is.na(m)){next}
                if(m>0 ){
                    fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                    }
            }

                converging_stat <- 0
                converging <- 0
                not_converging <- 0
                Intensifying <- 0
                gray_area <- 0
                for (i in 1:numcountries){
                    if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                            
                            if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                    theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                    var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                    conf95 <-  (var^0.5)*1.65 #one-tail  95%
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                         
                        if(fmod_fft$sign[2+5*(i-1)]==0){gray_area <- c(gray_area,(1+5*(i-1)):(5+5*(i-1))) }else{
                                not_converging <- c(not_converging,(1+5*(i-1)):(5+5*(i-1)))
                        
                        if(theta+conf95 < 0 ){
                            Intensifying <- c(Intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                        }}
                    }
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                            converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                            if(theta-conf95 >= 0 ){
                            converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                        }
                    }
                    
                    
                        
                        if((lastfreq+5*(i-1)) %notin% uncertain){
                        fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                        if((1+5*(i-1)) %notin% uncertain){
                            fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                    }
                    
                    
                filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                
                fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                fmod_fft$category <- "Undefined"
                fmod_fft$category[converging] <- "converging"
                fmod_fft$category[not_converging] <- "Not Converging"
                fmod_fft$category[converging_stat] <- "converging_stat"
                fmod_fft$category[Intensifying] <- "*Intensifying"
                fmod_fft$category[gray_area] <- "gray_area"

                table(fmod_fft$category)/5
                
                fmod_fft$Category <- "Undefined"
                fmod_fft$Category[converging] <- "converging"
                fmod_fft$Category[not_converging] <- "Not Converging"
                table(fmod_fft$Category)/5
                
                fmod_fft$CoeffDifferents <- 0
                fmod_fft$CoeffDifferents[Intensifying] <- 1
                fmod_fft$CoeffDifferents[converging_stat] <- 1
                table(fmod_fft$CoeffDifferents)/5
                
                table(fmod_fft$lowsignificant)/5
                table(fmod_fft$unfilteredsignificant)/5
                table(fmod_fft$category[which(fmod_fft$lowsignificant==1)])/5

                
                table(fmod_fft$category[which(fmod_fft$unfilteredsignificant==1)])/5
                table(fmod_fft$category[which(fmod_fft$lowsignificant==1)])/5
                table(fmod_fft$sign[which(fmod_fft$lowsignificant==1)])/5

                fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                library('rnaturalearth')
                world <- ne_countries(scale = "medium", returnclass = "sf")

                
                
                fmod_fft$iso_a3 <- fmod_fft$countrycode
                fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                

                #factor(fmod_fft$category)
                fmod_fft$category <- factor(fmod_fft$category, levels = c("*Intensifying", "Not Converging", "converging","converging_stat","gray_area"))
                df1 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                pal1 <- sequential_hcl(7, palette = "Oranges")
                pal2 <- sequential_hcl(7, palette = "Teal")
                df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                df1$lowsignificant <- factor(df1$lowsignificant)
                    
                hist_bf <- ggplot() + theme_bw() + 
                    geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                    #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                    stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                    #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"), 
                    name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                    scale_x_discrete(name="",
                        labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                        expression(paste("|",theta[f],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                df2 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                hist_bu <- ggplot() + theme_bw()+
                    geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                    #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                    scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                    #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"), 
                    name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                    scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                 
                fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(abs(fmod_fft$Estimate)>quantile(abs(fmod_fft$Estimate), na.rm=TRUE,0.99))],]
                #fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
               
               plot_fg <- ggplot(data=fmod_fft,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(category)))+
                    geom_line()+
                    scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],"dimgray",pal2[1])) +
                    theme_bw() + xlab("Minimum Periodicity after Filtering")+
                    geom_hline(yintercept=0,lty=2)+
                    #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    ylab("Estimated Effect of \n 1 additional milimeter of rain on Growth (pp)") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),legend.position="none") +
                    #panel.background = element_blank()) +
                    ggtitle("Precipitation")+
                        theme(plot.title = element_text(hjust = 0.5))
                    plot_fg      

                table(fmod_fft$sign)/5
                table(fmod_fft$sign)/5
                table(fmod_fft$category[fmod_fft$sign==1])/5
                table(fmod_fft$category[fmod_fft$sign==0])/5
                table(fmod_fft$unfilteredsignificant[fmod_fft$sign==0])/5
                
                plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"),name="",labels=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified"))
                leg_cat <- get_legend(plot_categoriescolors)                    
                
                fig_3 <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                fig_3
                e_wb <- fmod_fft
                #ggsave("Figures/Fig3.png",dpi=600)
                
                # map_categories <- ggplot(data=fmod_fft_map) +
                #             geom_sf(data = fmod_fft_map_na,aes(fill=NA))+
                #             theme_void()+
                #             geom_sf(aes(fill = category))+
                #             scale_fill_manual(name = "",
                #                             #labels=c("Converging","Not Converging","*Intensifying"),
                #                             values=c(pal1[1],pal2[3],"gray",pal1[3]))+
                #             theme(legend.position="none")
                #             #ggtitle("Behavior of coefficients")
                #             map_categories
                    
                    
                    #ggarrange(map_categories,leg_cat,ncol=2,nrow=1,widths=c(5,1))
                    
                #     #ggsave("Fig3_map.png",dpi=600)

                 fmod_fft_map_low <- fmod_fft_map[which(fmod_fft_map$filters=="15 years"),]
                 missingcountries15 <- fmod_fft_map_low$countrycode[is.na(fmod_fft_map_low$Estimate)]
                 fmod_fft_map_10 <- fmod_fft_map[which(fmod_fft_map$filters=="10 years"),]
                
                 fmod_fft_map_low$Estimate[which(fmod_fft_map_low$countrycode %in% missingcountries15)] <- fmod_fft_map_10$Estimate[which(fmod_fft_map_10$countrycode %in% missingcountries15)]
                    library("ggpattern")

                    fmod_fft_map_low$signlofreq[fmod_fft_map_low$signlofreq==0] <- "Negative"
                    fmod_fft_map_low$signlofreq[fmod_fft_map_low$signlofreq==1] <- "Positive"
                   map_cat_sign <-  ggplot(fmod_fft_map_low) +
                            geom_sf(data = fmod_fft_map_na,aes(fill=NA))+ 
                        geom_sf_pattern(data = fmod_fft_map_low, 
                        aes(fill = category, pattern = factor(signlofreq))) + # adding a mapping for the pattern
                        scale_fill_manual(name = "",
                                            #labels=c("Converging","Not Converging","*Intensifying"),
                                            values=c(pal1[1],pal2[3],"gray",pal1[3]),
                                            #values=c(pal1[1],pal2[3],pal2[1],"gray",pal1[3]),#,
                                            guide="none")+
                                            theme_void()+
                                            guides(pattern=guide_legend("Sign of Effect"))+
                                            theme(legend.position="bottom")
                                            #map_cat_sign
                        
                         #ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),map_cat_sign,ncol=1,nrow=3) 
                         #ggsave("Figures/Fig3_PRECIP_withMap_bothFiltered.png",dpi=400)

                         library("ggalluvial")


                         
                        fmod_fft2 <- fmod_fft
                        names(fmod_fft2)
                        fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                        fmod_fft2$signsign_filt <- paste(fmod_fft2$signlofreq,fmod_fft2$lowsignificant)
                        fmod_fft2$gl <- "Unclassified"
                        fmod_fft2$gl[fmod_fft2$lowsignificant ==1 ] <- "Presence of growth effects"         
                        fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                        
                        
                        alluvial_1 <- data.frame(fmod_fft2[,c(20,21,22)])
                        es_alluvial <- plyr::count(alluvial_1)
                        #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                        
                        glimpse(es_alluvial)
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 0"] <- "Negative not significant"
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 1"] <- "Negative significant"
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 0"] <- "Positive not significant"
                        es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 1"] <- "Positive significant"

                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                        es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                        


                        glimpse(es_alluvial)
                        es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial$signsign_filt <- factor(es_alluvial$signsign_filt, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial$gl <- factor(es_alluvial$gl, levels=c("Presence of growth effects","Only level effects detected","Unclassified"))
                        

                        pal1 <- sequential_hcl(7, palette = "Purples")
                            pal2 <- sequential_hcl(7, palette = "Greens")

                        library('ggalluvial')
                        sankeyplot_pval <- ggplot(as.data.frame(es_alluvial),
                        aes(y = freq/5, axis1 = factor(signsign_unf), axis2 = factor(signsign_filt))) +
                        geom_alluvium(aes(fill = factor(gl)), width = 1/12) +
                        geom_stratum(width = 1/12, aes(fill =  gl), color = "grey")+
                        #geom_stratum() +
                        geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                        scale_x_discrete(limits = c("Unfiltered", "15-year filtered estimate"), expand = c(0.1, .1))  +
                        ylab("Number of countries")+
                        ggtitle("Evidence from point estimates")+
                        scale_fill_manual(values=c("#d3818c","#7375a4","gray"),name="Type of Impact",breaks = levels(es_alluvial$gl))+
                        theme(plot.title = element_text(hjust = 0.5)) +
                        theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
                        guides(fill = guide_legend(title.position="top", title.hjust = 0.5))
                        ggsave("TempEffect_GvL.png",dpi=600)

                        fmod_fft3 <- fmod_fft
                        names(fmod_fft3)
                        #fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                        fmod_fft3$signsign_filt <- paste(fmod_fft3$signlofreq,fmod_fft3$lowsignificant)
                        fmod_fft3$gl <- "Unclassified"
                        fmod_fft3$gl[fmod_fft3$category =="*Intensifying" ] <- "Presence of growth effects detected"  
                        fmod_fft3$gl[fmod_fft3$category =="Not Converging" ] <- "Presence of growth effects detected"   
                        fmod_fft3$gl[fmod_fft3$category =="converging" ] <- "Only level effects detected"       
                        #fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                        
                        fmod_fft3$cat <- "Unclassified"
                        fmod_fft3$cat[fmod_fft3$category =="*Intensifying" ] <- "*Intensifying"  
                        fmod_fft3$cat[fmod_fft3$category =="Not Converging" ] <- "Not converging to zero"   
                        fmod_fft3$cat[fmod_fft3$category =="converging" ] <- "Converging to zero"       

                        alluvial_2 <- data.frame(fmod_fft3[,c(20,21,22)])
                        es_alluvial2 <- plyr::count(alluvial_2)
                        #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                        
                        glimpse(es_alluvial2)
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 0"] <- "Negative not significant"
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 1"] <- "Negative significant"
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 0"] <- "Positive not significant"
                        es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 1"] <- "Positive significant"


                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                        #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                        


                        glimpse(es_alluvial2)
                        #es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial2$signsign_filt <- factor(es_alluvial2$signsign_filt, 
                            levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                        es_alluvial2$gl <- factor(es_alluvial2$gl, 
                            levels=c("Presence of growth effects detected","Only level effects detected","Unclassified"))
                        
                        es_alluvial2$cat <- factor(es_alluvial2$cat, 
                            levels=c("*Intensifying","Not converging to zero","Converging to zero","Unclassified"))
                        factor(es_alluvial2$cat)

                        pal1 <- sequential_hcl(7, palette = "Oranges")
                        pal2 <- sequential_hcl(7, palette = "Teal")

                        library('ggalluvial')
                        glimpse(es_alluvial2)

                        sankeyplot_cat <- ggplot(as.data.frame(es_alluvial2),
                        aes(y = freq/5, axis1 = factor(signsign_filt), axis2 = factor(cat))) +
                        geom_alluvium(aes(fill = factor(cat)), width = 1/12) +
                        geom_stratum(width = 1/5, aes(fill = cat), color = "grey")+
                        #geom_stratum() +
                        geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                        scale_x_discrete(limits = c("15-year filtered estimate", "Convergence to zero category"), expand = c(0.1, .1))  +
                        ylab("")+
                        scale_fill_manual(values=c(pal2[3],pal1[3],"dimgray"),name="Type of Impact") +
                        ggtitle("Evidence from trends")+
                        theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")
                        #guides(fill = guide_legend(title.position="top", title.hjust = 0.5))

                        ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom", widths = c(2, 1) )

                        #ggsave("sankey_diag_3.pdf",device="pdf",dpi=300)

                        ggarrange(plot_fg,ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom", widths = c(1, 1) ),nrow=2,ncol=1)
                        ggsave("Figure3PRECIP_estimates_sankey_precipitation.png",dpi=300)

                           
        # Categorizing statistically different estimates (end)
    ## 3.3. categorizing - PRECIP (end) 

    ## 3.4. FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
        fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                origin = "iso3c",
                                destination = "continent")
        fmod_fft$invse <- 1/fmod_fft$StandardError
        fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
        f_unfiltpos <- fmod_fft[fmod_fft$signunfilt==1,]
        f_unfiltneg <- fmod_fft[fmod_fft$signunfilt==0,]
        
        #f_samesign <- f_samesign[f_samesign$signlofreq==1,]
        felm_1 <- felm(Estimate ~ factor(frequencies)|0|0|continent, data =f_unfiltpos,weights = f_unfiltpos$invse)
        felm_2 <- felm(Estimate ~ factor(frequencies)|0|0|continent, data =f_unfiltneg,weights = f_unfiltneg$invse)
        
        #stargazer(felm_2,felm_1, type="html", out="felm_1_2.html")
        stargazer(felm_2,felm_1,type="text")
        #write.csv(fmod_fft,'Growth_Estimates_Filters.csv')

        #felm_1_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft,weights = fmod_fft$invse)
        #felm_2_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft)
        #summary(felm_1_within)
        #stargazer(felm_2_within,felm_1_within,type="text")
        
        #stargazer(felm_2,felm_1,felm_2_within,felm_1_within, type="html", out="felm_1_2_within.html")
        
    ## 3.4. FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
    
    ## 3.5. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2 (start)
        # Adding socioeconomic variables (start)
            wb2 <- wb[!is.na(wb$gdppc),]
            glimpse(wb2)
            minyear <- aggregate(wb2$year, by = list(wb2$countrycode), FUN = min, na.rm = TRUE )
            names(minyear) <- c("countrycode","year")
            minyear$extra <- paste(minyear$countrycode,minyear$year,sep="")
            wb2$extra <- paste(wb2$countrycode,wb2$year,sep="")
            wb2 <- wb2[which(wb2$extra %in% minyear$extra),]
            wb2$loggdppc <- log(wb2$gdppc)
            wb2 <- wb2[,which(names(wb2) %in% c("countrycode","loggdppc"))]
            fmod_fft2 <- merge(fmod_fft, wb2,by = "countrycode")
            meangdppc <- aggregate(wb$gdppc, by = list(wb$countrycode), FUN = mean, na.rm = TRUE )
            names(meangdppc ) <-  c("countrycode", "meangdppc ")
            fmod_fft2 <- merge(fmod_fft2, meangdppc , by = "countrycode")
            pop_data <- wb_data("SP.POP.TOTL", start_date = 2019, end_date = 2019)
            names(pop_data)[2] <- "countrycode"
            fmod_fft2 <- merge(fmod_fft2, pop_data,by = "countrycode") 
            a <- wb_search("GDP.*PPP")
            a <- a[6,1]
            gdp_data <- wb_data(a, start_date = 2019, end_date = 2019)
            names(gdp_data)[2] <- "countrycode"
            fmod_fft2 <- merge(fmod_fft2, gdp_data,by = "countrycode")         
            a <- wb_search("GDP per capita")
            a <- a[10,1]
            gdp_data <- wb_data(a, start_date = 2019, end_date = 2019)
            names(gdp_data)[2] <- "countrycode"
            fmod_fft2 <- merge(fmod_fft2, gdp_data,by = "countrycode")
            meanT <- aggregate(wb$UDel_pop_temp, by = list(wb$countrycode), FUN = mean, na.rm = TRUE )
            names(meanT ) <-  c("countrycode", "meanT")
            fmod_fft2 <- merge(fmod_fft2, meanT,by = "countrycode")
            fmod_fft2 <- fmod_fft2[, !duplicated(colnames(fmod_fft))]
            #glimpse(fmod_fft2)
            #ggplot(data=fmod_fft2, aes(x=meanT,y= log(NY.GDP.PCAP.PP.KD), color=category))+geom_point()+theme_bw()
        # Adding socioeconomic variables (end)

        # weighted average (start)

            meanestimate <- data.frame(mean = rep(0,15),filter = rep(0,15), weight = rep(0,15), var = rep(0,15))
            frequency <- c("Unfiltered","3 years","5 years","10 years","15 years")
            fmod_fft15 <- fmod_fft2
            fmod_fft2[fmod_fft2$countrycode=="USA",]
            #fmod_fft15 <- fmod_fft15[!is.na(fmod_fft15$Estimate),]
            #fmod_fft15 <- fmod_fft15[!is.na(fmod_fft15$SP.POP.TOTL),]
            #fmod_fft15 <- fmod_fft15[!is.na(fmod_fft15$NY.GDP.MKTP.PP.KD),]
            c15 <- levels(factor(fmod_fft15$countrycode))
            fmod_fft15 <- fmod_fft2[(fmod_fft2$countrycode %in% c15),]
            weighted.var <- function(x, w, na.rm = FALSE) {
                if (na.rm) {
                    w <- w[i <- !is.na(x)]
                    x <- x[i]
                }
                sum.w <- sum(w)
                sum.w2 <- sum(w^2)
                mean.w <- sum(x * w) / sum(w)
                (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                na.rm)
                }
            for (i in 1:5){
                meanestimate$filter[(1 + ((i-1)*3)):(3 + ((i-1)*3))] <- frequency[i]
                meanestimate$mean[1 + ((i-1)*3)]=mean(fmod_fft15$Estimate[fmod_fft15$filters==frequency[i]],na.rm=TRUE)
                meanestimate$weight[1 + ((i-1)*3)]="Unweighted"
                meanestimate$var[1 + ((i-1)*3)] <- var(fmod_fft15$Estimate[fmod_fft15$filters==frequency[i]],na.rm=TRUE)


                x1 <- fmod_fft15[!is.na(fmod_fft15$SP.POP.TOTL),]
                x1 <- x1[!is.na(x1$Estimate),]
                x1 <- x1[x1$lowsignificant==1,]
                est <- x1$Estimate[x1$filters==frequency[i]]
                popw <- x1$SP.POP.TOTL[x1$filters==frequency[i]]
                gdpw <- x1$NY.GDP.MKTP.PP.KD[x1$filters==frequency[i]]
                countrypop <- x1$countrycode[x1$filters==frequency[i]]
                hist((est*popw)/sum(popw))
                datpop <- data.frame(estimate=est,population=popw,country=countrypop)
                ggplot(datpop,aes(x=estimate,y=population,label=country))+geom_label()
                meanestimate$mean[2 + ((i-1)*3)] <-sum((est*popw)/sum(popw),na.rm=TRUE)
                meanestimate$weight[2 + ((i-1)*3)] <- "Population in 2019"
                meanestimate$var[2 + ((i-1)*3)] <- weighted.var(x=est, w = popw,na.rm=TRUE)
                sum((est*(gdpw/popw))/sum(gdpw/popw,na.rm=TRUE),na.rm=TRUE)
                datpop <- data.frame(est=est,gdppopw=gdpw/popw,country=countrypop)
                #ggplot(datpop,aes(x=est,y=gdppopw,label=country))+geom_label()
                
                #x1 <- fmod_fft15[!is.na(fmod_fft15$Estimate) & !is.na(fmod_fft15$NY.GDP.PCAP.PP.KD),]
                x1 <- fmod_fft15[!is.na(fmod_fft15$Estimate) & !is.na(fmod_fft15$NY.GDP.MKTP.PP.KD),]
                #x1 <- x1[x1$lowsignificant==1,]
                gdpw <- x1$NY.GDP.MKTP.PP.KD[x1$filters==frequency[i]]
                #gdpw <- x1$NY.GDP.PCAP.PP.KD[x1$filters==frequency[i]]
                est <- x1$Estimate[x1$filters==frequency[i]]
                country <- x1$countrycode[x1$filters==frequency[i]]
                datpop <- data.frame(est=est,gdpw=gdpw,country=country)
                
                meanestimate$mean[3 + ((i-1)*3)] <- sum((est*gdpw)/sum(gdpw),na.rm=TRUE)
                meanestimate$weight[3 + ((i-1)*3)] <- "GDP, PPP constant 2017 USD"
                meanestimate$var[3 + ((i-1)*3)]  <- weighted.var(x=x1$Estimate[x1$filters==frequency[i]], w = x1$NY.GDP.MKTP.PP.KD[x1$filters==frequency[i]],na.rm=TRUE)

            }
            #See weighted means
            meanestimate
        # weighted average (end)

        # Socioeconomic characteristics across categories (start)
            
            # fmod_fft2$catsign <- paste(fmod_fft2$category,fmod_fft2$signlofreq,sep=".")
            # fmod_fft2 <- fmod_fft2[fmod_fft2$category!="other",]
            # fmod_fft2$logGDPpc <- log(fmod_fft2$NY.GDP.PCAP.PP.KD)
            # densityGDP6 <-  ggstatsplot::ggbetweenstats(
            #     data = fmod_fft2,
            #     x = catsign,
            #     y = logGDPpc
            #     )

            # densitymeanT6 <-  ggstatsplot::ggbetweenstats(
            # data = fmod_fft2,
            # x = catsign,
            # y = meanT
            # )

            # densityGDP3 <-  ggstatsplot::ggbetweenstats(
            # data = fmod_fft2,
            # x = category,
            # y = logGDPpc
            # )

            # densitymeanT3 <-  ggstatsplot::ggbetweenstats(
            # data = fmod_fft2,
            # x = category,
            # y = meanT
            # )

            # ggarrange(densityGDP3,densitymeanT3)
            # #ggsave("Category_byTandGDP.png",dpi=600)
        
            # cat <- data.frame(T = rep(0,6),minT=rep(0,6),maxT=rep(0,6),G = rep(0,6),minG=rep(0,6),maxG=rep(0,6),categories=c("Constant.negative","Constant.positive","Converging.negative","Converging.positive","not_converging.negative","not_converging.positive"))
            # cat$T <- aggregate(fmod_fft2$meanT, list(fmod_fft2$catsign), mean, na.rm=TRUE)[1:6,2]
            # cat$G <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$catsign), mean, na.rm=TRUE)[1:6,2]
            # cat$minT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$catsign), min, na.rm=TRUE)[1:6,2]
            # cat$maxT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$catsign),max, na.rm=TRUE)[1:6,2]
            # cat$minG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$catsign), min, na.rm=TRUE)[1:6,2]
            # cat$maxG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$catsign),max, na.rm=TRUE)[1:6,2]
            # cat$category <- c("Constant","Constant","Converging","Converging","Not Converging","Not Converging")
            # cat$signlofreq <- c("Negative","Positive","Negative","Positive","Negative","Positive")
            
            # cat6spread <- ggplot(data=cat, aes(x=T,y=G,xmin=minT,xmax=maxT,ymin=minG,ymax=maxG, shape=category, color = signlofreq,group=categories))+
            # geom_point()+geom_errorbar()+geom_errorbarh()+
            # scale_color_discrete(name="Direction of the effect")+
            # scale_color_brewer(name="Direction of the effect",palette="Dark2")+
            # scale_shape_discrete(name="Category")+
            # labs(x="Country's population-weighted temperature", y="log(GDP per capita)")+theme_bw()

            # cat6datapoints <- ggplot(data=fmod_fft2, aes(x=meanT,y=logGDPpc, shape=category, color = factor(signlofreq)))+
            # geom_point()+
            # scale_color_discrete(name="Direction of the effect")+
            # scale_color_brewer(name="Direction of the effect",palette="Dark2")+
            # scale_shape_discrete(name="Category")+
            # labs(x="Country's population-weighted temperature", y="log(GDP per capita)")+theme_bw()

            # ggarrange(cat6spread,cat6datapoints,ncol=1,common.legend=TRUE)
            # #ggsave("SupFig3_cat6.png",dpi=600)
            # cat <- data.frame(T = rep(0,3),minT=rep(0,3),maxT=rep(0,3),G = rep(0,3),minG=rep(0,3),maxG=rep(0,3),category=c("Constant","Converging","Not Converging"))
            # cat$T <- aggregate(fmod_fft2$meanT, list(fmod_fft2$category), mean, na.rm=TRUE)[,2]
            # cat$G <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$category), mean, na.rm=TRUE)[,2]
            # cat$minT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$category), min, na.rm=TRUE)[,2]
            # cat$maxT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$category),max, na.rm=TRUE)[,2]
            # cat$minG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$category), min, na.rm=TRUE)[,2]
            # cat$maxG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$category),max, na.rm=TRUE)[,2]
            # cat3spread <- ggplot(data=cat, aes(x=T,y=G,xmin=minT,xmax=maxT,ymin=minG,ymax=maxG,color=category))+
            # geom_point()+geom_errorbar()+geom_errorbarh()+
            # labs(x="Mean Temperature", y="logGDPpc")+theme_bw()

            # cat3datapoints <- ggplot(data=fmod_fft2, aes(x=meanT,y=logGDPpc, color=category))+
            # geom_point()+
            # labs(x="Mean Temperature", y="logGDPpc")+theme_bw()

            # ggarrange(densityGDP6,densitymeanT6,cat6spread,cat6datapoints,ncol=1,nrow=4)
            # #ggsave("soecioeconomic_6cat.png",dpi=600)
            # ggarrange(densityGDP3,densitymeanT3,cat3spread,cat3datapoints,ncol=1,nrow=4)
            # #ggsave("soecioeconomic_3cat.png",dpi=600)
            # head(fmod_fft2[fmod_fft2$filters=="Unfiltered",])
            glimpse(panel_data)
            glimpse(fmod_fft2)

            panel_data$countrycode <- as.factor(panel_data$countrycode)
            numobs <- aggregate( numobs ~countrycode,data=panel_data[panel_data$econdata=="wb",],FUN="mean")
            tempvar <- aggregate( tempsd ~countrycode,data=panel_data[panel_data$econdata=="wb",],FUN="mean")

            
            #panel_data$numobs[panel_data$countrycode=="AFG"]
            
            fmod_fft2_2 <- merge(fmod_fft2,numobs, by="countrycode")
            fmod_fft2_2 <- merge(fmod_fft2_2,tempvar, by="countrycode")

            library('ggrepel')
            glimpse(fmod_fft2_2)

            #a1=ggplot(fmod_fft2[fmod_fft2$filters=="Unfiltered",], 
            a1=ggplot(fmod_fft2_2[fmod_fft2$filters=="Unfiltered" & fmod_fft2$lowsignificant==1,], 
            aes(x=meanT,y=Estimate))+
            geom_text_repel(aes(label=iso_a3))+
            #geom_point()+
            geom_smooth(method = "lm",
                    mapping = aes(weight = (1/StandardError)),
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    labs(title="Temperature effect Filter = Unfiltered") +
                xlab("Mean temperature")
            

            b1=ggplot(fmod_fft2[fmod_fft2$filters=="Unfiltered"& fmod_fft2$unfilteredsignificant==1,],
            #aes(x=log(NY.GDP.PCAP.PP.KD),
            #b1=ggplot(fmod_fft2[fmod_fft2$filters=="Unfiltered",],
            aes(x=log(NY.GDP.PCAP.PP.KD),y=Estimate))+
            #geom_point()+
            geom_text_repel(aes(label=iso_a3))+#geom_errorbar()+
            geom_smooth(method = "lm", 
            mapping = aes(weight = (1/StandardError)),
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    xlab("log GDP per capita")+
                    labs(title="Temperature effect Filter = Unfiltered")
            
            
            c1 <- ggplot(fmod_fft2_2[fmod_fft2_2$filters=="Unfiltered",],
            aes(x=tempsd,y=abs(Estimate)/StandardError))+
            geom_point(aes(color=numobs, shape=factor(unfilteredsignificant)))+
            #geom_text_repel(aes(label=iso_a3))+#geom_errorbar()+
            geom_smooth(method = "lm", 
            #mapping = aes(weight = (1/StandardError)),
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) + 
            scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
                    guides(color=guide_legend("Number of observations"),shape=guide_legend("Statistically significant"))+
                    ylab("t-stat")+
                    xlab("SD of country-level temperature timeseries")+
                    ylim(0,4)+
                    labs(title="Temperature effect Filter = Unfiltered")

                 c <- ggplot(fmod_fft2_2[fmod_fft2_2$filters=="15 years",],
            aes(x=tempsd,y=abs(Estimate)/StandardError))+
            geom_point(aes(color=numobs, shape=factor(lowsignificant)))+
            #geom_text_repel(aes(label=iso_a3))+#geom_errorbar()+
            geom_smooth(method = "lm", 
            #mapping = aes(weight = (1/StandardError)),
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) + 
            scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
                    guides(color=guide_legend("Number of observations"),shape=guide_legend("Statistically significant"))+
                    ylab("t-stat")+
                    xlab("SD of country-level temperature timeseries")+
                    #xlim(c(0,4))+
                    labs(title="Temperature effect Filter = 15 years")
        

            a=ggplot(fmod_fft2[fmod_fft2$filters=="15 years" & fmod_fft2$lowsignificant==1 ,], 
            #a=ggplot(fmod_fft2[fmod_fft2$filters=="15 years" ,], 
            #a=ggplot(fmod_fft2[fmod_fft2$filters=="15 years"  & !is.na(fmod_fft2$NY.GDP.PCAP.PP.KD),], 
            #a=ggplot(fmod_fft2[fmod_fft2$filters=="15 years",], 
            aes(x=meanT,y=Estimate))+
            #geom_point()+ 
            geom_text_repel(aes(label=iso_a3))+
            xlab("Mean temperature")+
            #geom_errorbar(aes(ymin=low95,ymax=high95,width=0.1))+
            geom_smooth(method = "lm",
            #geom_smooth(method = "lm", formula=y ~ poly(x, 2),
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    labs(title="Temperature effect Filter = 15 year")
            

            b=ggplot(fmod_fft2[fmod_fft2$filters=="15 years"& fmod_fft2$lowsignificant==1,],aes(x=log(NY.GDP.PCAP.PP.KD),
            #b=ggplot(fmod_fft2[fmod_fft2$filters=="15 years",],aes(x=log(NY.GDP.PCAP.PP.KD),
            y=Estimate))+
            #geom_point()+
            #geom_errorbar()+
            geom_text_repel(aes(label=iso_a3))+
            geom_smooth(method = "lm", 
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    xlab("log GDP per capita")+
                    labs(title="Temperature effect Filter = 15 year")


          
            ggarrange(ggarrange(a1,a), ggarrange(b1, b), ggarrange(c1, c, common.legend=TRUE), ncol=1,nrow=3)
            

            glimpse(fmod_fft2)
            glimpse(df1)
            names(df1)
            fmod_fft2 <- df1
            names(fmod_fft2)
            fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
            fmod_fft2$signsign_filt <- paste(fmod_fft2$signlofreq,fmod_fft2$lowsignificant)
            fmod_fft2$gl <- "Gray Area"
            fmod_fft2$gl[fmod_fft2$lowsignificant ==1 ] <- "Growth Effect"         
            fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Levels Effect"
            
            
            alluvial_1 <- data.frame(fmod_fft2[,c(20,21,22)])
            es_alluvial <- plyr::count(alluvial_1)
            #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
            
            glimpse(es_alluvial)
            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 0"] <- "Negative not significant"
            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 1"] <- "Negative significant"
            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 0"] <- "Positive not significant"
            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 1"] <- "Positive significant"

            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

            


            glimpse(es_alluvial)
            es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
            es_alluvial$signsign_filt <- factor(es_alluvial$signsign_filt, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
            es_alluvial$gl <- factor(es_alluvial$gl, levels=c("Growth Effect","Levels Effect","Gray Area"))
            

            pal1 <- sequential_hcl(7, palette = "Purples")
                pal2 <- sequential_hcl(7, palette = "Greens")

            library('ggalluvial')
            ggplot(as.data.frame(es_alluvial),
            aes(y = freq, axis1 = factor(signsign_unf), axis2 = factor(signsign_filt))) +
            geom_alluvium(aes(fill = factor(gl)), width = 1/12) +
            geom_stratum(width = 1/12, aes(fill =  gl), color = "grey")+
            #geom_stratum() +
            geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
            scale_x_discrete(limits = c("Unfiltered", "15-year Filter"), expand = c(0.1, .1))  +
            ylab("Number of countries")+
            scale_fill_manual(values=c("#d3818c","#7375a4","gray"),name="Type of Impact",breaks = levels(es_alluvial$gl)) 
            ggsave("TempEffect_GvL.png",dpi=600)
       #+
       #ggtitle("Temperature effect country-level estimates")
        #guides(fill="none")
  
            #class(es_alluvial)
            #glimpse(es_alluvial)    
            
            table(fmod_fft2[,c(10,11)])
            
            fmod_fft2[,c(1,10,11)] %>% group_by(lowsignificant, unfilteredsignificant) %>%
            summarise(n = sum) -> tit2d


            glimpse(fmod_fft2)
            lodes <- to_lodes_form(as.data.frame(fmod_fft2[,c(1,10,11)]),
                           #axes = c(2,3),
                           id = "countrycode")
            glimpse(lodes)

            as.data.frame((fmod_fft2[,c(10,11)]))


             #ggsave("Estimates_meanTlogGDP.png",dpi=600)

             fmod_fft2[fmod_fft2$filters=="Unfiltered"& fmod_fft2$countrycode=="RUS",]
        # Socioeconomic characteristics across categories (end)
    ## 3.5. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2 (end)

    ## 3.6. Other datasets - Table 1, columns 3 to 6; Supp Fig 3 (start)
        # Barro
            ## 3.2. categorizing - Plotting Figure 3  (start)

                    # Categorizing statistically different estimates (start)
                        fmod_fft <- new_fft
                        fmod_fft$Variance <- fullmods_filter_var$Variance
                        
                        #fmod_fft <- fullmods_filter
                        fmod_fft <- fmod_fft[fmod_fft$econdata=="barro",]
                        fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                        fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                        fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                        numcountries <- length(unique(fmod_fft$countrycode))
                        #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                        fmod_fft$lowsignificant <- 0
                        fmod_fft$unfilteredsignificant <- 0
                        fmod_fft$signlofreq <- 0
                        fmod_fft$signunfilt <- 0
                                        
                        uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                        `%notin%` <- Negate(`%in%`)
                        fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                        
                        numcountries <- length(unique(fmod_fft$countrycode))
                        glimpse(fmod_fft)
                        for (i in 1:numcountries){
                            if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                            if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                lastfreq <- 5
                                }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                    lastfreq <- 4
                                }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                    lastfreq <- 3
                                }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                    lastfreq <- 2
                                }else{next}
                            if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                            if(fmod_fft$Estimate[1+5*(i-1)] >0){fmod_fft$signunfilt[(1+5*(i-1)):(5+5*(i-1))]=1}
                            m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                            if (is.na(m)){next}
                            if(m>0 ){
                                fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                                }
                        }

                            converging_stat <- 0
                            converging <- 0
                            not_converging <- 0
                            Intensifying <- 0
                            gray_area <- 0
                            for (i in 1:numcountries){
                                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                        lastfreq <- 5
                                        }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                            lastfreq <- 4
                                        }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                            lastfreq <- 3
                                        }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                            lastfreq <- 2
                                        }else{next}
                                        
                                        if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                                theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                                var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                                conf95 <-  (var^0.5)*1.65 #one-tail  95%
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                                    
                                    if(fmod_fft$sign[2+5*(i-1)]==0){gray_area <- c(gray_area,(1+5*(i-1)):(5+5*(i-1))) }else{
                                            not_converging <- c(not_converging,(1+5*(i-1)):(5+5*(i-1)))
                                    
                                    if(theta+conf95 < 0 ){
                                        Intensifying <- c(Intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                                    }}
                                }
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                                        converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                                        if(theta-conf95 >= 0 ){
                                        converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                    }
                                }
                                
                                
                                    
                                    if((lastfreq+5*(i-1)) %notin% uncertain){
                                    fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                    if((1+5*(i-1)) %notin% uncertain){
                                        fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                }
                                
                                
                            filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                            
                            fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                            fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                            fmod_fft$category <- "Undefined"
                            fmod_fft$category[converging] <- "converging"
                            fmod_fft$category[not_converging] <- "Not Converging"
                            fmod_fft$category[converging_stat] <- "converging_stat"
                            fmod_fft$category[Intensifying] <- "*Intensifying"
                            fmod_fft$category[gray_area] <- "gray_area"

                            table(fmod_fft$category)/5
                            
                            fmod_fft$Category <- "Undefined"
                            fmod_fft$Category[converging] <- "converging"
                            fmod_fft$Category[not_converging] <- "Not Converging"
                            
                            table(fmod_fft$Category)/5
                            
                            fmod_fft$CoeffDifferents <- 0
                            fmod_fft$CoeffDifferents[Intensifying] <- 1
                            fmod_fft$CoeffDifferents[converging_stat] <- 1
                            table(fmod_fft$CoeffDifferents)/5
                            
                            table(fmod_fft$lowsignificant)/5
                            table(fmod_fft$unfilteredsignificant)/5

                            fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                            library('rnaturalearth')
                            world <- ne_countries(scale = "medium", returnclass = "sf")

                            
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
                        
                            fmod_fft$iso_a3 <- fmod_fft$countrycode
                            fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                            fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                            fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                            

                            #factor(fmod_fft$category)
                            fmod_fft$category <- factor(fmod_fft$category, levels = c("*Intensifying", "Not Converging", "converging","converging_stat"))
                            df1 <- fmod_fft[which(fmod_fft$econdata=="barro" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                            pal1 <- sequential_hcl(7, palette = "Oranges")
                            pal2 <- sequential_hcl(7, palette = "Teal")
                            df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                            df1$lowsignificant <- factor(df1$lowsignificant)
                                
                            hist_bf_barro <- ggplot() + theme_bw() + 
                                geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                                #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                                stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Not Converging","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",
                                    labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                                    expression(paste("|",theta[f],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                            df2 <- fmod_fft[which(fmod_fft$econdata=="barro" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                            df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                            glimpse(df2)     
                            hist_bu_barro <- ggplot() + theme_bw()+
                                geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                                #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Not Converging","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                            
                            plot_fg_barro <- ggplot(data=fmod_fft,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(category)))+
                                geom_line()+
                                scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],pal2[1]),
                                labels=c("Converging","Not Converging","*Intensifying","Statistically Converging")) +
                                theme_bw() + xlab("Minimum Periodicity after Filtering")+
                                geom_hline(yintercept=0,lty=2)+
                                #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                                ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                                theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),legend.position="none") +
                                #panel.background = element_blank()) +
                                ggtitle("Barro-Ursua") 

                            table(fmod_fft$sign)/5
                            table(fmod_fft$sign)/5
                            table(fmod_fft$category[fmod_fft$sign==1])/5
                            table(fmod_fft$category[fmod_fft$sign==0])/5
                            table(fmod_fft$unfilteredsignificant[fmod_fft$sign==0])/5
                            
                            plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                            scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"),name="",labels=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified"))
                            leg_cat <- get_legend(plot_categoriescolors)                    
                            
                            fig_3_barro <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                            fig_3_barro
                            e_barro <- fmod_fft
                                
                            library("ggalluvial")

                            
                            fmod_fft2 <- fmod_fft
                            names(fmod_fft2)
                            fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                            fmod_fft2$signsign_filt <- paste(fmod_fft2$signlofreq,fmod_fft2$lowsignificant)
                            fmod_fft2$gl <- "Unclassified"
                            fmod_fft2$gl[fmod_fft2$lowsignificant ==1 ] <- "Presence of growth effects"         
                            fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                            
                            
                            alluvial_1 <- data.frame(fmod_fft2[,c(20,21,22)])
                            es_alluvial <- plyr::count(alluvial_1)
                            #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                            
                            glimpse(es_alluvial)
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 0"] <- "Negative not significant"
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 1"] <- "Negative significant"
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 0"] <- "Positive not significant"
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 1"] <- "Positive significant"

                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                            


                            glimpse(es_alluvial)
                            es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial$signsign_filt <- factor(es_alluvial$signsign_filt, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial$gl <- factor(es_alluvial$gl, levels=c("Presence of growth effects","Only level effects detected","Unclassified"))
                            

                            pal1 <- sequential_hcl(7, palette = "Purples")
                                pal2 <- sequential_hcl(7, palette = "Greens")

                            library('ggalluvial')
                            sankeyplot_pval <- ggplot(as.data.frame(es_alluvial),
                            aes(y = freq/5, axis1 = factor(signsign_unf), axis2 = factor(signsign_filt))) +
                            geom_alluvium(aes(fill = factor(gl)), width = 1/12) +
                            geom_stratum(width = 1/12, aes(fill =  gl), color = "grey")+
                            #geom_stratum() +
                            geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                            scale_x_discrete(limits = c("Unfiltered", "15-year filtered estimate"), expand = c(0.1, .1))  +
                            ylab("Number of countries")+
                            ggtitle("Evidence from trends")+
                            scale_fill_manual(values=c("#d3818c","#7375a4","gray"),name="Type of Impact",breaks = levels(es_alluvial$gl))+
                            theme(plot.title = element_text(hjust = 0.5)) +
                            theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
                            guides(fill = guide_legend(title.position="top", title.hjust = 0.5))
                            ggsave("TempEffect_GvL.png",dpi=600)

                            fmod_fft3 <- fmod_fft
                            names(fmod_fft3)
                            #fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                            fmod_fft3$signsign_filt <- paste(fmod_fft3$signlofreq,fmod_fft3$lowsignificant)
                            fmod_fft3$gl <- "Unclassified"
                            fmod_fft3$gl[fmod_fft3$category =="*Intensifying" ] <- "Presence of growth effects detected"  
                            fmod_fft3$gl[fmod_fft3$category =="Not Converging" ] <- "Presence of growth effects detected"   
                            fmod_fft3$gl[fmod_fft3$category =="converging" ] <- "Only level effects detected"       
                            #fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                            
                            fmod_fft3$cat <- "Unclassified"
                            fmod_fft3$cat[fmod_fft3$category =="*Intensifying" ] <- "*Intensifying"  
                            fmod_fft3$cat[fmod_fft3$category =="Not Converging" ] <- "Not converging to zero"   
                            fmod_fft3$cat[fmod_fft3$category =="converging" ] <- "Converging to zero"       

                            alluvial_2 <- data.frame(fmod_fft3[,c(20,21,22)])
                            es_alluvial2 <- plyr::count(alluvial_2)
                            #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                            
                            glimpse(es_alluvial2)
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 0"] <- "Negative not significant"
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 1"] <- "Negative significant"
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 0"] <- "Positive not significant"
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 1"] <- "Positive significant"


                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                            


                            glimpse(es_alluvial2)
                            #es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial2$signsign_filt <- factor(es_alluvial2$signsign_filt, 
                                levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial2$gl <- factor(es_alluvial2$gl, 
                                levels=c("Presence of growth effects detected","Only level effects detected","Unclassified"))
                            
                            es_alluvial2$cat <- factor(es_alluvial2$cat, 
                                levels=c("*Intensifying","Not converging to zero","Converging to zero","Unclassified"))
                            factor(es_alluvial2$cat)

                            pal1 <- sequential_hcl(7, palette = "Oranges")
                            pal2 <- sequential_hcl(7, palette = "Teal")

                            library('ggalluvial')
                            glimpse(es_alluvial2)

                            sankeyplot_cat <- ggplot(as.data.frame(es_alluvial2),
                            aes(y = freq/5, axis1 = factor(signsign_filt), axis2 = factor(cat))) +
                            geom_alluvium(aes(fill = factor(cat)), width = 1/12) +
                            geom_stratum(width = 1/5, aes(fill = cat), color = "grey")+
                            #geom_stratum() +
                            geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                            scale_x_discrete(limits = c("15-year filtered estimate", "Convergence to zero category"), expand = c(0.1, .1))  +
                            ylab("")+
                            scale_fill_manual(values=c(pal2[3],pal1[3],"dimgray"),name="Type of Impact") +
                            ggtitle("Evidence from point estimetaes")+
                            theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")
                            #guides(fill = guide_legend(title.position="top", title.hjust = 0.5))

                            ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom", widths = c(2, 1) )

                            #ggsave("sankey_diag_3.pdf",device="pdf",dpi=300)

                            ggarrange(plot_fg_barro,ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom" ),nrow=2,ncol=1)
                            ggsave("Figure3Barro_estimates_sankey.png",dpi=300)

                            #ggsave("Fig3_barro.png",dpi=600)
                                    
                    # Categorizing statistically different estimates (end)
            ## 3.2. categorizing - Plotting Figure 3 (end) 

            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
                fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                        origin = "iso3c",
                                        destination = "continent")
                fmod_fft$invse <- 1/fmod_fft$StandardError
                fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                f_unfiltpos <- fmod_fft[fmod_fft$signunfilt==1,]
                f_unfiltneg <- fmod_fft[fmod_fft$signunfilt==0,]
                felm_3 <- felm(Estimate ~ factor(frequencies)|0|0|continent, data =f_unfiltpos,weights = f_unfiltpos$invse)
                felm_4 <- felm(Estimate ~ factor(frequencies)|0|0|continent, data =f_unfiltneg,weights = f_unfiltneg$invse)
        
                #felm_3 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
                #felm_4 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (end)

            
            ## 3.3. LAGS categorizing - Plotting Figure 3   (start)

                # Categorizing statistically different estimates (start)
                
                    fmod_fft <- new_fft_l
                    #fmod_fft$Variance <- fullmods_filter_var$Variance
                    
                    #fmod_fft <- fullmods_filter
                    fmod_fft <- fmod_fft[fmod_fft$econdata=="barro",]
                    fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                    fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                    fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                    numcountries <- length(unique(fmod_fft$countrycode))
                    #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                    fmod_fft$lowsignificant <- 0
                    fmod_fft$unfilteredsignificant <- 0
                    fmod_fft$signlofreq <- 0
                    fmod_fft$signunfilt <- 0
                                    
                    uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                    `%notin%` <- Negate(`%in%`)
                    fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                    
                    numcountries <- length(unique(fmod_fft$countrycode))
                    glimpse(fmod_fft)
                    for (i in 1:numcountries){
                        if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                        if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                        if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                        if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                        if(fmod_fft$Estimate[1+5*(i-1)] >0){fmod_fft$signunfilt[(1+5*(i-1)):(5+5*(i-1))]=1}
                        m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                        if (is.na(m)){next}
                        if(m>0 ){
                            fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                            }
                    }

                        converging_stat <- 0
                        converging <- 0
                        not_converging <- 0
                        Intensifying <- 0
                        gray_area <- 0
                        for (i in 1:numcountries){
                            if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                    lastfreq <- 5
                                    }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                        lastfreq <- 4
                                    }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                        lastfreq <- 3
                                    }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                        lastfreq <- 2
                                    }else{next}
                                    
                                    if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                                    if(is.na(fmod_fft$StandardError[lastfreq+5*(i-1)])){next}
                            theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                            var <- fmod_fft$StandardError[lastfreq+5*(i-1)]^2 + fmod_fft$StandardError[1+5*(i-1)]^2
                            conf95 <-  (var^0.5)*1.65 #one-tail  95%
                            
                            if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                                
                                if(fmod_fft$sign[2+5*(i-1)]==0){gray_area <- c(gray_area,(1+5*(i-1)):(5+5*(i-1))) }else{
                                        not_converging <- c(not_converging,(1+5*(i-1)):(5+5*(i-1)))
                                
                                if(theta+conf95 < 0 ){
                                    Intensifying <- c(Intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                                }}
                            }
                            
                            if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                                    converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                                    if(theta-conf95 >= 0 ){
                                    converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                }
                            }
                            
                            
                                
                                if((lastfreq+5*(i-1)) %notin% uncertain){
                                fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                            }
                                if((1+5*(i-1)) %notin% uncertain){
                                    fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                            }
                            }
                            
                            
                        filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                        
                        #fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                        #fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                        fmod_fft$category <- "Undefined"
                        fmod_fft$category[converging] <- "converging"
                        fmod_fft$category[not_converging] <- "Not Converging"
                        fmod_fft$category[converging_stat] <- "converging_stat"
                        fmod_fft$category[Intensifying] <- "*Intensifying"
                        fmod_fft$category[gray_area] <- "gray_area"

                        table(fmod_fft$category)/5
                        
                        fmod_fft$Category <- "Undefined"
                        fmod_fft$Category[converging] <- "converging"
                        fmod_fft$Category[not_converging] <- "Not Converging"
                        table(fmod_fft$Category)/5
                        
                        fmod_fft$CoeffDifferents <- 0
                        fmod_fft$CoeffDifferents[Intensifying] <- 1
                        fmod_fft$CoeffDifferents[converging_stat] <- 1
                        table(fmod_fft$CoeffDifferents)/5
                        
                        table(fmod_fft$lowsignificant)/5
                        table(fmod_fft$unfilteredsignificant)/5

                        
                        table(fmod_fft$category[which(fmod_fft$unfilteredsignificant==1)])/5
                        table(fmod_fft$category[which(fmod_fft$lowsignificant==1)])/5
                        table(fmod_fft$category[which(fmod_fft$lowsignificant==0)])/5

                        #fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                            #fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                            #fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                        fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                        library('rnaturalearth')
                        world <- ne_countries(scale = "medium", returnclass = "sf")

                        
                        
                        fmod_fft$iso_a3 <- fmod_fft$countrycode
                        fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                        fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                        fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                        

                        #factor(fmod_fft$category)
                        fmod_fft$category <- factor(fmod_fft$category, levels = c("*Intensifying", "Not Converging", "converging","converging_stat","gray_area"))
                        df1 <- fmod_fft[which(fmod_fft$econdata=="barro" &fmod_fft$climdata=="UDel" & fmod_fft$Lags==0),]
                        pal1 <- sequential_hcl(7, palette = "Oranges")
                        pal2 <- sequential_hcl(7, palette = "Teal")
                        df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                        df1$lowsignificant <- factor(df1$lowsignificant)
                            
                        hist_bf <- ggplot() + theme_bw() + 
                            geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                            #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                            #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                            stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                            name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                            scale_x_discrete(name="",
                                labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                                expression(paste("|",theta[f],"| >0",sep=""))))+
                            scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                        df2 <- fmod_fft[which(fmod_fft$econdata=="barro" &fmod_fft$climdata=="UDel" & fmod_fft$Lags==0 & fmod_fft$unfilteredsignificant==1),]
                        df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                        hist_bu <- ggplot() + theme_bw()+
                            geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                            #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                            #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                            scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                            name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                            scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                            scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                        
                        fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(abs(fmod_fft$Estimate)>quantile(abs(fmod_fft$Estimate), na.rm=TRUE,0.99))],]
                        #fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
                    
                    plot_fg <- ggplot(data=fmod_fft,aes(x=Lags,y=Estimate*100, group = countrycode,color=factor(category)))+
                            geom_line()+
                            scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],"dimgray",pal2[1])) +
                            theme_bw() + xlab("Number of Lags")+
                            geom_hline(yintercept=0,lty=2)+
                            #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                            ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                            theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            panel.background = element_blank(),legend.position="none") +
                            #panel.background = element_blank()) +
                            ggtitle("")
                            plot_fg      

                        table(fmod_fft$sign)/5
                        table(fmod_fft$sign)/5
                        table(fmod_fft$category[fmod_fft$sign==1])/5
                        table(fmod_fft$category[fmod_fft$sign==0])/5
                        table(fmod_fft$unfilteredsignificant[fmod_fft$sign==0])/5
                        
                        plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                        scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"),name="",labels=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified"))
                        leg_cat <- get_legend(plot_categoriescolors)                    
                        
                        fig_3_b_l <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                        fig_3_b_l
                        e_barro_lags <- fmod_fft
                        ggsave("Figures/Fig3_barro_lags.png",dpi=600)
                        
                        # map_categories <- ggplot(data=fmod_fft_map) +
                        #             geom_sf(data = fmod_fft_map_na,aes(fill=NA))+
                        #             theme_minimal()+
                        #             geom_sf(aes(fill = category))+
                        #             scale_fill_manual(name = "",
                        #                             #labels=c("Converging","Not Converging","*Intensifying"),
                        #                             values=c(pal2[3],pal1[3],pal1[1]))+
                        #             theme(legend.position="none")+
                        #             ggtitle("Behavior of coefficients")
                        #     ggarrange(map_categories,leg_cat,ncol=2,nrow=1,widths=c(5,1))
                            
                        #     #ggsave("Fig3_map.png",dpi=600)

                        # fmod_fft_map_low <- fmod_fft_map[which(fmod_fft_map$filters=="15 years"),]
                        # missingcountries15 <- fmod_fft_map_low$countrycode[is.na(fmod_fft_map_low$Estimate)]
                        # fmod_fft_map_10 <- fmod_fft_map[which(fmod_fft_map$filters=="10 years"),]
                        
                        # fmod_fft_map_low$Estimate[which(fmod_fft_map_low$countrycode %in% missingcountries15)] <- fmod_fft_map_10$Estimate[which(fmod_fft_map_10$countrycode %in% missingcountries15)]
                        # ggplot(data = fmod_fft_map_low) +
                        #                 geom_sf(data=fmod_fft_map_na,fill=NA)+theme_minimal()+
                        #                 geom_sf(data=fmod_fft_map_low,aes(fill = Estimate*100))+
                        #                 scale_fill_gradient2(
                        #                     name = "Estimated impact \n (% per Degree)",
                        #                     low = "red",
                        #                     mid = "white",
                        #                     high = "#00BFC4",
                        #                     midpoint = 0,
                        #                     space = "Lab",
                        #                     na.value = "gray",
                        #                     guide = "colourbar",
                        #                     aesthetics = "fill")+
                        #                     ggtitle("Growth effects")+
                        #                 theme(legend.position="bottom")
                            #ggsave("Map_correcteddata.png",dpi=600)
                                
                # Categorizing statistically different estimates (end)
            ## 3.3. LAGS categorizing - Plotting Figure 3 (end) 

            ## 3.4. LAGS FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
                fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                        origin = "iso3c",
                                        destination = "continent")
                fmod_fft$invse <- 1/fmod_fft$StandardError
                fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                f_unfiltpos <- fmod_fft[fmod_fft$signunfilt==1,]
                f_unfiltneg <- fmod_fft[fmod_fft$signunfilt==0,]
                
                #f_samesign <- f_samesign[f_samesign$signlofreq==1,]
                felm_3_l <- felm(Estimate ~ factor(Lags)|0|0|continent, data =f_unfiltpos,weights = f_unfiltpos$invse)
                felm_4_l <- felm(Estimate ~ factor(Lags)|0|0|continent, data =f_unfiltneg,weights = f_unfiltneg$invse)
                
                #stargazer(felm_2,felm_1, type="html", out="felm_1_2.html")
                stargazer(felm_3_l,felm_4_l,type="text")
                stargazer(felm_2_l,type="text")
                #write.csv(fmod_fft,'Growth_Estimates_Filters.csv')

                #felm_1_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft,weights = fmod_fft$invse)
                #felm_2_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft)
                #summary(felm_1_within)
                #stargazer(felm_2_within,felm_1_within,type="text")
                
                #stargazer(felm_2,felm_1,felm_2_within,felm_1_within, type="html", out="felm_1_2_within.html")
                
            ## 3.4. LAGS FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
    
        # Barro  

        # Maddison
            ## 3.2. categorizing - Plotting Figure 3  (start)

                    # Categorizing statistically different estimates (start)
                        fmod_fft <- new_fft
                        fmod_fft$Variance <- fullmods_filter_var$Variance
                        
                        #fmod_fft <- fullmods_filter
                        fmod_fft <- fmod_fft[fmod_fft$econdata=="mad",]
                        fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                        fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.96
                        fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.96
                        numcountries <- length(unique(fmod_fft$countrycode))
                        #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                        fmod_fft$lowsignificant <- 0
                        fmod_fft$unfilteredsignificant <- 0
                        fmod_fft$signlofreq <- 0
                        fmod_fft$signunfilt <- 0
                                        
                        uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                        `%notin%` <- Negate(`%in%`)
                        fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                        
                        numcountries <- length(unique(fmod_fft$countrycode))
                        glimpse(fmod_fft)
                        for (i in 1:numcountries){
                            if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                            if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                lastfreq <- 5
                                }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                    lastfreq <- 4
                                }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                    lastfreq <- 3
                                }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                    lastfreq <- 2
                                }else{next}
                            if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                            if(fmod_fft$Estimate[1+5*(i-1)] >0){fmod_fft$signunfilt[(1+5*(i-1)):(5+5*(i-1))]=1}
                            m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                            if (is.na(m)){next}
                            if(m>0 ){
                                fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                                }
                        }

                            converging_stat <- 0
                            converging <- 0
                            not_converging <- 0
                            Intensifying <- 0
                            gray_area <- 0
                            for (i in 1:numcountries){
                                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                        lastfreq <- 5
                                        }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                            lastfreq <- 4
                                        }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                            lastfreq <- 3
                                        }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                            lastfreq <- 2
                                        }else{next}
                                        
                                        if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                                theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                                var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                                conf95 <-  (var^0.5)*1.65 #one-tail  95%
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                                    
                                    if(fmod_fft$sign[2+5*(i-1)]==0){gray_area <- c(gray_area,(1+5*(i-1)):(5+5*(i-1))) }else{
                                            not_converging <- c(not_converging,(1+5*(i-1)):(5+5*(i-1)))
                                    
                                    if(theta+conf95 < 0 ){
                                        Intensifying <- c(Intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                                    }}
                                }
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                                        converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                                        if(theta-conf95 >= 0 ){
                                        converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                    }
                                }
                                
                                
                                    
                                    if((lastfreq+5*(i-1)) %notin% uncertain){
                                    fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                    if((1+5*(i-1)) %notin% uncertain){
                                        fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                }
                                
                                
                            filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                            
                            fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                            fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                            fmod_fft$category <- "Undefined"
                            fmod_fft$category[converging] <- "converging"
                            fmod_fft$category[not_converging] <- "Not Converging"
                            fmod_fft$category[converging_stat] <- "converging_stat"
                            fmod_fft$category[Intensifying] <- "*Intensifying"
                            fmod_fft$category[gray_area] <- "gray_area"

                            table(fmod_fft$category)/5
                            
                            fmod_fft$Category <- "Undefined"
                            fmod_fft$Category[converging] <- "converging"
                            fmod_fft$Category[not_converging] <- "Not Converging"
                            
                            table(fmod_fft$Category)/5
                            
                            fmod_fft$CoeffDifferents <- 0
                            fmod_fft$CoeffDifferents[Intensifying] <- 1
                            fmod_fft$CoeffDifferents[converging_stat] <- 1
                            table(fmod_fft$CoeffDifferents)/5
                            
                            table(fmod_fft$lowsignificant)/5
                            table(fmod_fft$unfilteredsignificant)/5

                            fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                            library('rnaturalearth')
                            world <- ne_countries(scale = "medium", returnclass = "sf")

                            
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
                        
                            fmod_fft$iso_a3 <- fmod_fft$countrycode
                            fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                            fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                            fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                            

                            #factor(fmod_fft$category)
                            fmod_fft$category <- factor(fmod_fft$category, levels = c("*Intensifying", "Not Converging", "converging","converging_stat"))
                            df1 <- fmod_fft[which(fmod_fft$econdata=="mad" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                            pal1 <- sequential_hcl(7, palette = "Oranges")
                            pal2 <- sequential_hcl(7, palette = "Teal")
                            df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                            df1$lowsignificant <- factor(df1$lowsignificant)
                                
                            hist_bf_mad <- ggplot() + theme_bw() + 
                                geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                                #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                                stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Not Converging","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",
                                    labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                                    expression(paste("|",theta[f],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                            df2 <- fmod_fft[which(fmod_fft$econdata=="mad" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                            df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                            glimpse(df2)     
                            hist_bu_mad <- ggplot() + theme_bw()+
                                geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                                #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Not Converging","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                            
                            plot_fg_mad <- ggplot(data=fmod_fft,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(category)))+
                                geom_line()+
                                scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],pal2[1]),
                                labels=c("Converging","Not Converging","*Intensifying","Statistically Converging")) +
                                theme_bw() + xlab("Minimum Periodicity after Filtering")+
                                geom_hline(yintercept=0,lty=2)+
                                #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                                ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                                theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),legend.position="none") +
                                #panel.background = element_blank()) +
                                ggtitle("Maddison")

                            table(fmod_fft$sign)/5
                            table(fmod_fft$sign)/5
                            table(fmod_fft$category[fmod_fft$sign==1])/5
                            table(fmod_fft$category[fmod_fft$sign==0])/5
                            table(fmod_fft$unfilteredsignificant[fmod_fft$sign==0])/5
                            
                            plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                            scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"),name="",labels=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified"))
                            leg_cat <- get_legend(plot_categoriescolors)                    
                            
                            fig_3_mad <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                            fig_3_mad
                            e_mad <- fmod_fft
                            #ggarrange(ggarrange(plot_fg_mad,plot_fg_barro,ncol=2,nrow=1),ggarrange(hist_bu_mad,hist_bf_mad,leg_cat,hist_bu_barro,hist_bf_barro,ncol=5,nrow=1),ncol=1,nrow=2)
                            
                            library("ggalluvial")

                            
                            fmod_fft2 <- fmod_fft
                            names(fmod_fft2)
                            fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                            fmod_fft2$signsign_filt <- paste(fmod_fft2$signlofreq,fmod_fft2$lowsignificant)
                            fmod_fft2$gl <- "Unclassified"
                            fmod_fft2$gl[fmod_fft2$lowsignificant ==1 ] <- "Presence of growth effects"         
                            fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                            
                            
                            alluvial_1 <- data.frame(fmod_fft2[,c(20,21,22)])
                            es_alluvial <- plyr::count(alluvial_1)
                            #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                            
                            glimpse(es_alluvial)
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 0"] <- "Negative not significant"
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="0 1"] <- "Negative significant"
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 0"] <- "Positive not significant"
                            es_alluvial$signsign_filt[es_alluvial$signsign_filt=="1 1"] <- "Positive significant"

                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                            es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                            


                            glimpse(es_alluvial)
                            es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial$signsign_filt <- factor(es_alluvial$signsign_filt, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial$gl <- factor(es_alluvial$gl, levels=c("Presence of growth effects","Only level effects detected","Unclassified"))
                            

                            pal1 <- sequential_hcl(7, palette = "Purples")
                                pal2 <- sequential_hcl(7, palette = "Greens")

                            library('ggalluvial')
                            sankeyplot_pval <- ggplot(as.data.frame(es_alluvial),
                            aes(y = freq/5, axis1 = factor(signsign_unf), axis2 = factor(signsign_filt))) +
                            geom_alluvium(aes(fill = factor(gl)), width = 1/12) +
                            geom_stratum(width = 1/12, aes(fill =  gl), color = "grey")+
                            #geom_stratum() +
                            geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                            scale_x_discrete(limits = c("Unfiltered", "15-year filtered estimate"), expand = c(0.1, .1))  +
                            ylab("Number of countries")+
                            ggtitle("Evidence from trends")+
                            scale_fill_manual(values=c("#d3818c","#7375a4","gray"),name="Type of Impact",breaks = levels(es_alluvial$gl))+
                            theme(plot.title = element_text(hjust = 0.5)) +
                            theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
                            guides(fill = guide_legend(title.position="top", title.hjust = 0.5))
                            #ggsave("TempEffect_GvL.png",dpi=600)

                            fmod_fft3 <- fmod_fft
                            names(fmod_fft3)
                            #fmod_fft2$signsign_unf <- paste(fmod_fft2$signunfilt,fmod_fft2$unfilteredsignificant)
                            fmod_fft3$signsign_filt <- paste(fmod_fft3$signlofreq,fmod_fft3$lowsignificant)
                            fmod_fft3$gl <- "Unclassified"
                            fmod_fft3$gl[fmod_fft3$category =="*Intensifying" ] <- "Presence of growth effects detected"  
                            fmod_fft3$gl[fmod_fft3$category =="Not Converging" ] <- "Presence of growth effects detected"   
                            fmod_fft3$gl[fmod_fft3$category =="converging" ] <- "Only level effects detected"       
                            #fmod_fft2$gl[fmod_fft2$unfilteredsignificant == 1 & fmod_fft2$lowsignificant ==0] <- "Only level effects detected"
                            
                            fmod_fft3$cat <- "Unclassified"
                            fmod_fft3$cat[fmod_fft3$category =="*Intensifying" ] <- "Intensifying"  
                            fmod_fft3$cat[fmod_fft3$category =="Not Converging" ] <- "Not converging to zero"   
                            fmod_fft3$cat[fmod_fft3$category =="converging" ] <- "Converging to zero"       

                            alluvial_2 <- data.frame(fmod_fft3[,c(20,21,22)])
                            es_alluvial2 <- plyr::count(alluvial_2)
                            #names(es_alluvial) <- c("ES.Category","ES.General","ES.Specific","ES.Sub.Specific","Valuation.General","Freq")
                            
                            glimpse(es_alluvial2)
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 0"] <- "Negative not significant"
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="0 1"] <- "Negative significant"
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 0"] <- "Positive not significant"
                            es_alluvial2$signsign_filt[es_alluvial2$signsign_filt=="1 1"] <- "Positive significant"


                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 0"] <- "Negative not significant"
                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="0 1"] <- "Negative significant"
                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 0"] <- "Positive not significant"
                            #es_alluvial$signsign_unf[es_alluvial$signsign_unf=="1 1"] <- "Positive significant"

                            


                            glimpse(es_alluvial2)
                            #es_alluvial$signsign_unf <- factor(es_alluvial$signsign_unf, levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial2$signsign_filt <- factor(es_alluvial2$signsign_filt, 
                                levels=c("Positive significant","Positive not significant","Negative not significant","Negative significant"))
                            es_alluvial2$gl <- factor(es_alluvial2$gl, 
                                levels=c("Presence of growth effects detected","Only level effects detected","Unclassified"))
                            
                            es_alluvial2$cat <- factor(es_alluvial2$cat, 
                                levels=c("Intensifying","Not converging to zero","Converging to zero","Unclassified"))
                            factor(es_alluvial2$cat)

                            pal1 <- sequential_hcl(7, palette = "Oranges")
                            pal2 <- sequential_hcl(7, palette = "Teal")

                            library('ggalluvial')
                            glimpse(es_alluvial2)

                            sankeyplot_cat <- ggplot(as.data.frame(es_alluvial2),
                            aes(y = freq/5, axis1 = factor(signsign_filt), axis2 = factor(cat))) +
                            geom_alluvium(aes(fill = factor(cat)), width = 1/12) +
                            geom_stratum(width = 1/5, aes(fill = cat), color = "grey")+
                            #geom_stratum() +
                            geom_label(stat = "stratum", aes(label = after_stat(stratum))) + theme_minimal()+
                            scale_x_discrete(limits = c("15-year filtered estimate", "Convergence to zero category"), expand = c(0.1, .1))  +
                            ylab("")+
                            scale_fill_manual(values=c(pal2[3],pal1[1],pal1[3],"dimgray"),name="Type of Impact") +
                            ggtitle("Evidence from point estimetaes")+
                            theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")
                            #guides(fill = guide_legend(title.position="top", title.hjust = 0.5))

                            ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom", widths = c(2, 1) )

                            #ggsave("sankey_diag_3.pdf",device="pdf",dpi=300)

                            fig3_mad <- ggarrange(plot_fg_mad,ggarrange(sankeyplot_pval,sankeyplot_cat,common.legend=TRUE,legend="bottom" ),nrow=2,ncol=1)
                            ggsave("Figure3Mad_estimates_sankey.png",dpi=300)
                                    
                    # Categorizing statistically different estimates (end)
            ## 3.2. categorizing - Plotting Figure 3 (end) 

            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
                fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                        origin = "iso3c",
                                        destination = "continent")
                fmod_fft$invse <- 1/fmod_fft$StandardError
                fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                f_unfiltpos <- fmod_fft[fmod_fft$signunfilt==1,]
                f_unfiltneg <- fmod_fft[fmod_fft$signunfilt==0,]
                felm_5 <- felm(Estimate ~ factor(frequencies)|0|0|continent, data =f_unfiltpos,weights = f_unfiltpos$invse)
                felm_6 <- felm(Estimate ~ factor(frequencies)|0|0|continent, data =f_unfiltneg,weights = f_unfiltneg$invse)
                #felm_5 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
                #felm_6 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
                stargazer(felm_2,felm_1,felm_4,felm_3,felm_6,felm_5,type="html",out="Regression_on_estimates_negpos.html")
                
                dataset <- data.frame(estimate=double(),cil=double(),cih=double(),filter=factor(),dataset=factor(),sign=factor())
                felm_all <- felm(Estimate~signunfilt*factor(frequencies), data=fmod_fft)
                stargazer(felm_all,type="text")
                for(i in 1:3){
                for (j in 1:2){    
                    for (fi in 1:5){
                    if(i==1){
                        datasetname <- "World Bank"
                        if(j==1){
                            model <- felm_1
                            sign <- factor("positive")
                        } else {model <- felm_2
                        sign <- factor("negative")}
                    }
                    if(i==2){
                        datasetname <- "Barro-Ursua"
                        if(j==1){
                            model <- felm_3
                            sign <- factor("positive")
                        } else {model <- felm_4
                        sign <- factor("negative")}
                    }
                    if(i==3){
                        datasetname <- "Maddison"
                        if(j==1){
                            model <- felm_5
                            sign <- factor("positive")
                        } else {model <- felm_6
                        sign <- factor("negative")}
                    }
                    Sigma <- vcov(model)
                    coefT <- "(Intercept)"
                    if(fi==1){
                        sigma = Sigma[1,1]
                        beta.hat <- coef(model)[1]
                        xmat <-1
                        gestimated <- colSums(beta.hat*t(xmat)) 
                        cih <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        cil <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        Filter <- "Unfiltered"
                    }else{
                    if(fi==2){
                        coefT5 <- "factor(frequencies)3"
                        Filter <- "3-years"
                    }
                    if(fi==3){
                        coefT5 <- "factor(frequencies)5"
                        Filter <- "5-years"
                    }
                    if(fi==4){
                        coefT5 <- "factor(frequencies)10"
                        Filter <- "10-years"
                    }
                    if(fi==5){
                        coefT5 <- "factor(frequencies)15"
                        Filter <- "15-years"
                    }
                        start1 <- which(names(coef(model))==coefT)
                        end1 <- which(names(coef(model))==coefT5)
                        sigma = Sigma[c(start1,end1),c(start1,end1)]
                        beta.hat <- coef(model)[c(start1,end1)]
                        xmat <-cbind(1,1)
                        gestimated <- colSums(beta.hat*t(xmat)) 
                        cih <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        cil <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        }
                        dataset <- rbind(dataset,data.frame(estimate=gestimated,cil=cil,cih=cih,filter=factor(Filter),dataset=factor(datasetname),sign=factor(sign)))
                }
                }    
                }
                #dataset$filter <- factor(dataset$filter, levels = c("Unfiltered", "5 years", "10 years", "15 years"))
                
                
                cols <- c("#66c2a5","#fc8d62","#8da0cb")

                a <- ggplot(dataset[which(dataset$dataset=='World Bank'),],aes(x=filter,y=estimate*100,fill=dataset,color=dataset,group=interaction(sign,dataset)))+
                geom_line(color="#009E73")+theme_bw()+
                geom_hline(yintercept=0)+
                xlab("Filters")+ylab("Estimated pooled effect of \n 1 Degree Warming on Growth (pp)")+
                #scale_colour_manual(labels = c("a","b","c"),values = c(cols[1],cols[2],cols[3]))+
                geom_ribbon(aes(ymin=cih*100,ymax=cil*100),alpha=0.2,colour=NA,fill="#009E73")+
                ggtitle("World Bank")+
                ylim(-2.5,4.1)

                b <- ggplot(dataset[which(dataset$dataset=='Barro-Ursua'),],aes(x=filter,y=estimate*100,fill=dataset,color=dataset,group=interaction(sign,dataset)))+
                geom_line(color="#E69F00")+theme_bw()+
                geom_hline(yintercept=0)+
                xlab("Filters")+ylab("")+
                #scale_colour_manual(labels = c("a","b","c"),values = c(cols[1],cols[2],cols[3]))+
                geom_ribbon(aes(ymin=cih*100,ymax=cil*100),alpha=0.2,colour=NA,fill="#E69F00")+
                ggtitle("Barro Ursua")+
                ylim(-2.5,4.1)
                
                c <- ggplot(dataset[which(dataset$dataset=='Maddison'),],aes(x=filter,y=estimate*100,fill=dataset,color=dataset,group=interaction(sign,dataset)))+
                geom_line(color="#CC79A7")+theme_bw()+
                geom_hline(yintercept=0)+
                xlab("Filters")+ylab("")+
                #scale_colour_manual(labels = c("a","b","c"),values = c(cols[1],cols[2],cols[3]))+
                geom_ribbon(aes(ymin=cih*100,ymax=cil*100),alpha=0.2,colour=NA,fill="#CC79A7")+
                ggtitle("Maddison")+
                ylim(-2.5,4.1)
                


                ggarrange(a,b,c,ncol=3,nrow=1)
                #ggsave("Fig4_PooledEffect_NOWEIGHTS.png")

                #e_wb <- e_wb[which(e_wb$lowsignificant==1),]
                #table(e_wb$sign)
                
                #e_barro <- new_fft[new_fft$econdata=="barro",]
                #e_mad <- new_fft[new_fft$econdata=="mad",]

                #est <- merge(e_wb,e_barro, by=c("countrycode","frequencies"))
                #est <- est[which(est$Estimate.x>quantile(est$Estimate.x, na.rm=TRUE,0.01)),] 
                #est <- est[which(est$Estimate.x<quantile(est$Estimate.x, na.rm=TRUE,0.99)),] 
                #est <- est[which(est$Estimate.y>quantile(est$Estimate.y, na.rm=TRUE,0.01)),] 
                #est <- est[which(est$Estimate.y<quantile(est$Estimate.y, na.rm=TRUE,0.99)),] 
                    
                #ggplot(est, aes(x=Estimate.x,y=Estimate.y,color=factor(frequencies)))+
                #geom_point()+
                #facet_wrap(~frequencies)+
                #xlim(-0.25,0.25)+ylim(-0.25,0.25)
                #xlab("WB Estimates") + ylab("Maddison Estimates")

                #summary(lm(Estimate.x~0+Estimate.y,data=est))
            
            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (end)

            ## 3.3. LAGS categorizing - Plotting Figure 3   (start)

                # Categorizing statistically different estimates (start)
                
                    fmod_fft <- new_fft_l
                    #fmod_fft$Variance <- fullmods_filter_var$Variance
                    
                    #fmod_fft <- fullmods_filter
                    fmod_fft <- fmod_fft[fmod_fft$econdata=="mad",]
                    fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                    fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                    fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                    numcountries <- length(unique(fmod_fft$countrycode))
                    #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                    fmod_fft$lowsignificant <- 0
                    fmod_fft$unfilteredsignificant <- 0
                    fmod_fft$signlofreq <- 0
                    fmod_fft$signunfilt <- 0
                                    
                    uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                    `%notin%` <- Negate(`%in%`)
                    fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                    
                    numcountries <- length(unique(fmod_fft$countrycode))
                    glimpse(fmod_fft)
                    for (i in 1:numcountries){
                        if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                        if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                        if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                        if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                        if(fmod_fft$Estimate[1+5*(i-1)] >0){fmod_fft$signunfilt[(1+5*(i-1)):(5+5*(i-1))]=1}
                        m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                        if (is.na(m)){next}
                        if(m>0 ){
                            fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                            }
                    }

                        converging_stat <- 0
                        converging <- 0
                        not_converging <- 0
                        Intensifying <- 0
                        gray_area <- 0
                        for (i in 1:numcountries){
                            if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                    lastfreq <- 5
                                    }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                        lastfreq <- 4
                                    }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                        lastfreq <- 3
                                    }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                        lastfreq <- 2
                                    }else{next}
                                    
                                    if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                                    if(is.na(fmod_fft$StandardError[lastfreq+5*(i-1)])){next}
                            theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                            var <- fmod_fft$StandardError[lastfreq+5*(i-1)]^2 + fmod_fft$StandardError[1+5*(i-1)]^2
                            conf95 <-  (var^0.5)*1.65 #one-tail  95%
                            
                            if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                                
                                if(fmod_fft$sign[2+5*(i-1)]==0){gray_area <- c(gray_area,(1+5*(i-1)):(5+5*(i-1))) }else{
                                        not_converging <- c(not_converging,(1+5*(i-1)):(5+5*(i-1)))
                                
                                if(theta+conf95 < 0 ){
                                    Intensifying <- c(Intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                                }}
                            }
                            
                            if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                                    converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                                    if(theta-conf95 >= 0 ){
                                    converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                }
                            }
                            
                            
                                
                                if((lastfreq+5*(i-1)) %notin% uncertain){
                                fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                            }
                                if((1+5*(i-1)) %notin% uncertain){
                                    fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                            }
                            }
                            
                            
                        filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                        
                        #fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                        #fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                        fmod_fft$category <- "Undefined"
                        fmod_fft$category[converging] <- "converging"
                        fmod_fft$category[not_converging] <- "Not Converging"
                        fmod_fft$category[converging_stat] <- "converging_stat"
                        fmod_fft$category[Intensifying] <- "*Intensifying"
                        fmod_fft$category[gray_area] <- "gray_area"

                        table(fmod_fft$category)/5
                        
                        fmod_fft$Category <- "Undefined"
                        fmod_fft$Category[converging] <- "converging"
                        fmod_fft$Category[not_converging] <- "Not Converging"
                        table(fmod_fft$Category)/5
                        
                        fmod_fft$CoeffDifferents <- 0
                        fmod_fft$CoeffDifferents[Intensifying] <- 1
                        fmod_fft$CoeffDifferents[converging_stat] <- 1
                        table(fmod_fft$CoeffDifferents)/5
                        
                        table(fmod_fft$lowsignificant)/5
                        table(fmod_fft$unfilteredsignificant)/5

                        
                        table(fmod_fft$category[which(fmod_fft$unfilteredsignificant==1)])/5
                        table(fmod_fft$category[which(fmod_fft$lowsignificant==1)])/5
                        table(fmod_fft$category[which(fmod_fft$lowsignificant==0)])/5

                        #fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                            #fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                            #fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                        fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                        library('rnaturalearth')
                        world <- ne_countries(scale = "medium", returnclass = "sf")

                        
                        
                        fmod_fft$iso_a3 <- fmod_fft$countrycode
                        fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                        fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                        fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                        

                        #factor(fmod_fft$category)
                        fmod_fft$category <- factor(fmod_fft$category, levels = c("*Intensifying", "Not Converging", "converging","converging_stat","gray_area"))
                        df1 <- fmod_fft[which(fmod_fft$econdata=="mad" &fmod_fft$climdata=="UDel" & fmod_fft$Lags==0),]
                        pal1 <- sequential_hcl(7, palette = "Oranges")
                        pal2 <- sequential_hcl(7, palette = "Teal")
                        df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                        df1$lowsignificant <- factor(df1$lowsignificant)
                            
                        hist_bf <- ggplot() + theme_bw() + 
                            geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                            #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                            #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                            stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                            name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                            scale_x_discrete(name="",
                                labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                                expression(paste("|",theta[f],"| >0",sep=""))))+
                            scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                        df2 <- fmod_fft[which(fmod_fft$econdata=="mad" &fmod_fft$climdata=="UDel" & fmod_fft$Lags==0 & fmod_fft$unfilteredsignificant==1),]
                        df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                        hist_bu <- ggplot() + theme_bw()+
                            geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                            #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                            #name="",labels=c("*Intensifying","Not Converging","Converging","*Converging"))+
                            scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],"dimgray"), 
                            name="",labels=c("*Intensifying","Not Converging","Converging","Unclassified"))+ #No converging statistically significant found
                            scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                            scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                        
                        fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(abs(fmod_fft$Estimate)>quantile(abs(fmod_fft$Estimate), na.rm=TRUE,0.99))],]
                        #fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
                    
                    plot_fg <- ggplot(data=fmod_fft,aes(x=Lags,y=Estimate*100, group = countrycode,color=factor(category)))+
                            geom_line()+
                            scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],"dimgray",pal2[1])) +
                            theme_bw() + xlab("Number of Lags")+
                            geom_hline(yintercept=0,lty=2)+
                            #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                            ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                            theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            panel.background = element_blank(),legend.position="none") +
                            #panel.background = element_blank()) +
                            ggtitle("")
                            plot_fg      

                        table(fmod_fft$sign)/5
                        table(fmod_fft$sign)/5
                        table(fmod_fft$category[fmod_fft$sign==1])/5
                        table(fmod_fft$category[fmod_fft$sign==0])/5
                        table(fmod_fft$unfilteredsignificant[fmod_fft$sign==0])/5
                        
                        plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                        scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1],"dimgray"),name="",labels=c("*Intensifying","Not Converging","Converging","*Converging","Unclassified"))
                        leg_cat <- get_legend(plot_categoriescolors)                    
                        
                        fig_3_mad <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                        fig_3_mad
                        e_mad_lags <- fmod_fft
                        ggsave("Figures/Fig3_mad_lags.png",dpi=600)
                        
                        # map_categories <- ggplot(data=fmod_fft_map) +
                        #             geom_sf(data = fmod_fft_map_na,aes(fill=NA))+
                        #             theme_minimal()+
                        #             geom_sf(aes(fill = category))+
                        #             scale_fill_manual(name = "",
                        #                             #labels=c("Converging","Not Converging","*Intensifying"),
                        #                             values=c(pal2[3],pal1[3],pal1[1]))+
                        #             theme(legend.position="none")+
                        #             ggtitle("Behavior of coefficients")
                        #     ggarrange(map_categories,leg_cat,ncol=2,nrow=1,widths=c(5,1))
                            
                        #     #ggsave("Fig3_map.png",dpi=600)

                        # fmod_fft_map_low <- fmod_fft_map[which(fmod_fft_map$filters=="15 years"),]
                        # missingcountries15 <- fmod_fft_map_low$countrycode[is.na(fmod_fft_map_low$Estimate)]
                        # fmod_fft_map_10 <- fmod_fft_map[which(fmod_fft_map$filters=="10 years"),]
                        
                        # fmod_fft_map_low$Estimate[which(fmod_fft_map_low$countrycode %in% missingcountries15)] <- fmod_fft_map_10$Estimate[which(fmod_fft_map_10$countrycode %in% missingcountries15)]
                        # ggplot(data = fmod_fft_map_low) +
                        #                 geom_sf(data=fmod_fft_map_na,fill=NA)+theme_minimal()+
                        #                 geom_sf(data=fmod_fft_map_low,aes(fill = Estimate*100))+
                        #                 scale_fill_gradient2(
                        #                     name = "Estimated impact \n (% per Degree)",
                        #                     low = "red",
                        #                     mid = "white",
                        #                     high = "#00BFC4",
                        #                     midpoint = 0,
                        #                     space = "Lab",
                        #                     na.value = "gray",
                        #                     guide = "colourbar",
                        #                     aesthetics = "fill")+
                        #                     ggtitle("Growth effects")+
                        #                 theme(legend.position="bottom")
                            #ggsave("Map_correcteddata.png",dpi=600)
                                
                # Categorizing statistically different estimates (end)
            ## 3.3. LAGS categorizing - Plotting Figure 3 (end) 

            ## 3.4. LAGS FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
                fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                        origin = "iso3c",
                                        destination = "continent")
                fmod_fft$invse <- 1/fmod_fft$StandardError
                fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                f_unfiltpos <- fmod_fft[fmod_fft$signunfilt==1,]
                f_unfiltneg <- fmod_fft[fmod_fft$signunfilt==0,]
                
                felm_5_l <- felm(Estimate ~ factor(Lags)|0|0|continent, data =f_unfiltpos,weights = f_unfiltpos$invse)
                felm_6_l <- felm(Estimate ~ factor(Lags)|0|0|continent, data =f_unfiltneg,weights = f_unfiltneg$invse)
                
                stargazer(felm_5_l,felm_6_l,type="text")
                
            ## 3.4. LAGS FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
            
            ## Figure Lags
                dataset <- data.frame(estimate=double(),cil=double(),cih=double(),filter=factor(),dataset=factor(),sign=factor())
                felm_all_l <- felm(Estimate~signunfilt*factor(Lags), data=fmod_fft_l)
                stargazer(felm_all_l,type="text")
                for(i in 1:3){
                for (j in 1:2){    
                    for (fi in 1:5){
                    if(i==1){
                        datasetname <- "World Bank"
                        if(j==1){
                            model <- felm_1_l
                            sign <- factor("positive")
                        } else {model <- felm_2_l
                        sign <- factor("negative")}
                    }
                    if(i==2){
                        datasetname <- "Barro-Ursua"
                        if(j==1){
                            model <- felm_3_l
                            sign <- factor("positive")
                        } else {model <- felm_4_l
                        sign <- factor("negative")}
                    }
                    if(i==3){
                        datasetname <- "Maddison"
                        if(j==1){
                            model <- felm_5_l
                            sign <- factor("positive")
                        } else {model <- felm_6_l
                        sign <- factor("negative")}
                    }
                    Sigma <- vcov(model)
                    coefT <- "(Intercept)"
                    if(fi==1){
                        sigma = Sigma[1,1]
                        beta.hat <- coef(model)[1]
                        xmat <-1
                        gestimated <- colSums(beta.hat*t(xmat)) 
                        cih <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        cil <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        Filter <- "Unfiltered"
                    }else{
                    if(fi==2){
                        coefT5 <- "factor(Lags)3"
                        Filter <- "3-Lags"
                    }
                    if(fi==3){
                        coefT5 <- "factor(Lags)5"
                        Filter <- "5-Lags"
                    }
                    if(fi==4){
                        coefT5 <- "factor(Lags)10"
                        Filter <- "10-Lags"
                    }
                    if(fi==5){
                        coefT5 <- "factor(Lags)15"
                        Filter <- "15-Lags"
                    }
                        start1 <- which(names(coef(model))==coefT)
                        end1 <- which(names(coef(model))==coefT5)
                        sigma = Sigma[c(start1,end1),c(start1,end1)]
                        beta.hat <- coef(model)[c(start1,end1)]
                        xmat <-cbind(1,1)
                        gestimated <- colSums(beta.hat*t(xmat)) 
                        cih <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        cil <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                        }
                        dataset <- rbind(dataset,data.frame(estimate=gestimated,cil=cil,cih=cih,filter=factor(Filter),dataset=factor(datasetname),sign=factor(sign)))
                }
                }    
                }
                #dataset$filter <- factor(dataset$filter, levels = c("Unfiltered", "5 years", "10 years", "15 years"))
                
                
                cols <- c("#66c2a5","#fc8d62","#8da0cb")

                a <- ggplot(dataset[which(dataset$dataset=='World Bank'),],aes(x=filter,y=estimate*100,fill=dataset,color=dataset,group=interaction(sign,dataset)))+
                geom_line(color="#009E73")+theme_bw()+
                geom_hline(yintercept=0)+
                xlab("Lags")+ylab("Estimated pooled effect of \n 1 Degree Warming on Growth (pp)")+
                #scale_colour_manual(labels = c("a","b","c"),values = c(cols[1],cols[2],cols[3]))+
                geom_ribbon(aes(ymin=cih*100,ymax=cil*100),alpha=0.2,colour=NA,fill="#009E73")+
                ggtitle("World Bank")+
                ylim(-3.4,5.7)

                b <- ggplot(dataset[which(dataset$dataset=='Barro-Ursua'),],aes(x=filter,y=estimate*100,fill=dataset,color=dataset,group=interaction(sign,dataset)))+
                geom_line(color="#E69F00")+theme_bw()+
                geom_hline(yintercept=0)+
                xlab("Lags")+ylab("")+
                #scale_colour_manual(labels = c("a","b","c"),values = c(cols[1],cols[2],cols[3]))+
                geom_ribbon(aes(ymin=cih*100,ymax=cil*100),alpha=0.2,colour=NA,fill="#E69F00")+
                ggtitle("Barro Ursua")+
                ylim(-3.4,5.7)
                
                c <- ggplot(dataset[which(dataset$dataset=='Maddison'),],aes(x=filter,y=estimate*100,fill=dataset,color=dataset,group=interaction(sign,dataset)))+
                geom_line(color="#CC79A7")+theme_bw()+
                geom_hline(yintercept=0)+
                xlab("Lags")+ylab("")+
                #scale_colour_manual(labels = c("a","b","c"),values = c(cols[1],cols[2],cols[3]))+
                geom_ribbon(aes(ymin=cih*100,ymax=cil*100),alpha=0.2,colour=NA,fill="#CC79A7")+
                ggtitle("Maddison")+
                ylim(-3.4,5.7)
                


                ggarrange(a,b,c,ncol=3,nrow=1)
            ## Figure Lags
        # Maddison
    ## 3.5. Other datasets - Table 1, columns 3 to 6; Supp Fig 3 (end)

## 3. Empricial analysis - Figure 3 and Table 1  (end)

