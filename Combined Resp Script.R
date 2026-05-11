#Combined respirometry script
#Matthew Powers 
#powerm3@oregonstate.edu
#Combines RI, Pcrit, and alpha calculations for analysis of respirometry data
#Last edited 06-10-23 
#Built originally on R version 4.2.2 "Innocent and Trusting"
#Tested last on R version 4.3.3 "Angel Food Cake"

#load required packages
library(MESS) #data wrangling and modeling
library(stats) #data wrangling and modeling
library(respirometry) #pcrit, alpha, nlr (MUST USE VERSION 1.3.0, NEWER VERSIONS PRODUCE ERRORS)
library(ggpubr) #Also loads basic ggplot2
library(cowplot) #Pretty ggplots!
library(reshape2) #Data wrangling (melting)
library(dplyr) #Data wrangling
library(tidyverse) #Data wrangling
library(stringr) #data wrangling for splitting column characters
library(MetBrewer) #Pretty colors!
library(runner) #Do functions on sliding windows of data
library(presens) #O2 conversion using presens package


#Data formatting
#To start, the data file needs to be in .csv format
#The first column should be the time in seconds labeled time.sec

#The columns of copepod data should be labeled with individual ID as a single string with the group and replicate number
    #example BOB.1 BOB.2 SD.1 SD.2 BR.3 BR.5 and so on...

#The columns of blank data should be labeled to start with the exact string 'blank' followed by .1 through .n number of blanks
    #example blank.1 blank.2 up to blank.n

#Read in data of your choosing
#IMPORTANT: Data needs to be formatted and saved as a microsoft .csv (comma separated values) file


# O2 unit formatting: our system recorded oxygen units as mg_per_l. 
# Therefore, the data should be converted to kPa for comparison to other studies

############################################################################################

#Read in mass or length data for all your samples to normalize MO2 values
##IMPORTANT: You must read in this frame before running the for loop to calculate Pcrit, RI and other statistics if you want 
#mass specific MO2 values
#The mass or length MUST be in the second column for this to work as desired. The sample column MUST be titled "ID"
  datum = read.csv(file="Length data for R.csv", header=TRUE)


#Data read-in for respirometry data file (should pop up window to select .csv file)
  dat=read.csv(file.choose())
  
  str(dat) #View variable formats to make sure all data is integers (time.sec) or numbers (well data)
  head(dat) #View top of data frame to check contents



#If you recorded every 30 seconds, run this first before anything else
#dat = dat[seq(1, nrow(dat), 10), ] #Only keeps every 10th row which corresponds to every 300 seconds

#### Data trimming ####

#Choose the option that best suits you! Don't run all these back to back, that would rewrite and cause havoc in the data
#Just pick one burn in and trimming option. 
#I recommend dafault to 1 hour burn unless you don't have much data, then use half hour burn instead. 
#I recommend the individual, unique trimming at the end of the data using the for loops further down this script. 

#this package requires the time to be in minutes, so convert as needed:
  time.min = dat$time.sec/60
  dat = cbind(time.min, dat)
  head(dat)

#Below, you can choose to trim the data different ways to test the effect on the statistics
#Just reread the data in again as "dat" each time

##RECOMMENDED
#remove first 1 hours. Change number of rows to fit data file.
  dat=dat[-(1:13),]
  head(dat)
  dat$time.min=dat$time.min - 65   #Subtract from time.min column to artificially change time values to start at 0
  head(dat)
  dat=dat[,-(2)]
  head(dat)

#Other trimming options

# #remove first 30 min. Change number of rows to fit data file.
#   dat=dat[-(1:7),]
#   head(dat)
#   dat$time.min=dat$time.min - 35   #Subtract from time.min column to artificially change time values to start at 0
#   head(dat)
#   dat=dat[,-(2)]
#   head(dat)
# 
# #Test taking the first two hours off instead
# dat=dat[-(1:25),]
# head(dat)
# dat$time.min=dat$time.min - 125   #Subtract from time.min column to artificially change time values to start at 0
# head(dat)
# dat=dat[,-(2)]
# head(dat)
# 
# #Test taking five hours off the start 
  # dat=dat[-(1:61),]
  # head(dat)
  # dat$time.min=dat$time.min - 305   #Subtract from time.min column to artificially change time values to start at 0
  # head(dat)
  # dat=dat[,-(2)]
  # head(dat)
# 
# #Trim everything past 24 hours (1440 minutes) into another data frame called datb
# #use after trimming first hour off above and resetting rownames
  # dat=dat[-(290:516),]
  # tail(dat)



#Reset row numbers after trimming
  rownames(dat) <- NULL
  head(dat)


  
# Prior to calculating respiratory statistics, convert all values from mg/L O2 (what our system recorded) to kPa
# This will help compare to other studies which mostly report Pcrit values in kPa
# Uses the o2_unit_conv function from the 'presens' package
  for(i in 2:ncol(dat)) {
    dat[ , i] <- o2_unit_conv(o2 = dat[ , i], from = "mg_per_l", to = "kPa", salinity = 35,
                               temp = 20, air_pres = 1.013253)
  }
  
  
## Individual (unique) trimming of tail end of data

#Plot data without trimming
# Plot oxygen consumption graph for whole plate (use dat or datb to test trimming)
  dat_long <- melt(dat, id = "time.min")
  head(dat_long)
  levels(dat_long$variable)
  dat_long[c('group', 'replicate')] <- str_split_fixed(dat_long$variable, '[.]', 2)

#Check column formats
  str(dat_long)

#Reclassify group as factor
  dat_long$group = as.factor(dat_long$group)

#Relevel group labels
#dat_long$group <- factor(dat_long$group, levels = c("blank", "MBR", "FBR", "MSD", "FSD"))

#Check order of variables and put in order of peaks appearance
  levels(dat_long$group)



#Color pallete tweaking for ggplot (optional)
  br_pal <- met.brewer("Juarez")
  scales::show_col(br_pal)
  
  my_pal <- c("black", br_pal[c(1,2,3,5)])
  scales::show_col(my_pal)

#Make plot
  p1 <- ggplot(dat_long,            
               aes(x = time.min, y = value, group = variable, color=group))+
    theme_bw()+
    theme(legend.position = "bottom", legend.direction = "horizontal", 
          legend.key.size = unit(1, "lines"))+
    theme(panel.grid.minor = element_blank())+
    scale_x_continuous(n.breaks = 10, limits = c(0, 1600), expand= c(0.01,0),
                       sec.axis = sec_axis(~ . /60, name="Time in hours", breaks=c(0,5,10,15,20,25,30,35,40)))+
    scale_y_continuous(n.breaks=9, limits=c(-1,25))+
    ylab(expression(paste(PO[2], " (kPa ", O[2], ")")))+xlab("Time in minutes")+
    geom_line(linewidth=0.75)+
    scale_color_manual("Groups:", values=my_pal)
  
  p1 #View graph


#Export graph (run three lines as single block of code, change name to whatever you want to call it)
  # jpeg(file="test full.jpg", units="in", width=9, height=6, res=500)
  # p1
  # dev.off()




#Trim based on individual sample minimums


#Remove columns that start with the string 'blank' to leave only the copepod data
#Use dat or datb on first line to test trimming
  dat %>% 
    select(-contains("blank")) -> dat2
#Reset row numbers
  rownames(dat2) <- NULL

#This code will calculate the slope of the DO2 data vs time (i.e., MO2) 
#When MO2 reaches below a certain threshold (i.e., where the slope > -0.001, a.k.a when the data plateaus)
#it will mark this point in the data, add an hour, and then trim the data beyond that added hour for each individual sample

#create data frame with 0 rows and 4 columns to hold plateau point data for later
  plateau.dat <- data.frame(matrix(ncol = 4, nrow = 0)) 
#provide column names for empty data frame
  colnames(plateau.dat) <- c('variable', 'rowindex', 'time.min', 'value')
  
  for(i in 2:ncol(dat2)) {
    
    cutoff.df <- data.frame(a = dat2[,1], b = dat2[,i]) #Assign time and sample columns to data frame
    rowindex <- nrow(cutoff.df) #Get length of rows in dataframe
    
    #Using a sliding window of 45 data points (~3 hours worth), calculate the slope of the O2 values (i.e., MO2) and save as object
    #Use lag of 10 to make the window start 30 data points in just in case the first 50 minutes happen to have a slope too close to 0
    slidingslopes <- runner(x = cutoff.df, k = 45, lag = -30,
                   f = function(x) {
                     model <- lm(b ~ a, data = x)
                     coefficients(model)[]
                   }
    )
  
    slidingslopes2 <- as.data.frame(t(slidingslopes)) #assign slopes to data frame
    slidingslopes2$a <- format(slidingslopes2$a, scientific = F) #remove scientific notation
    suppressWarnings(slidingslopes2$a <- round(as.numeric(slidingslopes2$a), digits = 5)) #round vlaues to 5 decimal points
    
  #Apply threshold for slope. Default to slope less than 0.0001 in magnitude
  #Save line number of where slope becomes less than 0.001
  #If no slope is calculated to be close to 0 (i.e., flat), it defaults to NA.
    slopecutoff <- which(slidingslopes2$a > -0.0001) 
  
  #if statement to take care of instance where we never hit a slope of 0
  #Tests to see if slopecutoff is NA. If it is, it defaults slopecutoff to just be the last row of the dataset to use all the data
    slopecutoff <- ifelse(is.na(slopecutoff[1]), rowindex-12, slopecutoff) 
    
    cutoff.start <- slopecutoff[1] #Get row where plateau starts
    cutoff.start.plushour <- cutoff.start+12 #add an hour+1 past plateau start
    cutoff.start.plushour.andone <- cutoff.start.plushour+1 #Store next data point row number
  
  #As long as cutoff is less than last row of dataset, apply NA's to values after the cutoff
    if (cutoff.start.plushour < rowindex) {
      dat2[cutoff.start.plushour.andone:rowindex, i] <- NA
      } else {
      print(paste("All data used with", colnames(dat2[i])))
      }
    
    #Assign plateau point information to plotting data frame for making figure below
    plateau.dat[nrow(plateau.dat) + 1,] <- list(colnames(dat2[i]),cutoff.start.plushour, dat2[cutoff.start.plushour,1], dat2[cutoff.start.plushour,i])
    
  } 
  plateau.dat[c('group', 'replicate')] <- str_split_fixed(plateau.dat$variable, '[.]', 2) #split ID column into groups
#Need to sort frame based on plateau point order
  plateau.dat <- arrange(plateau.dat, time.min)
#create another variable with sequence of y-values for plotting so that the labels decrease high to low on the plot
  yvalues <- seq(from =17, to=9, length.out=nrow(plateau.dat))
  plateau.dat <- cbind(plateau.dat, yvalues)
  plateau.dat %>% replace_na(list(time.min = max(dat2$time.min), value = 0)) -> plateau.dat

#Relevel group labels, use line above again if you want to confirm it worked
#plateau.dat$group <- factor(plateau.dat$group, levels = c("blank", "MBOB", "FBOB", "MSD", "FSD"))
#Check order of variables and put in order of peaks appearance
#levels(plateau.dat$group)

#Make plot of individually trimmed data
# Plot oxygen consumption graph for whole plate (use dat or datb to test trimming)
  dat_long_unique <- melt(dat2, id = "time.min")
  head(dat_long_unique)
  dat_long_unique[c('group', 'replicate')] <- str_split_fixed(dat_long_unique$variable, '[.]', 2)

#Check column formats
  str(dat_long_unique)

#Append blank wells back onto the end of the cutoff data frame
  dat_long_unique <- rbind(dat_long_unique, dat_long[dat_long$group == "blank",])

#Reclassify group as factor
  dat_long_unique$group = as.factor(dat_long_unique$group)

#Check order of variables and put in order of peaks appearance
  levels(dat_long_unique$group)

#Relevel group labels, use line above again if you want to confirm it worked
#dat_long_unique$group <- factor(dat_long_unique$group, levels = c("blank", "MBR", "FBR", "MSD", "FSD"))

#Make plot2
  p2 <- ggplot(dat_long_unique,            
               aes(x = time.min, y = value, group = variable, color=group))+
    theme_bw()+
    theme(legend.position = "bottom", legend.direction = "horizontal", 
          legend.key.size = unit(1, "lines"))+
    theme(panel.grid.minor = element_blank())+
    scale_x_continuous(n.breaks = 10, limits = c(0, 1600), expand= c(0.01,0),
                       sec.axis = sec_axis(~ . /60, name="Time in hours", breaks=c(0,5,10,15,20,25,30,35,40)))+
    scale_y_continuous(n.breaks=9, limits=c(-1,25))+
    ylab(expression(paste(PO[2], " (kPa ", O[2], ")")))+xlab("Time in minutes")+
    geom_line(linewidth=0.75)+
    scale_color_manual("Groups:", values=my_pal)+
    geom_segment(data = plateau.dat, aes(time.min, value, xend = time.min, yend = yvalues))+
    geom_text(data = plateau.dat,
              aes(time.min-5, yvalues, label = paste(variable,":",time.min,"min.")), check_overlap = FALSE, size=2.5, hjust=-0.1)
  
  p2 #View graph

#Export graph (run three lines as single block of code, change name to whatever you want to call it)
  # jpeg(file="test UNIQUE.jpg", units="in", width=9, height=6, res=500)
  # p2
  # dev.off()



#### P-crit, alpha, RI, and plotting loop ####

###IMPORTANT: Prior to running for loop, need to create calcs dataframe each time if you want to save data
#To run as block, highlight code from "calcs <- data.frame(matrix(ncol = 7, nrow = 0))" to "dev.off()"

#Create empty data frame to record pcrit and alpha calculations
#create data frame with 0 rows and 2 columns
  calcs <- data.frame(matrix(ncol = 5, nrow = 0)) 

#provide column names
  colnames(calcs) <- c('ID', 'Alpha', 'Breakpoint', 'NLR', 'RI')


          #Or use other data frame naming conventions to test trimming
          # calcs2 <- data.frame(matrix(ncol = 5, nrow = 0)) #2 hour trim
          # calcs5 <- data.frame(matrix(ncol = 5, nrow = 0))#5 hour trim
          # calcs24 <- data.frame(matrix(ncol = 5, nrow = 0)) #cap at 24 hours trim
          # calcsunique <- data.frame(matrix(ncol = 5, nrow = 0)) #cap by sample trim          
          
          #Use above with:
          # colnames(calcs2) <- c('ID', 'Alpha', 'Breakpoint', 'NLR', 'RI')
          # colnames(calcs5) <- c('ID', 'Alpha', 'Breakpoint', 'NLR', 'RI')
          # colnames(calcs24) <- c('ID', 'Alpha', 'Breakpoint', 'NLR', 'RI')
          # colnames(calcsunique) <- c('ID', 'Alpha', 'Breakpoint', 'NLR', 'RI')
          
          

#define axes labels with superscript to add to plots
  Xlab_expression <- expression(paste(PO[2], " (kPa ", O[2], ")"))
  Ylab_expression <- expression(paste("Metabolic rate ", MO[2] , " (uMol ", O[2], " ", hr^-1, ")"))
  


#### Calculate pcrit and alpha with respirometry package ####

#To have all plots export to pdf, run pdf line at the start and then dev.off() at the end after for loop
pdf(file='test plots.pdf')
#Run for loop over data frame with oxygen data
  for(i in 2:ncol(dat2)) {

    
#make bins that account for changes of rates at different PO2.
  bins=make_bins(
    o2= dat2[ , i],
    duration=dat2$time.min,
    min_o2_width = 1/100,
    max_o2_width = 1/20,
    n_thresholds = 10
  )

#Calculate MO2 and mean 02 (and other statistics)
#MO2 is expressed as unit oxygen in umol (micromoles) consumed per hour but can convert below
  mo=calc_MO2(
    duration=dat2$time.min,
    o2=dat2[ , i],
    o2_unit = "kPa",
    bin_width = bins,
    vol = 0.00008,
    temp = 20,
    sal = 35,
    atm_pres = 1013.25,
    good_data = TRUE
  )

#Remove lines with NA's so the calculation below does not freeze
  #mo = na.omit(mo)
  mo=mo[!is.na(mo$O2_MEAN) & !is.na(mo$MO2),]
  
  # mo=droplevels(mo)
  # rownames(mo) <- NULL

#If you want to express MO2 in kPa O2 per min instead. Default is umol O2 per hour. 
  #mo$MO2 <- conv_resp_unit(value = mo$MO2, from = 'umol / hr', to = 'kPa / min', temp = 20)

#Calculate mass or length specific MO2 values for a given sample
    #Get sample ID of current loop iteration
    sampleID <- colnames(dat2[i])
    #Get row from mass or length from data frame
    sampleIDrow <- subset(datum, ID == sampleID)
    #Isolate just the length or mass value
    sampleIDvalue <- sampleIDrow[1,2]
    #Normalize all MO2 values by mass or length value and save as new column in 'mo' dataframe
    mo$MO2.norm <- mo$MO2/sampleIDvalue

#Calculate pcrit
  calc <- calc_pcrit(
    po2=mo$O2_MEAN,
    mo2=mo$MO2.norm, 
    method = "All",
    avg_top_n = 1, 
    level = 0.95,
    iqr = 1.5, 
    NLR_m = 0.065,
    MR = NULL,
    mo2_threshold = Inf,
    return_models = FALSE)

#Calculate alpha
  #Units will be in kPa O2 / mm / hr / kPa
  calc.a <- calc_alpha(po2=mo$O2_MEAN, 
             mo2=mo$MO2.norm, 
             avg_top_n = 1, 
             MR = NULL, 
             mo2_threshold = Inf
            )

#plot pcrit and alpha
  plot_pcrit(
    po2=mo$O2_MEAN,
    mo2=mo$MO2.norm,
    method = "All",
    avg_top_n = 1,
    level = 0.95,
    iqr = 1.5,
    NLR_m = 0.065,
    MR = NULL,
    mo2_threshold = Inf,
    ylim = c(min(mo$MO2.norm-0.0001), max(mo$MO2.norm+0.0001)),
    ylab=Ylab_expression,
    xlab=Xlab_expression,
  )
  title(main = colnames(dat2[i]), adj=0)
#text(x = max(mo$O2_MEAN-0.2), y = min(mo$MO2), colnames(dat2[i]), col="darkgreen", cex=1.2)

#RI calculation and plotting
  test.frame = as.data.frame(cbind(mo$O2_MEAN, mo$MO2.norm))
  test.frame <- rename(test.frame, DO=V1)
  test.frame <- rename(test.frame, MO2=V2)

#Get average oxygen consumption at early and late stages of trial, if desired
#From normoxia to 3.5 mg/L O2
#   upper.test.frame=subset(test.frame, DO >= 10 & MO2 != "NA")
#   upper.rate <- mean(upper.test.frame$MO2)
# #From 3.5 mg/L O2 to 0.1 mg/L O2
#   lower.test.frame=subset(test.frame, DO <= 10 & DO >= 0 & MO2 != "NA")
#   lower.rate <- mean(lower.test.frame$MO2)

  #For each of the below calculations, the spline, regline and conline, use the same span

  #Calculate area under curve for sample data using spline
      #Runs your linear model over a certain span
      lo <- loess(test.frame$MO2~test.frame$DO, span=0.50)
      #Gathers your sequence of dissolved oxygen values using min and max values in data file
      xl <- seq(min(test.frame$DO),max(test.frame$DO), (max(test.frame$DO) - min(test.frame$DO))/1000)
      #Make a new data frame using only complete cases with DO2 and MO2 values
      na <-test.frame[complete.cases(test.frame[, 2]), ]
      #Generate predicted MO2 values using model from above and DO2 values
      MO2<- predict(lo,na$DO)
      #Bind model predicted values and DO2 values together
      resp<-cbind(na$DO, MO2)
      resp <- as.data.frame(resp)
      colnames(resp) <- c("DO", "MO2")
      resp <- aggregate(MO2 ~ DO, data = resp, FUN = mean)
      resp <- resp[order(resp$DO), ]
      #Calculate area under the curve bounded by X axis using basic splice
      pattern <- auc(resp$DO, resp$MO2, type = "spline")
      pattern
      
      #Assign predicted MO2 values and DO2 values to data frame for plotting. 
      resp.plot <- resp  # already has DO and MO2 columns
  
  #Calculate auc for perfect regulatory line
      lo <- loess(test.frame$MO2~test.frame$DO, span=0.50)
      xl <- seq(min(test.frame$DO),max(test.frame$DO), (max(test.frame$DO) - min(test.frame$DO))/1000)
      na <-test.frame[complete.cases(test.frame[, 2]), ]
      n<-length(na$MO2)
      #Get max predicted MO2
      max_MO2<- max(predict(lo,xl), na.rm=T)
      #Get max MO2 at highest DO2 value
      #max_MO2 <- resp[1,2]
      #Repeat max MO2 value as many times over as the number of DO2 values you have
      max_MO2_rep<-rep(max_MO2, n)
      #Bind maxMO2 values to DO2 values and deduplicate
      max_MO2_rep <- as.data.frame(cbind(na$DO, max_MO2))
      colnames(max_MO2_rep) <- c("DO", "MO2")
      max_MO2_rep <- aggregate(MO2 ~ DO, data = max_MO2_rep, FUN = mean)
      max_MO2_rep <- max_MO2_rep[order(max_MO2_rep$DO), ]
      #Calculate auc for regulatory line
      max_auc <- auc(max_MO2_rep$DO, max_MO2_rep$MO2, type = "spline")
      
      #Add regulatory line values to plotting data frame
      resp.plot$regline <- max_MO2_rep$MO2
  
  
  #Calculate auc for perfect conformity line
      r = subset(test.frame, test.frame$MO2>0)
      lo <- loess(test.frame$MO2~test.frame$DO, span=0.50)
      xl <- seq(min(test.frame$DO),max(test.frame$DO), (max(test.frame$DO) - min(test.frame$DO))/1000)
      #Get max M02 from predicted data
      #max_MO2<- max(predict(lo,xl), na.rm=T)
      #Get max MO2 at highest DO2 value
      max_MO2 <- resp[1,2]
      #Get starting DO from resp frame, should match starting DO in test.frame (the original raw data)
      min.DO<-resp[1, 1]
      #Filter for complete cases
      na <-resp[complete.cases(resp[, 2]), ]
      #Get the slope between max_MO2 and startin DO2
      slope<-max_MO2/min.DO
      #Multiply slope down list of DO2 values to get MO2 y coordinates to plot
      resp_min<-slope*na[ ,1]
      #Bind y coordinates with DO2 values
      con_MO2<-cbind(na[ ,1], resp_min)
      con_MO2 <- as.data.frame(con_MO2)
      colnames(con_MO2) <- c("DO", "MO2")
      con_MO2 <- aggregate(MO2 ~ DO, data = con_MO2, FUN = mean)
      con_MO2 <- con_MO2[order(con_MO2$DO), ]
      #Calculate auc for conformity line
      con_auc <- auc(con_MO2$DO, con_MO2$MO2, type = "spline")
      
      #Add conformitory values to plotting data frame
      resp.plot$conline <- con_MO2$MO2
  
  
#Cobmine filename and auc's of data line and two perfect lines
  RIplot<-cbind(pattern, max_auc, con_auc)
  RIplot<-as.data.frame(RIplot)
  RIplot$pattern = as.numeric(RIplot$pattern)
  RIplot$max_auc = as.numeric(RIplot$max_auc)
  RIplot$con_auc = as.numeric(RIplot$con_auc)
  
#For RI>0 (RI_pos)
  RIplot$RI_pos<-(RIplot$pattern-RIplot$con_auc)/(RIplot$max_auc-RIplot$con_auc)
  
#For RI<0 (RI_neg)
  RIplot$RI_neg<-(-1)*(RIplot$con_auc-RIplot$pattern)/(RIplot$con_auc)
  
  RIplot$RI_final<-RIplot$RI_pos
  
#When RI$RI_pos is negative, use r$RI_neg
  RIplot<-within(RIplot, RIplot$RI_final[RIplot$RI_final<0] <- (RIplot$RI_neg[RIplot$RI_final<0]))
  
  RIplot$RIplot$RI_final <- round(RIplot$RIplot$RI_final, 4)
  
  
#Print name of column (i.e., name of well on microplate)
  print(colnames(dat2[i]))
#Print alpha and pcrit calculations
  print(calc[2:5])
  cat("Alpha =",calc.a[[1]])
#Print RI for that individual
  cat(" RI =",RIplot$RIplot$RI_final)
  #print("")
  
#Generate plot with all three lines
  print(
    ggplot(resp.plot, aes(x=DO)) +
    geom_line(aes(y=MO2), color="black")+
    geom_line(aes(y=regline), lty = "twodash")+
    geom_line(aes(y=conline), lty = "twodash")+
    geom_ribbon(data=subset(resp.plot, 0 <= DO & DO <= 100), 
                aes(ymin=conline,ymax=MO2), fill="red", alpha=0.5) +
    geom_ribbon(data=subset(resp.plot, 0 <= DO & DO <= 100), 
                aes(ymin=MO2,ymax=regline), fill="blue", alpha=0.5) +
    geom_point(aes(y=MO2))+
    theme_cowplot()+
    xlab(Xlab_expression)+ylab(Ylab_expression)+
      annotate("text", 
               x = max(resp.plot$DO) - 1, 
               y = min(resp.plot$MO2), 
               label = paste("RI =", RIplot$RIplot$RI_final),
               hjust = 1)+
      ggtitle(colnames(dat2[i]))
    )
  
  plot.dat=cbind(dat2$time.min, dat2[i])
  
  print(
    ggplot(plot.dat, aes(x = plot.dat[,1], y=plot.dat[,2]))+
      theme_bw()+
      theme(panel.grid.minor = element_blank())+
      scale_x_continuous(n.breaks = 10, limits = c(0, 1600), expand= c(0.01,0))+
      scale_y_continuous(n.breaks=8, limits=c(-0.5,23))+
      ylab(expression(paste(PO[2], " (kPa ", O[2], ")")))+xlab("Time in minutes")+
      geom_line(size=1)+
      ggtitle(colnames(dat2[i]))
  )
  
  #Add pcrit values to data frame for export later
  #Calcs is saved as a data frame so you can view the summary data for each well in tabular format
  #Change name of data frame (calcs, calcs2, etc) to add data to whichever frame you want matching names from above
  
  #IMPORTANT FOR ALPHA. 
  #To add the alpha value (not alpha pcrit) to the final frame calcs, you just need to replace 
  #the second item in the list (calc[[1]]) with calc.a[[1]]. This pulls the actual alpha value from the calc.a frame instead.
  
    calcs[nrow(calcs) + 1,] <- list(colnames(dat2[i]),calc.a[[1]], calc[[2]], calc[[4]], RIplot$RIplot$RI_final)

  } 
  dev.off() #End for loop



#export calcs data frame as excel .csv to use with other programs if desired
#Changed cencoring parameters and saved other statistics data frames for comparison
  write.csv(calcs,"Plate_Resp_Statistics_test.csv", row.names = FALSE)

print("Goodbye")
