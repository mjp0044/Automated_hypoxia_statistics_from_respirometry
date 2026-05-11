# Automated hypoxia statistics from respirometry data
An automated looping script designed to calculated hypoxia statistics from respirometry chambers with multiple organisms at once. 

The point of making this script was to avoid having to calculate statistics for one organism at a time. I wanted a fast way to take data from multiple chambers of a respirometer, and analyzing it all at once. For context, you can see the type of respirometer plate I use to collect data, along with a zoomed in look at the aquatic chamber where I put small marine zooplankton. 

<table>
  <tr>
    <td><img src="Figures/Plate.jpg" width="300"></td>
    <td><img src="Figures/Well.jpg" width="200"></td>
  </tr>
</table>



## Problem 1: I needed to be able to analyze data from dozens of individuals at once measured simultaneously on a plate respirometer. 
Per individual organism, I needed to be able to calculate multiple statistics that estimate tolerance to decreasing oxygen in the environment. 

These included: 
- Alpha (a.k.a. oxygen supply capacity)
- *P*<sub>crit</sub> 
- non-linear *P*<sub>crit</sub> 
- regulation index (RI)

For extra information on what these statistics are and represent, check the dropdown box below. 
<details>
<summary>Hypoxia statistics glossary</summary>
 
 - **Alpha**: Alternatively known as the oxygen supply capacity. This represents the supply of oxygen needed to support a given metabolic rate. 
 
 - ***P*<sub>crit</sub>**: Perhaps the most widely used statistic for estimating how tolerant an organism is to hypoxia. It is also known as critical oxygen tension. It represents the point at which an organism can no longer actively control its respiratory rate as oxygen is depleted. *P*<sub>crit</sub> and Alpha are on the same scale (units of oxygen partial pressure), where lower values are associated with increased tolerance. In fact, *P*<sub>crit</sub> is the breakpoint at which alpha reaches is maximum value. Oxygen levels below this breakdown becomes a limiting factor during respiration. 
 
 - **Non-linear *P*<sub>crit</sub>**: A form of calculating *P*<sub>crit</sub> using non-linear regression. Typically, *P*<sub>crit</sub> is calculated using linear equations, where the data ends up looking like a broken stick. Non-linear *P*<sub>crit</sub> identifies *P*<sub>crit</sub> at a predifined part of the slope of a curve instead. 

 - **Regulation index (RI)**: Captures an organism's abilityt to regulate oxygen across all oxygen values from normoxia to anoxia instead of estimating a single breakpoint, like *P*<sub>crit</sub> or Alpha.
   
</details>


I also needed to generate *multiple* graphs showing these calculations on a per individual basis in a way that I could quickly scroll through them. 

## Problem 2: I needed an automated way to standardize the respirometry data to help improve consistency in my calculations. 
Not all individuals consume oxygen at the same rate. The statistics I calculate use the entire curve of oxygen consumption, from normoxia to when the organisms stops actively respiring. So, I need to trim noisy data from the start of runs and find the point at which oxygen consumption stops on a per individual basis. 

### I generated a script that takes care of these issues with one click. 
The script called `Combined Resp Script.R ` is available in this repository, along with sample data called `example_resp_data.csv`. Once the user has updated the names of output files in the script, it will produce graphs, pdfs, and results tables in one go. 

Below, let's walk through how to set up data to use this script and what it does. 

First, we'll go through the required packages:

```r
#load required packages
library(MESS) #data wrangling and modeling
library(stats) #data wrangling and modeling
library(respirometry) #pcrit, alpha, nlr (MUST USE VERSION 1.3.0, NEWER VERSIONS PRODUCE ERRORS)
library(ggpubr) #Also loads basic ggplot2
library(cowplot) #Pretty ggplots!
library(reshape2) #Data wrangling
library(dplyr) #Data wrangling
library(tidyverse) #Data wrangling
library(stringr) #data wrangling for splitting column characters
library(MetBrewer) #Pretty colors!
library(runner) #Do functions on sliding windows of data
library(presens) #O2 conversion using presens package
```

Some of the packages that will be doing heavy lifting in our script are `respirometry`, which we need to calculate Alpha and *P*<sub>crit</sub>, `runner`, which we need to estimate the point at which organisms stop respiring, and `presens` which has automated functions for converting oxygen units. The rest are preffered packages I use for data wrangling and plotting. 

## Data formatting

### Respirometry data
The script is set up to handle respirometry data in a simple, but specific format. The data must be: 
- Organized by column, where each column represents one individual
- The column names must be in the format `group.1`, where group represents treatments (if present) and the number represents the replicate for that treatment
- There must be a `.` between the group and replicate values. We will use this separater to parse this info later
- Blanks must be denoted using the same pattern, in the format `blank.1`
- The first column is named `time.sec` and contains the respirometry run time in seconds (this can be modified if desired, as we wlll go over below)

Here is a snapshot of the first 10 columns of the example data in its proper format in the dropdown below:


 ```r
 > head(dat[,c(1:10)])
  time.sec MSH.1 FSH.1 FSH.2 blank.1 blank.2 blank.3 FSH.3 MSH.3 FSH.4
1        0 6.789 6.408 6.328   7.636   7.696   7.639 6.617 6.998 6.555
2      300 6.745 6.293 6.142   7.606   7.651   7.639 6.544 6.984 6.484
3      600 6.686 6.251 6.042   7.576   7.666   7.609 6.458 6.881 6.370
4      900 6.657 6.194 5.971   7.606   7.666   7.594 6.385 6.881 6.327
5     1200 6.555 6.109 5.844   7.576   7.636   7.624 6.270 6.852 6.256
6     1500 6.570 6.052 5.830   7.606   7.651   7.609 6.227 6.779 6.228
 ```

### Length or mass data for standardizing metabolic rate

The script is set up to automatically standardize metabolic rate data by any desired values such as the body length or mass of an individual. 

This data needs to be in a simple two column format where the first column is named `ID` and the second column can have whatever name you choose. 

The `ID` column should contain sample names in the format matching the column name format in the main respirometry data file. 

For simplicity, the script is set up to automatically query this list of data to find matching values to the respirometry data. This means that you can read in the data from one respirometry run at a time without having to reload the length or mass data frame. 

```r
> head(datum)
      ID Totlen
1 MBOB.1  1.494
2 MBOB.3  1.548
3 FBOB.3  1.640
4 FSD.21  1.736
5 MSD.22  1.667
6 FSD.22  1.463
```

## Reading in data

Reading the data into R is the only manual part of the script, but this could be automated if desired. 

The first data read in is the length or mass data. Here, we used body length. 

```r
#The mass or length MUST be in the second column for this to work as desired. The sample column MUST be titled "ID"
  datum = read.csv(file="Length data for R.csv", header=TRUE)
```

Next, we select the respiromtry file of our choosing using a blank `file.choose()` call so that you can select the desired file in your device's file browser. 
```r
#Data read-in for respirometry data file (should pop up window to select .csv file)
  dat=read.csv(file.choose())
```

NOTE: This could be modified to loop over many files in a directory, if desired. For my purposes, I needed to work with one respirometer plate run at a time to trouble shoot our respirometer or inspect individual plate data, so it made more sense to operate at the single file level. 

The respirometry system I used records the time in seconds during the run. However, the R package expects the time to be in minutes. 

So, we convert by:

```r
time.min = dat$time.sec/60
  dat = cbind(time.min, dat)
> head(dat[,c(1:10)])
   time.min time.sec MSH.1 FSH.1 FSH.2 blank.1 blank.2 blank.3 FSH.3 MSH.3
14       65     3900 6.308 5.618 5.339   7.606   7.636   7.624 5.843 6.547
15       70     4200 6.280 5.576 5.270   7.606   7.666   7.639 5.829 6.490
16       75     4500 6.265 5.521 5.298   7.651   7.681   7.670 5.815 6.446
17       80     4800 6.193 5.507 5.270   7.666   7.696   7.716 5.815 6.475
18       85     5100 6.193 5.507 5.229   7.726   7.726   7.731 5.801 6.461
19       90     5400 6.193 5.507 5.215   7.741   7.801   7.792 5.702 6.475
```

Next, the script trims the first hour of the run from the data to remove the initial temperature equilibration period and noise at the start of the run. 

```r
#remove first 1 hours. Change number of rows to fit data file.
  dat=dat[-(1:13),]
  dat$time.min=dat$time.min - 65   #Subtract from time.min column to artificially change time values to start at 0
  dat=dat[,-(2)]
 > head(dat[,c(1:10)])
   time.min MSH.1 FSH.1 FSH.2 blank.1 blank.2 blank.3 FSH.3 MSH.3 FSH.4
27       65 6.050 5.300 4.817   7.801   7.876   7.822 5.590 6.360 5.378
28       70 6.050 5.313 4.749   7.831   7.891   7.868 5.576 6.346 5.378
29       75 6.035 5.217 4.776   7.846   7.891   7.822 5.534 6.274 5.310
30       80 6.007 5.204 4.722   7.861   7.922   7.884 5.534 6.274 5.269
31       85 6.007 5.176 4.627   7.891   7.937   7.914 5.548 6.246 5.228
32       90 5.993 5.163 4.586   7.891   7.937   7.914 5.492 6.117 5.200
```

Next, we reset the row names to start the numbering back at row 1 after removing the first hour. 

```r
rownames(dat) <- NULL
> head(dat[,c(1:10)])
  time.min MSH.1 FSH.1 FSH.2 blank.1 blank.2 blank.3 FSH.3 MSH.3 FSH.4
1       65 6.050 5.300 4.817   7.801   7.876   7.822 5.590 6.360 5.378
2       70 6.050 5.313 4.749   7.831   7.891   7.868 5.576 6.346 5.378
3       75 6.035 5.217 4.776   7.846   7.891   7.822 5.534 6.274 5.310
4       80 6.007 5.204 4.722   7.861   7.922   7.884 5.534 6.274 5.269
5       85 6.007 5.176 4.627   7.891   7.937   7.914 5.548 6.246 5.228
6       90 5.993 5.163 4.586   7.891   7.937   7.914 5.492 6.117 5.200
```

Our system records oxygen values in mg oxygen per liter. However, the `respirometry` package calculations perform better when the units are converted to kPa. Plus most publications report *P*<sub>crit</sub> in kPa. 

To make this conversion, we take advantage of the `presens` package's built-in conversion function.

We will loop over the columns in the data frame and convert each to kPa in one go, making sure to set the salinity, temperature, and air pressure that match our conditions.

```r
for(i in 2:ncol(dat)) {
    dat[ , i] <- o2_unit_conv(o2 = dat[ , i], from = "mg_per_l", to = "kPa", salinity = 35,
                               temp = 20, air_pres = 1.013253)
  }
```

Before the script calculates the respirometry stats, I also wanted to generate plots that let us visualize the run in full. 

To do that, we will flip the data lengthwise using the `melt` function from `reshape2` where we organize the data by time. Below, you can see this creates a data frame with three columns, the first with the time values, the second with the sample ID and the last with the oxygen value in kPa. Using the `levels` function, we can see all the samples are now ordered lengthwise in the variable column. 

```r
# Plot oxygen consumption graph for whole plate (use dat or datb to test trimming)
  dat_long <- melt(dat, id = "time.min")
  >   head(dat_long)
  time.min variable    value
1       65    MSH.1 17.24216
2       70    MSH.1 17.24216
3       75    MSH.1 17.19942
4       80    MSH.1 17.11962
5       85    MSH.1 17.11962
6       90    MSH.1 17.07972

>   levels(dat_long$variable)
 [1] "MSH.1"   "FSH.1"   "FSH.2"   "blank.1" "blank.2" "blank.3" "FSH.3"   "MSH.3"   "FSH.4"   "MSH.5"   "FSH.5"   "FSH.6"   "MBR.39"  "MBR.40"  "MBR.53" 
[16] "FBR.53"  "MBR.54"  "FBR.54"  "MBR.56"  "FBR.47"  "MBR.47"  "FBR.49"  "MBR.49"  "FBR.50"
```

Next, we will split the variable column into two new columns that help us organize the data and plot by treatment group. I use the `str_split_fixed` function from `stringr` and reclassify the group variable as a factor in case I want to re-order the groups later. 

```r
dat_long[c('group', 'replicate')] <- str_split_fixed(dat_long$variable, '[.]', 2)

#Reclassify group as factor
  dat_long$group = as.factor(dat_long$group)

#Relevel group labels (OPTIONAL)
#dat_long$group <- factor(dat_long$group, levels = c("blank", "MBR", "FBR", "MSD", "FSD"))

#Check order of variables and put in order of peaks appearance
>   levels(dat_long$group)
[1] "blank" "FBR"   "FSH"   "MBR"   "MSH"
```

Now we can start plotting! 

Making beatiful plots is an art form (no seriously), so I like to integrate art into my color palettes. 

I love the `met.brewer` package by Blake Mills. It is a collection of palettes inspired by pieces at the metropolitan museum of art. 

I recommend checking the `met.brewer` GitHub page [here](https://github.com/BlakeRMills/MetBrewer): 

I chose to pull colors from the "Juarez" palette. Using the `show_col` function from `scales` we can preview it: 

```r
br_pal <- met.brewer("Juarez")
  scales::show_col(br_pal)
```
<img src="Figures/br_pal.png" width="300">

We will pull the first, second, third, and fifth colors since I am plotting four groups at a time in the example data. More can be take from other palettes if desired using the same method. Then, we combine with the color black for the blank chambers. 

```r
my_pal <- c("black", br_pal[c(1,2,3,5)])
  scales::show_col(my_pal)
```
<img src="Figures/my_pal.png" width="300">

Now, lets use our lengthwise data to generate a line graph showing the decline in oxygen over the run. 

We will color code the lines by the treatment groups and group the values by the sample ID (which still is just named 'variable'). 

We will also display a second x-axis on the top of the graph in hours to accomodate showing very long runs using the `sec.axis` option in `scale_x_continuous`. 
```r
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
```
<img src="Figures/02-28-23 722 full.jpg" width="600">
Now, we can see the blanks at the top and all the oxygen values falling in chambers containing organisms. 

Now, let's address problem 2. We need to standardize these oxygen curves. 

Right now, they plateau at different times and this could bias calculations on a per individual basis. 

So, I used the `runner` package to calculate slopes in a sliding window along the curve and truncate the data at one hour post plateau. If the data ends before the hour cut-off, the code automatically tells you that all the data was used. 

This code block runs a for-loop over the columns and replaces oxygen values past the one-hour mark with NAs. 

It also marks the cutoff point at the hour mark so we can label this point on the graph we will make. It generates two new data frames for plotting and for the rest of the script called "plateau.dat" and "dat_long_unique". 

It's a lot of code, so I won't break it down line by line. But it is fully annotated below:  

```r
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

```

The end product when we generate another line graph looks like this: 

```r
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
```

<img src="Figures/02-28-23 722 UNIQUE.jpg" width="600">

