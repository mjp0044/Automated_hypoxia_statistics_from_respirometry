# Automated_hypoxia_statistics_from_respirometry
An automated looping script designed to calculated hypoxia statistics from respirometry plates with multiple organisms per plate

## Problem 1: I needed to be able to analyze data from dozens of individuals at once measured simultaneously on a plate respirometer. 
Per individual organism, I needed to be able to calculate multiple statistics that estimate tolerance to decreasing oxygen in the environment. 
These included: Alpha, Pcrit (segmented broken stick), non-linear Pcrit (NLR), and regulation index (RI). 

I also needed to generate multiple graphs showing these calculations on a per individual basis in a way that I could quickly scroll through them. 

## Problem 2: I needed an automated way to standardize the respirometry data to help improve consistency in my calculations. 
Not all individuals consume oxygen at the same rate. 

The statistics I calculate use the entire curve of oxygen consumption, from normoxia to when the organisms stops actively respiring. 

We need to trim noisy data from the start of runs and find the point at which oxygen consumption stops on a per individual basis. 

### I generated a script that takes care of these issues with one click. 
It is available in this repository, along with sample data. 

Below, let's walk through how to set up data to use this script and what it does. 
