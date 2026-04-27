# Automated hypoxia statistics from respirometry data
An automated looping script designed to calculated hypoxia statistics from respirometry plates with multiple organisms per plate

<img src="Figures/02-28-23 722 UNIQUE.jpg" width="600">
A typical round of oxygen consumption on a plate respirometer. Each line represents a single organism in a single well on the plate.

## Problem 1: I needed to be able to analyze data from dozens of individuals at once measured simultaneously on a plate respirometer. 
Per individual organism, I needed to be able to calculate multiple statistics that estimate tolerance to decreasing oxygen in the environment. 

These included: 
- Alpha (oxygen supply capacity)
- *P*<sub>crit</sub> (segmented broken stick)
- non-linear *P*<sub>crit</sub> (NLR) 
- regulation index (RI)

For extra information on what these statistics are and represent, check the dropdown box below. 
<details>
<summary>Hypoxia statistics glossary</summary>
 
 - Alpha: Alternatively known as the oxygen supply capacity. This represents the supply of oxygen needed to support a given metabolic rate. 
 
 - *P*<sub>crit</sub>: Perhaps the most widely used statistic for estimating how tolerant an organism is to hypoxia. It is also known as critical oxygen tension. It represents the point at which an organism can no longer actively control its respiratory rate as oxygen is depleted. *P*<sub>crit</sub> and Alpha are on the same scale (units of oxygen partial pressure), where lower values are associated with increased tolerance. In fact, *P*<sub>crit</sub> is the breakpoint at which alpha reaches is maximum value. Oxygen levels below this breakdown becomes a limiting factor during respiration. 
 
 - Non-linear *P*<sub>crit</sub>: A form of calculating *P*<sub>crit</sub> using non-linear regression. Typically, *P*<sub>crit</sub> is calculated using linear equations, where the data ends up looking like a broken stick. Non-linear *P*<sub>crit</sub> identifies *P*<sub>crit</sub> at a predifined part of the slope of a curve instead. 

 - Regulation index (RI): Captures an organism's abilityt to regulate oxygen across all oxygen values from normoxia to anoxia instead of estimating a single breakpoint, like *P*<sub>crit</sub> or Alpha.
   
</details>


I also needed to generate *multiple* graphs showing these calculations on a per individual basis in a way that I could quickly scroll through them. 

## Problem 2: I needed an automated way to standardize the respirometry data to help improve consistency in my calculations. 
Not all individuals consume oxygen at the same rate. The statistics I calculate use the entire curve of oxygen consumption, from normoxia to when the organisms stops actively respiring. So, I need to trim noisy data from the start of runs and find the point at which oxygen consumption stops on a per individual basis. 

### I generated a script that takes care of these issues with one click. 
The script called `Combined Resp Script.R ` is available in this repository, along with sample data called `example_resp_data.csv`. 

Below, let's walk through how to set up data to use this script and what it does. 
