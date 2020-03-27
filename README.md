# woodlice-and-nursery-web-spiders
Data and analysis from three behavioural experiments examining how nursery web spiders and grasshoppers respond to the presence of woodlice.

The folder "Data" contains all the datasets used in the analysis and is loaded in largely during the assembly of data in the code "behaviour_March2019.R". There are two kinds of data: behavioural data and cage data. The former is from behavioural observations, while the latter is the end of summer cage harvesting data and survival counts. This code is sourced in the analysis file "GreenBrownAnalysis.R" because the code in "behaviour_March2019.R" only rearranges and merges the raw data for analysis. 

The code automatically produces plots found in the paper along with other diagnostic plots.

Data Analysis: GreenBrownAnalysis.R

The model 1 measures the effect of woodlice and temperature on the height of grasshoppers and spiders. It includes default priors and is run in STAN.

The model 1.1 adds priors for grasshopper and spider height at different based on data collected by Brandon Barton and presented in his 2009 MS. This is the model included in the MS, but the addition of priors does not change the qualitative conclusions of our analysis.

The model 1.2 only tests for a temperature effect. This is not included in the MS.

The model 2 tests for the relationship between grasshopper survival and predicted spider attack rate using a poisson distribution.

The model 2.1 is the model 2 without random effects.

The model 2.2 tests for the relationship between grasshopper survival and predicted spider attack rate using a binomial distribution. This was identified as the best model and is presented in the MS.

The final three sections of the code create SI graphics, Figure 1 graphics, and include some waste code with additional model examples.

Model Simulations: theoretical_model_analysis.R

The first section of code (1 -> 1.6) loads in the data and calculates the distribution of each animal. Isopods are gamma distrubuted, while the other species are normally distributed. Then (1.5) the probability of encounter is calculated by looking at each distrubution. Finally, this section contains some notes on the energy cost data. 

The second section of the code runs the signal detection theory model. The 'outcome' function runs the analysis and takes inputs of the number of attack attepts (Nin) and the height of the spiders (htmira). This function can then be optimized to find the best spider behaviour.

The third section runs a similar model, but assumes that spiders always attack. The function 'alwaysattack' can be optimized in the same way as the 'outcome' function above. Now the inputs are spider height (HTMIRA), attack success, encounters per day (en.per.day), and presence of woodlice.

The fourth section runs the random movement simulations based on the decision rules set out in 'simfunc.' Parameters are defined above the function. The results of these simulations are saved in the "SimRes" folder.

The fifth section produces the plots for the MS.

1 versus 3 dimensions: compare_1D_3D.R

Compares the 1-D separation and the 3-D separation and produces Figure S1

HOBO Temperature data: HOBOdataanalysis.R

Produces the temperature plots from HOBO logger data in the 2018 cages

