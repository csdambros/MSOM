# Multi Species Occupancy Models
Implement and exemplify the use of multi-species occupancy models (MSOM)

This project contains several builds of multi-species occupancy models and serves the purpose of
1 - Demonstrating the how MSOM works
2 - Showing how standard ecological data can be organized to run the models
3 - Showing how flexible the model is for several types of ecological data
4 - Creating a backbone for other users to improve the model (please create branches!)

# Structure

## Start

The file ExampleData1.csv contains a simple ecological dataframe that can be used as a model. The data is represented by a matrix with sites as rows and species as columns. The numbers filling the matrix represent species INCIDENCES (not abundance!). The number of incidences of a species in a site (number filling the matrix) could represent:

1 - number of pitfall traps in which a species was found in a site   
2 - number of observers that saw a bird in a site during a bird watching event   
3 - number of hours in a day in which a butterfly species was observed in a site   
4 - number of trees in a site where an epiphyte species was observed   
5 - number of nets in which a bat species was observed in a site   
6 - number of sections within a transect in which a plant or a termite species was observed   
7 - many other examples   

Note that in all cases, there is a cap on the maximum number of incidences in a site. For example, if a butterfly species was observed every hour during a 6 hour period, then the minimum observation of the species is zero and the maximum is 6. Therefore, the values filling the matrix MUST lay somewhere in the range from zero to 6. This does not mean that some species must have been observed in observed hours (i.e. the maximum value in the matrix can be lower than 6).

The last column of this data has the maximum number of incidences that could be obtained.

This data structure allows:
1 - Estimating the occupancy (true presence/absence) of a each species in each site   
2 - Estimating the detection (how easy it is to observe) of each species in each site   
3 - Estimating the variability among species in occupancy   
4 - Estimating the variability among species in detection   
5 - Estimating the number of undetected species   
6 - Estimating other community metrics (eg. similarity in species composition) when some species were not detected   
7 - Estimating how species occupancy and detection vary as a function of site covariate (if environmental covariates were measured in each site)   
8 - Estimating how species richness change with scale (from group of sites (region) to individual sites to individual replicates within sites).   

This data structure DOES NOT allows:

1 - Estimating how occupancy and detection of species respond to environmental co-variates measured WITHIN sites (eg. difference in forest cover from one pitfall to the next).


script Demonstration.R creates a simple ec
