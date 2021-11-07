# Random-variate-generation

A class facilitating random variate generation from Gamma, Chi square, Exponential, Pareto, Exponential Gamma, Beta, Fisher, Weibull, Gumbel, Extreme value, Geometric, Logistic, Cauchy, Pareto, Log Normal, Student T, Poisson, and Negative Binomial distributions.

The Gamma distribution is based on the method described by Marsaglia and Tsang:

    // A Simple Method for Generating Gamma Variables
    // ACM Transactions on Mathematical Software Vol. 26, No.3
    // September 2000, Pages 363-372

The Student T distribution is based on the paper:

    // Polar generation of random variates with the t-distribution
    // By Ralph W. Bailey, Mathematics of Computation, Volume 62, Number 206
    // April 1994, Pages 779-781
    
The Poisson distribution was implemented from the method described by Hörmann:

    // The transformed rejection method for generating Poisson random variables
    // Insurance: Mathematics and Economics 12, 1993, Pages 39-45
    
The Exponential distribution is based on table lookup. Briefly, the distribution is divided into n parts
all of equal probability. A variable can be generated with equal probability from each of these regions.
The variable is randomized further with a uniform distribution within each region.
The method can be sped up considerably by only using one random draw.

The remaining probability distributions are combinations of the above described distributions.

Java 17 LTS contains some very nice and effecient versions of the Exponential and Normal distributionm based on the paper:

    // Christopher D. McFarland (2015): A modified ziggurat algorithm for generating exponentially and normally distributed
    // pseudorandom numbers, Journal of Statistical Computation and Simulation, DOI: 10.1080/00949655.2015.1060234
    // To link to this article: http://dx.doi.org/10.1080/00949655.2015.1060234

DISCLAIMER: Although the distributions have been tested there will most likely be many bugs so use at your own risk.
            Some of the implemented distributions are based on the Gamma, Poisson, Student T and Exponential distributions
            and may therefore not be very efficient.
            
Usage:  Just include the ExtendedDistributions class and instantiate an object, i.e. ExtendedDistributions e = new ExtendedDistributions();
        Then you can just use e.nextGamma(), e.nextPoission() etc. See code for more info and for mandatory paramters to include.
