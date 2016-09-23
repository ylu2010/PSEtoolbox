# PSEtoolbox
## Parameter Space Exploration Toolbox

  This is software package for exploring the parameter space of a numerical model using Markov-Chain Monte-Carlo, Latin Hyper-Cube sampling, and others. 
I write these packages in my research on astrophysics. My research focuses on the modeling the formation of galaxies in the Universe. I build complex models for galaxies and use observational data to constrain my models. 
These models are often high-dimensional and uncertain. The best way to make process in improving the models is to use observational data to constrain the model and to study the constrained models in their high-dimensional paramter space. 

## MCMC


### compilation
1). enter the directory PSEtoolbox/mcmc/

2). type "python setup.py mvnormal ./project_mvn"

3). enter the directory ./project_mvn

4). type "make"

### run the program
1). enter the directory run/

2). type "../project_mvn/mcmc_test -f ./MCMC_PARAMETER -r 0" to run

3). type "../project_mvn/mcmc_test -f -f ./MCMC_PARAMETER_RESUME -r 10" to resume the run

### add your own likelihood function (project)
