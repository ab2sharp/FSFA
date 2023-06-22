# FSFA
All simulations for FSFA

# Directories

## Parameter Recovery

  The file `frankenclustAuto.R` contains the functions for fitting the model
  
  The file `parrec_funs.R` contains functions used specifically for parameter recovery, including parameter and data generation
  
  The file ` parrec_genparams.R` is used to generate the true parameter sets
  
  The file `genparameters.sh` runs `parrec_genparams.R` for the required specifications
  
  The files `trueparamsX.R` are the results of running the bash file
  
  Finally, `parrec_simfile.R` runs the parameter recovery simulation. 

Data files to be added shortly. 

## Pitch Clustering

This contains code related to unsupervised clustering done on baseball pitch trajectories.

  The file `scrape_clean_statcast.R` contains functions and code from getting and cleaning the data from statcast.

  The file `get_raw_trajectories.R` grabs all pitcher/pitchtype combinations that were seen at least 100 times in a season and collects the associated data.

More to come soon.
