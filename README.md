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

  The file `scrape_clean_statcast.R` contains functions and code for getting and cleaning the data from statcast.

  The file `get_raw_trajectories.R` grabs all pitcher/pitchtype combinations that were seen at least 100 times in a season and collects the associated data.

More to come soon.


## Comparison: Scenario B+

  The folder `R` contains two files, `comps_funs.R` and `frankenclustAuto.R`, which contain the functions used for running the simulation.
  The file `triangledata_ARI_sim.R` runs the simulation, and produces the file `sim_tridat_results.RData`
  There will be an additional file which runs the simulation for the functional k-means approach as well. This file produces `sim_tridata_kmeans.RData`
  The file `plot_triangledata_sim.R` produces the histogram seen in the paper
  
