# Multi-mode Double-Beamforming Technique


## Description
This code aims to extract phase velocities in different modes in Rayleigh wave via double-beamforming technique. 
1. To reference a webpage, you can use the following format:
   ```markdown
   [Title of the article](https://doi.org/10.1093/gji/ggae283)

## Workflow
1. Create source-receiver pair: python 01_create_src_rec_pair/get_double_beamforming_source_and_receiver_stations.py
2. Execute double-beamforming with selected pair and coarse/fine increments for slowness grid search: python 02_double_beamforming/double_beamforming.py M02 M26 3 0.05 0.005
3. Plot the output: sh 03_plot_DB_results/db_2d.sh M02_M26 3

# Install the required dependencies (if applicable)
1. pip install -r requirement
2. The plotting code requires the Global Mapping Tools (GMT) version above 5.