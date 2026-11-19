# Community diversity favors coexistence between bacteria and parasitic bacteriophages

Microbial communities often contain diverse bacteria and their specific viruses 
(bacteriophages). Lytic phages frequently select for host bacteria to evolve 
resistance, in-turn threatening phage extinction. However, in nature 
phage-bacteria interactions do not occur in isolation, but within complex 
communities containing other species. To test how community context can alter 
resistance evolution and phage persistence, we experimentally evolved 
Pseudomonas aeruginosa bacteria and a Pseudomonas-targeting phage with and 
without a community representing a cystic fibrosis disease pulmonary infection. 
Phages selected for rapid evolution of receptor-mediated resistance, regardless 
of community context. However, susceptible bacteria often persisted due to 
resistance-growth trade-offs. A diverse community favored increased frequencies 
of susceptibility, fueling higher phage densities and fewer phage extinctions. 
Our findings show that community complexity can constrain the evolution of 
phage resistance in bacteria, promoting host-parasite coexistence and fostering 
microbial diversity.

# Getting started

R is needed to re-run analysis files included (described below). All exact 
version and package dependencies are described in the renv.lock file (described 
below). In order to re-create the environment where analyses were run, download 
all files into a single directory. Install R 4.5.1, set the working directory 
to the directory containing files, and download renv 1.1.8 (as in the renv.lock 
file). Then run renv::restore()

The `clean_data` directory contains nicely-formatted files with all the data.
These files were generated from the files in `raw_data` by `Clean_files.R`,
and `Clean_sequencing_files.R`. To re-create all figures and analyses in the
manuscript, run `Analyze_data.R` and `Analyze_sequencing_data.R` in whatever
order (the sequencing analysis depends on `Define_colgroups_sequencing.R`).
Plots will be generated in the `plots` directory.