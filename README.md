# EWAS Twin Simulation

## Usage

### Run the script

Run script "Simulation_Main.R". The simulation uses lock file mechanism, so that multiple instances could run at same time as long as they have shared disk.

### Change parameters

Paramters could be changed in "Parameter_space.R" if you wish.

### Results

All results are located in "result" folder. File names are formated as:

*P<sub>K</sub> _istwin_h<sup>2</sup>_E(M)_R<sup>2</sup><sub>M,E</sub>_Var(M)_ρ<sub>ε</sub>*

We used lm column in the results. For ordinary twin design, this column is the power calculation from linear regression. For discordant twin design, it is calculated from linear regression on the pairwise difference. 
