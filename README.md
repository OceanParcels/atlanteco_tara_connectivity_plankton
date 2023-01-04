### Source code for article titled: Computing Marine Plankton Connectivity under Thermal Constraints


We used [Ocean Parcels](https://oceanparcels.org) to run the simulations with high resolution model data provided by CMCC. The code used for this study is organized in folder structures. 

```data/```: Data files for sample stations in the Atlantic Ocean and thermal constraints settings used in this study. 

###### Programs
Code for the analysis can be used in the following order-

1. ```processing/AtlanticReleaseLocationsH3.ipynb ```: To create the set of release locations and mapping of particles to H3 grids HexId.
2. ```simulations/Atlantic_TSSampling_Args.py```
3. ```processing/MonthlySparseTM.py``` and then ```processing/AnnualAverageTM.py```: to compute binary connectivity and temperature matrices from simulation outputs. 
4. ```analysis/AllStationsConnectivities.py```: compute connectivity between stations for non-constrained and thermally constrained scenarios
 Bootstrap
5. ```processin/Bootstrap_MonthlyTM.py```, ```processing/Bootstrap_Annual_TM.py``` and ```analysis\BootParticleEnsembles.py``` for the bootstrapping analysis.    

###### Figures
Plots in the manuscript have been generated using the following codes:
- Figure 1: ```visualizations/PlotAtlanticStations.py```
- Figure 2A: ```visualizations/ReleasePointsGridding.py```
- Figure 4A: ```visualizations/PlotAllStationConnectivities.py```
- Figure 4B, G-J and 7C-D: ```visualizations/StationsPairConnectivity.py```
- Figure 4C-F, 7A-B: ```visualizations/PairFractionalChange.py```
- Figure 5: ```visualizations/TimeMeanVariables.ipynb```
- Figure 6: ```visualizations/BootParticleSensitivity.py```
- Figure 8, 9: ```visualizations/PlotFractionalChange.py```
