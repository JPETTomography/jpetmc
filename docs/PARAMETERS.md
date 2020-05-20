# Parameters for MC simulations with Geant4 macros
Following options can be added to macro files, that are read by Geat4. Example files are in `scripts` folder, and also those files are copied to `bin` directory during program build.  

## Using different geometries  
* `/jpetmc/detector/loadJPetBasicGeom`  
 * 3 layers of scintillators (48, 48, 96)  
 * each scintillator: 1.9x0.7x50 cm^3 wrapped in kapton foil  
 * detector frame: loaded from CAD file  
 * standard setup for runs 1-7  
* `/jpetmc/detector/loadGeomForRun runNr`  
 * loads standard geometry `JPetBasicGeom`  
 * and the chamber specific for desired measurement: `runNr=3/5/6/7`  
* `/jpetmc/detector/loadJPetExtendedGeom`  
 * loads the geometry of the 4th layer - modular detector
* `/jpetmc/detector/loadOnlyScintillators`  
 * shows only scintillator wrapped in a kapton foil  
* `/jpetmc/detector/loadIdealGeom`  
 * not implemented yet  
 * should load n-layer fully packed with scintillators barrel  

## General parameters:  
* Hit merging time:  
 * define time range, between hits in scintillator, which will be classified as single hit:  
 * `/jpetmc/detector/hitMergingTime`  
* Adding date and time to the name of the output file, so multiple executions of the simulation does not overwrite the default file  
 * `/jpetmc/output/AddDatetime 1`  

## Additional parameters:
* simulate only oPs gamma:  
 `/jpetmc/material/threeGammaOnly`  
* simulate only 2 gamma:  
 `/jpetmc/material/twoGammaOnly`  
* simulate only pick-off process (2g & long lifetime):
 `/jpetmc/material/pickOffOnly`  
* adding oPs component for material given further by reload material command
 `/jpetmc/material/oPsComponent (mean lifetime in ns) (probability 0-100 %)`  
* setting pPs component for material given further by reload material command
 `/jpetmc/material/pPsComponent (mean lifetime in ns) (fraction of pPs to oPs)`  
* adding oPs component for material given further by reload material command
 `/jpetmc/material/directComponent (mean lifetime in ns) (probability 0-100 %)`  
* reloading lifetime parameters for material given by command material (xad4, kapton, aluminium, plexiglass, pa6)
 `/jpetmc/material/reloadMaterials (material)`  
* set center of the annihilation chamber:  
 `/jpetmc/run/setChamberCenter`  
* for run5: define range where we expect annihilation to occur:   
 `/jpetmc/run/setEffectivePositronRange`  
* save true(!) generated events based on multiplicity (0,2-10):  
 `/jpetmc/event/saveEvtsDetAcc`
* change lower value of saved multiplicity (0):
 `/jpetmc/event/minRegMulti`
  (above valid only with: /jpetmc/event/saveEvtsDetAcc)
* change upper value of saved multiplicity (10):
 `/jpetmc/event/maxRegMulti`
  (above valid only with: /jpetmc/event/saveEvtsDetAcc)
* change excluded value of multiplicity (1):
 `/jpetmc/event/excludedMulti`
  (above valid only with: /jpetmc/event/saveEvtsDetAcc)
* print how many events were generated  
 `/jpetmc/event/printEvtStat`  
* print out option during execution of the simulation - X in divisor (10^X) for number of printed events  
 `/jpetmc/event/printEvtFactor`  
* show generation progress (in %)  
 `/jpetmc/event/ShowProgress`  
* Give nema point number to simulate (1-6)  
 `/jpetmc/source/nema`  
* set parameters of gamma beam  
 `/jpetmc/source/gammaBeam/setEnergy setPosition setMomentum`  
* set parameters of the isotope source, shape, if cylinder: radius, z-lenght (half)  
 `/jpetmc/source/isotope/setShape cylinder`  
 `/jpetmc/source/isotope/setShape/cylinderRadius`  
 `/jpetmc/source/isotope/setShape/cylinderZ`  
 `/jpetmc/source/isotope/setPosition`  
* set number of gamma quanta to generate 1 / 2 / 3 by the isotope  
 `/jpetmc/source/isotope/setNGamma`  
* generate only cosmic muons  
 `/jpetmc/source/cosmiscOnly`  
