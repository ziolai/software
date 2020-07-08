Abstract
The objective of this humble endeavor is to describe how the celebrated Pitz-Daily incompressible tutorial (or any its variants) for simpleFoam can be extended with a Reynold stress tensor model for turbulence. This endeavor has been pursued by many (including the giants in this community) in the past and is by no means new. Our endeavor can in fact be stated as reproducing the tutorial [1,2] (among possible various others). Our motivation for stirring in this old pot once again is contribute to the documentation on how to set-up the case. We wish to complement existing tutorials with our insights. The Reynolds stress model solves for the six (in 2D four) components of the symmetric stress tensor R and the turbulent dissipation epsilon. The magnitude of the tensor R is equal to epsilon. Our biggest take-away lessons is that the over-relaxation of the equation (should this field?) for epsilon is essential to obtain convergence of the model. 

Joint work with Kundan Mishra (iamkundan@gmail.com)


1/ Introduction
We aim at extending the Pitz-Daily incompressible tutorial for simpleFoam with a Reynold Stress Tensor Model (RSM) for turbulent flow closure. Motivation for this endeavor is to be able to increase resolution in modeling flow with large recirculation. The Reynolds stress model is a six turbulence equation model for the six components of the symmetric Reynolds stress tensor. The Launders, Reece and Rodi variant of this model [3,4] is implemented as the LRR model in OpenFoam. Previous discussions on the use of the LRR model in this forum include [5,6]. Below we outline features of the the Pitz-Daily test case in Section 2, our understanding on how the model should be set up starting from the tutorial model in Section 3 and on how the model converges with the LRR RASModel activated. 

2/ Pitz-Daily out of the box in OpenFoam-v1906 
The Pitz-Daily test case as defined in OpenFoam-v1906 converges in 283 iterations (run simpleFoam) and meet the criterium for y+ on the mesh (run simpleFoam -postProcess -func yPlus -latestTime). The 283 iterations differs significantly from the tutorial [2]. This difference is likely due to different settings used.  

3/ Case Set-Up 
In this section we explain how the Pitz-Daily tutorial can be modified to run with the LRR RASModel in three steps. 

Step 1/3: Setting up the system-folder
In the first step we set up the system-folder by editting the system-folder/fvSolution-file and the system-folder/fvSchemes-file. We edit the system-folder/fvSolution-file and add settings to solve for R including relaxation as described in e.g. [5]. We  edit the system-folder/fvSchemes-file and add setting to discretize div(R) and div(phi,R) and settings for the wallDist dictionary. Sample settings for  div(R) and div(phi,R) can be copied from other tutorials. To find such setting, use e.g. for i in `find /opt/OpenFOAM/OpenFOAM-v1906/tutorials/ -name fvSchemes`; do echo $i; grep $i "div(R)" ; done . 

Step 2/3: Setting up the 0-folder 
In the second step we edit the 0-folder. We do *not* define the boundary conditions for the components of the tensor R for inlet, outlet and wall patches directly. Instead, we resort to a two-substep procedure in which in the first substep, a simulation using a two-equation turbulence model runs for a number of iterations. In the second substep, the stress tensor R is computed in post-processing stage for the latest iteration. This computed stress tensor is used as boundary value setting for the simulation using the LRR model.  

More specifically, in the first substep, we run simpleFoam on the out-of-the-box tutorial. At the last iteration (iteration 283 in our case) we compute the stress tensor R by post-processing computed results and running simpleFoam -postProcess -func R -latestTime. This command produce a file called turbulenceProperties: R in the latest time folder (folder 283). This file can be renamed to R and copied to the 0-folder. Observe that the settings for the turbulence equation used in the first run is used as input to compute R. One can thus *not* modify the constant-folder/turbulenceProperties-file prior to computing R. One can visualize components and magnitude of the tensor R in paraview. 

Step 3/3: Running the simulations
In the third third we (prepare to) run simpleFoam with the LRR model. We edit the constant-folder/turbulenceProperties-file to change RASModel from current value into LRR. We edit the system-folder/fvSolution file and choose sufficiently low relaxation parameters for R and epsilon. We can now run simpleFoam with LRR. Alternatively, we can avoid copying the 283-folder/R-file to the 0-folder and run simpleFoam with LRR and startTime equal to the latestTime of the previous simulation. 

4/ Discussion and Conclusions 
We discussed how the incompressible tutorial Pitz-Daily for simpleFoam can be modified to include the Reynolds stress tensor wonder. We wonder whether the preprocessing using a two-equation turbulence model can be avoided. 

The above description is likely to be incomplete and fragile as time evolves. We theferore wonder whether a tutorial using the LRR model can be made available as tutorial in the OpenFoam distribution. 

References 
[1]: Turbulence Steady State, Tutorial Six of TU Vienna: https://www.cfd.at/sites/default/files/tutorialsV7/6-ExampleSix.pdf: shows change divSchemes in fvSchemes to take R into account: no settings explained for fvSolution: decrease in turbulent viscosity is shown;
[2]: Introduction to Stationary Turbulence Modeling, Jozsef Nagy, https://www.youtube.com/watch?v=-46pgYQYER8. The LRR model is discussed from minute 18 of 26 onwards  
[3]: NASA Langley Turbulence Modeling Resource on Launders, Reece and Rodi model: https://turbmodels.larc.nasa.gov/rsm-ssglrr.html ; 
[4]: Ansys Fluent Manual online   
[5] Jasak on cfd.online.com on how to use the LRR Model: details of how to run the LRR model and showing convergence of the residual for the stress tensor components:  https://www.cfd-online.com/Forums/openfoam-solving/60032-problems-rsm-simplefoam.html , 2006. 
[6]: discussion regarding the LRR RASModel on cfdonline.com: settings for volSymmTensorField R in fvSolution are shown. No settings for fvSchemes are given unfortunately: https://www.cfd-online.com/Forums/openfoam/86577-usage-r-lrr.html , 2010. 
