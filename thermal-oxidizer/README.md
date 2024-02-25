# Thermal Oxidizer of Acid Gasses CFD Model

## Section 1: Introduction 
Here goes 
1. an introduction of the device: thermal processing of [hydrogen sulfide](https://en.wikipedia.org/wiki/Hydrogen_sulfide) (an example of an acid gas) in typical [treatment plants](https://f.hubspotusercontent10.net/hubfs/7724363/Zeeco-En_June_2020/pdf/Gas-Processing-Apps-TO-RGB.pdf); 
2. the need to simulate the acid gas removal process in order to optimize fuel efficiency, fuel versatility and throughput; 
3. embedding into previous literature on OpenFoam for chemical reactors and combustion systems; 
4. summary of the results obtained; 

## Section 2: Geometry and Mesh

See [notebook](./notebooks/geometry-mesh.ipynb).

## Section 4: Tutorial Models We Looked Into 

<b>compressible/rhoPimpleFoam/RAS/angledDuct</b> and <b>compressible/rhoPimpleFoam/RAS/angledDuctLTS</b>
1. why is the mixingLength of the inlet set to 0.005 = 5e-3 ? 
2. why is value on inlet patches set to internalField value? What is the importance of this setting?  

## Section 5: Alternative Choices for Model Setup 
1. Set internalField for U, k, epsilon and nut equal to zero;
2. Should relaxation factors be adapated;  
3. Should discretization schemes be adapted? Copy those from Pitz-Daily? Use the banana method to recover all avialable methods. Do settings change going from simpleFoam to rhoSimpleFoam? 

## Section 5: Incompressible s/pimpleFoam Models 

See [notebook](./notebooks/incompressible-simplefoam-models.ipynb).

## Section 6: Temperature-Dependent rhoS/PimpleFoam Models 

See [notebook](./notebooks/temperature-rhosimplefoam-models.ipynb).

## Section 7: reactingFoam Computations with Combustion Switched Off

## Section 8: reactingFoam Computations with Combustion Switched On 


```julia

```
