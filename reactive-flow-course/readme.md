# Conversation with Eric Daymo on Reactive Flow Course

##  Section 1: Introduction 

## Section 2: Course Description 

### Section 1.2: Target Audience 
Assume prior exposure to use of OpenFoam for 
1. (in)stationaryturbulent incompressible flow (s/pimleFoam) (well familiar; geometry definition, mesh generation, case set up, post-processing, residual control);  
2. buoyant heat transfer (bouyantS/PimpleFoam) (prevous exposure); 
3. radiative heat transfer (P1 and fvDOM) (limited exposure);
4. conjugate heat transfer (fluid/solid domain with interface coupling) (limited exposure); 

Course will use [Foundation Version 11](https://github.com/OpenFOAM/OpenFOAM-11) of OpenFoam. 

Virtually no coding assignment. Course limited to setting up cases, running cases, interpretation of results, discussing what-if scenarios. 

reactingFoam: describe what equations reactingFoam solves. describe how reactingFoam solves these equations? 

What other solvers does the course intend to use (aside from reactingFoam)?  

### Section 2.2: Course Content 
#### Day-1: Copy-Paste from Eric here
#### Day-2: Copy-Paste from Eric here
#### Day-3: Copy-Paste from Eric here

### Section 2.3: Course References 
1. wiki 
2. Ansys Fluent Manual 
3. cfd-online material 
4. Dacolt course notes; 
5. F. Contino course notes; 





## Section 3: reactingFoam Case-Setup Description 
Goal: explain case set-up by explaining dictionaries in 0-dir/files, system-dir/files and constant-dir/files. In particular 

### Section 1.3: reactingFoam case-dir/constant-dir  

1. turbulenceProperties: see e.g. [Turbulence](https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence.html#sec-turbulence-usage) and other sources; 

2. physicalProperties (formerly thermophysicalProperties): see e.g. [Heat transfer](https://www.openfoam.com/documentation/guides/latest/doc/guide-heat-transfer.html#sec-heat-transfer-usage) and other sources; 

3. thermo.compressibleGasGRI

4. radiationProperties: no good reference see. See [Radiation](https://www.openfoam.com/documentation/guides/latest/doc/guide-heat-transfer.html)

5. boundaryRadiationProperties

6. chemistryProperties: describe [chemkinToFoam](https://github.com/OpenFOAM/OpenFOAM-11/blob/master/applications/utilities/thermophysical/chemkinToFoam/chemkinToFoam.C) and suggest cfd-online forum for trouble shooting;   

7. reactionsGRI

8. combustionProperties

9. momentumTransport: allows Le <>1 and/or Sc <>1. No documentation aside from H-file; (or use Cantera/CoolProp/Duochem to compute coefficients); [example file](https://github.com/OpenFOAM/OpenFOAM-11/blob/master/tutorials/multicomponentFluid/SandiaD_LTS/constant/momentumTransport); 

#### Specifying Le<>1 in physicalProperties (or thermophysicalProperties?)
\begin{verbatim}

    RAS
    {
        model           FickianEddyDiffusivity;
 
        mixtureDiffusionCoefficients yes;
 
        Prt             0.85;
        Sct             0.7;
 
        Dm // [m^2/s]
        {
            O2 6.64443705E-05;
            H2 4.22575501E-04;
            N2 4.74274421E-05;
        }
 
    }
\end{verbatim}

### Section 2.3: reactingFoam case-dir/system-dir

Explain relevant parts of fvSolution (solvers.Yi{} and PIMPLE{Yref{}})  and of fvSchemes (divSchemes.div(phi,Yi) and divSchemes.div(phi,K)). Refer to other courses for parts not explained here. 

### Section 3.3: reactingFoam case-dir/0-dir 

Explain Ydefault and species-name. 

## Section 4: Laminar (non-premixed? hydrogen) ToroFlames:
See [laminarSMOKE](https://github.com/acuoci/laminarSMOKE) and 
[ToroFlames/F3](https://github.com/acuoci/laminarSMOKE/tree/master/run/validation/ToroFlames/F3) 

Validation by Eric using DUO. Need to include results here.

## Section 5: Turbulent Small-Scale Non-Premixed Combustion: SandiaD test case 

### Section 1.5: Test case description 
1. describe the case: non-premixed combustion on CH4; pilot-burner stabilized; reference results from measurements of TNF workshop and from literature; 
2. transient solution to reach steady state; 
3. geometry definition, mesh generation and mesh check;
4. describe patch settings: inlet, outlet, walls, symmetry settings;  
5. running the non-reactive and reactive case; (with and without LTS, with and without tabulation and chemistry reduction); 
6. post-processing results and residual; 

### Section 2.5: Numerical Considerations 

1. compare runs with one-step chemistry, two-step chemistry, various forms of four-step chemistry, GRI-2.0 and GRI-3.0;
2. compare thermal NO post-processing with NO computed from chemical mechanism; 
3. compare runs with and without local time stepping (LTS); is this reported in the literature? 
4. compare runs with and without TDAC chemistry tabulation and reduction: literature? 
5. compare runs on various meshes using remap from coarse to fine mesh; 
6. compare various wall temperature boundary conditions (insulating, radiating); 
7. extend to (D)DES and LES; 
8. extend to other flames (SandiaE and SandiaF) in TNF workshop series; 

### Section 3.5: Reference results
1. [report](borujeni-norouzi-SandiaD-reactingFoam.pdf) by Borujeni and Norouzi; 
2. [report](talice-SandiaD-reactingFoam.pdf) by Marco Talice using ESI v1906; 
3. [paper](http://web.stanford.edu/group/pitsch/publication/PitschPhysFluids_12_2000.pdf) by Pitsch and Steiner using LES and Lagrangian Flamelet Model; 

## Section 6: Turbulent Large-Scale Non-Premixed Combustion: DLR test case|

[DLRCJH](https://develop.openfoam.com/committees/hpc/-/tree/develop/combustion/XiFoam/DLRCJH) 

## Section 7: Turbulent Non-Premixed Combustion with Heat Transfer to Solid: reverseBurner test case 

[reverseBurner](https://mega.nz/fm/LJd1XQoS)


```julia

```


```julia

```
