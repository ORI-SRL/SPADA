# SPADA
SPADA (**S**oft **P**neumatic **A**ctuator **D**esign fr**A**mework) is a toolbox developed to facilitate the design, simulation and optimisation of bellow soft pneumatic actuctors (bellow-SPAs) for shape matching. It consists of two tabs: simulation and optimisation. 
- Simulation:
  - Design bellow-SPAs by tunning geometric parameters, stacking modules and defining material properties.
  - Use background COMSOL Finite Elment Method (FEM) to simulate the deformation of the designed bellow-SPA based on the input interal pressure. It also provides information on material stress and strain of modules.
  - Generate material-specified FEM datasets.
 
- Optimisation:
  - Train an Artificial Neural Network using the FEM dataset as training data to create a surrogate model.
  - Approximate a desired shape, which comprises sequential 3D coordinates of points, using constant curvature segments.
  - Optimise bellow-SPA design to match the shape based on the surrogate model.
  - Generate a ready-to-print CAD file.

## Intsallation Requirements
- Operating System: Windows 10/11
- COMSOL Multiphysics: Version 5.6 or later, with the following modules:
  - Nonlinear Structural Materials Module
  - Matlab Livelink
- Matlab: Version 2020b or later, with the following toolboxes installed:
  - Deep Learning toolbox
  - Global Optimization toolbox

### For the COMSOL-MATLAB integration, follow these steps:
1. Open COMSOL and go to Preferences.
2. Navigate to Livelink Connection and select LiveLink for MATLAB.
3. Choose the folder path for MATLAB 2020b or later.


## Setup
To set up and run SPADA, follow these steps:

0. Make sure your system meets the requirements
1. Clone and install this repository
2. Open COMSOL Multiphysics Matlab Livelink
3. Add the folder path of SPADA toolbox in Matlab
4. Execute the app.mlapp file to run the SPADA application.

## Explore the examples
To explore the examples provided in SPADA, you can refer to the following folders:
- FEM Datasets: The "DataSet_Agilus30.mat" file contains a dataset for the Agilus30 material.
- Input Shapes: There is a collection of example shapes available in ".mat" files. These files contain the sequential 3D coordinates of points that define the desired shapes. Additionally, you will find corresponding ".m" files that provide code to generate those shapes programmatically.
- CAD: A set of CAD files is available that correspond to the input shapes, provided in the ".stl" format.

## Demonstration
To view a demonstration of SPADA and some designed actuators in action, please refer to the video "SupplementalVideo.mp4" provided with the project.

## Citation
When citing SPADA, use this citation:

```bibtex
@misc{yao2023spada,
      title={SPADA: A Toolbox of Designing Soft Pneumatic Actuators for Shape Matching based on the Surrogate Model}, 
      author={Yao Yao and Liang He and Perla Maiolino},
      year={2023},
      eprint={2305.19509},
      archivePrefix={arXiv},
      primaryClass={cs.RO}
}
```


Additionally, the following work provides guidance on optimisation algorithms that complement the SPADA toolbox:
- Y. Yao, Y. Chen, L. He and P. Maiolino, "Design Optimization for Bellow Soft Pneumatic Actuators in Shape-Matching," 2023 IEEE International Conference on Soft Robotics (RoboSoft), Singapore, Singapore, 2023, pp. 1-7, doi: 10.1109/RoboSoft55895.2023.10121957.
