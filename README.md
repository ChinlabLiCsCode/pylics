# PyLiCs

PyLiCs will become a Python package containing Python tools for analysis and control of the Chin Lab Lithium-Cesium experiment. The code in this repository will initially be primarily a replacement for the existing MATLAB codebase (see [github.com/ChinLabLiCs/lics-codebase](https://github.com/ChinlabLiCsCode/lics-codebase)). However, the goal is to expand the capabilities of the codebase to include new features and tools for the analysis and control of the experiment.

## Open questions 

- How to handle daily params & other things that change over time? 
- How package things? Do all analysis notebooks, etc. go in the main repo? 

## Task list

- Port over MATLAB code from the existing codebase
- Build in tests for all functions and classes
- API 
- Extend labview sequence modification capabilities
    - Labview plotter interface
    - Modify 65500 parameters 
    - Copy subsequences 
- Make a requirements file