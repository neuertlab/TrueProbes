# TrueProbes
Quantitative RNA-FISH Probe Design Software



## Dependencies
[MATLAB 2023b (or higher)](https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml)
[MATLAB Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)
[MATLAB Bioinformatics Toolbox](https://www.mathworks.com/products/bioinfo.html)
[MATLAB Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html)
[MATLAB Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)
[MATLAB Curve Fitting Toolbox](https://www.mathworks.com/products/curvefitting.html)
[BLAST+ Support Package for Bioinformatics Toolbox](https://www.mathworks.com/products/curvefitting.html)

### Project Directories:

To run the `data` needs to be included in the main TrueProbes folder.
Directory: `src` contains the following: 
1. **Metadata**
   - Metadata annotation files for each organism
2.  **DatabaseData**
   - Database annotation files for each organism with Blast Databases, Reference Genome/Transcriptome, and Expression Data
   
Directory: `src` contains the following: 
1. **built-in-functions**
   - built in MATLAB functions used
2. **core:**
   - MATLAB codes for generating each step of the TrueProbes pipeline
   - steps: (Probe Generation, Probe Blasting, Gene Expression, Thermodynamic Information,
     Binding Site Mapping, ProbeDesignerStats, Probe Selection, Probe Metrics)
3. **modified-matlab-functions:**
   - Modified MATLAB functions to add additional outputs 

4. **third-party:**
   - NCBI Blast scripts, and other MATLAB software scripts and packages  

5. **util:**
   - general MATLAB functions utilized throughout software
   
6. **wrappers:**
   - general MATLAB wrapper functions


1. Check that you have permission to read and write the f