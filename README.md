# TrueProbes
Quantitative RNA-FISH Probe Design Software



## Dependencies    
[MATLAB 2022b (or higher)](https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml)    
[MATLAB Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)    
[MATLAB Bioinformatics Toolbox](https://www.mathworks.com/products/bioinfo.html)    
[MATLAB Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html)    
[MATLAB Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)    
[MATLAB Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)    
[MATLAB Curve Fitting Toolbox](https://www.mathworks.com/products/curvefitting.html)    
[NCBI-BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)    

### Project Directories:

To run the `data/DatabaseData` needs to be included in the main TrueProbes folder.
Directory: `DatabaseData` contains the following: 
1. **Metadata**
   - Metadata annotation files for each organism


Directory: `src` contains the following: 
1. **Metadata**
   - Metadata annotation files for each organism
2.  **DatabaseData**
   - Database annotation files for each organism with Blast Databases, Reference Genome/Transcriptome, and Expression Data
   
Directory: `src` contains the following: 
1. **blast**
   - installers for blast+
2. **core:**
   - MATLAB codes for generating each step of the TrueProbes pipeline
   - steps: (Probe Generation, Probe Blasting, Gene Expression, Thermodynamic Information,
     Binding Site Mapping, ProbeDesignerStats, Probe Selection, Probe Metrics) 
3. **third-party:**
   - Other MATLAB and linux software scripts and packages  
4. **util:**
   - general MATLAB functions utilized throughout software
5. **wrappers:**
   - general MATLAB wrapper functions

The TrueProbes software runs by performing eight steps sequentially.     
	1. Probe Generation. It generates all possible probes shared between a list of inclusion IDs and inclusion te xt files, but not in exclusion text files or exclusion IDs, within a set probe length range.   
	2. BLAST Alignment. All probes target hits in the reference genome and/or transcriptome are identified at least as long as the minimum homology length.    
	3. Binding Affinity Calculation. Binding affinities are calculated for all pairs of probe and target sequences homology matches in the BLAST results.   
	4. BLAST Target Gene/Transcript Expression. The gene expression and transcript expression levels are collected for reference expression databases specified in the TrueProbes settings for all targets in the BLAST results.   
	5. Probe-Target Binding Site Mapping. All blast hits and binding affinities are converted into a site-specific binding map to generate a formatted map by relative binding site position on each target gene, transcript, or chromosome.   
	6. Probe-Target Statistics. Generate statistics on blast hits, thermodynamics, and a comparison of probes sharing off-targets and relative trade-offs when quantifying off-targets by probe and comparing probes that bind them in a site-specific manner.    
	7. Probe Design. Sort probes by with/without expression data, the number of off-targets, and then by difference in on-target binding to off-target binding and secondary structure binding affinity to iteratively design probes. Print out the list in an Excel spreadsheet.   
	8. Model Evaluation. The final probe set and reference expression are combined to compute equilibrium probe binding and statistics, cumulative off-target binding, on-target binding, etc., when given reference values for cell size, probe concentration, and probe intensity.    

Probe design is run via the command line using A0_BKJH_ProbeDesign_Wrapper_cluster_V5(id,cluster).

The design file needs to be run with two inputs (id and cluster), with a table (inputs1) describing each target gene (inputs1), 
and a set of settings describing how the probe design will be run (settings).

Main Arguement:
	id: id is an integer and is the row of the input design table at the top of the script to run and design probes against. 
	cluster: Cluster is an integer and determines the software parallelization pool between local or remote servers when running the script via Slurm. 
		The only difference between running on a cluster is the number of cores in the Slurm file. cluster = 0 if run locally, and 1 if run remotely
	
	design_file: An excel spreadsheet file with the input files and parameter settings for the probe design
		Inputs1 Table: inputs1 is the table with the list of targets probes will be designed for. The table row has nine column entries. 
		Input table columns:
			1. Organism to design probes for.
			2. Included target accession IDs. Designs probes shared across all accession numbers 
			3. Excluded target accession IDs. Removes probes in exclusion accession numbers.			
			4. Text Sequence Files to Include (files). Default empty
			5. Text Sequence Files to Exclude (files). Default Empty

	
	Example: {'Human',{NM_000805.5'},{},{},{},;}...    

Settings.
	Below the inputs1 table is a list of primary and secondary settings which can be configured to determine how the probes are designed and evaluated.

Input Parameters:

**Main Probe Design Parameters:**    
- minProbeSize: min nt length of potential probes, default 20    
- maxProbeSize: max nt length of potential probes, default 20     
- MinProbeSpacing min spacing between probes, default 3]    
- MaxNumberOfProbes: Max number of probes to design, default 96    
- targetStrand: which strand to of target to design probes against, default 1 for ‘plus’, with 0 for ‘minus’]   
- MinHomologySearchTargetSize:Minimum homology length for BLAST alignments to be recorded and used in probe design and evaluation, default 15   
- BLASTrna: decide to blast RNA sequences, default 1  
- BLASTdna: decide to blast DNA sequences, default 0   
- ExpressionReferenceForProbeDesign:[which row across all expression reference files to use in probe design, with 0 meaning to not use expression data to design probes, default 0]    
	
Thermodynamic Settings:    
	Gibbs_model:[Which thermodynamic model to use for probe design and evaluation, default 4]    
		Model 1:Breslauer86. Breslauer K.J., Frank R., Blocker H., Marky L.A., (1986) Predicting DNA duplex stability from the base sequence Proc Natl Acad Sci U S A 83, 3746-3750    
		Model 2:SantaLucia96. SantaLucia, J., Allawi, H. T., and Seneviratne, P. A. (1996) Improved nearest-neighbor parameters for predicting DNA duplex stability Biochemistry 35, 3555-3562    
		Model 3: SantaLucia98. SantaLucia, J. (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics Proc Natl Acad Sci U S A 95, 1460-1465   
		Model 4:Sugimoto96. Sugimoto, N., Nakano, S., Yoneyama, M., and Honda, K. (1996) Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes Nucleic Acids Res 24, 4501-4505    
		Model 5:SantaLucia04. SantaLucia Jr, J., and Hicks, D. (2004) The thermodynamics of DNA structural motifs Annu Rev Biophys Biomol Struct 33, 415-440    
		Model 6: Allawi97. Allawi, H. T., and SantaLucia, J. (1997) Thermodynamics and NMR of internal G.T mismatches in DNA Biochemistry 36, 10581-10594   
		Model 7: Rejali21. Rejali, N. A., Ye, F. D., Zuiter, A. M., Keller, C. C., and Wittwer, C. T. (2021) Nearest-neighbour transition-state analysis for nucleic acid kinetics Nucleic Acids Res 49, 4574-4585   
		Model 8: Martins24. de Oliveira Martins, E., and Weber, G. (2024) Nearest-neighbour parametrization of DNA single, double and triple mismatches at low sodium concentration Biophys Chem 306, 107156   
	HybridizationTemperature:[Hybridization temperature, default 37C]   
	HeatCapacityReferenceTemperature:[Reference temperature for Cp measurement and Gibbs model, default 37C]   
	SaltConcentration:[Salt Concentration M, default 0.05]   
	PrimerConcentration:[Primer Concentration M, default 50e-6]   
	RemoveMisMatches:[Remove sequence mismatched base pairs before evaluating probe-target binding affinity, default 1,  including adds flanking sequences to alignments and uses Gibbs model 8 with mismatch base pair inclusion]   

Design Filtering Settings   
	RemoveProbesBindingOffTargetRibosomalHits:[Filter out probes with off-targets to ribosomal proteins, default 1]   
	packOptimal_ProbesWithNoOffTargets: [When designing probes without off-target use optimal packing to get as many probes with no off-targets as possible as opposed to normal sequential selection, default 1]   
	IncludeSelfHybridizationInProbeSelection:[When designing probes consider probe self-hybridization when ranking probes based on binding affinity, default 1]    

Parallelization Parameters: (Usually only changed when using longer genes)   
	Parallization_probeBatchSize: [number of probes to evaluate in a single batch when performing parallelized calculations, default 20]   
	Parallization_targetBatchSize: [number of targets to evaluate in a single batch when performing parallelized calculations, default 200]   
	ParsingPreference:[blast simultaneously in parallel(1) or blast probes sequentially (0), default 1]   
  
Gene Expression Parameters:	  
	UseGeneOverTranscLevelExpression:[use gene level (1) or transcript isoform level (0) gene expression values, default 0]   
	DoAllGenesHaveSameExpression:[decide to assume equal expression for all genes (1) or to use gene expression reference (0) , default 1]   
	UseRegularDNAExpression:[(0) use DNA expression from gene expression track in expression data, (1) set expression to 2 for DNA, default 1]   
	nullRNAcopynumber: [number of RNA copy when not using reference expression levels, default 100]   
	nullDNAcopynumber: [number of DNA copy number when not using reference expression levels, default 2]   
	TMM_LogRatioTrim:[when normalizing TPM expression data using TMM set log ratio trim threshold cutoff, default 0.3]   
	TMM_SumTrim:[when normalizing TPM expression data using TMM set sum trim threshold cutoff, default 0.05]   
	TMM_Acutoff:[when normalizing TPM expression data using TMM set A cutoff value, default -1e10]   
	TMM_doWeighting:[when normalizing TPM expression data weight terms using inverse of approximate asymptotic variance of the M-values to account for genes with higher read counts having lower variance on log scale and more reliable mean estimation, default 1]   
	
Make Blast Database Settings:      
	Parse_seqids: [when making blastdb in software decision to parse sequence ids from fasta files as blastdb ids, useful for blastdbcmd sequence retrieval, default true]   
	Hash_index:[when making blastdb make sequence hash indexes leads to faster exact match retrieval but less accurate range matches, default false]   

BLAST Parameters:    
	evalue:[Expectation value cutoff, default 1000]    
	Dust:[filter query sequences with DUST, default no]   
	Gapextend:[Cost to extend a gap (integer), default 2]   
	Gapopen:[cost to open a gap (integer), default 5]   
	Num_alignments:[number of database sequences to show num_alignments for, default 1000]   
	Penalty:[penalty for a nucleotide mismatch, default -3]   
	Reward:[reward for a nucleotide match, default 1]   
	Word_size:[word size for wordfinder algorithm, default 7]   

Model Simulation Settings    
	removeUndesiredIsoformsFromPrediction: [when evaluating main on/off-target binding should alternate isoforms of desired targets be removed and not included in off-target quantification, default 1]    
	ProbeConcentration_MicroMolar:[probe concentration in uM, default 5e-6]   
	CellRadius_Micron: [cell radius in microns for converting RNA molecule counts into concentration for solving binding equilibrium, default 10]   
	Dilution_Vector: [vector of probe dilutions to evaluate binding equilibrium and predictions at, default 1,1e-2,1e-4]   
	Gibbs_Model_Vector: [vector of gibbs thermodynamic models to use for evaluating binding equilibrium and predictions at,  default 1,2,3,4]    
	Temperature_Celsius_Model_Vector: [vector of hybridization temperatures in Celsius to use for evaluating binding equilibrium and predictions at, default 37,42,50,60]   
	InitialFreeSolutionGuessConcentration_MicroMolar:[initial solution guess for steady state free probe concentration in uM, default 1e-10]    
	SolutionErrorTolerance: [total tolerance level for error in final equilibrium solution, default 1]   
	MaxRecursiveEquilibriiumIterations: [max number of equilibrium equation iterations stopping at final calculated steady state, overrides solution error tolerance]   
	CellDiameter_Pixels: [cell pixel diameter for background predictions, default 50]   
	SpotRadius_Pixels: [ spot pixel radius for predictions, default 5]    
	NumberOfReferenceZStacks:[ number of z-stacks to spread background across for intensity predictions, default 67]   
	SignalStepSize: [step size for signal intensity bins in signal predictions, default 1e-1]    
	SignalMaxValue: [max intensity value for ranges in solution, default 3000]    
	AutoFluoresenceBackground_MEAN:[reference mean autofluorescence, default 278]   
	AutoFluoresenceBackground_STD: [reference autofluorescence standard deviation, default 33]   
	NumberOfProbesInReferenceSpots:[ number of probes in reference spot intensity for calibrating intensity predictions, default 48]   
	ReferenceSpotIntensity_MEAN: [mean reference spot intensity , default 827 ]   
	ReferenceSpotIntensity_STD: [reference spot intensity standard deviation, default 28]   

TrueProbe Design Software Output Files:    
	(GeneName)_AccessionID_probes_TrueProbes.mat [structure with probe sequences, location on on-target]   
	(GeneName)_AccessionID_hits_table_TrueProbes.mat [structure with information on BLAST hits]   
	(GeneName)_AccessionID_ExpressionInfo_TrueProbes.mat [Structure with expression data]   
	(GeneName)_AccessionID_dKbInfo_TrueProbes.mat [Structure with heat capacity for all target binding reactions]   
	(GeneName)_AccessionID_dCpInfo_TrueProbes.mat [Structure with heat capacity for all target binding reactions]   
	(GeneName)_AccessionID_dHInfo_TrueProbes.mat [Structure with enthalpy for all target bindings reactions]   
	(GeneName)_AccessionID_dSInfo_TrueProbes.mat [Structure with entropy for all target binding reactions]   
	(GeneName)_AccessionID_TmInfo_TrueProbes.mat [Structure with entropy for all target binding reactions]   
	(GeneName)_AccessionID_Tm[HybridizationTemperature]_OnOffThermoInfo_TrueProbes.mat [Structure with binding energy of all hits]   
	(GeneName)_binding_hits_map_TrueProbes [binding site map]   
	(GeneName)_AccessionID_BindingMatricies_TrueProbes.mat [Entropy, Enthalphy, and heat capacity in binding site map format for RNA]  
	(GeneName)_AccessionID_Tm[HybridizationTemperature]__BindingEnergyMatrix_TrueProbes.mat [Equilibrium Binding Energy in binding site map format]   
	(GeneName)_AccessionID_BindingMatricies2_TrueProbes.mat [Entropy, Enthalphy, and heat capacity in binding site map format for complementary DNA strand binding]   
	(GeneName)_AccessionID_Tm[HybridizationTemperature]__BindingEnergyMatrix2_TrueProbes.mat [Complementary DNA strand Equilibrium Binding Energy in binding site map format]   
	(GeneName)_AccessionID_Tm[HybridizationTemperature]_BasicDesignerStats_TrueProbes.mat [Index information on stats used for design probes]   
	(GeneName)_AccessionID_chosen.mat [List of chosen probe indexes]   
	(GeneName)_AccessionID_probes_final_MaxProbeNumber_max.xlsx [Excel spreadsheet with final probes, and some stats]   
	(GeneName)_AccessionID_Tm[HybridizationTemperature]_ModelMetrics_TrueProbes.mat [Structure with binding affinity calculations and probe design metrics]   
