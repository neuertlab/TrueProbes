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

The TrueProbes software runs by performing eight steps sequentially.
	1. Probe Generation. It generates all possible probes shared between a list of inclusion IDs and inclusion text files, but not in exclusion text files or exclusion IDs, within a set probe length range.
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
			1. Included target accession IDs. Designs probes shared across all accession numbers 
			2. Text Sequence Files to Include (files). Default empty
			3. Text Sequence Files to Exclude (files). Default Empty
			4. Organism to design probes for.
			5. Gene Name 1.
			6. Gene Name 2. Second potential gene name to use instead of the first.
			7. Chromosome. The chromosome number.
			8. Excluded target accession IDs. Removes probes in exclusion accession numbers.
			9. Strand. Which strand to design probes against for RNA default is ‘plus’.
	
	Example: {{NM_000805.5'},{},{}, 'Human','(GAST)','(GAST)','17',{},1  ;}...    

Settings.
	Below the inputs1 table is a list of primary and secondary settings which can be configured to determine how the probes are designed and evaluated.

Main Input Parameters:
	Locations
		SaveRoot: [Location where files are to be saved]
		customBlastDatabase_DNA:[Location of user-custom DNA blast database if user wants to use their own custom database, default N/A]
		customBlastDatabase_RNA:[Location of user-custom RNA blast database if user wants to use their own custom database, default N/A]

Main Probe Design Parameters:
	max_probes:	[Max number of probes to design, default 96]
	minProbeSize: [min nt length of potential probes, default 20]
	maxProbeSize: [max nt length of potential probes, default 20]
	MininumProbeSpacing: [min spacing between probes, default 3]
	BLASTrna: [Will BLAST include reference transcriptomic sequences, default 1]
	BLASTdna: [Will BLAST include reference genomic sequences, default 0]
	HybridizationTemperature: [Hybridization temperature used in probe evaluation default 37C]
	ExpressionReferenceForDesigningProbes: [Which expression from expression data is used for designing probes, default 0] no expression (0).
	
Secondary Parameters:
	Nmodel:	[Which thermodynamic model to use for probe design and evaluation, default 4]
	SaltConcentration: [Salt Concentration mM, default 0.05]
	RemoveRibosomalHits: [Filter out probes with targets to ribosomal proteins, default 1]
	MinHomologySearchTargetSize: [Minimum homology length for BLAST alignments to be recorded and used in probe design and evaluation, default 15]
	removeUndesiredIsos: [When getting probe target statistics used in designing probes omit from computation other isoforms of target transcript, default 1]
	packOptimal: [When designing probes without off-target use optimal packing to get as many probes with no off-targets as possible as opposed to normal sequential selection, default 1]
	UseSelfProb: [When designing probes consider probe self-hybridization when ranking probes based on binding affinity, default 1]
	RemoveMisMatches: [When desinging probes remove mismatches from nearest neighbor quantification of probe binding, default 1]
	RunOffline: [When getting probe sequences use reference genome in database data instead of going online to access genbank record, default 1]

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

BLAST Parameters:
	evalue:[Expectation value cutoff, default 1000]
	Dust:[filter query sequences with DUST, default no]
	Gapextend:[Cost to extend a gap (integer), default 2]
	Gapopen:[cost to open a gap (integer), default 5]
	Num_alignments:[number of database sequences to show num_alignments for, default 1000]
	Penalty:[penalty for a nucleotide mismatch, default -3]
	Reward:[reward for a nucleotide match, default 1]
	Word_size:[word size for wordfinder algorithm, default 7]

TrueProbe Design Software Output Files: 
	(GeneName)_TranscriptID_probes_TrueProbes.mat: [structure with probe sequences, location on on-target]
	(GeneName)_TranscriptID_hits_table_TrueProbes.mat:[structure with information on BLAST hits]
	(GeneName)_TranscriptID_ExpressionInfo_TrueProbes.mat:[Structure with expression in TCGA and GTEX]
	(GeneName)_TranscriptID_Tm37_OnOffThermoInfo_TrueProbes.mat:[Structure with binding energy of all hits]
	(GeneName) _TranscriptID _dCpInfo_TrueProbes.mat:[Structure with heat capacity for all target binding reactions]
	(GeneName)_TranscriptID_dHInfo_TrueProbes.mat:[Structure with enthalpy for all target bindings reactions]
	(GeneName)_TranscriptID_dSInfo_TrueProbes.mat:[Structure with entropy for all target binding reactions]	
	(GeneName)_binding_hits_map_TrueProbes	[binding site map]
	(GeneName)_TranscriptID_Tm37_BindingEnergyMatrix_TrueProbes.mat: [Equilibrium Binding Energy in binding site map format]
	(GeneName)_RefSeqID_BindingMatricies_TrueProbes.mat: [Entropy, Enthalphy, and heat capacity in binding site map format for RNA]
	(GeneName)_RefSeqID_BindingMatricies_TrueProbes.mat: [Entropy, Enthalphy, and heat capacity in binding site map format for complementary strand DNA binding]
	(GeneName)_RefSeqID_Tm37_BasicDesignerStats_TrueProbes.mat:[Index information on stats used for design probes]
	(GeneName)_RefSeqID_chosen.mat: [List of chosen probe indexes]
	(GeneName)_RefSeqID_probes_final_96max.xlsx:[Excel spreadsheet with final probes, and some stats]
	(GeneName)_RefSeqID_Tm37_ModelMetrics_TrueProbes.mat:[Structure with binding affinity calculations and probe design metrics]
