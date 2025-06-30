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
   	- Database annotation files for each organism with Blast Databases, Reference Genome/Transcriptome, Gene Annotation, and Expression Data
    - **Blast Databases:** TrueProbes/data/DatabaseData/Blast_Databases/Organism/
    - **Reference Genome/Transcriptome:** TrueProbes/data/DatabaseData/Blast_Databases/Organism
    - **GTF:** TrueProbes/data/DatabaseData/GTF_Databases/Organism
    - **GFF:** TrueProbes/data/DatabaseData/GFF3_Databases/Organism
    - **Gene Expression Data:** TrueProbes/data/DatabaseData/GeneExpressionData/Organism
   
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
   
## Basic Description of How TrueProbes Works
**The TrueProbes software runs by performing eight steps sequentially.**     
1. **Probe Generation.** It generates all possible probes shared between a list of inclusion IDs and inclusion te xt files, but not in exclusion text files or exclusion IDs, within a set probe length range.   
2. **BLAST Alignment.** All probes target hits in the reference genome and/or transcriptome are identified at least as long as the minimum homology length.    
3. **Binding Affinity Calculation.** Binding affinities are calculated for all pairs of probe and target sequences homology matches in the BLAST results.   
4. **BLAST Target Gene/Transcript Expression**. The gene expression and transcript expression levels are collected for reference expression databases specified in the TrueProbes settings for all targets in the BLAST results.   
5. **Probe-Target Binding Site Mapping**. All blast hits and binding affinities are converted into a site-specific binding map to generate a formatted map by relative binding site position on each target gene, transcript, or chromosome.   
6. **Probe-Target Statistics**. Generate statistics on blast hits, thermodynamics, and a comparison of probes sharing off-targets and relative trade-offs when quantifying off-targets by probe and comparing probes that bind them in a site-specific manner.    
7. **Probe Design**. Sort probes by with/without expression data, the number of off-targets, and then by difference in on-target binding to off-target binding and secondary structure binding affinity to iteratively design probes. Print out the list in an Excel spreadsheet.   
8. **Model Evaluation**. The final probe set and reference expression are combined to compute equilibrium probe binding and statistics, cumulative off-target binding, on-target binding, etc., when given reference values for cell size, probe concentration, and probe intensity.    

## Usage

TrueProbes Probe Design uses four main input files for specifying probe design
1. **TrueProbes_DesignTargets.csv:** CSV File with list of all design targets
2. **ProbeDesignSettings_Parameters.xml:** XML File where all design settings are specified
3. **DatabaseLocations.xml:** XML File with location of all database files needed in design when using NCBI or ENSEMBL reference genome, for any potential organism designed against
4. **GeneExpressionDataFileLocation.xml:** XML File with location of all gene expression files, schema, and sample label for all reference gene expression desired in design

Probe design is run via the command line using A0_BKJH_ProbeDesign_Wrapper_cluster_V5(id,cluster).

The design file needs to be run with two inputs (id and cluster), with a table (inputs1) describing each target gene (inputs1), 
and a set of settings describing how the probe design will be run (settings).

**Main Argument**      
	*id*: id is an integer and is the row of the input design table at the top of the script to run and design probes against.    
	*cluster*: Cluster is an integer and determines the software parallelization pool between local or remote servers when running the script via Slurm. 
		The only difference between running on a cluster is the number of cores in the Slurm file. cluster = 0 if run locally, and 1 if run remotely
	
	*design_file*: An excel spreadsheet file with the input files and parameter settings for the probe design
	
 ## Design Targets in TrueProbes_DesignTargets.csv
 
 Inputs1 Table: inputs1 is the table with the list of targets probes will be designed for. The table row has nine column entries. 
		Input table columns:
			1. Organism to design probes for.
			2. Included target accession IDs. Designs probes shared across all accession numbers 
			3. Excluded target accession IDs. Removes probes in exclusion accession numbers.			
			4. Text Sequence Files to Include (files). Default empty
			5. Text Sequence Files to Exclude (files). Default Empty

	
	Example: {'Human',{NM_000805.5'},{},{},{},;}...   

## Parameter Settings ProbeDesignSettings_Parameters.xml

Settings.
	Below the inputs1 table is a list of primary and 
 
Parameters are specified in TrueProbes_ParameterSettings.xml and are grouped into different categories which can be configured to determine how the probes are designed and evaluated.

**Main Probe Design Parameters:**    
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| minProbeSize | min nt length of potential probes | 20 | 
| maxProbeSize | max nt length of potential probes | 20 |    
| MinProbeSpacing | min spacing between probes | 3 |   
| MaxNumberOfProbes | Max number of probes to design | 96 |    
| targetStrand | which strand to of target to design probes against, 1 for ‘plus’ with 0 for ‘minus’ | 1 |   
| MinHomologySearchTargetSize | Minimum homology length for BLAST alignments to be recorded and used in probe design and evaluation | 15 |   
| BLASTrna | decide to blast RNA sequences | 1 |
| BLASTdna | decide to blast DNA sequences | 0 |   
| ExpressionReferenceForProbeDesign | which row across all expression reference files to use in probe design, with 0 meaning to not use expression data to design probes | 0 |    
	
**Thermodynamic Settings:**    
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| GibbsModel | Which thermodynamic model to use for probe design and evaluation | 4 |      
| HybridizationTemperature | Hybridization temperature in Celcius | 37 |  
| HeatCapacityReferenceTemperature | Reference temperature for Cp measurement and Gibbs model in Celcius | 37 |   
| SaltConcentration | Salt Concentration in M | 0.05 |
| PrimerConcentration | Primer Concentration M | 50e-6 |   
| RemoveMisMatches | Remove sequence mismatched base pairs before evaluating probe-target binding affinity, (1) including adds flanking sequences to alignments and uses Gibbs model 8 with mismatch base pair inclusion | 1 |   

**Gibbs Free Energy Models**
| Model Number  | Abbreviation | Reference |
| ----- | ----- | ----- |
| 1 | Breslauer86 | Breslauer K.J., Frank R., Blocker H., Marky L.A., (1986) Predicting DNA duplex stability from the base sequence Proc Natl Acad Sci U S A 83, 3746-3750 |
| 2 | SantaLucia96 | SantaLucia, J., Allawi, H. T., and Seneviratne, P. A. (1996) Improved nearest-neighbor parameters for predicting DNA duplex stability Biochemistry 35, 3555-3562 |   
| 3 | SantaLucia98 | SantaLucia, J. (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics Proc Natl Acad Sci U S A 95, 1460-1465 |
| 4 | Sugimoto96 | Sugimoto, N., Nakano, S., Yoneyama, M., and Honda, K. (1996) Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes Nucleic Acids Res 24, 4501-4505 |
| 5 | SantaLucia04 | SantaLucia Jr, J., and Hicks, D. (2004) The thermodynamics of DNA structural motifs Annu Rev Biophys Biomol Struct 33, 415-440 |    
| 6 | Allawi97 | Allawi, H. T., and SantaLucia, J. (1997) Thermodynamics and NMR of internal G.T mismatches in DNA Biochemistry 36, 10581-10594 |
| 7 | Rejali21 | Rejali, N. A., Ye, F. D., Zuiter, A. M., Keller, C. C., and Wittwer, C. T. (2021) Nearest-neighbour transition-state analysis for nucleic acid kinetics Nucleic Acids Res 49, 4574-4585 |
| 8 | Martins24 | de Oliveira Martins, E., and Weber, G. (2024) Nearest-neighbour parametrization of DNA single, double and triple mismatches at low sodium concentration Biophys Chem 306, 107156 |

**Design Filtering Settings:**   
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| RemoveProbesBindingOffTargetRibosomalHits | Filter out probes with off-targets to ribosomal proteins | 1 |   
| packOptimal_ProbesWithNoOffTargets | When designing probes without off-target use optimal packing to get as many probes with no off-targets as possible as opposed to normal sequential selection | 1 |   
| IncludeSelfHybridizationInProbeSelection | When designing probes consider probe self-hybridization when ranking probes based on binding affinity | 1 |    

**Parallelization Parameters:**    (Usually only changed when using longer genes)   
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| Parallization_probeBatchSize| number of probes to evaluate in a single batch when performing parallelized calculations | 20 |   
| Parallization_targetBatchSize | number of targets to evaluate in a single batch when performing parallelized calculations | 200 |   
| ParsingPreference | blast simultaneously in parallel(1) or blast probes sequentially (0) | 1 |   
  
**Gene Expression Parameters:**   	
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| UseGeneOverTranscLevelExpression | use gene level (1) or transcript isoform level (0) gene expression values | 0 |   
| DoAllGenesHaveSameExpression | decide to assume equal expression for all genes (1) or to use gene expression reference (0) |  0 |   
| UseRegularDNAExpression | (0) use DNA expression from gene expression track in expression data, (1) set expression to 2 for DNA | 1 |   
| nullRNAcopynumber | [number of RNA copy when not using reference expression levels | 100 |   
| nullDNAcopynumber | [number of DNA copy number when not using reference expression levels | 2 |   
| TMM_LogRatioTrim | when normalizing TPM expression data using TMM set log ratio trim threshold cutoff | 0.3 |   
| TMM_SumTrim | when normalizing TPM expression data using TMM set sum trim threshold cutoff | 0.05 |   
| TMM_Acutoff | when normalizing TPM expression data using TMM set A cutoff value | -1e10 |   
| TMM_doWeighting | when normalizing TPM expression data weight terms using inverse of approximate asymptotic variance of the M-values to account for genes with higher read counts having lower variance on log scale and more reliable mean estimation | 1 |   
	
**Make Blast Database Settings:**    
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| Parse_seqids | when making blastdb in software decision to parse sequence ids from fasta files as blastdb ids, useful for blastdbcmd sequence retrieval | true |   
| Hash_index | when making blastdb make sequence hash indexes leads to faster exact match retrieval but less accurate range matches | false |   

**BLAST Parameters:**     
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| evalue | expectation value cutoff | 1000 |    
| Dust | filter query sequences with DUST | no |   
| Gapextend | Cost to extend a gap (integer) | 2 |   
| Gapopen | cost to open a gap (integer) | 5 |   
| Num_alignments | number of database sequences to show num_alignments for | 1000 |   
| Penalty | penalty for a nucleotide mismatch | -3 |   
| Reward | reward for a nucleotide match | 1 |   
| Word_size | word size for wordfinder algorithm | 7 |   

**Model Simulation Settings:**       
| Name  | Description | Default Value |
| ----- | ----- | ----- |
| removeUndesiredIsoformsFromPrediction | when evaluating main on/off-target binding should alternate isoforms of desired targets be removed and not included in off-target quantification | 1 |  
| ProbeConcentration_MicroMolar | probe concentration in uM | 5e-6 |   
| CellRadius_Micron | cell radius in microns for converting RNA molecule counts into concentration for solving binding equilibrium | 10 |   
| Dilution_Vector | comma seperated vector of probe dilutions to evaluate binding equilibrium and predictions at | 1 |  
| Gibbs_Model_Vector | comma seperated vector of gibbs thermodynamic models to use for evaluating binding equilibrium and predictions at | 1,2,3,4 |    
| Temperature_Celsius_Model_Vector| comma seperated vector of hybridization temperatures in Celsius to use for evaluating binding equilibrium and predictions at | 37 |   
| InitialFreeSolutionGuessConcentration_MicroMolar | initial solution guess for steady state free probe concentration in uM | 1e-10 |    
| SolutionErrorTolerance | total tolerance level for error in final equilibrium solution | 1 |   
| MaxRecursiveEquilibriiumIterations | max number of equilibrium equation iterations stopping at final calculated steady state, overrides solution error tolerance | 40 |   
| CellDiameter_Pixels | cell pixel diameter for background predictions | 50 |   
| SpotRadius_Pixels | spot pixel radius for predictions | 5 |    
| NumberOfReferenceZStacks | number of z-stacks to spread background across for intensity predictions | 67 |   
| SignalStepSize | step size for signal intensity bins in signal predictions | 1e-1 |    
| SignalMaxValue | max intensity value for ranges in solution | 3000 |    
| AutoFluoresenceBackground_MEAN | reference mean autofluorescence | 278 |  
| AutoFluoresenceBackground_STD | reference autofluorescence standard deviation | 33 |   
| NumberOfProbesInReferenceSpots | number of probes in reference spot intensity for calibrating intensity predictions | 48 |   
| ReferenceSpotIntensity_MEAN | mean reference spot intensity | 827 |   
| ReferenceSpotIntensity_STD | reference spot intensity standard deviation |28 |   

## DatabaseLocations of Input Database files are put in DatabaseLocations.xml 

**EMBL_to_NCBI**    
Stores the location of files for mapping ENSEMBL gene and transcript accession numbers to NCBI RefSeq gene and accession numbers     
Each row of EMBL_to_NCBI is row for each organism with    
| Name  | Description |
| ----- | ----- |
| Organism | name of organism to search for in input table |
| StableIDs | location of file with EMBL and NCBI transcripts paired to one another for that organism |

**EMBL**     
Stores the location of files for using ENSEMBL annotation in probe design    
Each row of EMBL is row for each organism with     
| Name  | Description |
| ----- | ----- |
| Organism | name of organism to search for in input table |
| Root_FASTA | Location of folder with all ENSEMBL annotation dna and rna fasta files |
| BLASTDB_RNA | Location of ensembl blast genome DNA database files |
| BLASTDB_RNA | Location of ensembl blast transcript RNA database files | 
| GTF | Location of ensembl reference genome gtf file |
| GFF | Location of ensembl reference genome gff file |

**NCBI**     
Stores location of files for using NCBI RefSeq annotation in probe design     
Each row of NCBI is row for each organism with     
| Name  | Description |
| ----- | ----- |
| Organism |name of organism to search for in input table | 
| Root_FASTA | Location of folder with all refseq annotation dna and rna fasta files | 
| BLASTDB_RNA | Location of refseq blast genome DNA database files |
| BLASTDB_RNA | Location of refseq blast transcript RNA database files | 
| GTF | Location of refseq reference genome gtf file | 
| GFF | Location of refseq reference genome gff file |

## Gene Expression Data File Locations are in GeneExpressionDataLocations.xml

**Gene Expression File Locations**   
Stored for each organism location of all expression files with identifier name   
Row for each organism's reference gene or transcript level expression file   
| Name  | Description |
| ----- | ----- |
| Organism | name of organism to search for in input table |
| “Data Identifier Name” | Location of gene expression file with extension in TrueProbes Folder | 

**Gene Expression File Schema**     
Stored for each organism location of all expression schema files listing different cell sample names associated with each expression file    
Schema for each organism's reference gene or transcript level expression file   
| Name  | Description |
| ----- | ----- |
| Organism | name of organism to search for in input table |
| “Data Identifier Name” | Location of schema for gene expression file with extension |

**Gene Expression File Column Names**         
Stored also is the list of columns for each file output extension to use when finding expression data, and gene or target ids to map it to blast results    
tracks for output file types list of column information when reading the expression data    
| Name  | Description |
| ----- | ----- |
| “file extension” | column sorted list of variable names |

## TrueProbes Design Output Files
**TrueProbe Design Software Output Files:**       
| Name  | Description |
| ----- | ----- | 
| (GeneName)_AccessionID_probes_TrueProbes.mat | structure with probe sequences, location on on-target | 
| (GeneName)_AccessionID_hits_table_TrueProbes.mat | structure with information on BLAST hits |   
| (GeneName)_AccessionID_ExpressionInfo_TrueProbes.mat | Structure with expression data |   
| (GeneName)_AccessionID_dKbInfo_TrueProbes.mat | Structure with heat capacity for all target binding reactions |
| (GeneName)_AccessionID_dCpInfo_TrueProbes.mat | Structure with heat capacity for all target binding reactions |  
| (GeneName)_AccessionID_dHInfo_TrueProbes.mat | Structure with enthalpy for all target bindings reactions |   
| (GeneName)_AccessionID_dSInfo_TrueProbes.mat | Structure with entropy for all target binding reactions |   
| (GeneName)_AccessionID_TmInfo_TrueProbes.mat | Structure with entropy for all target binding reactions |  
| (GeneName)_AccessionID_TmHybridizationTemperature_OnOffThermoInfo_TrueProbes.mat | Structure with binding energy of all hits |
| (GeneName)_binding_hits_map_TrueProbes | binding site map |
| (GeneName)_AccessionID_BindingMatricies_TrueProbes.mat | Entropy, Enthalphy, and heat capacity in binding site map format for RNA |
| (GeneName)_AccessionID_TmHybridizationTemperature__BindingEnergyMatrix_TrueProbes.mat | Equilibrium Binding Energy in binding site map format |
| (GeneName)_AccessionID_BindingMatricies2_TrueProbes.mat | Entropy, Enthalphy, and heat capacity in binding site map format for complementary DNA strand binding |   
| (GeneName)_AccessionID_TmHybridizationTemperature__BindingEnergyMatrix2_TrueProbes.mat | Complementary DNA strand Equilibrium Binding Energy in binding site map format |
| (GeneName)_AccessionID_TmHybridizationTemperature_BasicDesignerStats_TrueProbes.mat | Index information on stats used for design probes |
| (GeneName)_AccessionID_chosen.mat | List of chosen probe indexes |
| (GeneName)_AccessionID_probes_final_MaxProbeNumber_max.xlsx | Excel spreadsheet with final probes, and some stats |
| (GeneName)_AccessionID_TmHybridizationTemperature_ModelMetrics_TrueProbes.mat | Structure with binding affinity calculations and probe design metrics |
