core_folders = genpath(strcat(pwd,filesep,'src',filesep,'core'));
modified_matlab_function_folders = genpath(strcat(pwd,filesep,'src',filesep,'modified_matlab_functions'));
util_folders = genpath(strcat(pwd,filesep,'src',filesep,'util'));
wrappers_folders = genpath(strcat(pwd,filesep,'src',filesep,'wrappers'));
mtp_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'modified_thirdparty'));
nd_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'ndSparse_G4_2025_01_19'));
multiprod_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'Multiprod_2009'));
xsum_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'XSum_2014_06_16'));
xml_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'xml2struct-master'));
trueSpot_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'TrueSpot-main'));
progbar_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'parwaitbar'));
VarGibbs_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'VarGibbs-4.1'));
hpf_folders = genpath(strcat(pwd,filesep,'src',filesep,'thirdparty',filesep,'HighPrecisionFloat'));
if ~(ismcc || isdeployed)
addpath(core_folders);
addpath(modified_matlab_function_folders);
addpath(util_folders);
addpath(wrappers_folders);
addpath(nd_folders);
addpath(mtp_folders);
addpath(multiprod_folders);
addpath(xsum_folders);
addpath(xml_folders);
addpath(trueSpot_folders);
addpath(progbar_folders);
addpath(VarGibbs_folders);
addpath(hpf_folders);
end
if (ispc)
path_info = strsplit(path,';');
end
if (ismac || isunix)
path_info = strsplit(path,':');
end
paths_add = path_info(contains(path_info,pwd));
src_code_scripts = CATnWrapper(cellfun(@(x) dir([x filesep '**/*.m']),paths_add,'Un',0),1);
scripts_locations = extractAfter(join([convertCharsToStrings({src_code_scripts.folder})' convertCharsToStrings({src_code_scripts.name})'],filesep),strcat(pwd,filesep))';
scripts_locations = unique(scripts_locations);
scripts_locations = scripts_locations(~arrayfun(@(x) double(sum(contains(codeIssues(scripts_locations(x)).Issues.Description,'Parse error'))>0),1:length(scripts_locations)));
scripts_locations = scripts_locations(~arrayfun(@(x) double(sum(contains(codeIssues(scripts_locations(x)).Issues.Description,'syntax error(s)'))>0),1:length(scripts_locations)));
scripts_locations = scripts_locations(~arrayfun(@(x) double(sum(contains(string(codeIssues(scripts_locations(x)).Issues.Severity),'error'))>0),1:length(scripts_locations)));


ProbeDesign_Software = 'A0_BKJH_ProbeDesign_Wrapper_cluster_V5.m';
appFile_ProbeDesign = fullfile(which(ProbeDesign_Software));
opts_ProbeDesign = compiler.build.StandaloneApplicationOptions(appFile_ProbeDesign);
outdir_ProbeDesign = opts_ProbeDesign.OutputDir;
outdir_ProbeDesign = extractBefore(outdir_ProbeDesign,'standalone');
if (ismac)
    outputDirectory_ProbeDesign = strcat(outdir_ProbeDesign,'_MacOS');
elseif (isunix)
    outputDirectory_ProbeDesign = strcat(outdir_ProbeDesign,'_Linux');
elseif (ispc)
    outputDirectory_ProbeDesign = strcat(outdir_ProbeDesign,'_Win');
end
opts_ProbeDesign.TreatInputsAsNumeric = 1;
opts_ProbeDesign.OutputDir = outputDirectory_ProbeDesign;
opts_ProbeDesign.AdditionalFiles = scripts_locations;
opts_ProbeDesign.Verbose = 'on';
opts_ProbeDesign.EmbedArchive = 'on';
ProbeDesign_builder = compiler.internal.build.builder.StandaloneApplication(opts_ProbeDesign);
build_info_ProbeDesign = ProbeDesign_builder.build;
ProbeDesign_input0 = ProbeDesign_builder.MccInfo.MccInput;
ProbeDesign_inputs1 = ['-R' '-softwareopengl' '-R' '-logfile,default_ProbeDesign_log' ProbeDesign_input0(1:6) '-A' 'all' ProbeDesign_input0(7:end)];
mcc(ProbeDesign_inputs1{:})
installer_opts_ProbeDesign = compiler.package.InstallerOptions(build_info_ProbeDesign);
installer_opts_ProbeDesign.InstallerName = strcat('TrueProbes_ProbeDesignInstaller_',computer('arch'));
installer_opts_ProbeDesign.OutputDir = outputDirectory_ProbeDesign;
compiler.package.installer(build_info_ProbeDesign,'Options',installer_opts_ProbeDesign);
% try
% if (isunix)
% docker_opts_ProbeDesign = compiler.package.DockerOptions(build_info_ProbeDesign,'ImageName','A0_BKJH_ProbeDesign_Wrapper_cluster_V5');
% compiler.package.docker(build_info_ProbeDesign,'Options',docker_opts_ProbeDesign)
% end
% compiler.runtime.download
% if (~isMATLABReleaseOlderThan("R2023b"))
% compiler.runtime.createDockerImage(build_info_ProbeDesign,'ImageName','A0_BKJH_ProbeDesign_Wrapper_cluster_V5');
% end
% catch
% end

ProbeComparison_Software = 'A0_BKJH_ProbeComparison_Wrapper_cluster_V5.m';
appFile_ProbeComparison = fullfile(which(ProbeComparison_Software));
opts_ProbeComparison = compiler.build.StandaloneApplicationOptions(appFile_ProbeComparison);
outdir_ProbeComparison = opts_ProbeComparison.OutputDir;
outdir_ProbeComparison = extractBefore(outdir_ProbeComparison,'standalone');
if (ismac)
    outputDirectory_ProbeComparison = strcat(outdir_ProbeComparison,'_MacOS');
elseif (isunix)
    outputDirectory_ProbeComparison = strcat(outdir_ProbeComparison,'_Linux');
elseif (ispc)
    outputDirectory_ProbeComparison = strcat(outdir_ProbeComparison,'_Win');
end
opts_ProbeComparison.TreatInputsAsNumeric = 1;
opts_ProbeComparison.OutputDir = outputDirectory_ProbeComparison;
opts_ProbeComparison.AdditionalFiles = scripts_locations;
opts_ProbeComparison.Verbose = 'on';
opts_ProbeComparison.EmbedArchive = 'on';
ProbeComparison_builder = compiler.internal.build.builder.StandaloneApplication(opts_ProbeComparison);
build_info_ProbeComparison = ProbeComparison_builder.build;
ProbeComparison_input0 = ProbeComparison_builder.MccInfo.MccInput;
ProbeComparison_inputs1 = ['-R' '-softwareopengl' '-R' '-logfile,default_ProbeDesign_log' ProbeComparison_input0(1:6) '-A' 'all' ProbeComparison_input0(7:end)];
mcc(ProbeComparison_inputs1{:})
installer_opts_ProbeComparison = compiler.package.InstallerOptions(build_info_ProbeComparison);
installer_opts_ProbeComparison.InstallerName = strcat('TrueProbes_ProbeComparisonInstaller_',computer('arch'));
installer_opts_ProbeComparison.OutputDir = outputDirectory_ProbeComparison;
compiler.package.installer(build_info_ProbeComparison,'Options',installer_opts_ProbeComparison);