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
% if (~isMATLABReleaseOlderThan("R2024a"))
% compiler.runtime.customInstaller(strcat('TrueProbes_ProbeDesignMinimalMCRInstaller_',computer('arch')),...
%     build_info_ProbeDesign,'OutputDir',installer_opts_ProbeDesign.OutputDir);%custom installer
% % compiler.runtime.customInstaller("matrixInstaller",[results1,results2],...
% % OutputDir="customInstallers",RuntimeDelivery="installer")
% end
if (isunix && ~ismac)
    [status, msg] = system('docker version');
    disp(msg);
    if ~status
        see_docker_images = 'docker images';
        [~, msg] = system(see_docker_images);
        if ~contains(msg,'ncbi/blast')
            get_ncbi_blast_image = 'docker pull ncbi/blast';
            [status_blast, msg] = system(get_ncbi_blast_image);
            disp(msg)
        else
            status_blast = 0;
        end
        if ~status_blast
            docker_blast_info = 'docker inspect ncbi/blast:latest';
            [status_blast_info, msg_blast_info] = system(docker_blast_info);
            blast_info = jsondecode(msg_blast_info);
            blast_env = blast_info.ContainerConfig.Env;
            blast_env_copies = cellfun(@(x) ['ENV ' x],blast_env,'Un',0);
            docker_blast_dependencies = 'docker sbom ncbi/blast:latest';
            [status_blast_dep, msg_blast_dep] = system(docker_blast_dependencies);
            blast_info = jsondecode(msg_blast_info);
            if (~isMATLABReleaseOlderThan("R2022b"))
                docker_MCR_baseline = compiler.runtime.createDockerImage(build_info_ProbeDesign,...
                    'DockerContext',strcat(pwd,filesep,'A0_BKJH_ProbeDesign_Wrapper_cluster_V5_Docker'),...
                    'ExecuteDockerBuild','on');
            else
                compiler.runtime.download
            end
            docker_opts_ProbeDesign = compiler.package.DockerOptions(build_info_ProbeDesign,'ImageName','a0_bkjh_probedesign_wrapper_cluster_v5');
            docker_opts_ProbeDesign.DockerContext = strcat(pwd,filesep,'A0_BKJH_ProbeDesign_Wrapper_cluster_V5_Docker');
            if (~isMATLABReleaseOlderThan("R2022b"))
            docker_opts_ProbeDesign.RuntimeImage = docker_MCR_baseline;
            end
            docker_opts_ProbeDesign.AdditionalInstructions = ...
                [{'RUN mkdir -p /data/databaseData/Blast_Databases'},...
                {'RUN mkdir -p /data/databaseData/ENSEMBL_NCBI_StableIDs'},...
                {'RUN mkdir -p /data/databaseData/Gene_Expression_Data'},...
                {'RUN mkdir -p /data/databaseData/GFF3_Databases'},...
                {'RUN mkdir -p /data/databaseData/GTF_Databases'},...
                {'COPY --from=ncbi/blast:latest /blast /'},...
                blast_env_copies(:)'];    
            compiler.package.docker(build_info_ProbeDesign,'Options',docker_opts_ProbeDesign)
        end
    end      
end
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
ProbeComparison_inputs1 = ['-R' '-softwareopengl' '-R' '-logfile,default_ProbeComparison_log' ProbeComparison_input0(1:6) '-A' 'all' ProbeComparison_input0(7:end)];
mcc(ProbeComparison_inputs1{:})
installer_opts_ProbeComparison = compiler.package.InstallerOptions(build_info_ProbeComparison);
installer_opts_ProbeComparison.InstallerName = strcat('TrueProbes_ProbeComparisonInstaller_',computer('arch'));
installer_opts_ProbeComparison.OutputDir = outputDirectory_ProbeComparison;
compiler.package.installer(build_info_ProbeComparison,'Options',installer_opts_ProbeComparison);
if (isunix && ~ismac)
    [status, msg] = system('docker version');
    disp(msg);
    if ~status
        see_docker_images = 'docker images';
        [~, msg] = system(see_docker_images);
        if ~contains(msg,'ncbi/blast')
            get_ncbi_blast_image = 'docker pull ncbi/blast';
            [status_blast, msg] = system(get_ncbi_blast_image);
            disp(msg)
        else
            status_blast = 0;
        end
        if ~status_blast
            docker_blast_info = 'docker inspect ncbi/blast:latest';
            [status_blast_info, msg_blast_info] = system(docker_blast_info);
            blast_info = jsondecode(msg_blast_info);
            blast_env = blast_info.ContainerConfig.Env;
            blast_env_copies = cellfun(@(x) ['ENV ' x],blast_env,'Un',0);
            docker_blast_dependencies = 'docker sbom ncbi/blast:latest';
            [status_blast_dep, msg_blast_dep] = system(docker_blast_dependencies);
            blast_info = jsondecode(msg_blast_info);
            if (~isMATLABReleaseOlderThan("R2022b"))
                docker_MCR_baseline = compiler.runtime.createDockerImage(build_info_ProbeComparison,...
                    'DockerContext',strcat(pwd,filesep,'A0_BKJH_ProbeComparison_Wrapper_cluster_V5_Docker'),...
                    'ExecuteDockerBuild','on');
            else
                compiler.runtime.download
            end
            docker_opts_ProbeComparison = compiler.package.DockerOptions(build_info_ProbeComparison,'ImageName','a0_bkjh_probecomparison_wrapper_cluster_v5');
            docker_opts_ProbeComparison.DockerContext = strcat(pwd,filesep,'A0_BKJH_ProbeComparison_Wrapper_cluster_V5_Docker');
            if (~isMATLABReleaseOlderThan("R2022b"))
            docker_opts_ProbeComparison.RuntimeImage = docker_MCR_baseline;
            end
            docker_opts_ProbeComparison.AdditionalInstructions = ...
                [{'RUN mkdir -p /data/databaseData/Blast_Databases'},...
                {'RUN mkdir -p /data/databaseData/ENSEMBL_NCBI_StableIDs'},...
                {'RUN mkdir -p /data/databaseData/Gene_Expression_Data'},...
                {'RUN mkdir -p /data/databaseData/GFF3_Databases'},...
                {'RUN mkdir -p /data/databaseData/GTF_Databases'},...
                {'COPY --from=ncbi/blast:latest /blast /'},...
                blast_env_copies(:)'];    
            compiler.package.docker(build_info_ProbeComparison,'Options',docker_opts_ProbeComparison)
        end
    end      
end


% imageToPull = 'mathworks/matlab:r2025a'; % Example: MATLAB R2025a image
%     command = ['docker pull ', imageToPull:2.16.0];
% 
%     [status, cmdout] = system(command);
% 
%     if status == 0
%         disp(['Successfully pulled Docker image: ', imageToPull]);
%         disp(cmdout); % Display output from the docker pull command
%     else
%         disp(['Error pulling Docker image: ', imageToPull]);
%         disp(cmdout); % Display error message
%     end
% end
