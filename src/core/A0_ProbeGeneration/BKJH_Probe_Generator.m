function [probes,tent_probes,seqs,gene_names,organisms] = BKJH_Probe_Generator(L,AccessionNumbers,TextInclude,TextExclude,AccessionExclude,strand,isOffline,SEQdbRoot)

%%% This function is meant to take acession numbers as input, use the matlab
%%% function to retrieve the gene sequences, and then return the designed
%%% probes that are common to all accession numbers as a cell array with
%%% strings for each probe. The first column is the index of the probe, and
%%% the second is the probe itself
%Xist Accession Numbers: NR_001463.3, NR_001570.2
% AccessionNumbers = {'NR_001463.3','NR_001570.2'}; %One of Xist accession numbers
%% Retrieve information from Accession Numbers
clear tent_probes seqs organisms gene_names all_tent_probes probes
probes = {};
tent_probes = {};   %Will store all initial probes determined from the sequences. 
all_tent_probes = {}; %Will store all initial probes determined from the sequences in one row
seqs = {};  %Will store the sequences for all the accession numbers
gene_names = {}; %Will store the gene name for all accession numbers
counter = 1;
probePos = {};
Ll = L(1);
Lh = L(2);
curr =0;
RV = @(x) (seqrcomplement(x));
if (isOffline)
   FNAFiles_FNA = dir([SEQdbRoot '**/*.fna']);
   FNAFiles_FA = dir([SEQdbRoot '**/*.fa']);
   if (~isempty(FNAFiles_FNA)&&~isempty(FNAFiles_FA))
        FNAFiles = cat(1,FNAFiles_FNA,FNAFiles_FA);
   elseif (~isempty(FNAFiles_FNA))
        FNAFiles = FNAFiles_FNA;
   elseif (~isempty(FNAFiles_FA))
        FNAFiles = FNAFiles_FA;
   end
   FNAFileLoc = cell(1,length(FNAFiles));
   BioIFobj = cell(1,length(FNAFiles));
   BioKeys = cell(1,length(FNAFiles));
   for k=1:length(FNAFiles)
        FNAFileLoc{k} = strcat(FNAFiles(k).folder,filesep,FNAFiles(k).name);
        BioIFobj{k} = BioIndexedFile('FASTA',FNAFileLoc{k});
        BioKeys{k}{:} = getKeys(BioIFobj{k});
   end        
end
organisms = cell(size(AccessionNumbers,2),1);
seqs = cell(size(AccessionNumbers,2),1);
gene_names = cell(size(AccessionNumbers,2),1);
if (strand==1)
for i=1:Lh-Ll+1
    Lp = Ll+i-2;
    if not(isempty(AccessionNumbers))
        for Acc_Num = 1:size(AccessionNumbers,2)
            if (isOffline)
                isInBioIFobj = zeros(1,length(BioIFobj));
                for k=1:length(BioIFobj)
                    isInBioIFobj(k) = sum(contains(BioKeys{k}{:},AccessionNumbers{Acc_Num}));
                end
                cObj = find(isInBioIFobj>0,1); 
                if (~isempty(cObj))
                    Data = read(BioIFobj{cObj},AccessionNumbers{Acc_Num});
                    clear isInBioIFobj cObj
                else
                    fprintf('\n')
                    fprintf("BioKey not found in reference genomes or transcriptomes.")
                    fprintf('\n')
                    fprintf(strcat("Trying to search for matches to ",extractBefore(AccessionNumbers{Acc_Num},'.')))
                    fprintf('\n')
                    %BioKey version not found
                      isInBioIFobj = zeros(1,length(BioIFobj));
                      isInBioIFobj_loc = zeros(1,length(BioIFobj));
                      for k=1:length(BioIFobj)
                            with_dot = find(contains(BioKeys{k}{:},'.'));
                            if (sum(strcmp(extractBefore(BioKeys{k}{1}(with_dot),'.'),extractBefore(AccessionNumbers{Acc_Num},'.')))>0)
                            isInBioIFobj_loc(k) = with_dot(find(strcmp(extractBefore(BioKeys{k}{1}(with_dot),'.'),extractBefore(AccessionNumbers{Acc_Num},'.'))));
                            isInBioIFobj(k) = 1;
                            end
                      end
                       cObj = find(isInBioIFobj>0,1); 
                       if (~isempty(cObj))
                             Data = read(BioIFobj{cObj},BioKeys{cObj}{1}{isInBioIFobj_loc(cObj)});
                             clear isInBioIFobj cObj isInBioIFobj_loc
                       else
                           msg = 'Error. Target Inclusion Accession Number not in reference genome or transcriptome';
                           error(msg)
                       end      
                end
            else
                Data = getgenbank(AccessionNumbers{Acc_Num});    %Retrieve accession number information
                pause(3);
            end
            
            temp_seq = Data.Sequence;
            seqs{Acc_Num,1} = Data.Sequence;
            if (isOffline)
                organisms{Acc_Num,1} = Data.Header;  
            else
                organisms{Acc_Num,1} = Data.Source;  
            end
            for base1 = 1:size(temp_seq,2)-Lp
                tent_probes{base1,Acc_Num} = temp_seq(base1:base1+Lp);
                all_tent_probes{1,counter} = temp_seq(base1:base1+Lp);
                ProbePos(counter) = base1;
                counter = counter+1; 
            end
            if (isOffline)
            gene_names{Acc_Num,1} = [];
            else  
            gene_names{Acc_Num,1} = Data.Definition;    %Store the gene name
            end
        end
    end

%% Add sequences from text files if applicable
    if not(isempty(TextInclude))
        for Txt_Num = 1:size(TextInclude,2)
        text1 = fileread(TextInclude{Txt_Num});
        temp_seq = text1;
        seqs{Txt_Num,1} = temp_seq;
        organisms{Txt_Num,1} = 'TextFile';
            for base1 = 1:size(temp_seq,2)-Lp
                tent_probes{base1,Txt_Num} = temp_seq(base1:base1+Lp);
                all_tent_probes{1,counter} = temp_seq(base1:base1+Lp);
                ProbePos(counter) = base1;
                counter = counter+1;
            end
            gene_names{Txt_Num,1} = 'TextFile';    %Store the gene name
        end
    end
end

%% Only pick probes common to all Sequences

    [C1, ia1, ic1] = unique(all_tent_probes);
% size(all_tent_probes)
% size(ic1)
    counter = 1; %Will be used for storing probes
    for num1 = 1:max(ic1)
        if sum(ic1 == num1) == size(AccessionNumbers,2)+size(TextInclude,2)
            index1 = find(ic1 == num1,1,'first');
            probes{counter,1} = index1; %Store the location of the original probe
%         probes{counter,2} = C1{num1}; %Store the original probe
            probes{counter,2} = all_tent_probes{index1};
            probes{counter,3} = ProbePos(index1);
            counter = counter+1;
        end
    end
    probes = sortrows(probes,1);    %Sort by the index in the original sequence

%% Remove Probes matching the input sequence  (such as to remove exons from intronic probes)
%%% Add probes from sequence
    if not(isempty(TextExclude))
        for exclude1 = 1:size(TextExclude,2)
        tentprobes1 = probes;  %will have the tentative probes plus probes to avoid
        clear probes
        probes = {};
        counter = 1;
        text1 = fileread(TextExclude{exclude1});
            for num1 = 1:size(tentprobes1,1)
                if isempty(strfind(text1,tentprobes1{num1,2}))
                    probes{counter,1} = tentprobes1{num1,1}; %Store the location of the original probe
            %         probes{counter,2} = C1{num1}; %Store the original probe
                    probes{counter,2} = tentprobes1{num1,2};
                    probes{counter,3} = tentprobes1{num1,3};
                    counter = counter+1;
                end
            end
        end
    end
%% Remove probes from Accession Numbers to Exclude
    if not(isempty(AccessionExclude))
        for exclude1 = 1:size(AccessionExclude,2)
        tentprobes1 = probes;  %will have the tentative probes plus probes to avoid
        clear probes
        probes = {};
        counter = 1;
            if (isOffline)
                for k=1:length(BioIFobj)
                    isInBioIFobj(k) = sum(contains(BioKeys{k}{:},AccessionExclude{exclude1}));
                end
                cObj = find(isInBioIFobj>0,1); 
                if (~isempty(cObj))
                    Data = read(BioIFobj{cObj},AccessionExclude{exclude1});
                    clear isInBioIFobj cObj
                else
                    clear isInBioIFobj cObj
                    continue;   
                end
            else
                Data = getgenbank(AccessionExclude{exclude1});    %Retrieve accession number information
                pause(3);
            end
        
        text1 = Data.Sequence;
            for num1 = 1:size(tentprobes1,1)
                if isempty(strfind(tentprobes1{num1,2},text1))
                    probes{counter,1} = tentprobes1{num1,1}; %Store the location of the original probe
            %         probes{counter,2} = C1{num1}; %Store the original probe
                    probes{counter,2} = tentprobes1{num1,2};
                    probes{counter,3} = tentprobes1{num1,3};
                    counter = counter+1;
                end
            end
        end
    end
    
else
for i=1:Lh-Ll+1
    Lp = Ll+i-2;
    if not(isempty(AccessionNumbers))
        for Acc_Num = 1:size(AccessionNumbers,2)
            if (isOffline)
                for k=1:length(BioIFobj)
                    isInBioIFobj(k) = sum(contains(BioKeys{k}{:},AccessionNumbers{Acc_Num}));
                end
                cObj = find(isInBioIFobj>0,1); 
                if (~isempty(cObj))
                    Data = read(BioIFobj{cObj},AccessionNumbers{Acc_Num});
                    clear isInBioIFobj cObj
                else
                    clear isInBioIFobj cObj
                    continue;   
                end
            else
                Data = getgenbank(AccessionNumbers{Acc_Num});    %Retrieve accession number information
                pause(3);
            end
            temp_seq = RV(Data.Sequence);
            seqs{Acc_Num,1} = RV(Data.Sequence);
            if (isOffline)
                organisms{Acc_Num,1} = Data.Header;  
            else
                organisms{Acc_Num,1} = Data.Source;
            end
            for base1 = 1:size(temp_seq,2)-Lp
                tent_probes{base1,Acc_Num} = temp_seq(base1:base1+Lp);
                all_tent_probes{1,counter} = temp_seq(base1:base1+Lp);
                ProbePos(counter) = base1;
                counter = counter+1; 
            end
            if (isOffline)
            gene_names{Acc_Num,1} = [];
            else  
            gene_names{Acc_Num,1} = Data.Definition;    %Store the gene name
            end
        end
    end

%% Add sequences from text files if applicable
    if not(isempty(TextInclude))
        for Txt_Num = 1:size(TextInclude,2)
        text1 = fileread(TextInclude{Txt_Num});
        temp_seq = RV(text1);
        seqs{Txt_Num,1} = RV(temp_seq);
        organisms{Txt_Num,1} = 'TextFile';
            for base1 = 1:size(temp_seq,2)-Lp
                tent_probes{base1,Txt_Num} = temp_seq(base1:base1+Lp);
                all_tent_probes{1,counter} = temp_seq(base1:base1+Lp);
                ProbePos(counter) = base1;
                counter = counter+1;
            end
            gene_names{Txt_Num,1} = 'TextFile';    %Store the gene name
        end
    end
end

%% Only pick probes common to all Sequences

    [C1, ia1, ic1] = unique(all_tent_probes);
% size(all_tent_probes)
% size(ic1)
    counter = 1; %Will be used for storing probes
    for num1 = 1:max(ic1)
        if sum(ic1 == num1) == size(AccessionNumbers,2)+size(TextInclude,2)
            index1 = find(ic1 == num1,1,'first');
            probes{counter,1} = index1; %Store the location of the original probe
%         probes{counter,2} = C1{num1}; %Store the original probe
            probes{counter,2} = all_tent_probes{index1};
            probes{counter,3} = ProbePos(index1);
            counter = counter+1;
        end
    end
    probes = sortrows(probes,1);    %Sort by the index in the original sequence

%% Remove Probes matching the input sequence  (such as to remove exons from intronic probes)
%%% Add probes from sequence
    if not(isempty(TextExclude))
        for exclude1 = 1:size(TextExclude,2)
        tentprobes1 = probes;  %will have the tentative probes plus probes to avoid
        clear probes
        probes = {};
        counter = 1;
        text1 = fileread(TextExclude{exclude1});
            for num1 = 1:size(tentprobes1,1)
                if isempty(strfind(RV(text1),tentprobes1{num1,2}))
                    probes{counter,1} = tentprobes1{num1,1}; %Store the location of the original probe
            %         probes{counter,2} = C1{num1}; %Store the original probe
                    probes{counter,2} = tentprobes1{num1,2};
                    probes{counter,3} = tentprobes1{num1,3};
                    counter = counter+1;
                end
            end
        end
    end
%% Remove probes from Accession Numbers to Exclude
    if not(isempty(AccessionExclude))
        for exclude1 = 1:size(AccessionExclude,2)
        tentprobes1 = probes;  %will have the tentative probes plus probes to avoid
        clear probes
        probes = {};
        counter = 1;
            if (isOffline)
                for k=1:length(BioIFobj)
                    isInBioIFobj(k) = sum(contains(BioKeys{k}{:},AccessionExclude{exclude1}));
                end
                cObj = find(isInBioIFobj>0,1); 
                if (~isempty(cObj))
                    Data = read(BioIFobj{cObj},AccessionExclude{exclude1});
                    clear isInBioIFobj cObj
                else
                    clear isInBioIFobj cObj
                    continue;   
                end
            else
                Data = getgenbank(AccessionExclude{exclude1});    %Retrieve accession number information
                pause(3);
            end
        text1 = RV(Data.Sequence);
            for num1 = 1:size(tentprobes1,1)
                if isempty(strfind(tentprobes1{num1,2},text1))
                    probes{counter,1} = tentprobes1{num1,1}; %Store the location of the original probe
            %         probes{counter,2} = C1{num1}; %Store the original probe
                    probes{counter,2} = tentprobes1{num1,2};
                    probes{counter,3} = tentprobes1{num1,3};
                    counter = counter+1;
                end
            end
        end
    end       
end
end
