function Seq = ReadChromosomeSeq(settings,Organism,Chr,Start,End)
    %FILE lOCATION of Reference Genome Assembly report with accession numbers of all chromosomes
    %and chromosome patches (edits or unlocalized DNA)
    %Reads Accession Number of chromosome and returns DNA sequence between
    %start and end point.
    %Try for multiple start and end regions
if (strcmp(settings.referenceType,'RefSeq'))
    if (strcmp(Organism,'Human'))
        S8 = readtable(settings.HumanGenomeAssemblyReportFile);
        AN_Chr = S8.RefSeq_Accn(1:639);%%chromosomes 1-22, X (23), Y(24), patchs/scaffolds 25-638, chrM MT (639)
        if (strcmpi(ChrNum,'X'))
            AN_ON_Chr = AN_Chr{23};
        elseif (strcmpi(ChrNum,'Y'))
            AN_ON_Chr = AN_Chr{24};
        elseif (strcmpi(ChrNum,'M'))
            AN_ON_Chr = 'NC_012920.1';
        else
            AN_ON_Chr = AN_Chr{str2double(ChrNum)};
        end
    elseif strcmp(Organism,'Mouse')
        S16 = readtable(settings.MouseGenomeAssemblyReportFile);
        AN_Chr = S16.RefSeq_Accn;%chromosomes 1-19, X (20), Y(21), patche scaffolds 22-60, chrM MT (61)
        if (strcmpi(ChrNum,'X'))
            AN_ON_Chr = AN_Chr{20};
        elseif (strcmpi(ChrNum,'Y'))
            AN_ON_Chr = AN_Chr{21};
        elseif (strcmpi(GeneChr,'MT'))
            AN_ON_Chr = 'NC_005089.1';
        else
            AN_ON_Chr = AN_Chr{str2double(ChrNum)};
        end
    elseif strcmp(Organism,'Yeast')
        S21 = readtable(settings.YeastGenomeAssemblyReportFile);
        AN_Chr = S21.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
        if (strcmpi(GeneChr,'MT'))
            AN_ON_Chr = AN_Chr{17};
        else
            AN_ON_Chr = AN_Chr{roman2num(GeneChr)};
        end
    else
        S22 = readtable(settings.CustomGenomeAssemblyReportFile);
        AN_Chr = S22.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
        if (strcmpi(ChrNum,'X'))
            AN_ON_Chr = AN_Chr{settings.Custom_X_ChromNumber};
        elseif (strcmpi(ChrNum,'Y'))
            AN_ON_Chr = AN_Chr{settings.Custom_Y_ChromNumber};
        elseif (strcmpi(GeneChr,'MT'))
            AN_ON_Chr = AN_Chr{settings.Custom_MT_ChromNumber};
        else
            AN_ON_Chr = AN_Chr{str2double(ChrNum)};
        end
    end
    Seq = getgenbank(AN_Chr{Chr},'PARTIALSEQ',[Start,End]).Sequence;
else
    Seq = [];
end
end