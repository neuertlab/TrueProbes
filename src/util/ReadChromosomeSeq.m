function Seq = ReadChromosomeSeq(Organism,Chr,Start,End)
    %FILE lOCATION of Reference Genome Assembly report with accession numbers of all chromosomes
    %and chromosome patches (edits or unlocalized DNA)
    %Reads Accession Number of chromosome and returns DNA sequence between
    %start and end point.
    
    %Try for multiple start and end regions
if (strcmp(Organism,'Human'))
    S8 = readtable('Metadata/Human/GCF_000001405.39_GRCh38.p13_assembly_report.txt');
    AN_DNA = S8.RefSeq_Accn(1:639);%%chromosomes 1-22, X (23), Y(24), patchs/scaffolds 25-638, chrM MT (639)
    Seq = getgenbank(AN_DNA{Chr},'PARTIALSEQ',[Start,End]).Sequence;
elseif strcmp(Organism,'Mouse')  
    S16 = readtable('Metadata/Mouse/GCA_000001635.9_GRCm39_assembly_report.txt'); 
    AN_DNA = S16.RefSeq_Accn;%chromosomes 1-19, X (20), Y(21), patche scaffolds 22-60, chrM MT (61)
    Seq = getgenbank(AN_DNA{Chr},'PARTIALSEQ',[Start,End]).Sequence;
elseif strcmp(Organism,'Yeast')
    S21 = readtable('Metadata/Yeast/S288C/GCF_000146045.2_R64_assembly_report.txt');
    AN_DNA = S21.RefSeq_Accn;%chromosome 1-16, chromosome M (Mitocondria) (17)
    Seq = getgenbank(AN_DNA{Chr},'PARTIALSEQ',[Start,End]).Sequence;
end
end