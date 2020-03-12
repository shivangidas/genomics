%Question 6
mito_gbk = getgenbank('NC_001643'); %chimp
mitochondria = mito_gbk.Sequence;
mitochondria_length = length(mitochondria)
first_300_bases = seqdisp(mitochondria(1:300))

figure
ntdensity(mitochondria)

bases = basecount(mitochondria)

compBases = basecount(seqrcomplement(mitochondria))

figure
basecount(mitochondria,'chart','pie');
title('Distribution of Nucleotide Bases for Chimpanzee Mitochondrial Genome');

%question 7
chimp = getgenbank('NC_001643','SequenceOnly',true);
whos chimp

codoncount(chimp)
orf = seqshoworfs(chimp,'MINIMUMLENGTH',1, 'geneticcode',2,'frames','all');
orf1=seqshoworfs(chimp(randperm(length(chimp))),'MINIMUMLENGTH',1,'geneticcode',2,'frames','all');
ORFLength=[];
for i=1:6
   for j=1:length(orf(i).Stop)
    ORFLength=[ORFLength; orf(i).Stop(j) - orf(i).Start(j)];
   end
end
length(ORFLength)

ORFLength1=[];
for i=1:6
   for j=1:length(orf1(i).Stop)
    ORFLength1=[ORFLength1; orf1(i).Stop(j) - orf1(i).Start(j)];
   end
end
length(ORFLength1)

max_threshold=max(ORFLength1)
n_max=length(find(ORFLength>=max_threshold))

orfNew=[];
for i=1:6
    for j=1:length(orf_threshold(i).Stop)
        orfNew=[orfNew; orf_threshold(i).Stop(j) - orf_threshold(i).Start(j)];
    end
end
length(orfNew)

protein = chimp(orf_threshold(1).Start(1):orf_threshold(1).Stop(1)+2)

