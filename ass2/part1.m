%% get the human cytb and cox1 genes for comparision
human_cytb = getgenpept('AFI52643','sequenceonly','true');
human_cox1 = getgenpept('AGW78696','sequenceonly','true');
%% Get the raw Pan troglodytes mitochondrion dna
chimp = getgenbank('NC_001643');
chimpSeq = chimp.Sequence;
%% get ORFs
orf = seqshoworfs(chimpSeq,'MINIMUMLENGTH',3, 'geneticcode',2,'frames','all','nodisplay','true');
%% get ORFs in a similar random permutaion
orf1 = seqshoworfs(chimpSeq(randperm(length(chimpSeq))),'MINIMUMLENGTH',3,'geneticcode',2,'frames','all','nodisplay','true');
%%
% collect all ORFs
ORFLength=[];
for i=1:6
   for j=1:length(orf(i).Stop)
    ORFLength=[ORFLength; orf(i).Stop(j)+2 - orf(i).Start(j)];
   end
end
disp("Number of ORFs in chimpanzee mitochondrion")
length(ORFLength)
%%
% thresholding using random permutaion sequence's max length ORF
ORFLength1=[];
for i=1:6
   for j=1:length(orf1(i).Stop)
    ORFLength1=[ORFLength1; orf1(i).Stop(j)+2 - orf1(i).Start(j)];
   end
end

%%
max_threshold=max(ORFLength1)

n_max=length(find(ORFLength>=max_threshold))
%%
orf_threshold = seqshoworfs(chimp,'MINIMUMLENGTH',max_threshold/3, 'geneticcode',2,'frames','all','nodisplay','true');

%% convert orf sequences to amino acids
codons=[];
aminoacid=[];
for i=1:6
    for j= 1:length(orf_threshold(i).Stop)
        proteinSeq = chimpSeq(orf_threshold(i).Start(j):orf_threshold(i).Stop(j)+2);
        codons= [codons; string(proteinSeq)];
        aminoacid = [aminoacid; string(nt2aa(proteinSeq, 'geneticcode',2)),orf_threshold(i).Start(j),orf_threshold(i).Stop(j)+2];
    end
end
%%

%% check for cytb protein
[startorf, endorf] = sigTest(aminoacid, human_cytb);
%%
chimp_cytb_aminoacid = nt2aa(chimpSeq(str2double(startorf):str2double(endorf)),'geneticcode',2);
chimp_cytb_aminoacid(1:50)
%% Check for cox1 protein
[startorf, endorf] = sigTest(aminoacid, human_cox1);
%%
chimp_cox1_aminoacid = nt2aa(chimpSeq(str2double(startorf):str2double(endorf)),'geneticcode',2);
chimp_cox1_aminoacid(1:50)
%%

%%
function [startorf, endorf] = sigTest(aminoacid, human_protein)
    scores = [];
    for j=1:length(aminoacid)
        scores = [scores; j,nwalign(char(human_protein), char(aminoacid(j,1)), 'scoringmatrix','blosum30','gapopen',5,'extendgap',5)];
    end
    startorf = 0;
    endorf = 0;
    for i=1:length(scores)
        c = scores(i,1);
        score = scores(i,2);
        chimp_protein = char(aminoacid(c,1));
        n = 1000;
        globalscores = zeros(n,1);
        chimpLen = length(chimp_protein);
        for j = 1:n
            perm = randperm(chimpLen);
            globalscores(j) = nwalign(human_protein,chimp_protein(perm),'scoringmatrix','blosum30','gapopen',5,'extendgap',5);
        end
        p_value = length(find(globalscores >= score)) / 1000;
        if (p_value < 0.01)
            figure
            buckets = ceil(n/5);
            hist(globalscores,buckets)
            hold on;
            stem(score,1,'k')
            xlabel('Score'); ylabel('Number of Sequences');
            hold off;
            score
            p_value
            startorf = aminoacid(c,2);
            endorf = aminoacid(c,3);
            strcat('Human protein is homologous to ORF in position (',aminoacid(c,2),':',aminoacid(c,3), ') in chimpanzee')
        end
    end
end
