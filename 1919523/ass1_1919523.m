%Question 6
mito_gbk = getgenbank('NC_001643'); %chimp
chimp = mito_gbk.Sequence;
mitochondria_length = length(chimp)
first_300_bases = seqdisp(chimp(1:300))

figure
ntdensity(chimp)

bases = basecount(chimp)

compBases = basecount(seqrcomplement(chimp))

figure
basecount(chimp,'chart','pie');
title('Distribution of Nucleotide Bases for Chimpanzee Mitochondrial Genome');

%question 7 and 8
chimp = getgenbank('NC_001643','SequenceOnly',true);
whos chimp
codoncount(chimp)

%get ORFs
orf = seqshoworfs(chimp,'MINIMUMLENGTH',3, 'geneticcode',2,'frames','all','nodisplay','true');
%get ORFs in a similar random permutaion
orf1 = seqshoworfs(chimp(randperm(length(chimp))),'MINIMUMLENGTH',3,'geneticcode',2,'frames','all','nodisplay','true');

%collect all ORFs
ORFLength=[];
for i=1:6
   for j=1:length(orf(i).Stop)
    ORFLength=[ORFLength; orf(i).Stop(j)+2 - orf(i).Start(j)];
   end
end
disp("Number of ORFs in chimpanzee mitochondrion")
length(ORFLength)

%thresholding using random permutaion sequence's max length ORF
ORFLength1=[];
for i=1:6
   for j=1:length(orf1(i).Stop)
    ORFLength1=[ORFLength1; orf1(i).Stop(j)+2 - orf1(i).Start(j)];
   end
end
length(ORFLength1)

% % threshold can also be done for multiple permutaions.
% ORFLength1=[];
% maxOrfLength = [];
% for k=1:100
%     orf1= seqshoworfs(chimp(randperm(length(chimp))),'MINIMUMLENGTH',1,'geneticcode',2,'frames','all','nodisplay','true');
%     for i=1:6
%        for j=1:length(orf1(i).Stop)
%         ORFLength1=[ORFLength1; orf1(i).Stop(j)+2 - orf1(i).Start(j)];
%        end
%     end
%     maxOrfLength=[maxOrfLength;max(ORFLength1)];
% end
% max_threshold=max(maxOrfLength)

max_threshold=max(ORFLength1)
disp("Number of significant ORFs in chimpanzee mitochondrion")
n_max=length(find(ORFLength>=max_threshold))

%get significant ORFs
orf_threshold = seqshoworfs(chimp,'MINIMUMLENGTH',max_threshold/3, 'geneticcode',2,'frames','all','nodisplay','true');
orfNew=[];
for i=1:6
    for j=1:length(orf_threshold(i).Stop)
        orfNew=[orfNew; orf_threshold(i).Stop(j)+2 - orf_threshold(i).Start(j)];
    end
end


%question9
protein = chimp(orf_threshold(1).Start(1):orf_threshold(1).Stop(1)+2) %first ORF selected
aminoacid1 = nt2aa(protein, 'geneticcode',2)
aminoacid1(1:50)
seq1 = getgenpept('Q9T9W3','SequenceOnly',true);
seq2 = getgenpept('AIV00479','SequenceOnly',true);
seq3 = getgenpept('AEQ36262','SequenceOnly',true);
seq4 = getgenpept('ANQ92411','SequenceOnly',true);
seq5 = getgenpept('AMB65312','SequenceOnly',true);
seq6 = getgenpept('AJO25286','SequenceOnly',true);
seqs = {seq1, seq2, seq3, seq4, seq5, seq6};
multialign(seqs)

%question10
human = getgenbank('NC_012920','SequenceOnly',true);

humanProtein = getgenpept('AFF91323','SequenceOnly',true);
chimpProtein = aminoacid1;
seqdotplot(humanProtein,chimpProtein)
xlabel('human');ylabel('chimpanzee');

[sc50,globAlig50] = nwalign(humanProtein,chimpProtein);
fprintf('Score = %g \n',sc50)
[sc30,globAlig30] = nwalign(humanProtein,chimpProtein,'scoringmatrix','blosum30');
fprintf('Score = %g \n',sc30)
showalignment(globAlig50);

%significance test
n = 1000;
globalscores = zeros(n,1);
chimpLen = length(chimpProtein);
for j = 1:n
    perm = randperm(chimpLen);
    globalscores(j) = nwalign(humanProtein,chimpProtein(perm),'scoringmatrix','blosum30','gapopen',5,'extendgap',5);
end
figure
buckets = ceil(n/5);
hist(globalscores,buckets)
hold on;
stem(sc50,1,'k')
xlabel('Score'); ylabel('Number of Sequences');
titleStr = strcat('Significance plot for human and chimpanzee NADH dehydrogenase subunit 1 (mitochondrion)');
title(titleStr);
hold off;
p_value50 = length(find(globalscores >= sc50)) / 1000
p_value30 = length(find(globalscores >= sc30)) / 1000

% question 11
h_orf = seqshoworfs(human,'MINIMUMLENGTH',1, 'geneticcode',2,'frames','all','nodisplay','true');
h_orf1= seqshoworfs(human(randperm(length(human))),'MINIMUMLENGTH',1,'geneticcode',2,'frames','all','nodisplay','true');
h_ORFLength=[];
for i=1:6
   for j=1:length(h_orf(i).Stop)
    h_ORFLength=[h_ORFLength; h_orf(i).Stop(j)+2 - h_orf(i).Start(j)];
   end
end
disp("Number of ORFs in human mitochondrion")
length(h_ORFLength)
h_ORFLength1=[];
for i=1:6
   for j=1:length(h_orf1(i).Stop)
    h_ORFLength1=[h_ORFLength1; h_orf1(i).Stop(j)+2 - h_orf1(i).Start(j)];
   end
end
length(h_ORFLength1)
h_max_threshold=max(h_ORFLength1)
disp("Number of significant ORFs in human mitochondrion")
h_n_max=length(find(h_ORFLength>=h_max_threshold))
h_orf_threshold = seqshoworfs(human,'MINIMUMLENGTH',h_max_threshold/3, 'geneticcode',2,'frames','all','nodisplay','true');

h_orfNew=[];
for i=1:6
    for j=1:length(h_orf_threshold(i).Stop)
        h_orfNew=[h_orfNew; h_orf_threshold(i).Stop(j)+2 - h_orf_threshold(i).Start(j)];
    end
end

%convert to amino acid and store start stop values
codons=[];
aminoacid=[];
for i=1:6
    for j= 1:length(orf_threshold(i).Stop)
        proteinSeq = chimp(orf_threshold(i).Start(j):orf_threshold(i).Stop(j)+2);
        codons= [codons; string(proteinSeq)];
        aminoacid = [aminoacid; string(nt2aa(proteinSeq, 'geneticcode',2)),orf_threshold(i).Start(j),orf_threshold(i).Stop(j)+2];
    end
end

% length(codons)
% length(aminoacid)
h_protein=[];
h_aminoacid=[];
for i=1:6
    for j= 1:length(h_orf_threshold(i).Stop)
        protein = human(h_orf_threshold(i).Start(j):h_orf_threshold(i).Stop(j)+2);
        h_protein= [h_protein; string(protein)];
        h_aminoacid = [h_aminoacid; string(nt2aa(protein, 'geneticcode',2)),h_orf_threshold(i).Start(j),h_orf_threshold(i).Stop(j)+2];
    end
end
% length(h_protein)
% length(h_aminoacid)


all_amino = [h_aminoacid;aminoacid];
% length(all_amino)
% exploring with multialign
ma = multialign(char(all_amino(:,1)))
seqalignviewer(ma)

%get alignment scores for all combinations and threshold at 100
scores = [];
for i=1: length(h_aminoacid)
    for j=1:length(aminoacid)
        scores = [scores; i,j,nwalign(char(h_aminoacid(i,1)), char(aminoacid(j,1)), 'scoringmatrix','blosum30','gapopen',5,'extendgap',5)];
    end
end

total_combinations = length(scores)

%test with significance plot
for i=1:total_combinations
    h = scores(i,1);
    c = scores(i,2);
    score = scores(i,3);
    human_protein = char(h_aminoacid(h,1));
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
         titleStr = strcat('Significance plot for ORF(' ,h_aminoacid(h,2),':',h_aminoacid(h,3), ...
             ') human and ORF(',aminoacid(c,2),':',aminoacid(c,3), ') chimpanzee');
        title(titleStr);
        hold off;
        [score1, alignment] = nwalign(human_protein, chimp_protein, ...
            'scoringmatrix','blosum30','gapopen',5,'extendgap',5);
        showalignment(alignment)
        score
        p_value
        strcat('ORF in position(' ,h_aminoacid(h,2),':',h_aminoacid(h,3), ...
              ') in human is homologous to ORF in position (',aminoacid(c,2),':',aminoacid(c,3), ') in chimpanzee')
    end
end

% % Uncomment to run question 12
% 
% 
% h_nt = h_protein;
% c_nt = codons;
% 
% scoresNT = [];
% for i=1: length(h_nt)
%     for j=1:length(c_nt)
%         scoresNT = [scoresNT; i,j,nwalign(char(h_nt(i,1)), char(c_nt(j,1)), 'scoringmatrix','blosum30','gapopen',5,'extendgap',5,'ALPHABET','NT')];
%     end
% end
% 
% total_combinations1 = length(scoresNT)
% 
% for i=1:total_combinations1
%     h = scoresNT(i,1);
%     c = scoresNT(i,2);
%     scoreNT = scoresNT(i,3);
%     human_NT = char(h_nt(h,1));
%     chimp_NT = char(c_nt(c,1));
%     n = 1000;
%     globalscores = zeros(n,1);
%     chimpLen = length(chimp_NT);
%     for j = 1:n
%         perm = randperm(chimpLen);
%         globalscores(j) = nwalign(human_NT,chimp_NT(perm),'scoringmatrix','blosum30','gapopen',5,'extendgap',5,'ALPHABET','NT');
%     end
%     p_value = length(find(globalscores >= scoreNT)) / 1000;
%     if (p_value < 0.01)
%         figure
%         buckets = ceil(n/5);
%         hist(globalscores,buckets)
%         hold on;
%         stem(scoreNT,1,'k')
%         xlabel('Score'); ylabel('Number of Sequences');
%          titleStr = strcat('Significance plot for ORF(' ,h_nt(h,2),':',h_nt(h,3), ...
%              ') human and ORF(',c_nt(c,2),':',c_nt(c,3), ') chimpanzee');
%         title(titleStr);
%         hold off;
%         [score1, alignment] = nwalign(human_NT, chimp_NT, ...
%             'scoringmatrix','blosum30','gapopen',5,'extendgap',5,'ALPHABET','NT');
%         showalignment(alignment)
%         scoreNT
%         p_value
%         strcat('ORF in position(' ,h_aminoacid(h,2),':',h_aminoacid(h,3), ...
%               ') in human is homologous to ORF in position (',aminoacid(c,2),':',aminoacid(c,3), ') in chimpanzee')
%     end
% end