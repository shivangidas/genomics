%% Task 1
% get the human cytb and cox1 genes for comparision
human_cytb = getgenpept('AFI52643','sequenceonly','true');
human_cox1 = getgenpept('AGW78696','sequenceonly','true');
% Get the raw Pan troglodytes mitochondrion dna
chimp = getgenbank('NC_001643');
chimpSeq = chimp.Sequence;
% get ORFs
orf = seqshoworfs(chimpSeq,'MINIMUMLENGTH',3, 'geneticcode',2,'frames','all','nodisplay','true');
% get ORFs in a similar random permutaion
orf1 = seqshoworfs(chimpSeq(randperm(length(chimpSeq))),'MINIMUMLENGTH',3,'geneticcode',2,'frames','all','nodisplay','true');
% collect all ORFs
ORFLength=[];
for i=1:6
   for j=1:length(orf(i).Stop)
    ORFLength=[ORFLength; orf(i).Stop(j)+2 - orf(i).Start(j)];
   end
end
disp("Number of ORFs in chimpanzee mitochondrion")
length(ORFLength)

% thresholding using random permutaion sequence's max length ORF
ORFLength1=[];
for i=1:6
   for j=1:length(orf1(i).Stop)
    ORFLength1=[ORFLength1; orf1(i).Stop(j)+2 - orf1(i).Start(j)];
   end
end

%
max_threshold=max(ORFLength1)

n_max=length(find(ORFLength>=max_threshold))
%
orf_threshold = seqshoworfs(chimp,'MINIMUMLENGTH',max_threshold/3, 'geneticcode',2,'frames','all','nodisplay','true');

% convert orf sequences to amino acids
codons=[];
aminoacid=[];
for i=1:6
    for j= 1:length(orf_threshold(i).Stop)
        proteinSeq = chimpSeq(orf_threshold(i).Start(j):orf_threshold(i).Stop(j)+2);
        codons= [codons; string(proteinSeq)];
        aminoacid = [aminoacid; string(nt2aa(proteinSeq, 'geneticcode',2)),orf_threshold(i).Start(j),orf_threshold(i).Stop(j)+2];
    end
end


% Check for cytb protein
[startorf, endorf] = sigTest(aminoacid, human_cytb);
chimp_cytb_aminoacid = nt2aa(chimpSeq(str2double(startorf):str2double(endorf)),'geneticcode',2);
chimp_cytb_aminoacid(1:50)

% Check for cox1 protein
[startorf, endorf] = sigTest(aminoacid, human_cox1);
chimp_cox1_aminoacid = nt2aa(chimpSeq(str2double(startorf):str2double(endorf)),'geneticcode',2);
chimp_cox1_aminoacid(1:50)
%% Task 3 Phylogenetic tree

data_cytb = {'Pan troglodytes'   'NP_008198';
        'Pan troglodytes verus' 'ACY65608';
        'Pan troglodytes ellioti' 'AMB65389';
        'Pan troglodytes troglodytes' 'ACY65530';
        'Pan troglodytes schweinfurthii' 'ACY65621';
        'Pan paniscus' 'ADA55535';
        'Papio cynocephalus' 'AFX82173';};

NumCytb=length(data_cytb);

for ind = 1:NumCytb
    cytb_data(ind).Header   = data_cytb{ind,1};
    cytb_data(ind).Sequence = getgenpept(data_cytb{ind,2},'sequenceonly','true');
end

distances = seqpdist(cytb_data,'method','jukes-cantor','indels','pairwise-delete','squareform',true);

% rooted tree
NJtree = seqneighjoin(distances,'equivar',cytb_data);
sel = getbyname(NJtree,'Papio cynocephalus');
tree1=reroot(NJtree,sel);
h = plot(tree1,'orient','left');
title('Neighbor-joining tree using Jukes-Cantor model');

%% Task 4
% full mtDNA accession numbers
 
data_mtDna = {'Pan troglodytes'   'NC_001643';
        'Pan troglodytes verus' 'JF727213';
        'Pan paniscus' 'NC_001644';
        'Papio cynocephalus' 'NC_020007';
        'Homo sapiens' 'NC_012920';
        'Gorilla gorilla gorilla' 'NC_011120';
        'Pongo abelii' 'NC_002083'};

% number of taxa
Num=length(data_mtDna);

% extract the cytb and cox1 genes

for ind = 1:Num
    cytb_nt_data(ind).Header   = data_mtDna{ind,1};
    cox1_nt_data(ind).Header = data_mtDna{ind,1};
    seqData = getgenbank(data_mtDna{ind,2});
    seq = seqData.Sequence;
    gene_coding_seq = seqData.CDS;
    for i=1:length(gene_coding_seq)
        if strcmp(gene_coding_seq(i).gene, 'CYTB')
            orfStart = gene_coding_seq(i).indices(1,1);
            orfEnd =  gene_coding_seq(i).indices(1,2);
            cytbSeq = seq(orfStart:orfEnd);
        elseif strcmp(gene_coding_seq(i).gene, 'COX1')
            orfStart = gene_coding_seq(i).indices(1,1);
            orfEnd =  gene_coding_seq(i).indices(1,2);
            cox1Seq = seq(orfStart:orfEnd);
        end
    end
    cytb_nt_data(ind).Sequence = cytbSeq;
    cox1_nt_data(ind).Sequence = cox1Seq;
end

% convert nucleotide sequences to amino acid

for ind = 1:Num
    cytb_aa_data(ind).Header = cytb_nt_data(ind).Header;
    cytb_aa = nt2aa(cytb_nt_data(ind).Sequence,'geneticcode',2);
    cytb_aa_data(ind).Sequence = cytb_aa;
    cox1_aa_data(ind).Header = cox1_nt_data(ind).Header;
    cox1_aa = nt2aa(cox1_nt_data(ind).Sequence,'geneticcode',2);
    cox1_aa_data(ind).Sequence = cox1_aa;
end

% genetic distances 

distances_cytb_nt = seqpdist(cytb_nt_data,'method','jukes-cantor','ALPHABET','NT','PAIRWISEALIGNMENT',true,'squareform',true)

distances_cytb_aa = seqpdist(cytb_aa_data,'method','jukes-cantor','PAIRWISEALIGNMENT',true,'squareform',true)

distances_cox1_nt = seqpdist(cox1_nt_data,'method','jukes-cantor','ALPHABET','NT','PAIRWISEALIGNMENT',true,'squareform',true)

distances_cox1_aa = seqpdist(cox1_aa_data,'method','jukes-cantor','PAIRWISEALIGNMENT',true,'squareform',true)


%% Task 5 Phylogenetic trees
% neighbour-joining tree for cytb nt
NJtree_cytb_nt = seqneighjoin(distances_cytb_nt,'equivar',cytb_nt_data);
sel = getbyname(NJtree_cytb_nt,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cytb_nt,sel);
tree = plot(rooted_tree,'orient','left');
set(tree.terminalNodeLabels,'FontSize',12);
title('CYTB NT');

% neighbour-joining tree for cytb aa
NJtree_cytb_aa = seqneighjoin(distances_cytb_aa,'equivar',cytb_aa_data);
sel = getbyname(NJtree_cytb_aa,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cytb_aa,sel);
plot(rooted_tree,'orient','left');
title('CYTB AA');

% neighbour-joining tree for cox1 nt
NJtree_cox1_nt = seqneighjoin(distances_cox1_nt,'equivar',cox1_nt_data);
sel = getbyname(NJtree_cox1_nt,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cox1_nt,sel);
plot(rooted_tree,'orient','left');
title('COX1 NT');

% neighbour-joining tree for cox1 aa
NJtree_cox1_aa = seqneighjoin(distances_cox1_aa,'equivar',cox1_aa_data);
sel = getbyname(NJtree_cox1_aa,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cox1_aa,sel);
plot(rooted_tree,'orient','left');
title('COX1 AA');
%% consensus tree using weighted average
wdistances_cytb_aa = seqpdist(cytb_aa_data,'method','jukes-cantor','PAIRWISEALIGNMENT',true);
wdistances_cox1_aa = seqpdist(cox1_aa_data,'method','jukes-cantor','PAIRWISEALIGNMENT',true);

weights = [sum(wdistances_cytb_aa) sum(wdistances_cox1_aa)];
weights = weights / sum(weights);
distances = wdistances_cytb_aa .* weights(1) + wdistances_cox1_aa .* weights(2);
NJtree = seqneighjoin(distances,'equivar',cox1_aa_data);
sel = getbyname(NJtree,'Papio cynocephalus');
rooted_tree = reroot(NJtree,sel);
plot(rooted_tree,'orient','left');
title('Consensus tree');

%%
wdistances_cytb_nt = seqpdist(cytb_nt_data,'method','jukes-cantor','ALPHABET','NT','PAIRWISEALIGNMENT',true);
wdistances_cox1_nt = seqpdist(cox1_nt_data,'method','jukes-cantor','ALPHABET','NT','PAIRWISEALIGNMENT',true);

weights = [sum(wdistances_cytb_nt) sum(wdistances_cox1_nt)];
weights = weights / sum(weights);
distances = wdistances_cytb_nt .* weights(1) + wdistances_cox1_nt .* weights(2);
NJtree = seqneighjoin(distances,'equivar',cox1_nt_data);
sel = getbyname(NJtree,'Papio cynocephalus');
rooted_tree = reroot(NJtree,sel);
plot(rooted_tree,'orient','left');
title('Consensus tree NT');
%% Not used, seqlinkage tree for cytb nt
% envtree = seqlinkage(distances_cytb_nt,'UPGMA',cytb_nt_data)
% %plot(envtree,'type','angular');
% sel = getbyname(envtree,'Papio cynocephalus');
% rooted_tree = reroot(envtree,sel);
% plot(rooted_tree,'type','angular');
% title('UPGMA CYTB NT')

%% Task 7 Multiple alignment
% cytb nt
ma_cytb_nt = multialign(cytb_nt_data,NJtree_cytb_nt);
seqalignviewer(ma_cytb_nt);
% cox1 nt
ma_cox1_nt = multialign(cox1_nt_data,NJtree_cox1_nt);
seqalignviewer(ma_cox1_nt);
% cytb aa
ma_cytb_aa = multialign(cytb_aa_data,NJtree_cytb_aa);
seqalignviewer(ma_cytb_aa);
% cox1 aa
ma_cox1_aa = multialign(cox1_aa_data,NJtree_cox1_aa);
seqalignviewer(ma_cox1_aa);
%% function for homologous genes
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
%%
