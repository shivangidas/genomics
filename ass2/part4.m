% full mtDNA accession number 
%% 
data_mtDna = {'Pan troglodytes'   'NC_001643';
        'Pan troglodytes verus' 'JF727213';
        'Pan paniscus' 'NC_001644';
        'Papio cynocephalus' 'NC_020007';
        'Homo sapiens' 'NC_012920';
        'Gorilla gorilla gorilla' 'NC_011120';
        'Pongo abelii' 'NC_002083'};
%% 
% number of taxa    
%% 
Num=length(data_mtDna);
%% 
%extract the cytb and cox1 genes
%% 
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
%%
% convert nucleotide sequences to amino acid
%% 
for ind = 1:Num
    cytb_aa_data(ind).Header = cytb_nt_data(ind).Header;
    cytb_aa = nt2aa(cytb_nt_data(ind).Sequence,'geneticcode',2);
    cytb_aa_data(ind).Sequence = cytb_aa;
    cox1_aa_data(ind).Header = cox1_nt_data(ind).Header;
    cox1_aa = nt2aa(cox1_nt_data(ind).Sequence,'geneticcode',2);
    cox1_aa_data(ind).Sequence = cox1_aa;
end
%% 

%genetic distances 
%% 
distances_cytb_nt = seqpdist(cytb_nt_data,'method','jukes-cantor','ALPHABET','NT','PAIRWISEALIGNMENT',true,'squareform',true)

distances_cytb_aa = seqpdist(cytb_aa_data,'method','jukes-cantor','PAIRWISEALIGNMENT',true,'squareform',true)

distances_cox1_nt = seqpdist(cox1_nt_data,'method','jukes-cantor','ALPHABET','NT','PAIRWISEALIGNMENT',true,'squareform',true)

distances_cox1_aa = seqpdist(cox1_aa_data,'method','jukes-cantor','PAIRWISEALIGNMENT',true,'squareform',true)
%% 

%trees

%% neighbour-joining tree for cytb nt
NJtree_cytb_nt = seqneighjoin(distances_cytb_nt,'equivar',cytb_nt_data);
sel = getbyname(NJtree_cytb_nt,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cytb_nt,sel);
tree = plot(rooted_tree,'orient','left');
set(tree.terminalNodeLabels,'FontSize',12);
title('CYTB NT');
%% seqlinkage tree for cytb nt
envtree = seqlinkage(distances_cytb_nt,'UPGMA',cytb_nt_data)
%plot(envtree,'type','angular');
sel = getbyname(envtree,'Papio cynocephalus');
rooted_tree = reroot(envtree,sel);
plot(rooted_tree,'type','angular');
title('UPGMA CYTB NT')
%% neighbour-joining tree for cytb aa
NJtree_cytb_aa = seqneighjoin(distances_cytb_aa,'equivar',cytb_aa_data);
%h = plot(NJtree_cytb_aa,'orient','left');
sel = getbyname(NJtree_cytb_aa,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cytb_aa,sel);
plot(rooted_tree,'orient','left');
title('CYTB AA');
%% neighbour-joining tree for cox1 nt
NJtree_cox1_nt = seqneighjoin(distances_cox1_nt,'equivar',cox1_nt_data);
%h = plot(NJtree_cox1_nt,'orient','left');
sel = getbyname(NJtree_cox1_nt,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cox1_nt,sel);
plot(rooted_tree,'orient','left');
title('COX1 NT');
%% neighbour-joining tree for cox1 aa
NJtree_cox1_aa = seqneighjoin(distances_cox1_aa,'equivar',cox1_aa_data);
%h = plot(NJtree_cox1_aa,'orient','left');
sel = getbyname(NJtree_cox1_aa,'Papio cynocephalus');
rooted_tree = reroot(NJtree_cox1_aa,sel);
plot(rooted_tree,'orient','left');
title('COX1 AA');
%% 
%consensus tree
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

%% Multiple alignement
%%
ma_cytb_nt = multialign(cytb_nt_data,NJtree_cytb_nt);
seqalignviewer(ma_cytb_nt);
%%
ma_cox1_nt = multialign(cox1_nt_data,NJtree_cox1_nt);
seqalignviewer(ma_cox1_nt);
%%
ma1 = multialign(cytb_aa_data,NJtree_cytb_aa);
seqalignviewer(ma1);
%%
ma2 = multialign(cox1_aa_data,NJtree_cox1_aa);
seqalignviewer(ma2);