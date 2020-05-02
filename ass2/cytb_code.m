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

distances = seqpdist(cytb_data,'method','jukes-cantor','indels','pairwise-delete','squareform',true)

NJtree = seqneighjoin(distances,'equivar',cytb_data);
h = plot(NJtree,'orient','left');
title('Neighbor-joining tree using Jukes-Cantor model');

tree1=reroot(NJtree,7);
title('Neighbor-joining tree using Jukes-Cantor model');
