function data = node_hit_q(data)

node  = data.node;
nmod  = length(node.bin1);
hit_q = zeros(nmod,1);

for inode = 1:nmod
    q = zeros(7,1);
    
    if node.bin1(inode)<5
        q(1) = node.bin1(inode);
    else
        q(1) = 5;
    end
    
    if node.bin2(inode)<5
        q(2) = node.bin2(inode);
    else
        q(2) = 5;
    end
    
    if node.bin3(inode)<5
        q(3) = node.bin3(inode);
    else
        q(3) = 5;
    end
    
    if node.bin4(inode)<5
        q(4) = node.bin4(inode);
    else
        q(4) = 5;
    end
    
    if node.core(inode)<5
        q(7) = node.core(inode);
    else
        q(7) = 5;
    end
    
    hit_q(inode) = sum(q)/25;
    
end

data.node.hit_q = hit_q;