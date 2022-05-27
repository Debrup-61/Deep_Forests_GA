function [in,n,h] = eval_block(Prey1,in,subnet,block_no)

%% This function evaluates one tree block given the input,subnet_no,block_no and the forest at hand named Prey1.
%% in:Input terminal set to be used in the block
%% n: number of nodes in the block
%% h: height of the block
%% subnet:subnet-number

%% Global variables
global T 

       
%% Check if prey is to be evaluated

if(Prey1.evaluate==1)
    if(block_no==1) 
        ter_set=T(subnet).terminals; 
        for i=1:size(in,2)
          eval([ter_set{i} '= in(:,i);']);       %example: ter_set{1}=X1,so X1=in(:,i) makes X1 a matrix 
        end     
    else
        eval(['L' '=in;']);                      %L is set to a column vector of outputs from previous block for m training examples
    end   
    
    block_tree=Prey1.subnet{subnet}{block_no};
    in=eval(cat(2,block_tree.tree{:}));          %in now has the output which is a column vector of m training examples
    
    h=length(dec2ari(max(block_tree.tree_index)));
    block_tree.tree_index(block_tree.T.index>0)=[];
    n=length(unique(block_tree.tree_index));
end
end   


