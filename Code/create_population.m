function [Prey,FVal] = create_population(pop_size)

%% CREATE_POPULATION Function to create a pop of trees 
%% Use looping to assert that the trees are not nan trees.
%% Use a Ramped Half and Half method to initialize trees.
%% Set the properties of the Prey like error,complexity,evaluate and so on.

%% General parameters defined as Global
global Setslog;
blocks=Setslog.n_blocks;                 % blocks has no of blocks in each subnet 
n_subnets=Setslog.n_subnets;             % No of subnets in the trees
max_bdepth=Setslog.max_block_depth;      % Max depth of a tree block
FVal=zeros(pop_size,2);                  % Defining FVal

%% Loop over the population

for i=1:pop_size
    for j=1:n_subnets
        for k=1:blocks(j)
            
            r=unidrnd(max_bdepth);       % Random height of the block by discrete uniform dist 
            r1=rand;                     % Random number to select full/grow method 
            pass=0;                      % pass var tells us whether to move to next block
            
            if(r1<0.5)    
               while(pass~=1)
                  tree_block=full_tree(r,'1',j,k);   % Use the full method
                  pass=test_block(tree_block,j,k);
                  if(pass==1)
                     Prey(i).subnet{j}{k}=tree_block;
                     
                  end    
               end
              
            else
                
                while(pass~=1)
                    tree_block=grow_tree(r,'1',j,k);   % Use the full method
                    pass=test_block(tree_block,j,k);
                    if(pass==1)
                      Prey(i).subnet{j}{k}=tree_block;
                    end    
                end
               
            end
            
        tree = Prey(i).subnet{j}{k}.tree_index;
        Prey(i).subnet{j}{k}.depth = length(dec2ari(max(tree)));
        tree(Prey(i).subnet{j}{k}.T.index > 0) = [];
        Prey(i).subnet{j}{k}.nodes = length(unique(tree));
            
            
        end
        
    end
    
    Prey(i).evaluate=1;                    % Need to evaluate tree
    
    [err,complexity,coeff,modelout,in_end,height_n,nodes,max_height]=DFOReval(i,Prey);
    FVal(i,:)=[err,complexity];            % Insert error and complexity in FVal
    Prey(i).error=err;
    Prey(i).complexity=complexity;
    Prey(i).coefficient=coeff;
    Prey(i).output=modelout;
    Prey(i).in_end=in_end;
    Prey(i).avg_subnet_ht=height_n;
    Prey(i).nodes=nodes;
    Prey(i).max_subnet_ht=max_height;
    Prey(i).evaluate=0;                    % Already evaluated tree
    
end

