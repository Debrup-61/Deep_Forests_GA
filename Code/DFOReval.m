function [err,complexity,coeff,modelout,in_end,height_n,nodes,max_height] = DFOReval(i,Prey)
%% DFOReVAL This function takes as input a forest individual and the terminal set T.It calculates the prediction made by the tree which 
%% is used to get the error.The complexity of the tree is also returned(size of the tree,etc can be used).We will also try and add a 
%% feature later so that very long blocks are replaced by smaller blocks possibly by splitting or other ways to counter the problem
%% of bloat.

%% Global variables
global Setslog run_no lambda F_bad

%% Defining variables
Prey1=Prey(i);
out = Setslog.dataset(run_no).out;                %True Output values from training set 
outmax = max(out);                                %Maximum value in output
outmin = min(out);                                %Minimum value in output
no_points = length(out);                          %No of training examples
num_subnets=Setslog.n_subnets;                    %No of Subnets
in_end=ones(no_points,num_subnets+1);             %Matrix of outputs from each subnet before LLSQ(first column is column of 1's for bias)      
nodes=0;                                          %Total Number of nodes in forest
n_blocks=Setslog.n_blocks;
max_blocks=max(n_blocks);
height=zeros(max_blocks,num_subnets);             %height of subnets in forest
%complexity=0;                                    %Overall complexity of the forest

%% Main loop on subnets
warning('off');
for subnet=1:num_subnets              %CAN EVALUATE ALL SUBNETS IN PARALLEL USING PARFOR
     in =Setslog.dataset(run_no).in;
     in = in(:,Setslog.Pop_str{subnet}{1});
     n_blocks=Setslog.n_blocks(subnet);
     
     for block_no=1:n_blocks
         
         [in,n,h]=eval_block(Prey1,in,subnet,block_no); %eval_block evaluates one tree block 
         height(subnet,block_no)=h;   % Store height of each tree block in matrix                        
         nodes=nodes+n;              % Store total number of nodes in structure
     end
    
%  if isnan(sum(in)) || isinf(sum(in))  % If sum(in) is nan or inf then there are undefined ops done
%    in=zeros(no_points,1);
%  end

in_end(:,subnet+1) =in;
end

height_n=mean(sum(height,2));         % height_n has the average height of the subnets
max_height=max(sum(height,2));        % max_height gives the height of a longest subnet

%% Replace the NaN or Inf value in in_end by 0

%in_end(isnan(in_end))=0;
%in_end(isinf(in_end))=0;

%% Assigning weights to each subnet and bias to the forest to solve LLSQ in_end*coeff=out

coeff = in_end\out;
modelout = in_end*coeff;
err = sqrt(sum(((out-modelout)/(outmax-outmin)).^2)/no_points);
if isnan(sum(err)) || isinf(sum(err))
    err=F_bad+eps;
    complexity=F_bad+eps;
else 
complexity=lambda*height_n+(1-lambda)*nodes;
end
warning('on');

end

