function pass= test_block(block_tree,subnet_no,block_no)
%% TEST_BLOCK tests if the tree block is independent of input or has NaN values(invalid) 
%% If block_no not equal to 1 then set L to a random column vector of arbitrary size and 
%% check if the value is Nan or Inf.

%% Definition of variables
   global T;
   set_size=200;
   
%% Check out the output for a random assignment of input variables  
    if(block_no==1) 
        ter_set=T(subnet_no).terminals;
        len=length(ter_set);
        in=rand(set_size,len);
        for i=1:size(in,2)
          eval([ter_set{i} '= in(:,i);']);       %example: ter_set{1}=X1,so X1=in(:,i) makes X1 a matrix 
        end     
    else
        in=rand(set_size,1);
        eval(['L' '=in;']);                      %L is set to a column vector of outputs from previous block for m training examples
    end   
    
   
    out=eval(cat(2,block_tree.tree{:}));         %out contains the output
    
%% Check if the output contains NaN or Inf

if isnan(sum(out)) || isinf(sum(out))
    pass=0;
else
    pass=1;
end

%% End of function
  
end

