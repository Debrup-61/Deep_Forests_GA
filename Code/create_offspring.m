function [Offspring,FVal_offsp] = create_offspring(parent1,parent2,Prey)
%Function Description
%{  
   This function returns the offspring and the fitness of the offspring
   given the index of parent1,parent2 and the Prey population.It should
   perform crossover between 2 forest individuals and then perform mutation
   on the offspring.We know that the forests can be pretty deep so we would
   want the crossover's between the tree blocks to be performed in parallel
   to speed up the computational time.
%}

%% Global variables definition
global Setslog
n_subnets=Setslog.n_subnets;            % Number of subnets in each forest
n_blocks=Setslog.n_blocks;              % Number of blocks in each forest

%% Assigning the parent forests
Prey1=Prey(1,parent1);                  %Prey 1
Prey2=Prey(1,parent2);                  %Prey 2

%% Loop over the subnets of Prey1 and Prey2
   %parfor_block performs crossover and mutation between the ith subnet of 2
   %parent forests

for i=1:n_subnets
    subnet=cell(2,n_blocks(i));
    subnet1=Prey1.subnet{i};
    subnet2=Prey2.subnet{i};
    parfor j=1:n_blocks(i)                                  %Parfor loop performs xover and mutation between all blocks in parallel                     
          arr=parfor_block([subnet1{j},subnet2{j}],i,j);    %arr is a cell array of size 2*1
          subnet(:,j)=arr;
    end
    Prey1.subnet{i}=subnet(1,:);                            %Assignment of new subnet i for Prey 1 and 2
    Prey2.subnet{i}=subnet(2,:);
end


Offspring(1,1)=Prey1;
Offspring(1,2)=Prey2;
Offspring(1,1).evaluate=1;
Offspring(1,2).evaluate=1;

FVal_offsp=zeros(2,2);
for i=1:2
    [err,complexity,endnet,out,in_end]=DFOReval(Offspring(1,i));
    FVal_offsp(i,:)=[err,complexity];
    Offspring(1,i).evaluate=0;
    Offspring(1,i).err = err;
    Offspring(1,i).complexity = complexity;
    Offspring(1,i).endnet = endnet;
    Offspring(1,i).out=out;
    Offspring(1,i).in_end=in_end;
end



