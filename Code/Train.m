function [Prey,FVal]=Train(Problem_name,parameters)
%TRAIN: 
%Training of a deep model for the dependent var using training data 
%with the use of subnets comprising tree structures.

%Save the data file containing training data in xlsx format to the 
%newmastercode directory before executing AutoRun.m

%===============================================================
global Setslog figure_handle 

figure_handle=[];                      % Array of the figure objects
RandStream('mt19937ar','seed', sum(100*clock));
Setslog = [];
in_index =parameters.in_index;         
out_index=parameters.out_index;        

%% Defining the filename and making the output folder
filename=[Problem_name '.xlsx'];
savedir = fullfile(pwd,'Output',Problem_name,'DFOR',parameters.name);
warning('off','all');
mkdir(savedir);           % Making an Output folder with a DFOR subfolder
warning('on','all');
%% Importing the data and defining Setslog
%[DataSet,paraname,DATA] = xlsread(filename);
DataSet=readmatrix(filename);
DATA=readcell(filename);
paraname=DATA(1,:);

Setslog.in_index = in_index;           %input variable cols in datafile
Setslog.out_index = out_index;         %output variable col in datafile
Setslog.DATA = DATA;                   %Cell Array of data in the xls file
Setslog.paraname= paraname;            %Cell Array of the column headings
Setslog.DataSet = DataSet;             %Matrix of the dataset 

%% Defining the function sets of the trees and storing in Setslog
F1_set={'sine' 'expo'};                      %Single Arity Functions
F2_set={'add' 'sub' 'mul' 'rdiv'};           %Double Arity Functions
F3_set={};                                   %Triple Arity Functions
F.set = {F1_set F2_set F3_set};
j = length(F.set);
temp ={}; F.sets = [];
for i = 1:j
    F.sets = [F.sets i*ones(1,length(F.set{i}))];
    temp = [temp F.set{i}];
end
F.set = temp;
Setslog.F = F;
clear count temp i j F1_set F2_set F3_set;

%% Defining terminal sets to be fed to each subnet
n_subnets=length(parameters.DFORtrain.Pop_str);
T=struct('terminals',cell(1,n_subnets));
for i=1:n_subnets
    T1_set=paraname(1,parameters.DFORtrain.Pop_str{i}{1});
    T2_set={'ERC'};                         %Ephemeral random floating point constants(ERC)
    T(i).terminals=[T1_set T2_set];    %T(i).terminals gets the terminals of 1st block of ith subnet
end    
Setslog.T=T;

%% Defining the parameters of DFOR to Setslog
Setslog.max_block_depth=parameters.DFORtrain.max_depth;
Setslog.prey_popsize=parameters.DFORtrain.Prey_popsize;
Setslog.no_Prey_preferred = parameters.DFORtrain.no_Prey_preferred;    %Desired popsize
Setslog.no_new_Prey=parameters.DFORtrain.no_new_Prey;                  %no of new prey every KillInterval 
Setslog.Predator_popsize = parameters.DFORtrain.Predator_popsize;      %Predator pop size
Setslog.generations = parameters.DFORtrain.generations;                %no of genearations
Setslog.ploton = parameters.DFORtrain.ploton;                          %Flag to plot
Setslog.sgen=parameters.DFORtrain.sgen;                                %No of Gen for Single OBJ optim
Setslog.max_rank=parameters.DFORtrain.max_rank;                        %max rank retained at KillInterval
Setslog.n_subnets=n_subnets;                                           %No of subnets
Setslog.KillInterval=parameters.DFORtrain.KillInterval;                %KillInterval
Setslog.no_x = parameters.DFORtrain.no_x;                              %lattice size (no of rows)
Setslog.no_y = parameters.DFORtrain.no_y;                              %lattice size (no of cols)
Setslog.tour_size=parameters.DFORtrain.tour_size;                      %Single Obj Function
Setslog.prey_add=parameters.DFORtrain.prey_add;
blocks=zeros(1,n_subnets);                         
for i=1:n_subnets
    blocks(i)= parameters.DFORtrain.Pop_str{i}{2}(2);
end    
Setslog.n_blocks=blocks;                                               %n_blocks has no of tree blocks in each subnet
Setslog.Pop_str=parameters.DFORtrain.Pop_str;                          %Subnet Info
Setslog.method=parameters.DFORtrain.method;

%% Defining the parameters associated with subsets of training data

set_size =length(DataSet(:,1)); Setslog.set_size = set_size;
subsets =parameters.DFORtrain.subsets;
overlap=parameters.DFORtrain.overlap;
subset_size = round((set_size + (subsets-1)*overlap) / subsets);
Setslog.subset_size = subset_size;

%% Setting the 4 figure handles
for out = out_index
Setslog.out_index = out;
if Setslog.ploton
    figure_handle = [figure(1) figure(2) figure(3) figure(4)];
    scrsz = get(0,'ScreenSize');
    set(figure_handle(1), 'OuterPosition', [0*scrsz(3) scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(1))
    set(figure_handle(2), 'OuterPosition', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(2))
    set(figure_handle(3), 'OuterPosition', [0*scrsz(3) 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(3))
    set(figure_handle(4), 'OuterPosition', [scrsz(3)/2 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(4))
    pause(0.1);
end
warning('off', 'all')
eval(['delete ' savedir '\Y' num2str(out-out_index(1)+1) '.mat']);
warning('on', 'all')

%% Making the different training datasets based on the no of subsets and subset size
if subsets == 1
    training_subsets = 1;
    no_run = 1; Setslog.no_run = no_run;
else
    training_subsets = [eye(subsets,subsets); ones(1,subsets)];%if for ex:subsets=2 then training_subsets:[1 0;0 1;1 1]
    no_run = 1+subsets; Setslog.no_run = no_run;               %We run the training for each subset seperately and one with all together
end

Setslog.chosen_training_sets = training_subsets;
data_index = (1:1:length(DataSet))';

%% Normalization of the X variable between eps and 1 
Xmin=eps;
Xmax=1;
Data_min=zeros(1,length(DataSet(1,:)));
Data_max=zeros(1,length(DataSet(1,:)));
DataSet_sc=zeros(size(DataSet));

for i=1:length(DataSet(1,:))
			Data_min(1,i) = min(DataSet(:,i));
			Data_max(1,i) = max(DataSet(:,i));
			DataSet_sc(:,i)=Xmin+(DataSet(:,i)-Data_min(1,i))/(Data_max(1,i)-Data_min(1,i))*(Xmax-Xmin);
end
Setslog.DataSet_sc=DataSet_sc;
Setslog.Xmin=Xmin;
Setslog.Xmax=Xmax;

%% Running the algorithm
for i = 1:no_run
    if(i==no_run)
        chosen_data_index=data_index;
        chosen_data_index = unique(chosen_data_index);
        Setslog.dataset(i).data_index = chosen_data_index;
        Setslog.dataset(i).in = Setslog.DataSet_sc(chosen_data_index,in_index);
        Setslog.dataset(i).out = Setslog.DataSet(chosen_data_index,Setslog.out_index);
        continue;
    end    
    
from = 1; to = from + subset_size - 1;
    for j = 1:subsets-1
        if training_subsets(i,j) == 1
            chosen_data_index=zeros(subset_size,1);
            chosen_data_index(:,1) = data_index(from:to,:);
            
        end
        from = to - overlap + 1; to = from + subset_size - 1;       % overlap:no of common examples between adjacent subsets
    end                                                             
    
    if training_subsets(i,subsets) == 1
        
        chosen_data_index=zeros(length(data_index(:,1))-from+1,1);
        chosen_data_index(:,1) =data_index(from:length(data_index(:,1)),:);
    end
    chosen_data_index = unique(chosen_data_index);
    Setslog.dataset(i).data_index = chosen_data_index;
    Setslog.dataset(i).in = Setslog.DataSet_sc(chosen_data_index,in_index);
    Setslog.dataset(i).out = Setslog.DataSet(chosen_data_index,Setslog.out_index);
    
end
Setslog
DataSet_sc
%% Training the model no_run times
for i = 1:no_run
    if i < no_run
        fprintf('\nTraining Dataset %d\n\n', i);
    else
        fprintf('\nTraining Whole Dataset\n\n');
    end
    [Prey,FVal]=PP_DFOR(i);
    %pause(0.1)
end
fprintf('\n\nTrainining Subsets\n');
disp(training_subsets);
end
 end
