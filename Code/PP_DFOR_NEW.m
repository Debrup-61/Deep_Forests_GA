function Prey = PP_DFOR_NEW(set_no)
%PP_DFOR_NEW Summary of this function goes here
%   Detailed explanation goes here
%% PPGA parameters
global generation;         % the current generation
global lattice no_x no_y;  % the lattice and the size of it
global Setslog run_no;     % Setslog and run_no

%% General GP parameters
global F T in_index out_index 
global max_block_depth max_arity min_arity
global F_bad                                   %Fitness assigned when nan or inf is encountered in tree fitness evaluation
global prey_popsize
global Predator_popsize 
global n_blocks
global lambda
global move_prey
global n_gen                                   % Number of generations

%% Tree Initialization parameter 
global choose_term                             %Probablity of choosing a terminal node during initialization using grow

%% Crossover parameters
global Px_stand Px_fair P_xover
global choose_xfunc
P_xover=1;
Px_stand = 0.5; Px_stand = P_xover*Px_stand;    %fraction of standard xover
Px_fair = P_xover - Px_stand;                   %fraction of Hieght fair xover

%% Mutation parameters
global Pm_stand Pm_small Pm_mono P_mut
P_mut = 0.35;                                 %prob of choosing a prey for mutation
Pm_stand = 0.45; Pm_stand = Pm_stand*P_mut;   %fraction of standard mutation
Pm_small = 0.15; Pm_small = Pm_small*P_mut;   %fraction of small mutation
Pm_mono = P_mut - Pm_stand - Pm_small;        %fraction of monoparental mutation

%% %% Defining the global parameters
choose_term=0.2;
run_no=set_no;
F_bad=1e8;
F=Setslog.F;
T=Setslog.T;
in_index=Setslog.in_index;
out_index=Setslog.out_index;
max_block_depth=Setslog.max_block_depth;
max_arity=max(F.sets);
min_arity=min(F.sets);
no_x=Setslog.no_x;
no_y=Setslog.no_y;
prey_popsize=Setslog.prey_popsize;
Predator_popsize = Setslog.Predator_popsize;    %Number of Predators
n_blocks=Setslog.n_blocks;
lambda=0.5;
Setslog.lambda=lambda;
move_prey=0.5;
n_gen=Setslog.generations;
choose_xfunc=0.9;

%% Initializing the lattice with zeros
FVal=[];
lattice = zeros(no_x+2,no_y+2);
%% Initializing the prey population and getting its fitness 
[Prey,FVal]=create_population(prey_popsize);

%% Placement of Prey Population
for i = 1:length(Prey)
    [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
    if ~isempty(emptyx > 0)
        j = ceil(rand*(length(emptyx)));
        lattice(emptyx(j)+1,emptyy(j)+1) = i;
    end
end
%% Initializing the predator population

Predators = (0:1/(Predator_popsize -1):1)';                          %value to weight F1- the error
%location of predators
for i = 1:length(Predators(:,1))
    [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
    if ~isempty(emptyx > 0)
        j = ceil(rand*(length(emptyx)));
        lattice(emptyx(j)+1,emptyy(j)+1) = -i;                       %Predator identified with neg values
    end
end

MovePrey(0, 0, 0);              % MovePrey(0,0,0) makes sure that the ends of the lattice are connected (left-right,top-bottom)                                                

%% Start of Generations
for generation=1:n_gen
    Pop_size = length(Prey);
    Prey_new = [];













end

