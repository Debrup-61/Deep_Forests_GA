function [Prey,FVal]=PP_DFOR(set_no)

% This function creates a population of deep forests and runs the predator 
% prey genetic algorithm on them to find the pareto optimal points
% corresponding to 2 objectives-complexity and accuracy of an individual 
% forest.

% set_no is the index of the run.

%% PPGA parameters
global generation;         % the current generation
global lattice no_x no_y;  % the lattice and the size of it
global Setslog run_no;     % Setslog and run_no
global max_rank KillInterval;

%% General GP parameters
global F T in_index out_index 
global max_block_depth max_arity min_arity
global F_bad         %Fitness assigned when nan or inf is encountered in tree fitness evaluation
global prey_popsize
global Predator_popsize 
global n_blocks
global lambda
global move_prey
global n_gen         % Number of generations
global sgen
global no_new_Prey
global no_Prey_preferred
global prey_add
%% Tree Initialization parameter 

global choose_term     %Probablity of choosing a terminal node during initialization using grow method
global choose_erc      %Given we pick a terminal probablity of choosing a erc

%% Single Objective Optimization Parameters

% Seperate Probablities during single and multi-obj crossovers since in
% single we want crossover to be lower in order to prevent loss of diversity in individuals

global Ps_xover Ps_xstand Ps_xfair
Ps_xover=0.60;
Ps_xstand = 0.5; Ps_xstand = Ps_xover*Ps_xstand;    %fraction of standard xover
Ps_xfair = Ps_xover - Ps_xstand;                    %fraction of Height fair xover

global Pms_stand Pms_small Pms_mono Ps_mut
Ps_mut = 0.15;                                      %prob of mutation of a Prey or Individual
Pms_stand = 0.45; Pms_stand = Pms_stand*Ps_mut;     %fraction of standard mutation
Pms_small = 0.15; Pms_small = Pms_small*Ps_mut;     %fraction of small mutation
Pms_mono = Ps_mut - Pms_stand - Pms_small;          %fraction of monoparental mutation

%% Crossover parameters
global Px_stand Px_fair P_xover 
global choose_xfunc
P_xover=0.9;
Px_stand = 0.5; Px_stand = P_xover*Px_stand;    %fraction of standard xover
Px_fair = P_xover - Px_stand;                   %fraction of Hieght fair xover

%% Mutation parameters                           

global Pm_stand Pm_small Pm_mono P_mut
P_mut = 0.35;                                 %prob of mutation of a Prey or Individual
Pm_stand = 0.45; Pm_stand = Pm_stand*P_mut;   %fraction of standard mutation
Pm_small = 0.15; Pm_small = Pm_small*P_mut;   %fraction of small mutation
Pm_mono = P_mut - Pm_stand - Pm_small;        %fraction of monoparental mutation

%% Plotting variables

avg_rms_err=zeros(1,n_gen);
avg_nodes=zeros(1,n_gen);
avg_height=zeros(1,n_gen);
avg_complexity=zeros(1,n_gen);
best_rms_err=zeros(1,n_gen);
best_complexity=zeros(1,n_gen);
CR=ones(1,n_gen);
n_pop=zeros(1,n_gen);


%% Defining the global parameters
choose_term=0.15;         % Probability of choosing a terminal node during initialization
choose_erc=0.10;
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
Predator_popsize = Setslog.Predator_popsize;  %Number of Predators
n_blocks=Setslog.n_blocks;
lambda=0.5;
Setslog.lambda=lambda;
move_prey=0.5;
max_rank=Setslog.max_rank;
KillInterval=Setslog.KillInterval;
n_gen=Setslog.generations;
no_Prey_preferred=Setslog.no_Prey_preferred ;
no_new_Prey =Setslog.no_new_Prey;
choose_xfunc=0.9;                     %Probability of choosing a functional node for crossover
sgen=Setslog.sgen;
prey_add=Setslog.prey_add;

h_dist=zeros(prey_popsize,n_gen);

%% Initializing the lattice with zeros
FVal=[];
lattice = zeros(no_x+2,no_y+2);
fitness=zeros(prey_popsize,n_gen);  % For fitness distribution
nodes=zeros(prey_popsize,n_gen);    % To look at nodes distribution

%% Initializing the prey population and getting its fitness 
[Prey,FVal]=create_population(prey_popsize);
n_pop(1,1)= prey_popsize;

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

MovePrey(0, 0, 0);   % MovePrey(0,0,0) makes sure that the ends of the lattice are connected (left-right,top-bottom)                                                

% Start of Generations
 for generation=1:n_gen
   avg_rms_err(generation)=mean([Prey.error]);
   avg_nodes(generation)=mean([Prey.nodes]);
   avg_height(generation)=mean([Prey.avg_subnet_ht]);
   avg_complexity(generation)=mean([Prey.complexity]);
   best_rms_err(generation)=min(FVal(:,1));
   best_complexity(generation)=min([Prey.complexity]);
   n_pop(1,generation)=length(unique([Prey.error]));
   fitness(:,generation)=sort(FVal(:,1),1);
   nodes(:,generation)=(sort(([Prey.nodes]),2))';
   h_dist(:,generation)=[Prey.avg_subnet_ht]';
   
   if(generation>=2)
    CR(1,generation)=best_rms_err(generation)/best_rms_err(generation-1);
    disp(CR);
   end 
   %Pop_size = length(Prey);
   
   %Prey_new = [];
  if(generation<=sgen)          %%Evolution strategy and parental and survivor selection
    
    [Prey,FVal]=Evol_Elite(Prey,FVal,Setslog.method);  % Performs xovers and mutations
    [Prey,FVal]=tournament_select(Prey,FVal);          % Tournament selection
    [fonrank, ~] = NONDOM_SORT([FVal]); 
    %disp(Prey);
    %disp(FVal);
    PlotPareto(FVal,fonrank,generation);
    
%     if(rem(generation,10)==0)
%         [Prey_add,FVal_add]=create_population(prey_add);
%         x=[Prey.err];
%         [~,pos]=maxk(x,prey_add);   % Replacing worst prey_add preys with new population
%         Prey(1,pos)=Prey_add;
%         FVal(pos,:)=FVal_add;
%     end    
% %  else    
% %  % Movement of preys one step based on (rand<move_prey)- given 10 chances
% %    for Pop_index = 1 : Pop_size
% %         if rand < move_prey
% %             %Identify location
% %             [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == Pop_index);
% %             xpos = xpos+1; ypos = ypos+1;
% %             for j = 1:10     %10 trials for free spot
% %                 dx = round(rand*(3-eps)-1.5);                                     %-1 to the left, 0 no move, 1 to the right
% %                 dy = round(rand*(3-eps)-1.5);                                     %-1 down, 0 no move, 1 up
% %                 if lattice(xpos+dx,ypos+dy) == 0
% %                     lattice(xpos,ypos) = 0;                                       %removing prey i
% %                     MovePrey(xpos+dx, ypos+dy, Pop_index);                        %Move prey to new position
% %                     break
% %                 end
% %             end
% %         end
% %     end
% % 
% % %  Breeding of the preys
% %   
% %     for Pop_index = 1 : Pop_size
% %         [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == Pop_index);
% %         xpos = xpos+1; ypos = ypos+1;
% %         moore = lattice(xpos-1:xpos+1,ypos-1:ypos+1);                             % Moore's neighbourhood of a prey
% %         %Remove the i-prey
% %         moore(2,2) = 0;
% %         [matex,matey] = find(moore >= 1 & moore <= Pop_size);
% %         if ~isempty(matex)
% %             comp=[];
% %             parent2 = ceil(rand*length(matex));
% %             for i=1:length(matex)
% %                 parent2 = lattice(xpos-2+matex(i),ypos-2+matey(i));       % For Crossover select the best prey in the moore's
% %                 e=Prey(parent2).err;c=Prey(parent2).complexity;           % neighbourhood. Form of Elitism-xover of good preys
% %                 comp=[comp;parent2,e,c];                                  % error given 1st priority then complexity
% %             end    
% %             comp=sortrows(comp,[2 3]);
% %             parent2=comp(1,1);
% %             %Crossover between two forest individuals
% %             [Offsprng, FVal_offsp] = create_offspring(Pop_index,parent2,Prey); 
% %             
% %             %Random placement of offsprings /10 trials
% %             for l = 1:2
% %                 for j=1:10 %10 trial for free spot
% %                     xpos = round(rand(1,1)*(no_x-1)+1)+1;
% %                     ypos = round(rand(1,1)*(no_y-1)+1)+1;
% %                     if lattice(xpos,ypos)==0
% %                         Prey = [Prey Offsprng(l)];                               %New offspring given ten chances to place in
% %                         FVal = [FVal ; FVal_offsp(l,:)];                         %any position of the lattice.Offspring added at
% %                         lattice(xpos,ypos) = length(Prey);                       %bottom of old prey structure.
% %                         break
% %                     end
% %                 end
% %             end
% %          end
% %         MovePrey(0, 0, 0);
% %    end
% %   
% %    [fonrank, front] = NONDOM_SORT([FVal]);                                         %Both fonrank and front are of size of pop
% %                                                                                    %have the fonseca rank and front no ,resp.
% %                                                                                 
% % % Removing Preys every KillInterval with (Fonseca rank> max_rank)
% % 
% %     if generation/KillInterval == round(generation/KillInterval)
% %         if generation < n_gen
% %              %Generate new prey to balance the population
% %              [Prey_new, FVal_new] = create_population(no_new_Prey);
% %          end
% %         indfr = find(fonrank > max_rank);
% %         FVal(indfr,:) = F_bad+eps;
% %         [Prey ,FVal] = KillBadPrey(Prey, FVal);  % Erases the prey's with fitness equal to F_bad
% %         [fonrank,front] = NONDOM_SORT([FVal]);   % Perform ranking after KillInterval
% %     
% %     end
% %     
% % % Killing of preys by predators
% % 
% %   crodit = CROW_SORT([FVal], front); 
% %   FValKill = [];
% %   FValKill(:,1) = (FVal(:,1) - min(FVal(:,1)))./(max(FVal(:,1)) - min(FVal(:,1)));
% %   FValKill(:,2) = (FVal(:,2) - min(FVal(:,2)))./(max(FVal(:,2)) - min(FVal(:,2)));
% % 
% %   PredMoves = floor((length(Prey) - no_Prey_preferred)/Predator_popsize);
% %   fprintf('\nGeneration %i: Predatorpop %i PredMoves %i\n',generation,length(Predators(:,1)),PredMoves);
% %   fprintf('Preypop before: %i; Preypop after: ', length(Prey))
% % 
% %     for i = 1:Predator_popsize
% %         for k = 1:PredMoves
% %             [xpos, ypos] = find(lattice(2:no_x+1,2:no_y+1) == -i);
% %             xpos = xpos+1; ypos = ypos+1;
% %             [matex,matey] = find(lattice(xpos-1:xpos+1,ypos-1:ypos+1) > 0);
% %             if length(matex) > 1 %prey available
% %                 rows=[];
% %                 for t = 1:length(matex)
% %                     rows = [rows; lattice(matex(t)+xpos-2,matey(t)+ypos-2)];
% %                 end
% %                 %Calculating Fitness Value for all Prey near Predator i, with Elitism
% %                 f = (FValKill(rows,:)*[Predators(i);1-Predators(i)]).*(front(rows)-1);
% % 
% %                 %Using Crowding in case all Prey have same fitness
% %                 if length(unique(front(rows))) == 1
% %                     f = (f+1)./(1+crodit(rows));
% %                 end
% %                 
% %                 %Killing Prey
% %                 [~,pos] = max(f);
% %                 j = rows(pos(1)); lattice(xpos,ypos) = 0; %removing predator i
% %                 [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == j);
% %                 xpos = xpos+1; ypos = ypos+1;
% %                 Prey = MovePredator(Prey, xpos, ypos, i, j);
% %                 FValKill(j,:) = []; FVal(j,:) = [];
% %                 front(j) = []; crodit(j) = []; 
% %                 fonrank(j) = [];
% %           else    %Only move
% %                   for j=1:10 %10 trial for free spot
% %                      dx=round(rand*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
% %                      dy=round(rand*(3-eps)-1.5); %-1 down, 0 no move, 1 up
% %                      if lattice(xpos+dx,ypos+dy)==0
% %                          lattice(xpos,ypos)=0; %removing predator i
% %                          Prey = MovePredator(Prey, xpos+dx, ypos+dy, i, inf);
% %                          break
% %                      end
% %                   end
% %             end
% %             
% %         end
% %         
% %     end
% % %Plotting Pareto Front and adding Prey_new if it is KillInterval 
% %     fprintf('%i\n' , length(Prey))
% %     PlotLattice 
% %     PlotPareto(FVal, fonrank,generation)  
% % 
% % %Random placement of New Prey /10 trials
% % 
% %     for i = 1:length(Prey_new)
% %         for k = 1:10 %10 trial for free spot
% %             [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
% %             if ~isempty(emptyx > 0)
% %                 j = ceil(rand*length(emptyx));
% %                 lattice(emptyx(j)+1,emptyy(j)+1) = length(Prey) + 1;
% %                 Prey = [Prey Prey_new(i)];
% %                 FVal = [FVal; FVal_new(i,:)];
% %                 break
% %             end
% %         end
% %     end
% % 
% %   
 end
n_pop(generation+1)=length(Prey);
 end   

%% PLOTTING THE TREND LINE OF PREDICTED DATA VS EXPERIMENTAL DATA
figure(5);
[~,i]=min([Prey.error]);
res(:,1)=Setslog.dataset(set_no).out;
res(:,2)=Prey(i).output;
res=sortrows(res,1);
x=res(:,1);
y=res(:,2);
c=polyfit(x,y,1);
y_fit=polyval(c,x);
plot(x,y_fit,"r--");
xlabel('Actual Value')
ylabel('Predicted Value')
title(["Trend line of Predicted Data vs Actual Data" "Slope of Fitted line:" num2str(c(1)) "Y-intercept" num2str(c(2))]);
%text(100,100,["Slope of Fitted Line:" num2str(c(1))],'fontsize', 30);
hold off;

%% Plotting the Rate of Convergence over generations
figure(6);
plot(1:n_gen,CR,"b-o");
xlabel('generation no');
ylabel('Convergence Rate');
title('Convergence Rate Vs Generations');

%% Have a visualisation of the distribution of fitness over generations


% Saving y_fit
save plotDFOR_vars.mat h_dist avg_height avg_rms_err avg_nodes  avg_complexity best_rms_err best_complexity CR n_pop fitness nodes y_fit
end
