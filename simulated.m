% Code to analyse the simulatation experiment from Cornish & Markowetz (2014),
% Supplementary code for the paper:
%
% Yuyanyuan Chen, Mingmin Xu, Xiaodan Fan, Cong Pian(2021),
% 'Identifying functional modules using energy minimization
%  with graph cuts'

% We performed 500 simulation runs and obtained the Knode and BioNet 
% output using the R code provided by Cornish & Markowetz (2014), PLOS 
% Computational Biology. The output was exported and saved as a .csv file. 
% Below, the output is read into MATLAB and the NEST, 'Degree', Limma,
% NePhe, Knode, BioNet, MRF, EMGC and EMLC hit lists are calculated.

% Note that running this code for all 500 runs takes ~1926 seconds using a
% Intel(R) Core(TM) with a 1.9 GHz Intel Core i5 processor (quad-core)
% and 8 GB of RAM.

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Load the simulated data

overall_adj = xlsread('overall_adj.csv');
overall_ctr = xlsread('overall_ctr.csv');
overall_exp = xlsread('overall_exp.csv');
overall_truehit = xlsread('overall_hit.csv');
overall_bionetresult = xlsread('overall_bionetresult.csv');
overall_knoderesult = xlsread('overall_knoderesult.csv');
overall_logfc = xlsread('overall_logfc.csv');
overall_pvalue = xlsread('overall_pvalue.csv');
overall_initial = xlsread('overall_initial.csv');


% Preallocate the arrays for the porportion of correctly identified hits
prop_correct_deg = NaN*ones(1000,1);
prop_correct_limma = NaN*ones(1000,1);
prop_correct_nephe = NaN*ones(1000,1);
prop_correct_nest = NaN*ones(1000,1);
prop_correct_knode = NaN*ones(1000,1);
prop_correct_bionet = NaN*ones(1000,1);
prop_correct_MRF = NaN*ones(1000,1);
prop_correct_EMGC = NaN*ones(1000,1);
prop_correct_EMLC = NaN*ones(1000,1);



%% For each of the 500 simulation runs
tic
for run = 1:1000
    
    if (mod(run,10) == 0)
        disp(['processed ',num2str(run), ' simulated runs out of 500']);
    end
    
    % construct the sparse adjacency matrix (Upper triangular matrix)
    m = 3*run-2;
    n = 3*run;
    adj = overall_adj(:,m:n);
    Adj = sparse(adj(:,1),adj(:,2),adj(:,3),1000,1000);
    
    % the vertex number of a simulated network 
    nVertices = length(Adj);
    
    % the vertex degree of a simulated network
    degree = full(sum(Adj + Adj'));
    
    % Allocate the output from the current simulation run 
    exp_mean = overall_exp(:,run);
    pvalue = overall_pvalue(:,run);
    knoderesults = overall_knoderesult(:,run);
    bionetresults = overall_bionetresult(:,run);
    truehits = overall_truehit(:,run);
  
    top = 30;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 
    % Calculate the 'degree', 'p-value', Knode and BioNet hits
    
    % Calculate the hits for 'degree'
    hitsdegree = zeros(1,nVertices);
    [~,ind] = sort(degree,'descend');
    hitsdegree(ind(1:top)) = 1;
    
    % Calculate the hits for 'p-value'
    hitslimma = zeros(1,nVertices);
    [~,ind] = sort(pvalue,'ascend');
    hitslimma(ind(1:top)) = 1;
    
    % Calculate the hits for 'Knode'
    hitsknode = zeros(1,nVertices);
    [~,ind] = sort(knoderesults,'descend');
    hitsknode(ind(1:top)) = 1;
    
    % Calculate the hits for 'BioNet'
    hitsbionet = bionetresults;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 
    % Calculate the NePhe scores and hits
    
    % We consider the shortest paths, average summation NePhe score (the
    % best performing pair of options) as proposed in Wang et al. (2009),
    % BMC Genomics.
    
    % The similarity matrix is 'shortest paths'
    similarity_matrix = exp( -graphallshortestpaths(Adj) );
    similarity_matrix(logical(eye(size(similarity_matrix)))) = 0; 
    
    % Since the NePhe method suggests binary scores as input we take the
    % p-value hits as the inital hits. Note that this actually gives better
    % results than using the scores themselves.
    
    % The summation is simply the average
    nephescore = (hitslimma*similarity_matrix)./sum(similarity_matrix);
    
    % Calculate the hit list for the NePhe method
    hitsnephe = zeros(1,nVertices);
    [~,ind] = sort(nephescore,'descend');
    hitsnephe(ind(1:top)) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 
    % Calculate the NEST scores and hits
    % We also consider the NEST score as proposed by Jiang et al. (2015),
    % Genome Biology.
    % Note that the NEST scores are just one possible combination within
    % the more general NePhe framework, however the difference is that
    % rather than inputting binary scores we input the scores themselves.
    % For both the NEST and NePhe methods using the binary scores actually
    % gives better output.
    
    nestscore = icdf('Normal',pvalue,0,1)'*Adj;
    
    % Calculate the hit list for the NEST method
    hitsnest = zeros(1,nVertices);
    [~,ind] = sort(nestscore,'ascend');
    hitsnest(ind(1:top)) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the MRF score and hits
    % We also consider the MRF score as proposed by Robinson, et al., (2017),
    % Bioinformatics.
    
    % For the MRF method, we used a general exponential
    % density with parameter ¦Ë = 0.1684 that
    % tercepts the standard uniform density at 0.3.
    Numlables = 2;
    totalsizes =length(exp_mean);
    
    % Calculate the unary potentials
    DataCost = NaN*ones(Numlables,totalsizes);
    
    % Define the distributions for each label
    pd_hits = makedist('Exponential','mu',0.1684); 
    pd_nonhits = makedist('Uniform','lower',0,'upper',1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Label 1: the 'hits'
    label = 1;
    
    y = pdf(pd_hits,pvalue);
    DataCost(label,:) = -log(y)';
    % Get rid of any Inf
    max_DataCost = max(DataCost(label, DataCost(label,:) < Inf));
    DataCost(label,DataCost(label,:) == Inf) = max_DataCost;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Label 2: the 'non hits'
    label = 2;
    
    y = pdf(pd_nonhits,pvalue);
    DataCost(label,:) = -log(y)';
    % Get rid of any Inf
    max_DataCost = max(DataCost(label, DataCost(label,:) < Inf));
    DataCost(label,DataCost(label,:) == Inf) = max_DataCost;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Minimise the energy
    
    % Set the resolution for beta
    n = 1000;
    run_beta = linspace(0,50,n+1);
    
    % Pre-allocate the array for the minimum energy labels
    min_energy_labels = NaN*ones(totalsizes,length(run_beta));
    
    for k = 1:length(run_beta)
        
        % Calculate the pairwise potentials
        Neighbours = run_beta(k)*Adj;
        
        % As gco requires integers
        scale = 1000;
        
        % Additionally supress related warnings to increase speed
        warning('off','all')
        
        % Create the gco object
        input_GCO = GCO_Create(totalsizes,Numlables);
        
        % Input the unary potenitals
        GCO_SetDataCost(input_GCO,scale*DataCost);
        
        % Input the form of the pairwise potentials
        SmoothCost = double(~eye(Numlables));
        GCO_SetSmoothCost(input_GCO,SmoothCost);
        
        % Input the pairwise potentials
        GCO_SetNeighbors(input_GCO,scale*Neighbours);
        
        % Minimise the energy using the alpha-expansion algorithm
        GCO_Expansion(input_GCO);
        
        % Obtain the minimum energy labels
        min_energy_labels(:,k) = GCO_GetLabeling(input_GCO);
        
        % Clear the gco object from memory
        GCO_Delete(input_GCO);
        
    end
    
    % Determine the instances before all nodes have the dominant label = 2
    beforedom = (sum(min_energy_labels) < 2*totalsizes);
    
    % Calculate the MRF scores
    MRFscore = sum(min_energy_labels(:,beforedom) == 1,2).*degree';
    
    % Calculate the hit list for the MRF method
    hitsMRF = zeros(1,totalsizes);
    [~,ind] = sort(MRFscore,'descend');
    hitsMRF(ind(1:top)) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

   %%
   % Calculate the MRF score and hits using EMGC method
   
   % Calculate the unary potentials (i.e.datacost) for EMGC method
   ctrs = overall_ctr(:,run);
   datacost1 = zeros(totalsizes,Numlables);
   for i = 1:totalsizes
       for j = 1:Numlables
            datacost1(i,j) = (exp_mean(i) - ctrs(j)).^2;  %¼ÆËãdatacostÏî
       end
   end
    datacost1 = datacost1'; 
    
    
   %%    
   % Set the resolution for beta
   n = 1000;
   upper = 50;
   run_beta = linspace(0,upper,n+1);
   
   % Pre-allocate the array for the minimum energy labels
   min_energy_labels_EMGC = NaN*ones(totalsizes,length(run_beta));
   %%
   for k = 1:length(run_beta)
       
       % Create the gco object
       h = GCO_Create(totalsizes,Numlables);
     
       % Additionally supress related warnings to increase speed
       warning('off','all');
       
       % Input the unary potenitals (datacost1)
       GCO_SetDataCost(h,1000*datacost1);
       
       % Input the form of the pairwise potentials
       SmoothCost = double(~eye(Numlables));
       GCO_SetSmoothCost(h,SmoothCost);
       
       % Calculate the pairwise potentials
       Neighbours = run_beta(k)*Adj; 
       GCO_SetNeighbors(h,1000*Neighbours); 
      
       % Minimise the energy using the alpha-expansion algorithm
       GCO_Expansion(h);
       
       % Obtain the minimum energy labels
       min_energy_labels_EMGC(:,k) = GCO_GetLabeling(h);
       
       % Clear the gco object from memory
       GCO_Delete(h);
   end 


     % Determine the instances before all nodes have the dominant label = 2
     beforedom_EMGC = (sum(min_energy_labels_EMGC) < 2*nVertices);   
     
     % Calculate the MRF scores for the EMGC method
     MRFscore_EMGC = sum(min_energy_labels_EMGC(:,beforedom_EMGC) == 1,2).*degree';
    
     % Calculate the hit list for the EMGC method
     hitsMRF_EMGC = zeros(1,totalsizes);
     [~,ind] = sort(MRFscore_EMGC,'descend');
     hitsMRF_EMGC(ind(1:top)) = 1;
    
     
     
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%
   % Calculate the MRF score and hits using EMLC method

     upper = 50;
     run_beta = linspace(0,upper,n+1);
     min_energy_labels_EMLC = NaN*ones(totalsizes,length(run_beta));
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   for k = 1:length(run_beta)
       
        % Reset previously switched off warnings
        warning('off','all')
        
        % Create the gco object
        h = GCO_Create(totalsizes,Numlables);
        
        % Input the unary potenitals (datacost1)
        GCO_SetDataCost(h,1000*datacost1); 
        
        SmoothCost = double(~eye(Numlables));
        GCO_SetSmoothCost(h,SmoothCost);
        
        % Input the form of the pairwise potentials
        Neighbours =run_beta(k)*Adj;   
        GCO_SetNeighbors(h,1000*Neighbours);
        
        % Input the labelcost 
        GCO_SetLabelCost(h,[1000 30]);
        
        % Minimise the energy using the alpha-expansion algorithm
        GCO_Expansion(h);
        
        % Obtain the minimum energy labels
        min_energy_labels_EMLC(:,k) = GCO_GetLabeling(h);
        
        % clear the GCO object
        GCO_Delete(h);
   end 
  
     beforedom_EMLC = (sum(min_energy_labels_EMLC) < 2*nVertices);   
     % Calculate the MRF scores for the EMLC method
     MRFscore_EMLC = sum(min_energy_labels_EMLC(:,beforedom_EMLC) == 1,2).*degree';
    %Calculate the hit list for the EMLC method
     hitsMRF_EMLC = zeros(1,nVertices);
     [~,ind] = sort(MRFscore_EMLC,'descend');
     hitsMRF_EMLC(ind(1:top)) = 1;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   %%
    % Work out the proportion of hits correctly identified for each method
    
    prop_correct_deg(run) = sum(hitsdegree.*truehits')/top;
    prop_correct_limma(run) = sum(hitslimma.*truehits')/top;
    prop_correct_nephe(run) = sum(hitsnephe.*truehits')/top;
    prop_correct_nest(run) = sum(hitsnest.*truehits')/top;
    prop_correct_knode(run) = sum(hitsknode.*truehits')/top;
    prop_correct_bionet(run) = sum(hitsbionet.*truehits)/top;
    prop_correct_MRF(run) = sum(hitsMRF.*truehits')/top;
    prop_correct_EMGC(run) = sum(hitsMRF_EMGC.*truehits')/top;
    prop_correct_EMLC(run) = sum(hitsMRF_EMLC.*truehits')/top;

    
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Plot the box and line plots of proportion of hits identified

figure
boxplot([prop_correct_nest,prop_correct_deg,prop_correct_limma,prop_correct_nephe,...
    prop_correct_knode,prop_correct_bionet,prop_correct_MRF,prop_correct_EMGC,prop_correct_EMLC],...
    'color','k', 'symbol', 'k*')

set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')

set(gca,'xtick',1:10)
set(gca,'xticklabel',{'NEST', 'Degree', 'Limma', 'NePhe', 'Knode', 'BioNet', 'MRF','EMGC','EMLC'})

set(gca,'ygrid','on')

ylabel('Proportion of nodes identified')

figure
plot([prop_correct_nest,prop_correct_deg,prop_correct_limma,prop_correct_nephe,...
    prop_correct_knode,prop_correct_bionet,prop_correct_MRF,prop_correct_EMGC,prop_correct_EMLC]',...
    '-*')

set(gca,'xtick',1:10)
set(gca,'xticklabel',{'NEST', 'Degree', 'Limma', 'NePhe', 'Knode', 'BioNet','MRF', 'EMGC','EMLC'})

set(gca,'ygrid','on')

ylabel('Proportion of nodes identified')