%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to output the toy example as given in Fig. 2.
%
% Supplementary code for the paper:
%
% Yuyanyuan Chen, Mingmin Xu, Xiaodan Fan, Cong Pian(2021),
% 'Identifying functional modules using energy minimization
%  with graph cuts'

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set the seed for the simulated observations (seed = 211 for paper)
seed = 211;
rng('default')
rng(seed)

% Total number of vertices
nVertices = 20;

% Build the base of the adjacency matrix
adj = double(~eye(nVertices));

adj(1:19,20) = 0;
adj(20,1:19) = 0;

% Add an edge between the bridge vertices (10 and 11)
adj(10,20) = 1;
adj(20,10) = 1;

figure
spy(adj)
title('Adjacency matrix')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the true labels and simulate the observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nLabels = 2;

% Vertices 1,2,...,19 have true 'blue' label = 1 and
% vertices 20 have true 'red' label = 2
true_labels = ones(1,nVertices);
true_labels(20) = 2;

% The standard deviation for the distributions
sd = 2;

% The distribution for the 'blue' label (centred at -1)
pd_blue = makedist('Normal','mu',-1,'sigma',sd);

% The distribution for the 'red' label (centred at 1)
pd_red = makedist('Normal','mu',1,'sigma',sd);

% Plot the densities for each label
figure
hold on
x = linspace(-6,6,100);
plot(x, pdf(pd_blue,x), 'b-')
plot(x, pdf(pd_red,x), 'r-')
plot([0 0], ylim,'k--')
title('Densities corresponding to each label')

% Simulate the observations
obs_values = NaN*ones(1,nVertices);

r = random(pd_blue,nVertices,1);
obs_values(true_labels == 1) = r(true_labels == 1);

r = random(pd_red,nVertices,1);
obs_values(true_labels == 2) = r(true_labels == 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the unary potentials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DataCost = NaN*ones(nLabels,nVertices);

% Label 1: the 'blue' vertices
label = 1;
center1 = mean(obs_values(true_labels == 1));
DataCost(label,:) = (obs_values-center1).^2;

% Label 2: the 'red' vertices
label = 2;
center2 = mean(obs_values(true_labels == 2));
DataCost(label,:) = (obs_values-center2).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimise the energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the resolution for beta
n = 10;
run_beta = linspace(0,1.35,n+1);

% Pre-allocate the array for the minimum energy labels
min_energy_labels = NaN*ones(nVertices,length(run_beta));
H = NaN*ones(1,length(run_beta));
for k = 1:length(run_beta);
    % Calculate the pairwise potentials
    Neighbours = run_beta(k)*adj;
    
    % As gco requires integers
    scale = 1000;
    % Additionally supress related warnings to increase speed
    warning('off','all')
    
    % Create the gco object
    input_GCO = GCO_Create(nVertices,nLabels);
    
    % Input the unary potenitals
    GCO_SetDataCost(input_GCO,scale*DataCost);
    
    % Input the form of the pairwise potentials
    SmoothCost = double(~eye(nLabels));
    GCO_SetSmoothCost(input_GCO,SmoothCost);
    
    % Input the pairwise potentials
    GCO_SetNeighbors(input_GCO,scale*Neighbours);
    
    % Input the label cost
    GCO_SetLabelCost(input_GCO,[1 scale*19]);

    % Minimise the energy using the alpha-expansion algorithm
    GCO_Expansion(input_GCO);
    
    % Obtain the minimum energy labels
    min_energy_labels(:,k) = GCO_GetLabeling(input_GCO);
    [E D S L] = GCO_ComputeEnergy(input_GCO);
    H(1,k) =E;
    
    % Clear the gco object from memory
    GCO_Delete(input_GCO);  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the minimum energy labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
imagesc(min_energy_labels)
title('Minimum energy labels')
xlabel('beta')
ylabel('Vertex and observed value')
colormap jet
set(gca,'XTick',(0:10)+0.5,'XTickLabel',run_beta)
yticklabel = cell(1,nVertices);
for k = 1:20
    yticklabel(k) = {[num2str(k),': ',num2str(round(obs_values(k)*100)/100)]};
end
set(gca,'YTick',(0:nVertices-1)+0.5,'YTickLabel',yticklabel)
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the MRF scores and save for cytoscape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

degree = full(sum(adj));

bluescore = sum(min_energy_labels == 1,2).*degree';
redscore = sum(min_energy_labels == 2,2).*degree';

MRFscore = redscore - bluescore;

% Additionally scale the scores to be between -1 and 1
MRFscore = (MRFscore - min(MRFscore))/(max(MRFscore) - min(MRFscore))*2 - 1;

% Write the edge information to file for input into cytoscape
fileID = fopen('toyexampleedges4cytoscape.txt','w');
[I,J] = find(triu(full(adj)));
fprintf(fileID,'%s;%s;%s\n','Node1','Node2','Edge');
for k = 1:length(I)
    fprintf(fileID,'%i;%i;%i\n',I(k),J(k),1);
end
fclose(fileID);

% Write the vertex information to file for input into cytoscape
fileID = fopen('toyexamplevertices4cytoscape.txt','w');
fprintf(fileID,'%s;%s;%s\n','Vertex','Observed value','MRF score');
for k = 1:nVertices
    fprintf(fileID,'%i;%i;%i\n',k,obs_values(k),MRFscore(k));
end
fclose(fileID);


