%%%% Loop through LSCMutBiasGeneral

clear
close all

%% Import values

k = 5;
        
% Import parameters sorted by mutational probability

tbl = readtable('MaxitinibParametersNoT315I_111919.csv','TreatAsEmpty','NA');
biasR = tbl.prob(2:end)';
alpha = tbl.alpha;

alpha_sen = alpha(1);
alpha_res = alpha(2:end);

n = length(biasR);

%% Loop through simulations

tic

for i = 1:(n-k+1)

    sen_idx = i:(i+k-1);

    alpha_res_k = alpha_res;
    alpha_res_k(sen_idx) = alpha_sen;

    alphas = [alpha_sen alpha_res_k'];

    output = Maxitinib_SL111919(alphas,biasR);

    filename = strcat('MaxitinibK',num2str(i,'%02.f'),'Full012719.csv');

    csvwrite(filename,output);

end 
toc

% Run simulations for original alpha values
alphas = alpha';

output = Maxitinib_SL111919(alphas,biasR);

filename = strcat('MaxitinibNoT315I_K',num2str(0,'%02.f'),'Full012719.csv');
csvwrite(filename,output);

