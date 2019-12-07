%%% Loop through ATSims
%%% Scott Leighow - 02/25/19

clear

pars = csvread('ATParameters113019.csv');

log_treat = 10;
log_remain = 6;

bias_vec = [0 1];

for i = 1:2 
    
    bias = bias_vec(i);
    if bias
        bias_str = '_WBias';
        pars_i = pars;
    else
        bias_str = '_WoBias';
        pars_i = [repmat(1/size(pars,1),size(pars,1),1) pars(:,2)];
    end
    
    simout = ATSims112519(10^log_treat,10^log_remain,pars_i);

    filename = strcat('ATSims',bias_str,'_113019.csv');
    csvwrite(filename,simout);
end