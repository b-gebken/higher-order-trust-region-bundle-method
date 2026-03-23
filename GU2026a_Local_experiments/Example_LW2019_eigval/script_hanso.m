% Script that applies HANSO to the eigenvalue problem. To run it, download
% HANSO from https://cs.nyu.edu/~overton/software/hanso/ and extract it so
% that the folder "hanso3_0" is in the same folder as this script.

% Flag for saving the result of HANSO
save_flag = 0;

%% Define the problem

addpath('LW2019_eigval');
n = 50;
m = 25;

rng("default")

A_tensor = zeros(n+1,m,m);
for i = 1:n+1
    tmp = 2*rand(m) - 1;
    A_tensor(i,:,:) = (1/2)*(tmp + tmp');
end

%% Prepare HANSO input

clear options

pars.nvar = n;
pars.f = @(x) lw2019_eigval_f(x,A_tensor);
pars.grad_f = @(x) lw2019_eigval_grad(x,A_tensor);
pars.fgname = 'hanso_objective';

options.x0 = zeros(n,1);

addpath("hanso3_0");
reps = 10;

disp(['Running HANSO ',num2str(reps),' times for the average...'])
total_time = 0;
for i = 1:reps
    start_tic = tic;
    [x_hanso,f_hanso,~,~,~,~,~,~,~,fevalrec_hanso,xrec_hanso] = bfgs(pars,options);
    total_time = total_time + toc(start_tic);
end
runtime_hanso = total_time/reps;
disp(['Average runtime = ',num2str(runtime_hanso)]);

%% Process and save results

f_hanso = min([f_hanso,[fevalrec_hanso{:}]]);

calls_per_iter = cellfun(@(in) size(in,2),fevalrec_hanso);

cumul_oracle_calls_hanso = zeros(1,size(xrec_hanso,2));
for j = 1:size(xrec_hanso,2)
    cumul_oracle_calls_hanso(j) = sum(calls_per_iter(1:j));
end

figure
plot(cumul_oracle_calls_hanso,log10(vecnorm(xrec_hanso - x_hanso,2,1)),'k.-')

if(save_flag)
    cDT = datetime('now');
    cDT.Format = 'dd-MMM-yyyy_HH_mm_ss';
    filename = strcat('hanso_sol_',string(cDT));

    save(filename,'x_hanso','f_hanso','fevalrec_hanso','xrec_hanso','cumul_oracle_calls_hanso','runtime_hanso');
end