%% Plot results of experiments
% Analyze and plot the output of the specified experiments for presentation
% and project report.

% turn on save_figs to save plots to file

%% Load the four experiments

au = load('data_scaling_uniform.mat');  % dAta scaling, Uniform
an = load('data_scaling_nonuniform.mat');  % dAta scaling, Nonuniform
iu = load('dim_scaling_uniform.mat');  % dImension scaling, Uniform
in = load('dim_scaling_nonuniform.mat');  % dImension scaling, Nonuniform

%% Export graphics for report/presentation?
save_figs = false;

%% Plot factorization + log-det times for uniform methodologies
fig = figure;
loglog(...
    iu.problem_sizes,iu.default_times(:,1),'-o',...
    iu.problem_sizes,iu.toeplitz_times(:,1),'-o',...
    iu.problem_sizes,iu.rskelf_times(:,1),'-o',...
    iu.problem_sizes,iu.hodlr_times(:,1),'-o',...
    iu.problem_sizes,iu.toeplitz_hodlr_times(:,1),'-o',...
    'linewidth',3);
legend('Default','Toeplitz','RS','HODLR','T-HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Dimension Scaling, Uniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'iu_factortimes.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    au.problem_sizes,au.default_times(:,1),'-o',...
    au.problem_sizes,au.toeplitz_times(:,1),'-o',...
    au.problem_sizes,au.rskelf_times(:,1),'-o',...
    au.problem_sizes,au.hodlr_times(:,1),'-o',...
    au.problem_sizes,au.toeplitz_hodlr_times(:,1),'-o',...
    'linewidth',3);
legend('Default','Toeplitz','RS','HODLR','T-HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Data Scaling, Uniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'au_factortimes.pdf','ContentType','vector');
end

%% Plot factorization + log-det times for nonuniform methodologies

fig = figure;
loglog(...
    in.problem_sizes,in.default_times(:,1),'-o',...
    in.problem_sizes,in.rskelf_times(:,1),'-o',...
    in.problem_sizes,in.hodlr_times(:,1),'-o',...
    'linewidth',3);
legend('Default','RS','HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Dimension Scaling, Nonuniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'in_factortimes.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    an.problem_sizes,an.default_times(:,1),'-o',...
    an.problem_sizes,an.rskelf_times(:,1),'-o',...
    an.problem_sizes,an.hodlr_times(:,1),'-o',...
    'linewidth',3);
legend('Default','RS','HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Data Scaling, Nonuniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'an_factortimes.pdf','ContentType','vector');
end

%% Plot time per CG iteration for uniform methodologies

fig = figure;
loglog(...
    iu.problem_sizes,iu.default_times(:,3)./iu.default_iters,'-o',...
    iu.problem_sizes,iu.toeplitz_times(:,3)./iu.toeplitz_iters,'-o',...
    iu.problem_sizes,iu.rskelf_times(:,2)./iu.rskelf_iters,'-o',...
    iu.problem_sizes,iu.hodlr_times(:,2)./iu.hodlr_iters,'-o',...
    iu.problem_sizes,iu.toeplitz_hodlr_times(:,2)./iu.toeplitz_hodlr_iters,'-o',...
    'linewidth',3);
legend('Default','Toeplitz','RS','HODLR','T-HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Dimension Scaling, Uniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Time per CG Iteration')

if save_figs
    exportgraphics(fig, 'iu_itertimes.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    au.problem_sizes,au.default_times(:,3)./au.default_iters,'-o',...
    au.problem_sizes,au.toeplitz_times(:,3)./au.toeplitz_iters,'-o',...
    au.problem_sizes,au.rskelf_times(:,2)./au.rskelf_iters,'-o',...
    au.problem_sizes,au.hodlr_times(:,2)./au.hodlr_iters,'-o',...
    au.problem_sizes,au.toeplitz_hodlr_times(:,2)./au.toeplitz_hodlr_iters,'-o',...
    'linewidth',3);
legend('Default','Toeplitz','RS','HODLR','T-HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Data Scaling, Uniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Time per CG Iteration')

if save_figs
    exportgraphics(fig, 'au_itertimes.pdf','ContentType','vector');
end

%% Plot time per CG iteration for nonuniform methodologies

fig = figure;
loglog(...
    in.problem_sizes,in.default_times(:,3)./in.default_iters,'-o',...
    in.problem_sizes,in.rskelf_times(:,2)./in.rskelf_iters,'-o',...
    in.problem_sizes,in.hodlr_times(:,2)./in.hodlr_iters,'-o',...
    'linewidth',3);
legend('Default','RS','HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Dimension Scaling, Nonuniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Time per CG Iteration')

if save_figs
    exportgraphics(fig, 'in_itertimes.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    an.problem_sizes,an.default_times(:,3)./an.default_iters,'-o',...
    an.problem_sizes,an.rskelf_times(:,2)./an.rskelf_iters,'-o',...
    an.problem_sizes,an.hodlr_times(:,2)./an.hodlr_iters,'-o',...
    'linewidth',3);
legend('Default','RS','HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Data Scaling, Nonuniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Time per CG Iteration')

if save_figs
    exportgraphics(fig, 'an_itertimes.pdf','ContentType','vector');
end

%% Plot comparisons between uniform/stationary and nonuniform/nonstationary experiments

fig = figure;
loglog(...
    iu.problem_sizes,iu.default_times(:,1),'-o',...
    au.problem_sizes,au.default_times(:,1),'-o',...
    in.problem_sizes,in.default_times(:,1),'-o',...
    an.problem_sizes,an.default_times(:,1),'-o',...
    'linewidth',3);
legend('DimUnif','DataUnif','DimNon','DataNon','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Default SKI Comparison','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'default_compare.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    iu.problem_sizes,iu.rskelf_times(:,1),'-o',...
    au.problem_sizes,au.rskelf_times(:,1),'-o',...
    in.problem_sizes,in.rskelf_times(:,1),'-o',...
    an.problem_sizes,an.rskelf_times(:,1),'-o',...
    'linewidth',3);
legend('DimUnif','DataUnif','DimNon','DataNon','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('RS SKI Comparison','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'rskelf_compare.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    iu.problem_sizes,iu.hodlr_times(:,1),'-o',...
    au.problem_sizes,au.hodlr_times(:,1),'-o',...
    in.problem_sizes,in.hodlr_times(:,1),'-o',...
    an.problem_sizes,an.hodlr_times(:,1),'-o',...
    'linewidth',3);
legend('DimUnif','DataUnif','DimNon','DataNon','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('HODLR SKI Comparison','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'hodlr_compare.pdf','ContentType','vector');
end
