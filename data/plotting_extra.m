%% Plot results of experiments
% Analyze and plot the output of the auxiliary experiments for presentation
% and project report

% turn on save_figs to save plots to file

%% Load the four experiments

u = load('dim_scaling_uniform_extra.mat');
n = load('dim_scaling_nonuniform_extra.mat');

%% Export graphics for report/presentation?
save_figs = false;

%% Plot factorization + log-det times for uniform methodologies
fig = figure;
loglog(...
    u.problem_sizes,u.default_times,'-o',...
    u.problem_sizes,u.toeplitz_times,'-o',...
    u.problem_sizes,u.rskelf_times,'-o',...
    u.problem_sizes,u.hodlr_times,'-o',...
    u.problem_sizes,u.toeplitz_hodlr_times,'-o',...
    'linewidth',3);
legend('Default','Toeplitz','RS','HODLR','T-HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Dimension Scaling, Uniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'u_factortimes_extra.pdf','ContentType','vector');
end

%% Plot factorization + log-det times for nonuniform methodologies

fig = figure;
loglog(...
    n.problem_sizes,n.default_times,'-o',...
    n.problem_sizes,n.rskelf_times,'-o',...
    n.problem_sizes,n.hodlr_times,'-o',...
    'linewidth',3);
legend('Default','RS','HODLR','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Dimension Scaling, Nonuniform','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'n_factortimes_extra.pdf','ContentType','vector');
end

%% Plot comparisons between uniform/stationary and nonuniform/nonstationary experiments

fig = figure;
loglog(...
    u.problem_sizes,u.default_times,'-o',...
    n.problem_sizes,n.default_times,'-o',...
    'linewidth',3);
legend('Uniform','Nonuniform','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('Default SKI Comparison','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'default_compare_extra.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    u.problem_sizes,u.rskelf_times,'-o',...
    n.problem_sizes,n.rskelf_times,'-o',...
    'linewidth',3);
legend('Uniform','Nonuniform','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('RS SKI Comparison','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'rskelf_compare_extra.pdf','ContentType','vector');
end

fig = figure;
loglog(...
    u.problem_sizes,u.hodlr_times,'-o',...
    n.problem_sizes,n.hodlr_times,'-o',...
    'linewidth',3);
legend('Uniform','Nonuniform','Location','north','Fontsize',12);
ax = gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;
title('HODLR SKI Comparison','Fontsize',14)
xlabel('Problem Size')
ylabel('Factor/Logdet Time')

if save_figs
    exportgraphics(fig, 'hodlr_compare_extra.pdf','ContentType','vector');
end
