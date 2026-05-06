clear; clc; close all;

% Yeji Han
% Stress relaxation test > Prony Series parameters (loop over multiple N)

%% Constants
epsilon0 = 0.3;  % Constant strain during test
% orders = [2, 3, 5, 9];  % Orders of Prony Series to loop over
orders = [2 5];  % Orders of Prony Series to loop over

area = pi*14.3^2; % mm^2, sample area

%% Load data
% time, strain, stress (kPa)
d1 = table2array(readtable('BBDINO30_CompRelaxation (1).csv'));
d2 = table2array(readtable('BBDINO30_CompRelaxation (2).csv'));
d3 = table2array(readtable('BBDINO30_CompRelaxation (3).csv'));

% Output folder for figures
outdir = 'plots';
if ~exist(outdir,'dir'); mkdir(outdir); end

min_idx = min([length(d1) length(d2) length(d3)]);

d1 = d1((1:end-(length(d1)-min_idx)),:);
d2 = d2((1:end-(length(d2)-min_idx)),:);
d3 = d3((1:end-(length(d3)-min_idx)),:);

%% average
dis_raw_1 = d1(:,2); dis_raw_2 = d2(:,2); dis_raw_3 = d3(:,2);
t_raw_1 = d1(:,1); t_raw_2 = d2(:,1); t_raw_3 = d3(:,1);
sigma_raw_1 = d1(:,3)/area*1000; sigma_raw_2 = d2(:,3)/area*1000; sigma_raw_3 = d3(:,3)/area*1000;

t_raw_avg = (t_raw_1 + t_raw_2 + t_raw_3)/3;
sigma_raw_avg = (sigma_raw_1 + sigma_raw_2 + sigma_raw_3)/3;
dis_raw_avg = (dis_raw_1 + dis_raw_2 + dis_raw_3)/3;

index_1 = find( abs(dis_raw_1 - dis_raw_1(end)) == 0, 3);
index_1 = index_1(1);
index_2 = find( abs(dis_raw_2 - dis_raw_2(end)) == 0, 3);
index_2 = index_2(1);
index_3 = find( abs(dis_raw_3 - dis_raw_3(end)) == 0, 3);
index_3 = index_3(1);
index_avg = find( abs(dis_raw_avg - dis_raw_avg(end)) == 0, 1);


t_1 = t_raw_1(index_1:end); t_2 = t_raw_2(index_2:end); t_3 = t_raw_3(index_3:end);
dis_raw_1 = dis_raw_1(index_1:end); dis_raw_2 = dis_raw_2(index_2:end); dis_raw_3 = dis_raw_3(index_3:end);
t_avg = t_raw_avg(index_avg:end);
t_avg = t_avg - t_avg(1);
dis_raw_avg = dis_raw_avg(index_avg:end);
sigma_1 = sigma_raw_1(index_1:end); sigma_2 = sigma_raw_2(index_2:end); sigma_3 = sigma_raw_3(index_3:end);
sigma_avg = sigma_raw_avg(index_avg:end);
Et_1 = sigma_1 / epsilon0; Et_2 = sigma_2 / epsilon0; Et_3 = sigma_3 / epsilon0;
Et_avg = sigma_avg / epsilon0;

Et_avg_set = zeros(length(Et_avg),2,length(orders));

sigma_diff = sigma_avg(1) - sigma_avg(end);
sigma_diff_per = sigma_diff/sigma_avg(1)*100;
fprintf("loss: %.3f kPa, loss in percentage: %.2f %%\n",sigma_diff, sigma_diff_per);

%% figure - avg
figure();
plot(t_avg, sigma_avg, 'r-', 'LineWidth', 1.5);

xlabel('Time [s]', 'FontName', 'Times new roman', 'FontSize', 25)
ylabel('Stress [kPa]', 'FontName', 'Times new roman', 'FontSize', 25)

xlim([0,t_avg(end)])
ylim([350,380])

title('Average curve', 'FontName', 'Times new roman','FontSize', 25);

ax = gca;
set(ax, 'FontSize', 25, 'FontName', 'Times new roman', 'LineWidth', 2);
set(gcf, 'Position', [500 50 800 600]);
box on;

saveas(gcf, fullfile(outdir, '2. comprelaxation_avg.png'));
saveas(gcf, fullfile(outdir,'2. comprelaxation_avg.fig'))

%% Loop over Prony orders
iter = 3;
for idx = 1:length(orders)
    N = orders(idx);  % current order

    % Initial guesses
    E_e_guess = Et_avg(end);  % or Et_avg(1) * 0.1
    E_guess = linspace(Et_avg(1), Et_avg(end), N);    
    tau_guess = logspace(log10(0.01), log10(100), N); 

    % Initial parameter vector: [E_e, E1, tau1, E2, tau2, ..., EN, tauN]
    p0 = [E_e_guess, reshape([E_guess; tau_guess], 1, [])];

    % Prony series
    prony = @(p, t) pronyfun(p,t,N);

    % Fit options
    options = optimoptions('lsqcurvefit', 'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 1e4, 'Display', 'off'); 
    % Lower bounds (all parameters ≥ 0)
    lb = zeros(size(p0));

    % Curve fit
    p_fit = lsqcurvefit(prony, p0, t_avg, Et_avg, lb, [], options);

    % save the variables
    Et_avg_set(:,1,idx) = t_avg;
    Et_avg_set(:,2,idx) = prony(p_fit, t_avg);

    % Plot (E-t)
    figure;
    % subplot(1,2,1)
    plot(t_avg, Et_avg*epsilon0, 'r-','Linewidth',2, 'DisplayName', 'Data');
    hold on;
    plot(t_avg, prony(p_fit, t_avg)*epsilon0, 'b--', 'LineWidth', 2, 'DisplayName', 'Prony Fit');
    hold off

    xlabel('Time [s]','FontName', 'Times new roman', 'FontSize', 25)
    ylabel('Compressive Stress [kPa]','FontName', 'Times new roman', 'FontSize', 25);

    lgd = legend('test data', 'curve fit','location','northeast');
    set(lgd,'FontSize', 25, 'FontName', 'Times new roman');
    legend boxon

    xlim([0,t_avg(end)]);
    ylim([350,380])

    % title(sprintf('Prony Series Fit (N = %d)', N), 'FontName', 'Times new roman','FontSize', 25);
    grid off;
    ax = gca;
    set(ax, 'FontSize', 25, 'FontName', 'Times new roman', 'LineWidth', 2);
    set(gcf, 'Position', [500 50 800 600]);

    filename = sprintf('%d.comprelaxation_Prony_%d.png', iter, N);
    filename2 = sprintf('%d.comprelaxation_Prony_%d.fig', iter, N);
    saveas(gcf, fullfile(outdir, filename));
    saveas(gcf, fullfile(outdir, filename2));

    iter = iter+1;

    % Extract fitted parameters
    E_e = p_fit(1);
    E_i = p_fit(2:2:end);
    tau_i = p_fit(3:2:end);
    E_i = E_i / Et_avg(1);

    % Display
    fprintf('\n====== Prony Order N = %d ======\n', N);
    fprintf('E_e = %.3f kPa\n', E_e);
    T = table((1:N)', E_i', E_i', tau_i', ...
        'VariableNames', {'Term', 'E_i [kPa]', 'E_i', 'tau_i [s]'});
    disp(T);

    E_model = pronyfun(p_fit, t_avg, N);

    residual = Et_avg - E_model;
    RMSE = sqrt(mean(residual.^2));
    MAE = mean(abs(residual));
    R2 = 1 - sum(residual.^2) / sum((Et_avg - mean(Et_avg)).^2);

    fprintf('\nFit Quality Metrics (N = %d)\n', N);
    fprintf('RMSE: %.3e Pa\n', RMSE);
    fprintf('MAE : %.3e Pa\n', MAE);
    fprintf('R²  : %.4f\n', R2);
end

filename = sprintf('Et_%g', epsilon0*100);
save(filename, 'Et_avg_set');

function E = pronyfun(p, t, N)
    E = p(1) * ones(size(t));  % E_e (equilibrium modulus)

    for i = 1:N
        E_i = p(2*i);
        tau_i = p(2*i + 1);
        E = E + E_i * exp(-t / tau_i);
    end
end