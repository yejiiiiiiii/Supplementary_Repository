% Yeji Han
% Curvefit for 1st order Mooney-Rivlin model parameters -- Uniaxial tension, compression, and planar shear test

clear; clc; close all force

%% ========================== Constants ==========================
colors = lines(6);

% --- Uniaxial tension geometry ---
area_uni = 2 * 6;       % mm^2
gauge_uni = 30;         % mm

% --- Planar shear geometry ---
area_planar = 1.5 * 215;  % mm^2
gauge_planar = 10;        % mm

% --- Compression geometry ---
diameter_comp = 28.6;        % mm
area_comp = pi * (diameter_comp/2)^2;  % circular cross-section area [mm^2]
gauge_comp = 12.5; % mm, height of the sample

% Output folder for figures
outdir = 'plots_MR';
if ~exist(outdir,'dir'); mkdir(outdir); end

outdir2 = 'Abaqus_testdata';
if ~exist(outdir2,'dir'); mkdir(outdir2); end

%% ========================== Import data ==========================
% Uniaxial tension
% input: line file, instron file, gauge length, area, sensitivity1 (line detector), sensitivity2 (instron)
% output: t, dis, lambda, eps, sigma, F
uni_1 = import_and_match('uniaxial_line1.csv', 'uniaxial_ins1.csv', gauge_uni, area_uni, 0.1, 0.1);
uni_2 = import_and_match('uniaxial_line2.csv', 'uniaxial_ins2.csv', gauge_uni, area_uni, 0.5, 0.5);
uni_3 = import_and_match('uniaxial_line3.csv', 'uniaxial_ins3.csv', gauge_uni, area_uni, 0.205, 0.205);
uni_abaqus = table2array(readtable("BBDINO_Abaqus_Tension_MR.csv"));

uni_1 = zero_data(uni_1);
uni_2 = zero_data(uni_2);
uni_3 = zero_data(uni_3);

% Planar shear
plan_1 = import_and_match('plan_line1.csv', 'plan_ins1.csv', gauge_planar, area_planar, 0.3, 0.2);
plan_2 = import_and_match('plan_line2.csv', 'plan_ins2.csv', gauge_planar, area_planar, 0.3, 0.2);
plan_3 = import_and_match('plan_line3.csv', 'plan_ins3.csv', gauge_planar, area_planar, 0.5, 0.2);

plan_1 = zero_data(plan_1);
plan_2 = zero_data(plan_2);
plan_3 = zero_data(plan_3);
plan_abaqus = table2array(readtable("BBDINO_Abaqus_PlanarShear_MR.csv"));

% Compression
% output: t, dis, lambda, eps, sigma, F
comp_1 = import_compression('comp_ins1.csv', area_comp, gauge_comp);
comp_2 = import_compression('comp_ins2.csv', area_comp, gauge_comp);
comp_3 = import_compression('comp_ins3.csv', area_comp, gauge_comp);
comp_abaqus = table2array(readtable("BBDINO_Abaqus_Compression_MR.csv"));

comp_1 = zero_data(comp_1);
comp_2 = zero_data(comp_2);
comp_3 = zero_data(comp_3);


%% ========================== Plot (stress-strain) ==========================
% %% Uniaxial
% figure();
% plot(uni_1.eps, uni_1.sigma, 'Color', colors(1,:), 'LineWidth', 2);
% hold on
% plot(uni_2.eps, uni_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
% plot(uni_3.eps, uni_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);
% hold off
% 
% axis tight
% xlabel('Strain', 'FontName', 'Arial', 'FontSize', 15);
% ylabel('Stress [MPa]', 'FontName', 'Arial', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3'}, 'Location', 'best', 'FontSize', 15);
% set(gca, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 2);
% grid off; set(gcf, 'Position', [500 50 800 600]); box on;
% saveas(gcf, fullfile(outdir,'1.stress_strain_uniaxial.png'));
% 
% %% Planar shear
% figure();
% plot(plan_1.eps, plan_1.sigma, 'Color', colors(1,:), 'LineWidth', 2);
% hold on
% plot(plan_2.eps, plan_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
% plot(plan_3.eps, plan_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);
% hold off
% 
% axis tight
% xlabel('Strain', 'FontName', 'Arial', 'FontSize', 15);
% ylabel('Stress [MPa]', 'FontName', 'Arial', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3'}, 'Location', 'best', 'FontSize', 15);
% set(gca, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 2);
% grid off; set(gcf, 'Position', [500 50 800 600]); box on;
% saveas(gcf, fullfile(outdir,'2.stress_strain_planarshear.png'));
% 
% %% Compression
% figure();
% plot(comp_1.eps, comp_1.sigma, 'Color', colors(1,:), 'LineWidth', 2);
% hold on
% plot(comp_2.eps, comp_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
% plot(comp_3.eps, comp_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);
% hold off
% 
% axis tight
% xlabel('Strain', 'FontName', 'Arial', 'FontSize', 15);
% ylabel('Stress [MPa]', 'FontName', 'Arial', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3'}, 'Location', 'best', 'FontSize', 15);
% set(gca, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 2);
% grid off; set(gcf, 'Position', [500 50 800 600]); box on;
% saveas(gcf, fullfile(outdir,'3.stress_strain_compression.png'));

%% ========================== CURVE FITTING ==========================
data.uni = {uni_1, uni_2, uni_3};
data.plan = {plan_1, plan_2, plan_3};
data.comp = {comp_1, comp_2, comp_3};

% Weights for each test type
weights = [1, 1, 1]; % [uniaxial, planar, compression]

% Initial guess for [C10, C01]
% params0 = [0.1 0.1];
% params0 = [0.5 0.5];
params0 = [1 1];
options = optimoptions('lsqnonlin',...
'Display','iter',...
'Algorithm','levenberg-marquardt',...
'MaxIterations',5000,...
'MaxFunctionEvaluations',50000);

% Perform nonlinear least squares fitting
% 3rd and 4th -> lower and upper bounds on the parameters.
params_fit = lsqnonlin(@(p) cost_fun_MR_all(p, data, weights), params0, [], [], options);
C10 = params_fit(1); 
C01 = params_fit(2);

fprintf('\nFitted Parameters (using ALL data points):\n');
fprintf('C10 = %.6f MPa\n', C10);
fprintf('C01 = %.6f MPa\n', C01);

%% ========================== MODEL COMPARISON ==========================
%% Uniaxial
figure();

plot(uni_1.eps, uni_1.sigma, 'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(uni_2.eps, uni_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
plot(uni_3.eps, uni_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);

plot(uni_abaqus(:,1), uni_abaqus(:,2), 'Color', colors(4,:), 'LineWidth', 4);

eps_model = linspace(min(uni_1.eps), max(uni_1.eps), 200);
lambda_model = exp(eps_model); % since ε_true = ln(λ)
sigma_model = mooney_rivlin(params_fit, lambda_model, 'uniaxial');
plot(eps_model, sigma_model, 'k-', 'LineWidth', 2.5);

hold off; axis tight;
xlabel('True Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
ylabel('True Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
legend({'sample 1', 'sample 2', 'sample 3', 'simulation', 'model'}, 'Location', 'northwest', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3', 'Model'}, 'Location', 'northwest', 'FontSize', 15);
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);

xlim([0,uni_3.eps(end)]);
grid off; set(gcf, 'Position', [500 50 800 600]); box on;

saveas(gcf, fullfile(outdir, '1.model_fit_comparison_uniaxial.png'));

%% Planar shear
figure();
plot(plan_1.eps, plan_1.sigma, 'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(plan_2.eps, plan_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
plot(plan_3.eps, plan_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);

plot(plan_abaqus(:,1), plan_abaqus(:,2), 'Color', colors(4,:), 'LineWidth', 4);

eps_model = linspace(min(plan_1.eps), max(plan_1.eps), 200);
lambda_model = exp(eps_model);
sigma_model = mooney_rivlin(params_fit, lambda_model, 'planar');
plot(eps_model, sigma_model, 'k-', 'LineWidth', 2.5);

hold off; axis tight;
xlabel('True Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
ylabel('True Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
legend({'sample 1', 'sample 2', 'sample 3', 'simulation', 'model'}, 'Location', 'northwest', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3', 'Model'}, 'Location', 'northwest', 'FontSize', 15);
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);

xlim([0,plan_3.eps(end)]);
grid off; set(gcf, 'Position', [500 50 800 600]); box on;

saveas(gcf, fullfile(outdir, '2.model_fit_comparison_planarshear.png'));

%% Compression
figure();
plot(comp_1.eps, comp_1.sigma, 'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(comp_2.eps, comp_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
plot(comp_3.eps, comp_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);

plot(comp_abaqus(:,1), comp_abaqus(:,2), 'Color', colors(4,:), 'LineWidth', 4);

eps_model = linspace(min(comp_1.eps), max(comp_1.eps), 200);
lambda_model = exp(eps_model);
sigma_model = mooney_rivlin(params_fit, lambda_model, 'compression');
plot(eps_model, sigma_model, 'k-', 'LineWidth', 2.5);


hold off; axis tight;
xlabel('True Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
ylabel('True Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
legend({'sample 1', 'sample 2', 'sample 3', 'simulation', 'model'}, 'Location', 'southeast', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3', 'Model'}, 'Location', 'southeast', 'FontSize', 15);
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);
grid off; set(gcf, 'Position', [500 50 800 600]); box on;
xlim([comp_1.eps(end),0]);
saveas(gcf, fullfile(outdir, '3.model_fit_comparison_compression.png'));

%% ========================== ENGINEERING COMPARISON PLOTS ==========================
% fprintf('\n--- Plotting Engineering Stress–Strain Comparisons ---\n');
% 
% % --- Helper to plot engineering stress-strain curves ---
% function plot_engineering_comparison(dat1, dat2, dat3, params_fit, mode, outdir, outname)
%     colors = lines(6);
%     figure(); hold on;
% 
%     % --- Plot experimental data (engineering) ---
%     samples = {dat1, dat2, dat3};
%     for k = 1:3
%         eps_eng = exp(samples{k}.eps) - 1;
%         sig_eng = samples{k}.sigma ./ exp(samples{k}.eps);
%         plot(eps_eng, sig_eng, 'Color', colors(k,:), 'LineWidth', 2);
%     end
% 
%     % --- Plot model curve (engineering) ---
%     eps_model = linspace(min(exp(dat1.eps)-1), max(exp(dat1.eps)-1), 200);
%     lambda_model = 1 + eps_model;
%     sigma_true = mooney_rivlin(params_fit, lambda_model, mode);
%     sigma_eng = sigma_true ./ lambda_model;
% 
%     plot(eps_model, sigma_eng, 'k-', 'LineWidth', 2.5);
% 
%     % --- Format plot ---
%     xlabel('Engineering Strain','FontName','Arial','FontSize',15);
%     ylabel('Engineering Stress [MPa]','FontName','Arial','FontSize',15);
%     legend({'sample 1','sample 2','sample 3','model'},'Location','best','FontSize',15);
%     set(gca,'FontSize',15,'FontName','Arial','LineWidth',2);
%     grid off; box on; axis tight;
%     set(gcf,'Position',[500 50 800 600]);
% 
%     saveas(gcf, fullfile(outdir, outname));
% end
% 
% % --- Uniaxial ---
% plot_engineering_comparison(uni_1, uni_2, uni_3, params_fit, 'uniaxial', outdir, '4.model_fit_comparison_uniaxial_engineering.png');
% 
% % --- Planar shear ---
% plot_engineering_comparison(plan_1, plan_2, plan_3, params_fit, 'planar', outdir, '5.model_fit_comparison_planarshear_engineering.png');
% 
% % --- Compression ---
% plot_engineering_comparison(comp_1, comp_2, comp_3, params_fit, 'compression', outdir, '6.model_fit_comparison_compression_engineering.png');
% 

%% ========================== EXPORT NOMINAL DATA FOR ABAQUS ==========================
% Abaqus expects: Nominal Stress [MPa], Nominal Strain [unitless ratio]

% -------------------------------------------------------------------------
% 1. Uniaxial tension
eps_nom_uni = exp(uni_1.eps) - 1;      % convert from true → nominal
sig_nom_uni = uni_1.sigma ./ (1 + eps_nom_uni);  % convert from true → nominal

uniaxial_nominal = [eps_nom_uni(:), sig_nom_uni(:)];  % [strain, stress]
writematrix(uniaxial_nominal, fullfile(outdir2, 'uniaxial_nominal.csv'));

% -------------------------------------------------------------------------
% 2. Compression
eps_nom_comp = exp(comp_1.eps) - 1;    % true → nominal
sig_nom_comp = comp_1.sigma ./ (1 + eps_nom_comp);

compression_nominal = [eps_nom_comp(:), sig_nom_comp(:)];
writematrix(compression_nominal, fullfile(outdir2, 'compression_nominal.csv'));

% -------------------------------------------------------------------------
% 3. Planar shear
eps_nom_plan = exp(plan_1.eps) - 1;
sig_nom_plan = plan_1.sigma ./ (1 + eps_nom_plan);

planarshear_nominal = [eps_nom_plan(:), sig_nom_plan(:)];
writematrix(planarshear_nominal, fullfile(outdir2, 'planarshear_nominal.csv'));

% -------------------------------------------------------------------------
% 4. Combined compression + tension (for Abaqus "Uniaxial Test Data")
uniaxial_combined = [eps_nom_comp(:), sig_nom_comp(:); ...
                     eps_nom_uni(:),  sig_nom_uni(:)];

% Sort by strain (ascending) to ensure continuity in Abaqus and plots
uniaxial_combined = sortrows(uniaxial_combined, 1);

% Export combined file
writematrix(uniaxial_combined, fullfile(outdir2, 'Uniaxial_comp_tensile_nominal.csv'));

% -------------------------------------------------------------------------
fprintf('\n=== Export Complete ===\n');
fprintf('Nominal data saved in folder: %s\n', outdir2);
fprintf('Files:\n');
fprintf(' - uniaxial_nominal.csv\n');
fprintf(' - compression_nominal.csv\n');
fprintf(' - planarshear_nominal.csv\n');
fprintf(' - Uniaxial_comp_tensile_nominal.csv (combined + sorted)\n');
fprintf('Columns: [Nominal Strain (ratio), Nominal Stress (MPa)]\n');

%% Import functions
% ================= UNAXIAL TENSION & PLANAR SHEAR =================
function data = import_and_match(line_file, ins_file, gauge, area, c1, c2)
    % Columns: time (s), distance (mm), distance (px), y1 (px), y2 (px)
    d1 = table2array(readtable(line_file));
    t1 = d1(:,1);
    dis_line = d1(:,2) - d1(1,2);      % displacement (mm)
    eps_eng_line = dis_line / gauge;   % engineering strain
    lambda = 1 + eps_eng_line;         % stretch

    % --- Load Instron data ---
    % Columns: time (s), displacement (mm), Force (N)
    d2 = table2array(readtable(ins_file));
    t2 = d2(:,1);
    dis_ins = d2(:,2) - d2(1,2);
    F = d2(:,3);

    % --- Velocity-based start detection for time alignment ---
    vel_line = gradient(dis_line);
    vel_ins = gradient(dis_ins);
    th1 = c1 * max(vel_line);
    th2 = c2 * max(vel_ins);
    idx_start_line = find(vel_line > th1, 1, 'first');
    idx_start_ins = find(vel_ins > th2, 1, 'first');

    % Align start times
    t1 = t1 - t1(idx_start_line);
    t2 = t2 - t2(idx_start_ins);

    % Trim data before motion start
    t1 = t1(idx_start_line:end);
    dis_line = dis_line(idx_start_line:end);
    lambda = lambda(idx_start_line:end);

    t2 = t2(idx_start_ins:end);
    F = F(idx_start_ins:end);
    %dis_ins = dis_ins(idx_start_ins:end);

    % --- Trim to peak deformation (first loading segment) ---
    [~, idx_peak] = max(dis_line);
    t1 = t1(1:idx_peak);
    dis_line = dis_line(1:idx_peak);
    lambda = lambda(1:idx_peak);

    % --- Interpolate both datasets to a common timeline ---
    N_target = min(length(t1), length(t2));
    t_common = linspace(0, min(t1(end), t2(end)), N_target);

    dis_line = interp1(t1, dis_line, t_common, 'linear', 'extrap');
    lambda = interp1(t1, lambda, t_common, 'linear', 'extrap');
    F = interp1(t2, F, t_common, 'linear', 'extrap');
    %dis_ins = interp1(t2, dis_ins, t_common, 'linear', 'extrap');

    % --- Compute true quantities (AFTER time matching) ---
    eps_true = log(lambda);            % true strain (ln λ)
    sigma_true = (F .* lambda) / area; % true stress (Cauchy) = Fλ/A₀

    % --- Output structured data ---
    data.t = t_common(:);
    data.dis = dis_line(:);
    %data.dis_ins = dis_ins(:);
    data.lambda = lambda(:);
    data.eps = eps_true(:);
    data.sigma = sigma_true(:);
    data.F = F(:);
end


% ========================== COMPRESSION ============================
function data = import_compression(ins_file, area_comp, gauge_comp)
    d = table2array(readtable(ins_file));
    t = d(:,1);
    dis = -d(:,2); % Convert compression to tensile-positive
    F = -d(:,3);   % Flip sign so compression stress is positive (tension-style)

    lambda = 1 + dis / gauge_comp;
    eps_true = log(lambda);
    sigma_true = (F .* lambda) / area_comp;  % True (Cauchy) stress = Fλ/A₀

    data.t = t(:);
    data.dis = dis(:);
    data.lambda = lambda(:);
    data.eps = eps_true(:);
    data.sigma = sigma_true(:);
    data.F = F(:);
end

%% ================== ZEROING FUNCTION ==================
function data = zero_data(data)
    % Force first point to be zero for all major quantities
    if isfield(data,'t'),        data.t = data.t - data.t(1); end
    if isfield(data,'dis'),      data.dis = data.dis - data.dis(1); end
    if isfield(data,'lambda'),   data.lambda = data.lambda - data.lambda(1) + 1; end
    if isfield(data,'eps'),      data.eps = data.eps - data.eps(1); end
    if isfield(data,'sigma'),    data.sigma = data.sigma - data.sigma(1); end
    if isfield(data,'F'),        data.F = data.F - data.F(1); end
end

%% ======================= COST FUNCTION =============================
function E_all = cost_fun_MR_all(params, data, weights)
    % Implements the cost function from the document:
    % E = Σ_i w_i e_i,  e_i = (1/n_i) Σ_j [(σ_exp - σ_model)/(σ_max - σ_min)]²

    modes = {'uni', 'plan', 'comp'};
    E_all = [];

    for i = 1:3
        mode = modes{i};
        w = weights(i);

        % Loop over all samples in this mode
        for k = 1:length(data.(mode))
            sample = data.(mode){k};
            lambda = sample.lambda;
            sigma_exp = sample.sigma;
            sigma_model = mooney_rivlin(params, lambda, mode);

            % Normalization for this mode
            ni = length(sigma_exp);
            s_range = max(sigma_exp) - min(sigma_exp);
            if s_range == 0, s_range = 1; end

            % Normalized error (dimensionless)
            % lsqnonlin minimizes sum(f.^2), so each element of 'e' must be the sqrt of the weighted, normalized residual.
            e = sqrt(w/ni) * ((sigma_exp - sigma_model) ./ s_range);

            % Append all residuals
            E_all = [E_all; e(:)];
        end
    end
end

%% ======================= MOONEY-RIVLIN MODEL =======================
function sigma = mooney_rivlin(params, lambda, mode)
    C10 = params(1); 
    C01 = params(2);

    switch lower(mode)
        case {'uni','uniaxial'}
            % σ11 = 2C10(λ^2 - λ^(-1)) + 2C01(λ - λ^(-2))
            sigma = 2*C10.*(lambda.^2 - lambda.^(-1)) + ...
                    2*C01.*(lambda - lambda.^(-2));

        case {'plan','planar'}
            % σ22 = 0
            % σ11 = 2C10(λ^2 - 1) - 2C01(λ^(-2) - 1)
            %sigma = 2*C10.*(lambda.^2 - 1) - 2*C01.*(lambda.^(-2)-1);

            % σ33 = 0
            % σ11 = 2(C10 + C01)(λ^2 - λ^(-2))
            sigma = 2*(C10 + C01).*(lambda.^2 - lambda.^(-2));

        case {'comp','compression'}
            % σ11 = 2C10(λ^2 - λ^(-1)) + 2C01(λ - λ^(-2))
            sigma = 2*C10.*(lambda.^2 - lambda.^(-1)) + ...
                    2*C01.*(lambda - lambda.^(-2));

        otherwise
            error('Unknown mode: %s', mode);
    end
end

%% ================== Interactive Visualization: C10, C01 True ==================
% %% --- Figure layout and styling ---
% f = uifigure('Name','Mooney–Rivlin Parameter Explorer (True Stress)', ...
%     'Position',[500 50 2400 650]);
% 
% % === Leave some bottom space for sliders ===
% t = tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
% t.Position = [0.05 0.25 0.9 0.7];
% 
% colors = lines(6);
% modes = {'uniaxial','planar','compression'};
% data_all = {data.uni, data.plan, data.comp};
% 
% % === Fixed axis ranges (adjust if needed) ===
% lims = {
%     [0, 0.8; 0, 0.8];   % uniaxial
%     [0, 1.0; 0, 1.5];   % planar
%     [-0.5, 0; -0.8, 0]; % compression
% };
% 
% % --- Axes initialization ---
% ax = gobjects(1,3);
% modelLine = gobjects(1,3);
% 
% for i = 1:3
%     ax(i) = nexttile(t,i);
%     hold(ax(i),'on');
%     dat = data_all{i};
%     for k = 1:length(dat)
%         plot(ax(i), dat{k}.eps, dat{k}.sigma, ...
%             'Color', colors(k,:), 'LineWidth', 2);
%     end
%     eps_model = linspace(min(dat{1}.eps), max(dat{1}.eps), 300);
%     modelLine(i) = plot(ax(i), eps_model, zeros(size(eps_model)), ...
%         'k-', 'LineWidth', 2.5);
% 
%     xlabel(ax(i),'True Strain','FontName','Arial','FontSize',15);
%     ylabel(ax(i),'True Stress [MPa]','FontName','Arial','FontSize',15);
%     title(ax(i), capitalizeFirst(modes{i}), 'FontSize',15,'FontWeight','bold');
%     legend(ax(i),{'sample 1','sample 2','sample 3','model'}, ...
%         'Location','southeast','FontSize',13);
%     set(ax(i),'FontSize',15,'FontName','Arial','LineWidth',2);
%     box(ax(i),'on'); grid(ax(i),'off');
%     xlim(ax(i), lims{i}(1,:));
%     ylim(ax(i), lims{i}(2,:));
%     hold(ax(i),'off');
% end
% 
% 
% 
% %% ================== SLIDERS + TEXT INPUTS ==================
% slider_ypos = 70;
% tick_vals = -1:0.2:1;
% C10_init = C10;
% C01_init = C01;
% 
% % --- C10 slider + numeric field ---
% uilabel(f, 'Text','C10', 'Position',[600 slider_ypos-20 40 22], 'FontSize',14);
% C10_slider = uislider(f, 'Limits',[-1,1], 'Value',C10_init, ...
%     'MajorTicks',tick_vals, 'MinorTicks',[], 'Position',[650 slider_ypos 400 20]);
% C10_edit = uieditfield(f,'numeric','Limits',[-1,1],'Value',C10_init, ...
%     'Position',[1070 slider_ypos-12 60 28],'FontSize',14,'ValueDisplayFormat','%.4g');
% 
% % --- C01 slider + numeric field ---
% uilabel(f, 'Text','C01', 'Position',[1300 slider_ypos-20 40 22], 'FontSize',14);
% C01_slider = uislider(f, 'Limits',[-1,1], 'Value',C01_init, ...
%     'MajorTicks',tick_vals, 'MinorTicks',[], 'Position',[1350 slider_ypos 400 20]);
% C01_edit = uieditfield(f,'numeric','Limits',[-1,1],'Value',C01_init, ...
%     'Position',[1770 slider_ypos-12 60 28],'FontSize',14,'ValueDisplayFormat','%.4g');
% 
% % --- Listeners for real-time updates ---
% addlistener(C10_slider,'ValueChanging', ...
%     @(src,evt) sync_update(evt.Value, C01_slider.Value, C10_edit, C01_edit, modes, data_all, modelLine));
% addlistener(C01_slider,'ValueChanging', ...
%     @(src,evt) sync_update(C10_slider.Value, evt.Value, C10_edit, C01_edit, modes, data_all, modelLine));
% 
% C10_edit.ValueChangedFcn = @(src,evt) manual_update(src.Value, C01_edit.Value, C10_slider, C01_slider, modes, data_all, modelLine);
% C01_edit.ValueChangedFcn = @(src,evt) manual_update(C10_edit.Value, src.Value, C10_slider, C01_slider, modes, data_all, modelLine);
% 
% update_plots(C10_init, C01_init, modes, data_all, modelLine);
% 
% %% --- Local subfunctions ---
% function strOut = capitalizeFirst(strIn)
%     strOut = lower(strIn);
%     strOut(1) = upper(strOut(1));
% end
% 
% function sync_update(C10,C01,C10_edit,C01_edit,modes,data_all,modelLine)
%     C10_edit.Value = C10;
%     C01_edit.Value = C01;
%     update_plots(C10,C01,modes,data_all,modelLine);
% end
% 
% function manual_update(C10,C01,C10_slider,C01_slider,modes,data_all,modelLine)
%     C10_slider.Value = C10;
%     C01_slider.Value = C01;
%     update_plots(C10,C01,modes,data_all,modelLine);
% end
% 
% function update_plots(C10,C01,modes,data_all,modelLine)
%     params = [C10, C01];
%     for i = 1:3
%         dat = data_all{i};
%         eps_model = linspace(min(dat{1}.eps), max(dat{1}.eps), 300);
%         lambda_model = exp(eps_model);
%         sigma_model = mooney_rivlin(params, lambda_model, modes{i});
%         modelLine(i).XData = eps_model;
%         modelLine(i).YData = sigma_model;
%     end
%     drawnow limitrate nocallbacks;
% end

% %% ================== INTERACTIVE: ENGINEERING STRESS–STRAIN ==================
% f2 = uifigure('Name','Mooney–Rivlin (Engineering) Parameter Explorer', ...
%     'Position',[200 50 2400 650]);
% 
% t2 = tiledlayout(f2,1,3,'Padding','compact','TileSpacing','compact');
% t2.Position = [0.05 0.25 0.9 0.7];
% 
% colors = lines(6);
% ax2 = gobjects(1,3);
% modelLine2 = gobjects(1,3);
% 
% for i = 1:3
%     ax2(i) = nexttile(t2,i);
%     hold(ax2(i),'on');
%     dat = data_all{i};
%     all_eps = []; all_sig = [];
% 
%     % --- Experimental engineering data ---
%     for k = 1:length(dat)
%         eps_eng = exp(dat{k}.eps) - 1;
%         sig_eng = dat{k}.sigma ./ exp(dat{k}.eps);
%         all_eps = [all_eps; eps_eng(:)];
%         all_sig = [all_sig; sig_eng(:)];
%         plot(ax2(i), eps_eng, sig_eng, 'Color', colors(k,:), 'LineWidth', 2);
%     end
% 
%     eps_model_eng = linspace(min(all_eps), max(all_eps), 300);
%     modelLine2(i) = plot(ax2(i), eps_model_eng, zeros(size(eps_model_eng)), ...
%         'k-', 'LineWidth', 2.5);
% 
%     xlabel(ax2(i),'Engineering Strain','FontName','Arial','FontSize',15);
%     ylabel(ax2(i),'Engineering Stress [MPa]','FontName','Arial','FontSize',15);
%     title(ax2(i), capitalizeFirst(modes{i}), 'FontSize',15,'FontWeight','bold');
%     legend(ax2(i),{'sample 1','sample 2','sample 3','model'}, ...
%         'FontSize',13,'Location','southeast');
%     set(ax2(i),'FontSize',15,'FontName','Arial','LineWidth',2);
%     box(ax2(i),'on'); grid(ax2(i),'off');
% 
%     % --- Auto axis fit ---
%     xpad = (max(all_eps)-min(all_eps)) * 0.05;
%     ypad = (max(all_sig)-min(all_sig)) * 0.1;
%     xlim(ax2(i), [min(all_eps)-xpad, max(all_eps)+xpad]);
%     ylim(ax2(i), [min(all_sig)-ypad, max(all_sig)+ypad]);
%     hold(ax2(i),'off');
% end
% 
% %% ================== SLIDERS (Engineering) ==================
% slider_ypos = 70;
% tick_vals = -1:0.2:1;
% C10_init = C10;
% C01_init = C01;
% 
% uilabel(f2,'Text','C10','Position',[500 slider_ypos-20 40 22],'FontSize',14);
% C10s = uislider(f2,'Limits',[-1 1],'Value',C10_init, ...
%     'MajorTicks',tick_vals,'MinorTicks',[], 'Position',[550 slider_ypos 400 20]);
% C10e = uieditfield(f2,'numeric','Limits',[-1 1],'Value',C10_init, ...
%     'Position',[970 slider_ypos-12 80 28],'FontSize',14,'ValueDisplayFormat','%.4g');
% 
% uilabel(f2,'Text','C01','Position',[1100 slider_ypos-20 40 22],'FontSize',14);
% C01s = uislider(f2,'Limits',[-1 1],'Value',C01_init, ...
%     'MajorTicks',tick_vals,'MinorTicks',[], 'Position',[1150 slider_ypos 400 20]);
% C01e = uieditfield(f2,'numeric','Limits',[-1 1],'Value',C01_init, ...
%     'Position',[1570 slider_ypos-12 80 28],'FontSize',14,'ValueDisplayFormat','%.4g');
% 
% % ================== LIVE UPDATES ==================
% addlistener(C10s,'ValueChanging', ...
%     @(src,evt) sync_update_eng(evt.Value,C01s.Value,C10e,C01e,modes,data_all,modelLine2));
% addlistener(C01s,'ValueChanging', ...
%     @(src,evt) sync_update_eng(C10s.Value,evt.Value,C10e,C01e,modes,data_all,modelLine2));
% 
% C10e.ValueChangedFcn = @(src,evt) manual_update_eng(src.Value,C01e.Value,C10s,C01s,modes,data_all,modelLine2);
% C01e.ValueChangedFcn = @(src,evt) manual_update_eng(C10e.Value,src.Value,C10s,C01s,modes,data_all,modelLine2);
% 
% update_plots_eng(C10_init,C01_init,modes,data_all,modelLine2);

% %% ================== SUBFUNCTIONS ==================
% function sync_update_eng(C10,C01,C10e,C01e,modes,data_all,modelLine2)
%     C10e.Value = C10; C01e.Value = C01;
%     update_plots_eng(C10,C01,modes,data_all,modelLine2);
% end
% 
% function manual_update_eng(C10,C01,C10s,C01s,modes,data_all,modelLine2)
%     C10s.Value = C10; C01s.Value = C01;
%     update_plots_eng(C10,C01,modes,data_all,modelLine2);
% end
% 
% function update_plots_eng(C10,C01,modes,data_all,modelLine2)
%     params = [C10, C01];
%     for i = 1:3
%         dat = data_all{i};
%         eps_model = linspace(min(dat{1}.eps), max(dat{1}.eps), 300);
%         lambda_model = exp(eps_model);
%         sigma_true = mooney_rivlin(params, lambda_model, modes{i});
%         eps_eng = lambda_model - 1;
%         sig_eng = sigma_true ./ lambda_model;
%         modelLine2(i).XData = eps_eng;
%         modelLine2(i).YData = sig_eng;
%     end
%     drawnow limitrate nocallbacks;
% end
