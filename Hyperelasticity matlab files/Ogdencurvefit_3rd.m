% Yeji Han
% Curvefit for Ogden model parameters -- Uniaxial tension, compression, and planar shear test

clear; clc; close all force

colors = lines(6);

% --- Uniaxial tension geometry ---
area_uni = 2 * 6;       % mm^2
gauge_uni = 30;         % mm

% --- Planar shear geometry ---
area_planar = 1.5 * 215;  % mm^2
gauge_planar = 10;        % mm

% --- Compression geometry ---
diameter_comp = 28.6;        % mm  (example)
area_comp = pi * (diameter_comp/2)^2;  % circular cross-section area [mm^2]
gauge_comp = 12.5; % mm

% Output folder for figures
outdir = 'plots_ogden_3rd';
if ~exist(outdir,'dir'); mkdir(outdir); end

% outdir2 = 'Abaqus_testdata';
% if ~exist(outdir2,'dir'); mkdir(outdir2); end

%% ========================== Import data ==========================
% Uniaxial tension
% input: line file, instron file, gauge length, area, sensitivity1, sensitivity2
% output: t, dis, lambda, eps, sigma, F
uni_1 = import_and_match('uniaxial_line1.csv', 'uniaxial_ins1.csv', gauge_uni, area_uni, 0.1, 0.1);
uni_2 = import_and_match('uniaxial_line2.csv', 'uniaxial_ins2.csv', gauge_uni, area_uni, 0.5, 0.5);
uni_3 = import_and_match('uniaxial_line3.csv', 'uniaxial_ins3.csv', gauge_uni, area_uni, 0.205, 0.205);
uni_abaqus = table2array(readtable("BBDINO_Abaqus_Tension_Ogden3.csv"));

uni_1 = zero_data(uni_1);
uni_2 = zero_data(uni_2);
uni_3 = zero_data(uni_3);

% Planar shear
plan_1 = import_and_match('plan_line1.csv', 'plan_ins1.csv', gauge_planar, area_planar, 0.3, 0.2);
plan_2 = import_and_match('plan_line2.csv', 'plan_ins2.csv', gauge_planar, area_planar, 0.3, 0.2);
plan_3 = import_and_match('plan_line3.csv', 'plan_ins3.csv', gauge_planar, area_planar, 0.5, 0.2);
plan_abaqus = table2array(readtable("BBDINO_Abaqus_PlanarShear_Ogden3.csv"));

plan_1 = zero_data(plan_1);
plan_2 = zero_data(plan_2);
plan_3 = zero_data(plan_3);

% Compression
% output: t, dis, lambda, eps, sigma, F
comp_1 = import_compression('comp_ins1.csv', area_comp, gauge_comp);
comp_2 = import_compression('comp_ins2.csv', area_comp, gauge_comp);
comp_3 = import_compression('comp_ins3.csv', area_comp, gauge_comp);
comp_abaqus = table2array(readtable("BBDINO_Abaqus_Compression_Ogden3.csv"));

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
% xlabel('Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
% ylabel('Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3'}, 'Location', 'best', 'FontSize', 15);
% set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);
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
% xlabel('Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
% ylabel('Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3'}, 'Location', 'best', 'FontSize', 15);
% set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);
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
% xlabel('Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
% ylabel('Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
% legend({'sample 1', 'sample 2', 'sample 3'}, 'Location', 'best', 'FontSize', 15);
% set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);
% grid off; set(gcf, 'Position', [500 50 800 600]); box on;
% saveas(gcf, fullfile(outdir,'3.stress_strain_compression.png'));

%% ========================== CURVE FITTING ==========================
data.uni  = {uni_1,uni_2,uni_3};
data.plan = {plan_1,plan_2,plan_3};
data.comp = {comp_1,comp_2,comp_3};
weights = [1 1 1];

% [mu1 mu2 mu3 alpha1 alpha2 alpha3]
params0 = [0.071846 0.071850 0.071848 1.953396 1.953507 1.953396];

options = optimoptions('lsqnonlin',...
'Display','iter',...
'Algorithm','levenberg-marquardt',...
'MaxIterations',5000,...
'MaxFunctionEvaluations',50000);

lb = [zeros(1,3), -10*ones(1,3)];
ub = [ones(1,3)*10,  10*ones(1,3)];

params_fit = lsqnonlin(@(p) cost_fun_Ogden_all(p,data,weights),params0,lb,ub,options);

mu1=params_fit(1); mu2=params_fit(2); mu3=params_fit(3);
a1=params_fit(4);  a2=params_fit(5);  a3=params_fit(6);

fprintf('\nFitted Ogden Parameters:\n');
fprintf('mu1=%.6f, mu2=%.6f, mu3=%.6f\n',mu1,mu2,mu3);
fprintf('a1 =%.6f, a2 =%.6f, a3 =%.6f\n',a1,a2,a3);

fprintf('\nFitted Ogden Parameters for Abaqus:\n');
fprintf('mu1=%.6f, mu2=%.6f, mu3=%.6f\n',mu1*a1/2,mu2*a2/2,mu3*a3/2);
fprintf('a1 =%.6f, a2 =%.6f, a3 =%.6f\n',a1,a2,a3);

%% ========================== MODEL COMPARISON ==========================
%% Uniaxial 
figure();
plot(uni_1.eps, uni_1.sigma, 'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(uni_2.eps, uni_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
plot(uni_3.eps, uni_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);

plot(uni_abaqus(:,1), uni_abaqus(:,2), 'Color', colors(4,:), 'LineWidth', 4);

% Generate model prediction
eps_model = linspace(min(uni_1.eps), max(uni_1.eps), 200);
lambda_model = exp(eps_model); % since ε_true = ln(λ)
sigma_model = ogden_model(params_fit, lambda_model, 'uniaxial');

% Plot model curve
plot(eps_model, sigma_model, 'k-', 'LineWidth', 2.5);

hold off; axis tight;
xlabel('True Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
ylabel('True Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
legend({'sample 1', 'sample 2', 'sample 3', 'simulation', 'model'}, 'Location', 'northwest', 'FontSize', 15);
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);
xlim([0,uni_3.eps(end)]);
grid off; set(gcf, 'Position', [500 50 800 600]); box on;

saveas(gcf, fullfile(outdir, '1.model_fit_comparison_uniaxial_true_Ogden.png'));

%% Planar shear
figure();
plot(plan_1.eps, plan_1.sigma, 'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(plan_2.eps, plan_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
plot(plan_3.eps, plan_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);

plot(plan_abaqus(:,1), plan_abaqus(:,2), 'Color', colors(4,:), 'LineWidth', 4);

eps_model = linspace(min(plan_1.eps), max(plan_1.eps), 200);
lambda_model = exp(eps_model);
sigma_model = ogden_model(params_fit, lambda_model, 'planar');
plot(eps_model, sigma_model, 'k-', 'LineWidth', 2.5);

hold off; axis tight;
xlabel('True Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
ylabel('True Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
legend({'sample 1', 'sample 2', 'sample 3', 'simulation', 'model'}, 'Location', 'northwest', 'FontSize', 15);
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);
xlim([0,plan_3.eps(end)]);
grid off; set(gcf, 'Position', [500 50 800 600]); box on;

saveas(gcf, fullfile(outdir, '2.model_fit_comparison_planarshear_true_Ogden.png'));

%% Compression
figure();
plot(comp_1.eps, comp_1.sigma, 'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(comp_2.eps, comp_2.sigma, 'Color', colors(2,:), 'LineWidth', 2);
plot(comp_3.eps, comp_3.sigma, 'Color', colors(3,:), 'LineWidth', 2);

plot(comp_abaqus(:,1), comp_abaqus(:,2), 'Color', colors(4,:), 'LineWidth', 4);

eps_model = linspace(min(comp_1.eps), max(comp_1.eps), 200);
lambda_model = exp(eps_model);
sigma_model = ogden_model(params_fit, lambda_model, 'compression');
plot(eps_model, sigma_model, 'k-', 'LineWidth', 2.5);

hold off; axis tight;
xlabel('True Strain', 'FontName', 'Times New Roman', 'FontSize', 15);
ylabel('True Stress [MPa]', 'FontName', 'Times New Roman', 'FontSize', 15);
legend({'sample 1', 'sample 2', 'sample 3', 'simulation', 'model'}, 'Location', 'southeast', 'FontSize', 15);
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 2);
xlim([comp_1.eps(end),0]);
grid off; set(gcf, 'Position', [500 50 800 600]); box on;

saveas(gcf, fullfile(outdir, '3.model_fit_comparison_compression_true_Ogden.png'));

%% ========================== MODEL COMPARISON (Engineering Stress–Strain) ==========================
% % Conversion: sigma_eng = sigma_true / lambda, eps_eng = lambda - 1
% 
% % ==== Uniaxial ====
% figure();
% plot(exp(uni_1.eps)-1, uni_1.sigma./exp(uni_1.eps), 'Color', colors(1,:), 'LineWidth', 2); hold on;
% plot(exp(uni_2.eps)-1, uni_2.sigma./exp(uni_2.eps), 'Color', colors(2,:), 'LineWidth', 2);
% plot(exp(uni_3.eps)-1, uni_3.sigma./exp(uni_3.eps), 'Color', colors(3,:), 'LineWidth', 2);
% 
% eps_model_eng = linspace(min(exp(uni_1.eps)-1), max(exp(uni_1.eps)-1), 200);
% lambda_model = 1 + eps_model_eng;
% sigma_true_model = ogden_model(params_fit, lambda_model, 'uniaxial');
% sigma_eng_model = sigma_true_model ./ lambda_model;
% plot(eps_model_eng, sigma_eng_model, 'k-', 'LineWidth', 2.5);
% axis tight;
% 
% xlabel('Engineering Strain', 'FontName','Times New Roman','FontSize',15);
% ylabel('Engineering Stress [MPa]', 'FontName','Times New Roman','FontSize',15);
% legend({'sample 1','sample 2','sample 3','Ogden model'},'FontSize',15,'Location','best');
% set(gca,'FontSize',15,'FontName','Times New Roman','LineWidth',2);
% grid off; box on;
% set(gcf,'Position',[400 50 800 600]);
% saveas(gcf, fullfile(outdir,'4.model_fit_comparison_uniaxial_eng_Ogden.png'));
% 
% % ==== Planar shear ====
% figure();
% plot(exp(plan_1.eps)-1, plan_1.sigma./exp(plan_1.eps), 'Color', colors(1,:), 'LineWidth', 2); hold on;
% plot(exp(plan_2.eps)-1, plan_2.sigma./exp(plan_2.eps), 'Color', colors(2,:), 'LineWidth', 2);
% plot(exp(plan_3.eps)-1, plan_3.sigma./exp(plan_3.eps), 'Color', colors(3,:), 'LineWidth', 2);
% 
% eps_model_eng = linspace(min(exp(plan_1.eps)-1), max(exp(plan_1.eps)-1), 200);
% lambda_model = 1 + eps_model_eng;
% sigma_true_model = ogden_model(params_fit, lambda_model, 'planar');
% sigma_eng_model = sigma_true_model ./ lambda_model;
% plot(eps_model_eng, sigma_eng_model, 'k-', 'LineWidth', 2.5);
% axis tight;
% 
% xlabel('Engineering Strain', 'FontName','Times New Roman','FontSize',15);
% ylabel('Engineering Stress [MPa]', 'FontName','Times New Roman','FontSize',15);
% legend({'sample 1','sample 2','sample 3','Ogden model'},'FontSize',15,'Location','best');
% set(gca,'FontSize',15,'FontName','Times New Roman','LineWidth',2);
% grid off; box on;
% set(gcf,'Position',[400 50 800 600]);
% saveas(gcf, fullfile(outdir,'5.model_fit_comparison_planarshear_eng_Ogden.png'));
% 
% % ==== Compression ====
% figure();
% plot(exp(comp_1.eps)-1, comp_1.sigma./exp(comp_1.eps), 'Color', colors(1,:), 'LineWidth', 2); hold on;
% plot(exp(comp_2.eps)-1, comp_2.sigma./exp(comp_2.eps), 'Color', colors(2,:), 'LineWidth', 2);
% plot(exp(comp_3.eps)-1, comp_3.sigma./exp(comp_3.eps), 'Color', colors(3,:), 'LineWidth', 2);
% axis tight;
% 
% eps_model_eng = linspace(min(exp(comp_1.eps)-1), max(exp(comp_1.eps)-1), 200);
% lambda_model = 1 + eps_model_eng;
% sigma_true_model = ogden_model(params_fit, lambda_model, 'compression');
% sigma_eng_model = sigma_true_model ./ lambda_model;
% plot(eps_model_eng, sigma_eng_model, 'k-', 'LineWidth', 2.5);
% 
% xlabel('Engineering Strain', 'FontName','Times New Roman','FontSize',15);
% ylabel('Engineering Stress [MPa]', 'FontName','Times New Roman','FontSize',15);
% legend({'sample 1','sample 2','sample 3','Ogden model'},'FontSize',15,'Location','best');
% set(gca,'FontSize',15,'FontName','Times New Roman','LineWidth',2);
% grid off; box on;
% set(gcf,'Position',[400 50 800 600]);
% saveas(gcf, fullfile(outdir,'6.model_fit_comparison_compression_eng_Ogden.png'));

%% ========================== EXPORT NOMINAL DATA FOR ABAQUS ==========================
% % Abaqus expects: Nominal Stress [MPa], Nominal Strain [unitless ratio]
% 
% % -------------------------------------------------------------------------
% % 1. Uniaxial tension
% eps_nom_uni = exp(uni_1.eps) - 1;      % convert from true → nominal
% sig_nom_uni = uni_1.sigma ./ (1 + eps_nom_uni);  % convert from true → nominal
% 
% uniaxial_nominal = [eps_nom_uni(:), sig_nom_uni(:)];  % [strain, stress]
% writematrix(uniaxial_nominal, fullfile(outdir2, 'uniaxial_nominal.csv'));
% 
% % -------------------------------------------------------------------------
% % 2. Compression
% eps_nom_comp = exp(comp_1.eps) - 1;    % true → nominal
% sig_nom_comp = comp_1.sigma ./ (1 + eps_nom_comp);
% 
% compression_nominal = [eps_nom_comp(:), sig_nom_comp(:)];
% writematrix(compression_nominal, fullfile(outdir2, 'compression_nominal.csv'));
% 
% % -------------------------------------------------------------------------
% % 3. Planar shear
% eps_nom_plan = exp(plan_1.eps) - 1;
% sig_nom_plan = plan_1.sigma ./ (1 + eps_nom_plan);
% 
% planarshear_nominal = [eps_nom_plan(:), sig_nom_plan(:)];
% writematrix(planarshear_nominal, fullfile(outdir2, 'planarshear_nominal.csv'));
% 
% % -------------------------------------------------------------------------
% % 4. Combined compression + tension (for Abaqus "Uniaxial Test Data")
% uniaxial_combined = [eps_nom_comp(:), sig_nom_comp(:); ...
%                      eps_nom_uni(:),  sig_nom_uni(:)];
% 
% % Sort by strain (ascending) to ensure continuity in Abaqus and plots
% uniaxial_combined = sortrows(uniaxial_combined, 1);
% 
% % Export combined file
% writematrix(uniaxial_combined, fullfile(outdir2, 'Uniaxial_comp_tensile_nominal.csv'));
% 
% % -------------------------------------------------------------------------
% fprintf('\n=== Export Complete ===\n');
% fprintf('Nominal data saved in folder: %s\n', outdir2);
% fprintf('Files:\n');
% fprintf(' - uniaxial_nominal.csv\n');
% fprintf(' - compression_nominal.csv\n');
% fprintf(' - planarshear_nominal.csv\n');
% fprintf(' - Uniaxial_comp_tensile_nominal.csv (combined + sorted)\n');
% fprintf('Columns: [Nominal Strain (ratio), Nominal Stress (MPa)]\n');

%% Import functions
% ================= UNAXIAL TENSION & PLANAR SHEAR =================
function data = import_and_match(line_file, ins_file, gauge, area, c1, c2)
    % Columns: time (s), distance (mm), distance (px), y1 (px), y2 (px)
    d1 = table2array(readtable(line_file));
    t1 = d1(:,1);
    dis_line = d1(:,2) - d1(1,2);      % displacement (mm)
    eps_eng_line = dis_line / gauge;   % engineering strain
    lambda = 1 + eps_eng_line;         % stretch ratio

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
function E_all = cost_fun_Ogden_all(params, data, weights)
    % params = [mu1 mu2 mu3 alpha1 alpha2 alpha3]
    % data: struct with fields data.uni, data.plan, data.comp
    % weights: weighting for each mode
    modes = {'uni', 'plan', 'comp'};
    E_all = [];

    for i = 1:3
        mode = modes{i};
        w = weights(i);

        for k = 1:length(data.(mode))
            sample = data.(mode){k};
            lambda = sample.lambda;
            sigma_exp = sample.sigma;
            sigma_model = ogden_model(params, lambda, mode);

            ni = length(sigma_exp);
            s_range = max(sigma_exp) - min(sigma_exp);
            if s_range == 0, s_range = 1; end

            e = sqrt(w/ni) * ((sigma_exp - sigma_model) ./ s_range);
            E_all = [E_all; e(:)];
        end
    end
end

%% ======================= OGDEN MODEL =======================
function sigma = ogden_model(params, lambda, mode)
    % params = [mu1 mu2 mu3 alpha1 alpha2 alpha3]
    mu = params(1:3);
    alpha = params(4:6);

    sigma = zeros(size(lambda));

    switch lower(mode)
        case {'uni','uniaxial'}
            % σ11 = Σ μp (λ^αp − λ^(−αp/2))
            for p = 1:3
                sigma = sigma + mu(p) .* (lambda.^alpha(p) - lambda.^(-alpha(p)/2));
            end

        case {'plan','planar'}
            % σ11 = Σ μp (λ^αp − λ^(−αp))
            for p = 1:3
                sigma = sigma + mu(p) .* (lambda.^alpha(p) - lambda.^(-alpha(p)));
            end

        case {'comp','compression'}
            % Compression follows same form as uniaxial (λ < 1)
            for p = 1:3
                sigma = sigma + mu(p) .* (lambda.^alpha(p) - lambda.^(-alpha(p)/2));
            end

        otherwise
            error('Unknown mode: %s', mode);
    end
end