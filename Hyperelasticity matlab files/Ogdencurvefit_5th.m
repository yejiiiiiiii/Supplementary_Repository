% Yeji Han
% Curvefit for Ogden (5th-order) model parameters — Uniaxial, Planar, Compression
clear; clc; close all force

colors = lines(6);

% --- Geometry -------------------------------------------------------------
area_uni    = 2 * 6;        gauge_uni    = 30;      % Uniaxial
area_planar = 1.5 * 215;    gauge_planar = 10;      % Planar shear
diam_comp   = 28.6;         area_comp    = pi*(diam_comp/2)^2;  gauge_comp = 12.5; % Compression

outdir = 'plots_ogden_5th';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ========================== Import data =================================
% Uniaxial
uni_1 = import_and_match('uniaxial_line1.csv','uniaxial_ins1.csv',gauge_uni,area_uni,0.1,0.1);
uni_2 = import_and_match('uniaxial_line2.csv','uniaxial_ins2.csv',gauge_uni,area_uni,0.5,0.5);
uni_3 = import_and_match('uniaxial_line3.csv','uniaxial_ins3.csv',gauge_uni,area_uni,0.205,0.205);
uni_1 = zero_data(uni_1); uni_2 = zero_data(uni_2); uni_3 = zero_data(uni_3);
uni_abaqus = table2array(readtable("BBDINO_Abaqus_Tension_Ogden5.csv"));

% Planar shear
plan_1 = import_and_match('plan_line1.csv','plan_ins1.csv',gauge_planar,area_planar,0.3,0.2);
plan_2 = import_and_match('plan_line2.csv','plan_ins2.csv',gauge_planar,area_planar,0.3,0.2);
plan_3 = import_and_match('plan_line3.csv','plan_ins3.csv',gauge_planar,area_planar,0.5,0.2);
plan_1 = zero_data(plan_1); plan_2 = zero_data(plan_2); plan_3 = zero_data(plan_3);
plan_abaqus = table2array(readtable("BBDINO_Abaqus_PlanarShear_Ogden5.csv"));

% Compression
comp_1 = import_compression('comp_ins1.csv',area_comp,gauge_comp);
comp_2 = import_compression('comp_ins2.csv',area_comp,gauge_comp);
comp_3 = import_compression('comp_ins3.csv',area_comp,gauge_comp);
comp_1 = zero_data(comp_1); comp_2 = zero_data(comp_2); comp_3 = zero_data(comp_3);
comp_abaqus = table2array(readtable("BBDINO_Abaqus_Compression_Ogden5.csv"));

%% ========================== Curve fitting ===============================
data.uni  = {uni_1,uni_2,uni_3};
data.plan = {plan_1,plan_2,plan_3};
data.comp = {comp_1,comp_2,comp_3};
weights   = [1 1 1];

% Initial guess: [mu1 mu2 mu3 mu4 alpha1 alpha2 alpha3 alpha4]
lb = [zeros(1,5), -10*ones(1,5)];
ub = [ones(1,5)*10,  10*ones(1,5)];

params0 = [0.005373 0.005659 0.000008 0.067388 0.137987 1.668473 1.664906 1.807965 1.892932 1.995314];

options = optimoptions('lsqnonlin',...
'Display','iter',...
'Algorithm','levenberg-marquardt',...
'MaxIterations',5000,...
'MaxFunctionEvaluations',50000);

params_fit = lsqnonlin(@(p) cost_fun_Ogden_all(p,data,weights),params0,[],[],options);

N = numel(params_fit)/2;
mu    = params_fit(1:N);
alpha = params_fit(N+1:end);
fprintf('\nFitted Ogden-5 Parameters:\n');
for i=1:N
    fprintf('μ%d = %.6f\t α%d = %.6f\n',i,mu(i),i,alpha(i));
end

fprintf('\nFitted Ogden-5 Parameters for Abaqus:\n');
for i=1:N
    fprintf('μ%d = %.6f\t α%d = %.6f\n',i,mu(i)*alpha(i)/2,i,alpha(i));
end

%% ========================== Model Comparison (True Stress) ==============
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
% legend({'sample 1', 'sample 2', 'sample 3', 'model'}, 'Location', 'best', 'FontSize', 15);
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
% legend({'sample 1', 'sample 2', 'sample 3', 'model'}, 'Location', 'best', 'FontSize', 15);
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

%% ========================================================================
% --------------------------- Helper Functions ---------------------------
function E_all = cost_fun_Ogden_all(params,data,weights)
    modes = {'uni','plan','comp'};
    E_all = [];
    for i = 1:3
        mode = modes{i}; w = weights(i);
        for k = 1:length(data.(mode))
            d = data.(mode){k};
            sig_exp = d.sigma;
            sig_model = ogden_model(params,d.lambda,mode);
            ni = length(d.sigma);
            s_range = max(sig_exp)-min(sig_exp); if s_range==0, s_range=1; end
            e = sqrt(w/ni) * ((sig_exp - sig_model) ./ s_range);
            E_all = [E_all; e(:)];
        end
    end
end

function sigma = ogden_model(params,lambda,mode)
    N = numel(params)/2;
    mu = params(1:N); alpha = params(N+1:end);
    sigma = zeros(size(lambda));
    switch lower(mode)
        case {'uni','uniaxial'}
            for p=1:N
                sigma = sigma + mu(p).*(lambda.^alpha(p) - lambda.^(-alpha(p)/2));
            end
        case {'plan','planar'}
            for p=1:N
                sigma = sigma + mu(p).*(lambda.^alpha(p) - lambda.^(-alpha(p)));
            end
        case {'comp','compression'}
            for p=1:N
                sigma = sigma + mu(p).*(lambda.^alpha(p) - lambda.^(-alpha(p)/2));
            end
        otherwise
            error('Unknown mode: %s',mode);
    end
end

function data = import_and_match(line_file,ins_file,gauge,area,c1,c2)
    d1 = table2array(readtable(line_file));
    t1 = d1(:,1); dis_line = d1(:,2)-d1(1,2);
    eps_eng_line = dis_line/gauge; lambda = 1+eps_eng_line;

    d2 = table2array(readtable(ins_file));
    t2 = d2(:,1); dis_ins = d2(:,2)-d2(1,2); F = d2(:,3);

    vel_line = gradient(dis_line); vel_ins = gradient(dis_ins);
    th1 = c1*max(vel_line); th2 = c2*max(vel_ins);
    idx1 = find(vel_line>th1,1); idx2 = find(vel_ins>th2,1);
    t1 = t1 - t1(idx1); t2 = t2 - t2(idx2);
    t1 = t1(idx1:end); dis_line = dis_line(idx1:end); lambda=lambda(idx1:end);
    t2 = t2(idx2:end); F = F(idx2:end);
    [~,idxp]=max(dis_line); t1=t1(1:idxp); dis_line=dis_line(1:idxp); lambda=lambda(1:idxp);
    N_target=min(length(t1),length(t2)); t_common=linspace(0,min(t1(end),t2(end)),N_target);
    dis_line=interp1(t1,dis_line,t_common,'linear','extrap');
    lambda=interp1(t1,lambda,t_common,'linear','extrap');
    F=interp1(t2,F,t_common,'linear','extrap');
    eps_true=log(lambda); sigma_true=(F.*lambda)/area;
    data=struct('t',t_common(:),'dis',dis_line(:),'lambda',lambda(:), ...
                'eps',eps_true(:),'sigma',sigma_true(:),'F',F(:));
end

function data = import_compression(ins_file,area_comp,gauge_comp)
    d=table2array(readtable(ins_file));
    t=d(:,1); dis=-d(:,2); F=-d(:,3);
    lambda=1+dis/gauge_comp; eps_true=log(lambda);
    sigma_true=(F.*lambda)/area_comp;
    data=struct('t',t(:),'dis',dis(:),'lambda',lambda(:), ...
                'eps',eps_true(:),'sigma',sigma_true(:),'F',F(:));
end

function data = zero_data(data)
    if isfield(data,'t'),data.t=data.t-data.t(1);end
    if isfield(data,'dis'),data.dis=data.dis-data.dis(1);end
    if isfield(data,'lambda'),data.lambda=data.lambda-data.lambda(1)+1;end
    if isfield(data,'eps'),data.eps=data.eps-data.eps(1);end
    if isfield(data,'sigma'),data.sigma=data.sigma-data.sigma(1);end
    if isfield(data,'F'),data.F=data.F-data.F(1);end
end

function s = capitalizeFirst(str)
    s = lower(str); s(1) = upper(s(1));
end
