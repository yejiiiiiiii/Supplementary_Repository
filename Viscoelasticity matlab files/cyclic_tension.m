% Yeji Han

clear
clc
close all

% Output folder for figures
outdir = 'plots';
if ~exist(outdir,'dir'); mkdir(outdir); end

%% Load data
% Time(sec) Displacement(mm) Force(N)
d1 = table2array(readtable('BBDINO30_CyclicTension (1).csv')); % swap to another dataset (2 or 3) if needed
d2 = table2array(readtable('BBDINO30_CyclicTension (2).csv'));
d3 = table2array(readtable('BBDINO30_CyclicTension (3).csv'));

%% Extract Top Peaks and Curve Fit

% Extract peak values from the force signal (column 3)
time = d1(:,1);
force = d1(:,3);

% Find local maxima (top peaks)
[peak_vals, peak_locs] = findpeaks(force, 'MinPeakHeight', 0.3, 'MinPeakDistance', 10);

peak_times = time(peak_locs);

%% Sampling rate
dt = mean(diff(time));
fs = 1/dt;
fprintf('Sampling rate: %.2f Hz\n', fs);

T_approx = 2.5;  
min_dist_samples = round(0.7 * T_approx / dt);  

%% Peaks
[peak_vals, peak_locs] = findpeaks(force, ...
    'MinPeakHeight',    0.35, ...          
    'MinPeakDistance',  min_dist_samples); 

peak_times = time(peak_locs);

fprintf('Number of peaks detected: %d\n', length(peak_vals));

%% Curve fitting
ft   = fittype('a*exp(b*x) + c', 'independent', 'x');
opts = fitoptions(ft);

opts.StartPoint = [0.1, -0.05, 0.35]; 
opts.Lower      = [-Inf, -Inf,  0  ];
opts.Upper      = [ Inf,  0,    Inf];   % b must be negative (decay)

[f_exp, gof_exp] = fit(peak_times, peak_vals, ft, opts);

%% Results
first_peak_time  = peak_times(1);
first_peak_value = peak_vals(1);
last_peak_time   = peak_times(end);
last_peak_value  = peak_vals(end);

fprintf('\n--- Results ---\n');
fprintf('First Peak : %.4f N  at t = %.2f sec\n', first_peak_value, first_peak_time);
fprintf('Last  Peak : %.4f N  at t = %.2f sec\n', last_peak_value,  last_peak_time);
fprintf('Exp fit    : a=%.4f  b=%.6f  c=%.4f  R²=%.4f\n', ...
    f_exp.a, f_exp.b, f_exp.c, gof_exp.rsquare);

% Predicted value at t=0 from exp fit
fprintf('Exp fit @ t=0  : %.4f N\n', f_exp.a + f_exp.c);
fprintf('Exp fit @ t=end: %.4f N\n', f_exp(peak_times(end)));

%% Plot
figure()
plot(time, force, 'r-', 'LineWidth', 1, 'DisplayName','Experiment');

xlabel('Time [sec]', 'FontSize', 25, 'FontName', 'Times new roman');
ylabel('Force [N]',  'FontSize', 25, 'FontName', 'Times new roman');

xlim([0 d1(end,1)]);  
ylim([-1, 5]);
ax = gca; set(ax, 'FontSize', 25, 'LineWidth', 1.5, 'FontName', 'Times new roman');
set(gcf, 'Position', [500 0 800 600]);

saveas(gcf, fullfile(outdir, '1. cyclictension_forcetime.png'));
saveas(gcf, fullfile(outdir, '1. cyclictension_forcetime.fig'));

%% Plot - curvefit
figure()
plot(time, force, 'r-', 'LineWidth', 1, 'DisplayName','Experiment');
hold on
plot(peak_times, peak_vals,          'bo',  'MarkerSize', 6,  'LineWidth', 1.5, 'DisplayName', 'Peaks');
plot(peak_times, f_exp(peak_times),  'g-',  'LineWidth', 1, 'DisplayName', sprintf('Curvefit'));

plot(first_peak_time, first_peak_value, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
    'DisplayName', sprintf('First: %.3f N', first_peak_value));
plot(last_peak_time, last_peak_value,   'ms', 'MarkerSize', 12, 'MarkerFaceColor', 'm', ...
    'DisplayName', sprintf('Last:  %.3f N', last_peak_value));
hold off

xlabel('Time [sec]', 'FontSize', 25, 'FontName', 'Times new roman');
ylabel('Force [N]',  'FontSize', 25, 'FontName', 'Times new roman');

legend('Location', 'northeast', 'FontSize', 25, 'FontName', 'Times new roman');
xlim([0 d1(end,1)]);  
ylim([3, 4.5]);
yticks(3:0.25:4.5)
ax = gca; set(ax, 'FontSize', 25, 'LineWidth', 1.5, 'FontName', 'Times new roman');
set(gcf, 'Position', [500 0 800 600]);

saveas(gcf, fullfile(outdir, '2. cyclictension_forcetime_curvefit.png'));
saveas(gcf, fullfile(outdir, '2. cyclictension_forcetime_curvefit.fig'));