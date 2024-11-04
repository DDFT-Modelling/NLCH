%% Errors for different values of ε and factors (grids in MS) for maximal sensible subdivision
% Warning: If you run the following code from scratch, it will take a
% couple of hours to finish. Average and total computing times are reported
% for each experiment.

% Extremely recomend to store working space: 'autosave(10,'workspace.mat')'
% Command to be stopped with 'autosave stop' followed by 'autosave delete'

%% ----------------------------------------- n = 10, n_i = 10 ----------------------------------------- %
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 24.9 s, mean 1.25 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 10, 10, 1);
    T(k)= toc;  % pair 1: toc
    
    Abs_10_10_a{k} = errors;
    Rel_10_10_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_10_10_b{k} = experiment;
    Rel_10_10_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_10_10_b{k} = smaller;

    Abs_10_10_c{k} = full_thing;
    Rel_10_10_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
    Exact_Full_10_10{k} = all_sq;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(10, 10, 20, Range_Epsilons, Abs_10_10_a, Abs_10_10_b, Abs_10_10_c, Exact_10_10_b, 1, true)

%% ----------------------------------------- n = 10, n_i = 20 ----------------------------------------- %
%[logspace(-7,-1, 100)*5.85/2.01, linspace(0.585/2,0.585,100)];   %logspace(-10,-1, 100)*5.85; %linspace(1e-10,0.585,200);
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 2 min 8.25 s, mean 6.41 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 10, 10, 2);
    T(k)= toc;  % pair 1: toc
    
    Abs_10_20_a{k} = errors;
    Rel_10_20_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_10_20_b{k} = experiment;
    Rel_10_20_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_10_20_b{k} = smaller;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(10, 10, 20, Range_Epsilons, Abs_10_20_a, Abs_10_20_b, [], Exact_10_20_b, 2, true)
close all

%% ----------------------------------------- n = 10, n_i = 40 ----------------------------------------- %
% 0.59 is like, the limit
%[logspace(-7,-1, 100)*5.85/2.01, linspace(0.585/2,0.585,100)];   %logspace(-10,-1, 100)*5.85; %linspace(1e-10,0.585,200);
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 5 min 55 s, mean 17.74 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 10, 10, 4);
    T(k)= toc;  % pair 1: toc
    
    Abs_10_40_a{k} = errors;
    Rel_10_40_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_10_40_b{k} = experiment;
    Rel_10_40_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_10_40_b{k} = smaller;

    Abs_10_40_c{k} = full_thing;
    Rel_10_40_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
    Exact_Full_10_40{k} = all_sq;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(10, 10, 20, Range_Epsilons, Abs_10_40_a, Abs_10_40_b, Abs_10_40_c, Exact_10_40_b, 4, true)
close all

%% ----------------------------------------- n = 10, n_i = 60 ----------------------------------------- %

k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 12 min 24.43 s, mean 37.22 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 10, 10, 6);
    T(k)= toc;  % pair 1: toc
    
    Abs_10_60_a{k} = errors;
    Rel_10_60_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_10_60_b{k} = experiment;
    Rel_10_60_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_10_60_b{k} = smaller;

    Abs_10_60_c{k} = full_thing;
    Rel_10_60_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
    Exact_Full_10_60{k} = all_sq;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(10, 10, 20, Range_Epsilons, Abs_10_60_a, Abs_10_60_b, Abs_10_60_c, Exact_10_60_b, 6, true)
close all

%% ----------------------------------------- n = 10, n_i = 80 ----------------------------------------- %
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 22 min 39 s, mean 1 min 8 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 10, 10, 8);
    T(k)= toc;  % pair 1: toc
    
    Abs_10_80_a{k} = errors;
    Rel_10_80_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_10_80_b{k} = experiment;
    Rel_10_80_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_10_80_b{k} = smaller;

    Abs_10_80_c{k} = full_thing;
    Rel_10_80_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
    Exact_Full_10_80{k} = all_sq;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(10, 10, 20, Range_Epsilons, Abs_10_80_a, Abs_10_80_b, Abs_10_80_c, Exact_10_80_b, 8, true)
close all

%% ----------------------------------------- n = 10, n_i = 100 ---------------------------------------- %
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 39 min 57.59 s, mean 1 min 59.88 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 10, 10, 10);
    T(k)= toc;  % pair 1: toc
    
    Abs_10_100_a{k} = errors;
    Rel_10_100_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_10_100_b{k} = experiment;
    Rel_10_100_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_10_100_b{k} = smaller;

    Abs_10_100_c{k} = full_thing;
    Rel_10_100_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
    Exact_Full_10_100{k} = all_sq;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(10, 10, 20, Range_Epsilons, Abs_10_100_a, Abs_10_100_b, Abs_10_100_c, Exact_10_100_b, 10, true)
close all

%% ----------------------------------------- n = 10, n_i = 160 ---------------------------------------- %
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 1 h 15 min 40.03 s, mean 3 min 47.0017 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 10, 10, 16);
    T(k)= toc;  % pair 1: toc
    
    Abs_10_160_a{k} = errors;
    Rel_10_160_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_10_160_b{k} = experiment;
    Rel_10_160_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_10_160_b{k} = smaller;

    Abs_10_160_c{k} = full_thing;
    Rel_10_160_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
    Exact_Full_10_160{k} = all_sq;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(10, 10, 20, Range_Epsilons, Abs_10_160_a, Abs_10_160_b, Abs_10_160_c, Exact_10_160_b, 16, true)
close all

%% Plot max errors

% Obtain max error for each series
All_Errs = zeros([20,7]);
for i = 1:20
    All_Errs(i,1) = max(Abs_10_10_a{i});
    All_Errs(i,2) = max(Abs_10_20_a{i});
    All_Errs(i,3) = max(Abs_10_40_a{i});
    All_Errs(i,4) = max(Abs_10_60_a{i});
    All_Errs(i,5) = max(Abs_10_80_a{i});
    All_Errs(i,6) = max(Abs_10_100_a{i});
    All_Errs(i,7) = max(Abs_10_160_a{i});
end
[MinA, MaxA] = bounds(All_Errs, 'all');

% Obtain truncation error
Exact_v = zeros([1,20]);
for i = 1:20
    Exact_v(i) = max(abs(Exact_10_80_b{i}));
end

% Colour selection
colour_A = [255, 225, 168]/255;
colour_B = [114, 61, 70]/255;
GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
G_colours = GRADIENT_flexible((1:7),7);

% Plot things
figure(1)
plot(Range_Epsilons, All_Errs, 'LineWidth',2);    hold on
yscale log
xscale log
ylim([MinA, MaxA].* [0.9, 1.1] );
xlim([Range_Epsilons(1), Range_Epsilons(end)]);
colororder(G_colours)
set(gca, 'FontName', 'CMR10')
xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
ylabel('Absolute errors','Interpreter','latex')
set(gca, 'YTick', 10.^(-20:1:0))

% Truncation error
Bound_G = @(x) x.^2 .* abs( -6 + log(4) + pi + 8 * log(x) ) / (2*pi);       % [0.5,0.5]
plot(Range_Epsilons, max( Exact_v, Bound_G(Range_Epsilons)), 'black','LineWidth',2,'LineStyle','--')

% Add some nice lines
for i = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    if i < 1e-2
        plot([Range_Epsilons(1),i],[Bound_G(i),Bound_G(i)], 'Color', [1,1,1]/1.5);
    end
    plot([i,i], [MinA, Bound_G(i)], 'Color', [1,1,1]/1.5);
end

lgd = legend({'1','2','4','6','8', '10','16'},'Location','SouthWest');
lgd.Title.String = 'Factors';

exportgraphics(figure(1), 'Errors_Max_Factor[10].pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300);
hold off

% Dividing the factor 2 steps:
% All_Errs(1,[1,2,3,5])./All_Errs(1,[2,3,5,7])
% > 23.3098   21.3960   21.0828   21.7605
% So each factor of 2 is around 20 times more accurate!

%% Plot modes of errors
% Obtain max error for each series
All_Errs = zeros([20,6]);
for i = 1:20
    All_Errs(i,1) = mode(Abs_10_10_a{i});
    All_Errs(i,2) = mode(Abs_10_20_a{i});
    All_Errs(i,3) = mode(Abs_10_40_a{i});
    All_Errs(i,4) = mode(Abs_10_60_a{i});
    All_Errs(i,5) = mode(Abs_10_80_a{i});
    All_Errs(i,6) = mode(Abs_10_100_a{i});
    All_Errs(i,7) = mode(Abs_10_160_a{i});
end
[MinA, MaxA] = bounds(All_Errs(All_Errs>0), 'all');

% Obtain truncation error
Exact_v = zeros([1,20]);
for i = 1:20
    Exact_v(i) = mode(abs(Exact_10_80_b{i}));
end

% Plot things
figure(2)
plot(Range_Epsilons, All_Errs, 'LineWidth',2);    hold on
yscale log
xscale log
ylim([MinA, MaxA].* [0.9, 1.1] );
xlim([Range_Epsilons(1), Range_Epsilons(end)]);
colororder(G_colours)
set(gca, 'FontName', 'CMR10')
xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
ylabel('Absolute errors','Interpreter','latex')
set(gca, 'YTick', 10.^(-20:1:0))

plot(Range_Epsilons, Exact_v, 'black','LineWidth',2,'LineStyle','--')
plot(Range_Epsilons, Bound_G(Range_Epsilons), 'Color',[1,1,1]/1.25,'LineWidth',1.5,'LineStyle','--')

% Add some nice lines
for i = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    if i < 1e-2
        plot([Range_Epsilons(1),i],[Bound_G(i),Bound_G(i)], 'Color', [1,1,1]/1.5);
    end
    plot([i,i], [MinA, Bound_G(i)], 'Color', [1,1,1]/1.5);
end

lgd = legend({'1','2','4','6','8', '10','16'},'Location','SouthWest');
lgd.Title.String = 'Factors';

exportgraphics(figure(2), 'Errors_Mode_Factor[10].pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300);
hold off

%%
%% ----------------------------------------- n = 20, n_i = 20 ----------------------------------------- %
% Behaves like 10-20
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 7 min 7 s, mean 21.4 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 20, 20, 1);
    T(k)= toc;  % pair 1: toc
    
    Abs_20_20_a{k} = errors;
    Rel_20_20_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_20_20_b{k} = experiment;
    Rel_20_20_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_20_20_b{k} = smaller;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(20, 20, 20, Range_Epsilons, Abs_20_20_a, Abs_20_20_b, [], Exact_20_20_b, 1, true)
close all

%% ----------------------------------------- n = 20, n_i = 40 ----------------------------------------- %
% Behaves like 10-40
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 27 min 22.5 s, mean 1 min 22.1 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 20, 20, 2);
    T(k)= toc;  % pair 1: toc
    
    Abs_20_40_a{k} = errors;
    Rel_20_40_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_20_40_b{k} = experiment;
    Rel_20_40_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_20_40{k} = smaller;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(20, 20, 20, Range_Epsilons, Abs_20_40_a, Abs_20_40_b, [], Exact_20_40, 2, true)
close all

%% ----------------------------------------- n = 20, n_i = 60 ----------------------------------------- %
% Behaves like 10-60
k = 1;
T = zeros(1,20);
Range_Epsilons = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 57 min 28.6 s, mean 3 min 52 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 20, 20, 3);
    T(k)= toc;  % pair 1: toc
    
    Abs_20_60_a{k} = errors;
    Rel_20_60_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_20_60_b{k} = experiment;
    Rel_20_60_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_20_60{k} = smaller;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(20, 20, 20, Range_Epsilons, Abs_20_60_a, Abs_20_60_b, [], Exact_20_60, 3, true)
close all

%% ----------------------------------------- n = 20, n_i = 80 ----------------------------------------- %

k = 1;
T = zeros(1,10);
Range_Epsilons = [logspace(-7,-3, 5), logspace(-3,-1, 5)*5.85/2.01];    %[logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 55.6 min, mean 5 min 45 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 20, 20, 4);
    T(k)= toc;  % pair 1: toc
    
    Abs_20_80_a{k} = errors;
    Rel_20_80_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_20_80_b{k} = experiment;
    Rel_20_80_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_20_80_b{k} = smaller;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(20, 20, 10, Range_Epsilons, Abs_20_80_a, Abs_20_80_b, [], Exact_20_80_b, 4, true)
close all

%% ----------------------------------------- n = 20, n_i = 160 ---------------------------------------- %

k = 1;
T = zeros(1,10);
Range_Epsilons = [logspace(-7,-3, 5), logspace(-3,-1, 5)*5.85/2.01];    %[logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
for i = Range_Epsilons % 3 h 14 min 15 s, mean 19 min 25.48 s
    tic;
    [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(i, 20, 20, 8);
    T(k)= toc;  % pair 1: toc
    
    Abs_20_160_a{k} = errors;
    Rel_20_160_a{k} = errors ./ max( abs(t_Full), 1e-6);

    Abs_20_160_b{k} = experiment;
    Rel_20_160_b{k} = experiment ./ max( abs(smaller), 1e-6);
    Exact_20_160_b{k} = smaller;

    k = k + 1;
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_eps(20, 20, 10, Range_Epsilons, Abs_20_160_a, Abs_20_160_b, [], Exact_20_160_b, 8, true)

%% Plot max errors

Range_Epsilons_a = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
% Obtain max error for each series
All_Errs_a = zeros([20,3]);
All_Errs_b = zeros([10,2]);
for i = 1:20
    All_Errs_a(i,1) = max(Abs_20_20_a{i});
    All_Errs_a(i,2) = max(Abs_20_40_a{i});
    All_Errs_a(i,3) = max(Abs_20_60_a{i});
end
for i = 1:10
    All_Errs_b(i,1) = max(Abs_20_80_a{i});
    All_Errs_b(i,2) = max(Abs_20_160_a{i});
end
[MinA_a, MaxA_a] = bounds(All_Errs_a, 'all');
[MinA_b, MaxA_b] = bounds(All_Errs_b, 'all');       % Quick eval says no need

% Obtain truncation error
Exact_v = zeros([1,20]);
for i = 1:20
    Exact_v(i) = max(abs(Exact_20_60{i}));
end

% Colour selection
colour_A = [255, 225, 168]/255;
colour_B = [114, 61, 70]/255;
GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
G_colours = GRADIENT_flexible((1:7),7);


% Plot things
figure(1)
plot(Range_Epsilons_a, All_Errs_a, 'LineWidth',2);    hold on
plot(Range_Epsilons, All_Errs_b, 'LineWidth',2);  % Add two additional lines
yscale log
xscale log
ylim([MinA_a, MaxA_a].* [0.9, 1.1] );
xlim([Range_Epsilons_a(1), Range_Epsilons_a(end)]);
colororder(G_colours([2,3,4,6,7],:))
set(gca, 'FontName', 'CMR10')
xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
ylabel('Absolute errors','Interpreter','latex')
set(gca, 'YTick', 10.^(-20:1:0))

% Truncation error
Bound_G = @(x) x.^2 .* abs( -6 + log(4) + pi + 8 * log(x) ) / (2*pi);       % [0.5,0.5]
plot(Range_Epsilons_a, max( Exact_v, Bound_G(Range_Epsilons_a)), 'black','LineWidth',2,'LineStyle','--')

% Add some nice lines
for i = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    if i < 1e-2
        plot([Range_Epsilons_a(1),i],[Bound_G(i),Bound_G(i)], 'Color', [1,1,1]/1.5);
    end
    plot([i,i], [MinA_a, Bound_G(i)], 'Color', [1,1,1]/1.5);
end

lgd = legend({'1','2','3','4','8'},'Location','SouthWest');
lgd.Title.String = 'Factors';

exportgraphics(figure(1), 'Errors_Max_Factor[20].pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300);
hold off

%% Plot mode errors

Range_Epsilons_a = [logspace(-7,-2, 12), logspace(-2,-1, 8)*5.85/2.01];
% Obtain max error for each series
All_Errs_a = zeros([20,3]);
All_Errs_b = zeros([10,2]);
for i = 1:20
    All_Errs_a(i,1) = mode(Abs_20_20_a{i});
    All_Errs_a(i,2) = mode(Abs_20_40_a{i});
    All_Errs_a(i,3) = mode(Abs_20_60_a{i});
end
for i = 1:10
    All_Errs_b(i,1) = mode(Abs_20_80_a{i});
    All_Errs_b(i,2) = mode(Abs_20_160_a{i});
end
[MinA_a, MaxA_a] = bounds(All_Errs_a, 'all');
[MinA_b, MaxA_b] = bounds(All_Errs_b, 'all');       % Quick eval says no need

% Obtain truncation error
Exact_v = zeros([1,20]);
for i = 1:20
    Exact_v(i) = max(abs(Exact_20_60{i}));
end

% Colour selection
colour_A = [255, 225, 168]/255;
colour_B = [114, 61, 70]/255;
GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
G_colours = GRADIENT_flexible((1:7),7);


% Plot things
figure(2)
plot(Range_Epsilons_a, All_Errs_a, 'LineWidth',2);    hold on
plot(Range_Epsilons, All_Errs_b, 'LineWidth',2);  % Add two additional lines
yscale log
xscale log
ylim([MinA_a, MaxA_a].* [0.9, 1.1] );
xlim([Range_Epsilons_a(1), Range_Epsilons_a(end)]);
colororder(G_colours([2,3,4,6,7],:))
set(gca, 'FontName', 'CMR10')
xlabel('Neighbourhood radius $\varepsilon$','Interpreter','latex')
ylabel('Absolute errors','Interpreter','latex')
set(gca, 'YTick', 10.^(-20:1:0))

% Truncation error
Bound_G = @(x) x.^2 .* abs( -6 + log(4) + pi + 8 * log(x) ) / (2*pi);       % [0.5,0.5]
plot(Range_Epsilons_a, max( Exact_v, Bound_G(Range_Epsilons_a)), 'black','LineWidth',2,'LineStyle','--')

% Add some nice lines
for i = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    if i < 1e-2
        plot([Range_Epsilons_a(1),i],[Bound_G(i),Bound_G(i)], 'Color', [1,1,1]/1.5);
    end
    plot([i,i], [MinA_a, Bound_G(i)], 'Color', [1,1,1]/1.5);
end

lgd = legend({'1','2','3','4','8'},'Location','SouthWest');
lgd.Title.String = 'Factors';

exportgraphics(figure(2), 'Errors_Mode_Factor[20].pdf', 'BackgroundColor','none', 'ContentType', 'vector', 'Resolution', 300);
hold off

%% It seems that times do not follow a 2^fact pattern:
plot([1,2,3,4,5,6],[3.4*60, 14.3*60, 31 * 60, 55.6*60, 3600 + 18.9933*60, 3600 + 59*60 + 6], '-x')
hold on
plot(1:0.1:6,(14.3*60) * 2.^(-1:0.1:4), '--')
yscale log
% if we fit, it is predicting 219.2799  * 7^1.9373 / 3600 = 2 h 38 min 31 s
% and 3 h 25 min 18 s for 8

% Maybe let's rather compute these fits by mean:
plot(1:7, 60*[0,0,1,2,3,5,7] + [10,42.99, 33.1, 47,57, 57, 23])
hold on
plot(1:0.1:7,(42.99) * 2.^(-1:0.1:5), '--')
yscale log
% 12.3185 * 8^1.8410 / 60 ~ 9 min 26 s [close]

%%
%% Fixed radius evaluation
% Overall the plots showcase a fair approximation for ε ≥ 1e-2. We observed
% modal values of around 1e-16 from not very complicated geometries (in
% terms of computational time). Let us recover such values. If we want to
% run smaller radii, we would have to sacrifice the error and increase the
% number of subdivisions in each MS. We will test this with some numerical 
% experiments.

%% 10 x 10 base grid
% Although we could certainly increase the refinement for small ε, since we
% want to compare against the 20 x 20 grid, we will use the same
% refinements

% The following takes 11 min 30.69 s
k = 1;
T = zeros(1,4);
Factors = [2,10,16];
Range_Epsilons = [1e-2, 1e-5];
for eps = Range_Epsilons 
    eps
    for i = Factors
        i
        tic;
        [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(eps, 10, 10, i);
        T(k)= toc;  % pair 1: toc
        
        Abs_10_a{k} = errors;
        Rel_10_a{k} = errors ./ max( abs(t_Full), 1e-6);
    
        Abs_10_b{k} = experiment;
        Rel_10_b{k} = experiment ./ max( abs(smaller), 1e-6);
        Exact_10_b{k} = smaller;
    
        Abs_10_c{k} = full_thing;
        Rel_10_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
        Exact_Full_10{k} = all_sq;
        
        Conv_Out_10{k} = Conv;
        k = k + 1;
    end
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_Facts(10, 10, 6, Range_Epsilons, Abs_10_a, Abs_10_b, Abs_10_c, Exact_10_b, Factors, true)

%% 20 x 20 base grid

% The following takes 56 min 21.06 s
k = 1;
T = zeros(1,4);
Factors = [1,5,8];
Range_Epsilons = [1e-2, 1e-5];
for eps = Range_Epsilons 
    eps
    for i = Factors
        i
        tic;
        [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(eps, 20, 20, i);
        T(k)= toc;  % pair 1: toc
        
        Abs_20_a{k} = errors;
        Rel_20_a{k} = errors ./ max( abs(t_Full), 1e-6);
    
        Abs_20_b{k} = experiment;
        Rel_20_b{k} = experiment ./ max( abs(smaller), 1e-6);
        Exact_20_b{k} = smaller;
    
        Abs_20_c{k} = full_thing;
        Rel_20_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
        Exact_Full_20{k} = all_sq;
        
        Conv_Out_20{k} = Conv;
        k = k + 1;
    end
end
[sum(T), mean(T)]

Plot_Errors_Singular_Conv_Facts(20, 20, 6, Range_Epsilons, Abs_20_a, Abs_20_b, Abs_20_c, Exact_20_b, Factors, true)

% The full integral will be approximated by [Conv * f] + [f ○ smaller]. As
% a result, we will store these two structures 

%% Store struct with outputs

%% 20 x 20 boxes
% 20 A: 1e-2
epsilon = Range_Epsilons(1);
Newtonians = Newtonian([aBox.Pts.y1_kv, aBox.Pts.y2_kv], epsilon);
Newtonians.eps = epsilon;  % This is ε
Newtonians.box = aBox;  % This is the domain
%Newtonians.NI; % is the convolution outside [0,1]^2 ∩ N(ε)
%Newtonians.NJ; % is the convolution [0,1]^2
%Newtonians.NG; % is the convolution in N(ε), this is `smaller`
Newtonians.N = 20;      % Number of points per dimension

% Now we add a substruct
Newtonians.Level = {};
Newtonians.Level.Available = Factors;     % Number of subdivisions per dimension
% Store Convolution matrices
Newtonians.Level.n2  = Conv_Out_20{1};
Newtonians.Level.n5 = Conv_Out_20{2};
Newtonians.Level.n8 = Conv_Out_20{3};

save('Singular_Kernels_Subs_20_epsA.mat', 'Newtonians');
disp('Contents of Singular_Kernels_Subs_20_epsA.mat:')
whos('-file', 'Singular_Kernels_Subs_20_epsA.mat')

clear('Newtonians')
%load('Singular_Kernels_Subs_20_epsA.mat', 'Newtonians')

% 20 B: 1e-5
epsilon = Range_Epsilons(2);
Newtonians = Newtonian([aBox.Pts.y1_kv, aBox.Pts.y2_kv], epsilon);
Newtonians.eps = epsilon;  % This is ε
Newtonians.box = aBox;  % This is the domain
%Newtonians.NI; % is the convolution outside [0,1]^2 ∩ N(ε)
%Newtonians.NJ; % is the convolution [0,1]^2
%Newtonians.NG; % is the convolution in N(ε), this is `smaller`
Newtonians.N = 20;      % Number of points per dimension

% Now we add a substruct
Newtonians.Level = {};
Newtonians.Level.Available = Factors;     % Number of subdivisions per dimension
% Store Convolution matrices
Newtonians.Level.n2  = Conv_Out_20{4};
Newtonians.Level.n5 = Conv_Out_20{5};
Newtonians.Level.n8 = Conv_Out_20{6};

save('Singular_Kernels_Subs_20_epsB.mat', 'Newtonians');
disp('Contents of Singular_Kernels_Subs_20_epsB.mat:')
whos('-file', 'Singular_Kernels_Subs_20_epsB.mat')

clear('Newtonians')
%load('Singular_Kernels_Subs_20_epsB.mat', 'Newtonians')



%% 10 x 10 boxes
geom.N = [10;10];
geom.y1Min = 0; geom.y1Max = 1; geom.y2Min = 0; geom.y2Max = 1;
aBox = Box(geom);

% 10 A: 1e-2
epsilon = Range_Epsilons(1);
Newtonians = Newtonian([aBox.Pts.y1_kv, aBox.Pts.y2_kv], epsilon);
Newtonians.eps = epsilon;  % This is ε
Newtonians.box = aBox;  % This is the domain
%Newtonians.NI; % is the convolution outside [0,1]^2 ∩ N(ε)
%Newtonians.NJ; % is the convolution [0,1]^2
%Newtonians.NG; % is the convolution in N(ε), this is `smaller`
Newtonians.N = geom.N;      % Number of points per dimension

% Now we add a substruct
Newtonians.Level = {};
Newtonians.Level.Available = Factors;     % Number of subdivisions per dimension
% Store Convolution matrices
Newtonians.Level.n2  = Conv_Out_10{1};
Newtonians.Level.n10 = Conv_Out_10{2};
Newtonians.Level.n16 = Conv_Out_10{3};

save('Singular_Kernels_Subs_10_epsA.mat', 'Newtonians');
disp('Contents of Singular_Kernels_Subs_10_epsA.mat:')
whos('-file', 'Singular_Kernels_Subs_10_epsA.mat')

clear('Newtonians')
%load('Singular_Kernels_Subs_10_epsA.mat', 'Newtonians')

% 10 B: 1e-5
epsilon = Range_Epsilons(2);
Newtonians = Newtonian([aBox.Pts.y1_kv, aBox.Pts.y2_kv], epsilon);
Newtonians.eps = epsilon;  % This is ε
Newtonians.box = aBox;  % This is the domain
%Newtonians.NI; % is the convolution outside [0,1]^2 ∩ N(ε)
%Newtonians.NJ; % is the convolution [0,1]^2
%Newtonians.NG; % is the convolution in N(ε), this is `smaller`
Newtonians.N = geom.N;      % Number of points per dimension

% Now we add a substruct
Newtonians.Level = {};
Newtonians.Level.Available = Factors;     % Number of subdivisions per dimension
% Store Convolution matrices
Newtonians.Level.n2  = Conv_Out_10{4};
Newtonians.Level.n10 = Conv_Out_10{5};
Newtonians.Level.n16 = Conv_Out_10{6};

save('Singular_Kernels_Subs_10_epsB.mat', 'Newtonians');
disp('Contents of Singular_Kernels_Subs_10_epsB.mat:')
whos('-file', 'Singular_Kernels_Subs_10_epsB.mat')

clear('Newtonians')
load('Singular_Kernels_Subs_10_epsB.mat', 'Newtonians')

%% 3D fun
% figure(1)
% [MX, MY] = meshgrid( linspace(0,1,100), linspace(0,1,100) );
% s = surf(MX, MY, rescale(-Newtonians.Level.n10, 0,1),'FaceAlpha',0.2 );
% s.EdgeColor = 'none';
% stlwrite('test.stl',MX,MY,rescale(-Newtonians.Level.n10, 0,1),'mode','ascii')


%% Plot errors for a fixed ε: 10

% 1, 4, 12, end
Plot_Errors_Singular_Conv_Domain(Range_Epsilons(4), Abs_10_100_a{4}, aBox.Pts, true); Range_Epsilons(4)


%% Compute angles without an ε

geom.N = [20;20];
geom.y1Min = 0;    geom.y1Max = 1;
geom.y2Min = 0;    geom.y2Max = 1;
% Create an abstract box, we will move this along the domain
aBox = Box(geom);

Pts = aBox.Pts;
Angles = zeros( Pts.N1 * Pts.N2, 8);
for i = 1:( Pts.N1 * Pts.N2 )
    if (0 < Pts.y1_kv(i) && Pts.y1_kv(i) < 1) && (0 < Pts.y2_kv(i) && Pts.y2_kv(i) < 1)
        Angles(i,1) = (1-Pts.y2_kv(i)) / Pts.y1_kv(i);
        Angles(i,2) =    Pts.y2_kv(i) / Pts.y1_kv(i);
        Angles(i,3) =    Pts.y2_kv(i) / (1-Pts.y1_kv(i));
        Angles(i,4) = (1-Pts.y2_kv(i)) / (1-Pts.y1_kv(i));

        Angles(i,5) = (1-Pts.y1_kv(i)) / Pts.y2_kv(i);
        Angles(i,6) =    Pts.y1_kv(i) / Pts.y2_kv(i);
        Angles(i,7) =    Pts.y1_kv(i) / (1-Pts.y2_kv(i));
        Angles(i,8) = (1-Pts.y1_kv(i)) / (1-Pts.y2_kv(i));

        %[i, Pts.y1_kv(i), Pts.y2_kv(i), rad2deg(atan(Angles(i,:))) ] 183
    end
end
%Worst_Angles = atan(max(Angles, [], 2));
Worst_Angles = max(atan(Angles), [], 2);
Worst_Angles( max(Angles, [], 2)==0) = 0;


Plot_Angles_Domain(Worst_Angles, aBox.Pts, true)













