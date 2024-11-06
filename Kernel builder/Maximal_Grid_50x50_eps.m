%% Approximate convolution for different values of ε and factors (grids in MS) for maximal sensible subdivision and large N
% Warning: If you run the following code from scratch, it will take a
% couple of hours to finish. Average and total computing times are reported
% for each experiment.

% Extremely recomend to store working space: 'autosave(10,'workspace.mat')'
% Command to be stopped with 'autosave stop' followed by 'autosave delete'

%% 50 x 50 base grid
% Factor 1     takes     19 min 40.66 s    /        21 min 10.44 s
% Factor 2     takes 1 h 14 min 20.90 s    /    1 h 19 min 53.69 s
% Factor 3.2   takes 4 h 24 min  1.62 s    /    4 h 53 min 23.34 s

% The following takes around 5 h 50 min
k = 1;
T = zeros(1,3);
Factors = [1,2,3.2];
Range_Epsilons = [1e-2];
for eps = Range_Epsilons 
    eps
    for i = Factors
        i
        tic;
        [aBox,MS, errors, t_Full, Conv, experiment, smaller, full_thing, all_sq] = SingularKernelIntergralMS_Conv_Maximal_Q(eps, 50, 50, i);
        T(k)= toc;  % pair 1: toc
        
        Abs_40_a{k} = errors;
        Rel_40_a{k} = errors ./ max( abs(t_Full), 1e-6);
    
        Abs_40_b{k} = experiment;
        Rel_40_b{k} = experiment ./ max( abs(smaller), 1e-6);
        Exact_40_b{k} = smaller;
    
        Abs_40_c{k} = full_thing;
        Rel_40_c{k} = full_thing ./ max( abs(all_sq), 1e-6);
        Exact_Full_40{k} = all_sq;
        
        Conv_Out_40{k} = Conv;
        k = k + 1;
    end
end
[sum(T), mean(T)]


%% Plot
Nx = 50;
Ny = 50;
de = 3;
colores = [255, 200, 87; 186, 45, 11; 86, 22, 67]/255;
%[186, 45, 11; 0, 105, 137; 255, 200, 87]/255; %; 46, 196, 182; 86, 22, 67]/255;
N = Nx * Ny;

% Define a colour gradient
colour_A = [255, 225, 168]/255;
colour_B = [114, 61, 70]/255;
GRADIENT_flexible = @(n,nn) interp1([1/nn 1],[colour_A;colour_B],n/nn);
G_colours = GRADIENT_flexible((1:de),de);

%% Swarms of each error set 
%%%%%%%%%%%
figure(2)
a(1) = plot(nan, nan, 'LineWidth', 2, 'LineStyle', '-', 'Color', 'black', 'DisplayName', 'Median'); hold on
a(2) = plot(nan, nan, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'black', 'DisplayName', 'Mode');
a(3) = plot(nan, nan, 'LineWidth', 1, 'LineStyle', '--', 'Color', 'black', 'DisplayName', 'Mean');
% 10 points
swarmchart(1.0*ones(N,1), Abs_40_a{1}, 10, colores(2,:), 'XJitterWidth', 0.3)
line([1-0.15,1+0.15],[median(Abs_40_a{1}),median(Abs_40_a{1})], 'LineWidth', 2, 'LineStyle', '-', 'Color', 'black') % colores( mod(1+4,5) + 1,:))
line([1-0.15,1+0.15],[mode(Abs_40_a{1}),mode(Abs_40_a{1})], 'LineWidth', 2, 'LineStyle', ':', 'Color', 'black') % colores( mod(1+4,5) + 1,:))
line([1-0.15,1+0.15],[mean(Abs_40_a{1}),mean(Abs_40_a{1})], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'black') % colores( mod(1+4,5) + 1,:))
yscale log
%hold on
for i = 2:de
    loc = (i+1)/2;
    swarmchart(loc*ones(N,1), Abs_40_a{i}, 10, colores( mod(i,3) + 1,:), 'XJitterWidth', 0.3)
    line([loc-0.15,loc+0.15],[median(Abs_40_a{i}),median(Abs_40_a{i})], ...
        'LineWidth', 2, 'LineStyle', '-', 'Color', 'black') %colores( mod(i+4,5) + 1,:))
    line([loc-0.15,loc+0.15],[mode(Abs_40_a{i}),mode(Abs_40_a{i})], ...
        'LineWidth', 2, 'LineStyle', ':', 'Color', 'black') %colores( mod(i+4,5) + 1,:))
    line([loc-0.15,loc+0.15],[mean(Abs_40_a{i}),mean(Abs_40_a{i})], ...
        'LineWidth', 1, 'LineStyle', '--', 'Color', 'black') %colores( mod(i+4,5) + 1,:))
end
xlim([0.75,(de+1)/2 + 0.25]);

MinMax = [1e+10,0];
for i = 1:de
    [mina, maxa] = bounds(Abs_40_a{i}(Abs_40_a{i}>0),'all');
    if mina < MinMax(1)
        MinMax(1) = mina;
    end
    if maxa > MinMax(2)
        MinMax(2) = maxa;
    end
end
%ylim([3e-18, 2e-5]);
ylim(MinMax .* [0.5, 2.5] );

% Additional decoration
lgd = legend(a, 'Interpreter','latex','Location','southwest');
title(lgd,'Descriptors','Interpreter','latex')
% 
% 
% 
xticks([1, 1.5, 2.0])
% 
divisions = {'50/1','50/2.5','50/4'};
xlabelArray = [ divisions; {' ',' 10^{-2}',' '} ]; 
xtickLabels = strtrim(sprintf('%s\\newline%s\n', xlabelArray{:}));
xticklabels(xtickLabels)

set(gca, 'FontName', 'CMR10')
xlabel('Grid resolution and detail factor','Interpreter','latex')
ylabel('Absolute errors','Interpreter','latex')
set(gca, 'YTick', 10.^(-20:1:-1)) %set(gca, 'YTick', logspace(-18,-1, 18))

file_name = strcat('Test_Singular_Convolution_Fixed_Swarm[', num2str(Nx), ',', num2str(Factors,'%.1f,'), num2str(de), '].pdf');
exportgraphics(figure(2), file_name, 'BackgroundColor','none','ContentType','vector')


%% Store struct with outputs
geom.N = [50;50];
geom.y1Min = 0; geom.y1Max = 1; geom.y2Min = 0; geom.y2Max = 1;
aBox = Box(geom);

epsilon = 1e-2;
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
Newtonians.Level.n1   = Conv_Out_40{1};
Newtonians.Level.n2 = Conv_Out_40{2};
Newtonians.Level.n3_2   = Conv_Out_40{3};

save('Singular_Kernels_Subs_50_epsA.mat', 'Newtonians');
disp('Contents of Singular_Kernels_Subs_10_epsB.mat:')
whos('-file', 'Singular_Kernels_Subs_10_epsB.mat')

clear('Newtonians')
load('Singular_Kernels_Subs_50_epsA.mat', 'Newtonians')

%% 3D fun
% figure(1)
% [MX, MY] = meshgrid( linspace(0,1,40^2), linspace(0,1,40^2) );
% s = surf(MX, MY, rescale(-Newtonians.Level.n4, 0,1),'FaceAlpha',0.2 );
% s.EdgeColor = 'none';
% stlwrite('test.stl',MX,MY,rescale(-Newtonians.Level.n4, 0,1),'mode','ascii')












