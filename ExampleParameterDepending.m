%% Set main parameters
clear; close all; addpath src; addpath TMATROM_OBJECT_ORIENTED_CORE

k = 2.39; % wavenumber
h = 0.02; % mesh size
p = 20; % number of plane wave directions

nv = 6; % number of vertices
for j=1:nv
    scatt.vertices(j,:) = 0.05.*[cos(j*2*pi/nv), sin(j*2*pi/nv)];
end
scatt.type='dir'; 

%% Compute T-matrix
[tmat,solver]=ComputeTMatrix(k,h,p,scatt);

%% Cicle on radius length
uinc = point_source(2+0i,k); %circular wave

Ri = 0.8; Rf = 1.4; step = 0.01; nexp = ceil((Rf-Ri)/step)+1;
l=1; r=0.5; L2_norm = zeros(nexp,1); RR = zeros(nexp,1);
for R=Ri:step:Rf

    % Define scatterers orientation and position
    NScat=30; % number of scatterers
    pos=zeros(NScat,1); rot=zeros(NScat,1);
    for j=1:NScat
        pos(j)=R*(cos(j*2*pi/NScat)+1i*sin(j*2*pi/NScat));
        rot(j)=pi/2+j*2*pi/NScat;
    end
    
    % define scatterer shapes
    ScatShape=ones(NScat,1); % they all have the same shape
    
    ScatArr.pos=pos; ScatArr.shape=ScatShape; ScatArr.rot=rot;
    
    % Solve and compute norm
    TmatCell = cell(1,1); TmatCell{1}=tmat;
    SolverCell = cell(1,1); SolverCell{1}=solver;

    disp(['Experiment ',num2str(l), ' out of ', num2str(nexp)])
    L2_norm(l) = MultiTmatEvaluate(k,uinc,TmatCell,SolverCell,ScatArr,r);
    RR(l) = R;
    l=l+1;
end

%% Plot norm
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


figure()
plot(RR,L2_norm,'.-r','LineWidth',1.5); grid
xlim([0.8 1.4]); 
LL = legend('$L^2$ norm','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Radius','FontSize',18, 'Interpreter','latex')
ylabel('Norm','FontSize',18, 'Interpreter','latex')


%% Plot solution
k = 2.39; % wavenumber
R = 1.05;

% Define scatterers orientation and position
NScat=30; % number of scatterers
pos=zeros(NScat,1); rot=zeros(NScat,1);
for j=1:NScat
    pos(j)=R*(cos(j*2*pi/NScat)+1i*sin(j*2*pi/NScat));
    rot(j)=pi/2+j*2*pi/NScat;
end

% define scatterer shapes
ScatShape=ones(NScat,1); % they all have the same shape

ScatArr.pos=pos; ScatArr.shape=ScatShape; ScatArr.rot=rot;

uinc = point_source(2+0i,k); %circular wave

% Solve and plot
PlotPar.inside = false; PlotPar.limX=[-1.5,2.5]; PlotPar.limY=[-2,2];
MultiTmatSolve(k,uinc,TmatCell,SolverCell,ScatArr,PlotPar);