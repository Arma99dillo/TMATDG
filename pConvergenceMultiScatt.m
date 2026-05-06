%% Set main parameters
clear; close all; addpath src; addpath TMATROM_OBJECT_ORIENTED_CORE

k = 10; % wavenumber
h = 0.5; % mesh size
p_ref = 25; % number of plane wave directions for the refined solution

% define the shape vertices as a Nx2 matrix
NShape = 3; ScatShape=cell(NShape,1);
ScatShape{1}.vertices = [1/3, 1/3; 1/3, 1; -1/3, 1; -1/3, 1/3; -1, 1/3;
    -1, -1/3; -1/3, -1/3; -1/3, -1; 1/3, -1; 1/3, -1/3; 1, -1/3; 1, 1/3]; % cross
ScatShape{1}.type = 'dir';
ScatShape{2}.vertices = [0, 1; -sqrt(3)/2, -1/2; sqrt(3)/2, -1/2]; % triangle
ScatShape{2}.type = 'trans'; ScatShape{2}.n_in = 2.5;
ScatShape{3}.vertices = [1, 0; 0, 1; -1, 0; 0,-1]; % square
ScatShape{3}.type = 'dir';

%define the scatterer type, position and rotation angle
ScatArr.shape = [1; 1; 2; 3; 2]; % scatterer shape
ScatArr.pos = [-4-4i; 4-3.5i; 0; -3+4i; 3.5+3i]; % scatterer center position
ScatArr.rot = [-pi/4; 0; 0; 0; pi]; % rotation angle

%% Compute refined T-matrices and solvers
[tmat_ref,solver_ref] = MultiTmat(k,h,p_ref,ScatShape);

%% p-convergence
theta = 3*pi/4; % incident angle
uinc = plane_wave(theta,k); % incident plane wave function
PlotPar.limX=[-7,7]; PlotPar.limY=[-7,7];
ni=p_ref-1; nf=p_ref-1; %first and last value of number of directions
L2Error=zeros(nf-ni+1,1); TmatErr=zeros(nf-ni+2,NShape); v=1;
for p=ni:nf

    disp(['Computing with p=',num2str(p), ' plane wave directions'])

    % Compute T-matrices and solvers for the corresponding value of p
    [tmat,solver] = MultiTmat(k,h,p,ScatShape);
    
    % Compute error of the T-matrices
    for l=1:NShape
        TmatErr(v,l)=max(max(abs(tmat{l}.matrix + tmat{l}.matrix' ...
                + 2 * tmat{l}.matrix' * tmat{l}.matrix)))/max(max(abs(tmat{l}.matrix)));
    end

    % Solve multiple scattering problems and compute error
    L2Error(v) = MultiTmatError(k,uinc,tmat_ref,solver_ref,tmat,solver,ScatArr,PlotPar);

    v=v+1;
end
for l=1:NShape
    TmatErr(v,l)=tmat_ref{l}.error();
end


%% Plot TDG error convergence
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


figure()
semilogy(ni:nf,L2Error,'*-r','LineWidth',1.2); grid
xlim([ni nf]);
LL = legend('$L^2$ error','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Number of local plane wave basis functions','FontSize',18, 'Interpreter','latex')
ylabel('Error','FontSize',18, 'Interpreter','latex')

%% Plot T-matrix error convergence
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


figure()
semilogy(ni:p_ref,TmatErr(:,1),'x-',ni:p_ref,TmatErr(:,2),'^-',...
    ni:p_ref,TmatErr(:,3),'s-','LineWidth',1.2); grid
xlim([ni p_ref]);
LL = legend('Cross','Triangle','Square','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Number of local plane wave basis functions','FontSize',18, 'Interpreter','latex')
ylabel('T-matrix error','FontSize',18, 'Interpreter','latex')