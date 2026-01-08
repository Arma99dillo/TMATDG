function PlotSolution(tmat,solver,uinc,PlotPar)

% Solve and plot T-matrix solution

% get parameters
nmax = tmat.order; radius = tmat.radius;
center = tmat.origin; vertices = solver.param.vertices;
type_plot = PlotPar.type_plot; inside = PlotPar.inside;
limX = PlotPar.limX; limY = PlotPar.limY;

disp('Computing wave function expansion of the scattered field')
tic()
a = regularwavefunctionexpansion(nmax,center,uinc); % create wave function expansion of incident plane wave
b = tmat * a; % compute wave function expansion of scattered wave using T-matrix
toc()

% Distinguish between total and scattered field
if strcmp(type_plot,'tot')
    disp('Evaluating and visualizing total field outside of the scatterer')
    tic()
    
    % setup a grid
    figure()
    t1=linspace(limX(1),limX(2),500); t2=linspace(limY(1),limY(2),500);
    [x,y]=meshgrid(t1,t2);
    z = x+y*1i;
    mask = abs(z-center) < 1.2*radius; % get a mask for the scatterer
    surf(x,y,real(b.evaluate(z,~mask)+uinc.evaluate(z,~mask)));
    view(2);
    shading interp;
    grid off;
    axis square
    colorbar
    toc()
elseif strcmp(type_plot,'scat')
    disp('Evaluating and visualizing scattered field outside of the scatterer')
    tic()

    figure()
    % setup a grid
    t1=linspace(limX(1),limX(2),500); t2=linspace(limY(1),limY(2),500);
    [x,y]=meshgrid(t1,t2);
    z = x+y*1i;
    mask = abs(z-center) < 1.2*radius; % get a mask for the scatterer
    surf(x,y,real(b.evaluate(z,~mask)));
    view(2);
    shading interp;
    grid off;
    axis square
    colorbar
    toc()
else
    error('Plot type invalid!');
end
hold on

% plot scattere profile
if strcmp(solver.type, 'dir')
    v=vertices;
    n=size(v,1); v(n+1,:)=v(1,:); v=v+[real(center),imag(center)];
    plot3(v(:,1),v(:,2),10*ones(size(v)),'k-',LineWidth=1.5)
    hold on
else
    v=vertices;
    n=size(v,1); v(n+1,:)=v(1,:); v=v+[real(center),imag(center)];
    if inside
        plot3(v(:,1),v(:,2),10*ones(size(v)),'w:',LineWidth=1.5)
        hold on
    else
        plot3(v(:,1),v(:,2),10*ones(size(v)),'k:',LineWidth=1.5)
        hold on
    end
end

% plot near the scatterer
if inside 
    solver.Visualize(uinc,b,center,type_plot)
end

axis equal; xlim(limX); ylim(limY);

end