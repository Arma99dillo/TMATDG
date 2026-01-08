function MultiScatt(kwave,h,nd,uinc,ScatShape,ScatDisp,PlotPar,SavePath)

% Computes T-matrices for every shape in ScatShape and solves the multiple
% scattering problem with disposion in ScatDisp
% Then plots the solution in the limits limX and limY

% get parameters
pos=ScatDisp.pos; rot=ScatDisp.rot; shape=ScatDisp.shape;
Nscat=size(shape,1); Ntype=size(ScatShape,1);
InsidePlot = PlotPar.inside; limX = PlotPar.limX; limY = PlotPar.limY;
SaveFile = SavePath.choice; FileName = SavePath.file;

% compute expansion orders for every shape and save the maximum
orders = zeros(Ntype,1);
for j=1:Ntype
    points = ScatShape{j}.vertices;
    points = points - mean(points, 1); % center in zero to easy compute distance
    distances = sqrt(sum(points.^2, 2)); % distance from the origin
    radius = max(distances);
    orders(j) = suggestedorder(kwave,radius);
end
nmax = max(orders);

% for each type of shape setup solver and compute T-matrix
parfor j=1:Ntype
    [tmat{j},solver{j}]=ComputeTMatrix(kwave,h,nd,ScatShape{j},nmax);
end

% save matrices and solvers in the given file
if SaveFile
    save(FileName, 'solver', 'tmat');
    disp(['Saved solvers and T-matrices in ', FileName, ' file.'])
end

% rotate matrices and solvers if necessary
v=0;
for j=1:Nscat
    if rot(j) ~= 0
        v=v+1;
        [rotTmat,rotSolver] = RotateTmat(tmat{shape(j)},solver{shape(j)},rot(j));
        shape(j) = Ntype+v;
        tmat{Ntype+v} = rotTmat;
        solver{Ntype+v} = rotSolver;
    end
end
Ntype = Ntype+v;

% Solve iteratively multiple particle scattering problem
disp(['Solving multiple scattering problem with ', num2str(Nscat), ' scatterers'])
tic()
% create wave function expansions of plane wave at the centers of
% each of the scatterers
parfor j=1:Nscat
    a{j} = regularwavefunctionexpansion(nmax,pos(j),uinc);
end

% apply the T-matrices to the incident field
parfor j=1:Nscat

    % we temporarily set the origin for the T-matrix object to the position
    % of the jth scatterer... this allows the T-matrix to interact with
    % wavefunction expansions with the same origin
    tmat{shape(j)}.setOrigin(pos(j));

    % apply the T-matrix to the incident field
    b{j} = tmat{shape(j)} * a{j};

end

% Note: GMRES works with vectors so we need to convert our cell array of
% wavefunction expansions into a vector...
% setup an array right hand side coefficients (to pass to GMRES)
rhs = pack(b);

% set number of GMRES iterations
nitns = min(100,floor(Nscat*(2*nmax+1)/2));

% Solve the linear system using GMRES
[x,~,~,~,~] = gmres(@matrix_product,rhs,nitns,1e-8,1);

% convert coefficients into wavefunction expansions
c = unpack(x);

toc()


% Visualize the total field
disp('Evaluating and visualizing total field outside of the scatterers')
tic()

% setup a grid
t=linspace(limY(1),limY(2),500); s=linspace(limX(1),limX(2),1000);
[x,y]=meshgrid(s,t);
z = x+y*1i;

% get the largest radius of the scatterers
maxrad = tmat{1}.radius;
for j=2:Ntype
    maxrad=max(maxrad,tmat{j}.radius);
end

% get a mask for the scatterers
mask = abs(z-pos(1)) < 1.2 * maxrad;
for j=1:Nscat
    mask = mask | abs(z-pos(j)) < 1.2 * maxrad;
end

% get the scattered field... this is just the sum of the radiating fields
% from each scatterer
scatfield = c{1}.evaluate(z,~mask);
for j=2:Nscat
    scatfield = scatfield + c{j}.evaluate(z,~mask);
end

figure()

% get the total field
totalfield = scatfield + uinc.evaluate(z,~mask);

% plot the total field
surf(x,y,real(totalfield))
view([0 90]);
shading interp;
colorbar; grid off
hold on
% plot the scatterer profile
for j=1:Nscat
    if strcmp(solver{shape(j)}.scatterer,'dir')
        v=solver{shape(j)}.param.vertices;
        n=size(v,1); v(n+1,:)=v(1,:); v=v+[real(pos(j)),imag(pos(j))];
        plot3(v(:,1),v(:,2),10*ones(size(v)),'k-',LineWidth=1.5)
        hold on
    else
        v=solver{shape(j)}.param.vertices;
        n=size(v,1); v(n+1,:)=v(1,:); v=v+[real(pos(j)),imag(pos(j))];
        if InsidePlot
            plot3(v(:,1),v(:,2),10*ones(size(v)),'w:',LineWidth=1.5)
            hold on
        else
            plot3(v(:,1),v(:,2),10*ones(size(v)),'k:',LineWidth=1.5)
            hold on
        end
    end
end
axis equal; xlim(limX); ylim(limY);
toc()

% plot near the scatterers 
if InsidePlot
    for j=1:Nscat
        disp(['Plotting near scatterer ', num2str(j)])
        tic()
        solver{shape(j)}.VisualizeMult(uinc,c,pos(j),Nscat);
        toc()
    end
end

%-----------------------------------------
% indented function to implement the matrix
% product in GMRES
%-----------------------------------------

    function y = matrix_product(x)
        
        % convert vector of coefficients into wavefunction expansions
        c = unpack(x);
        
        % apply matrix
        for j=1:Nscat
            
            % we temporarily set the origin for the T-matrix object to the position
            % of the jth scatterer... this allows the T-matrix to interact with
            % wavefunction expansions with the same origin
            tmat{shape(j)}.setOrigin(pos(j));
            
            % initialize sum
            csum = regularzero(nmax,pos(j),kwave);

            % sum contributions from the other scatterers
            for i=1:Nscat
                
                if i~=j
                    
                    % get the expansion of c{i} at pos{j}
                    csum = csum + regularwavefunctionexpansion(c{i},pos(j));
                    
                end
                
            end
            
            % apply the T-matrix to the sum
            d{j} = c{j} - tmat{shape(j)} * csum;
            
        end
        
        % convert coefficients into a vector
        y = pack(d);
        
    end

%-----------------------------------------
% indented function to pack wavefunction
% expansion coefficients into a vector
%-----------------------------------------

    function vec = pack(a)
        
        % create an array to hold the coefficients
        vec = zeros(2*nmax+1,Nscat);

        % copy the coefficients into the array
        for j=1:Nscat
            vec(:,j) = a{j}.getCoefficients();
        end
        
        % reshape the array into a vector
        vec = reshape(vec,[],1);
        
    end

%-----------------------------------------
% indented function to extract wavefunction
% expansion coefficients from a vector
%-----------------------------------------

    function a = unpack(vec)
        
        % reshape the vector into an array
        vec = reshape(vec,2*nmax+1,Nscat);

        % create radiating wavefunction expansions from the columns of the
        % array
        for j=1:Nscat
            a{j} = radiatingwavefunctionexpansion(nmax,pos(j),kwave,vec(:,j));
        end
        
    end

end