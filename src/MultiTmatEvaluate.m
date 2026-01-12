function L2_norm = MultiTmatEvaluate(kwave,uinc,tmat,solver,ScatDisp,r)

% Solve multiple scattering problem with tmat, solver and disposition given
% and evaluate L2 norm on the ball of radius r

% get parameters
pos=ScatDisp.pos; rot=ScatDisp.rot; shape=ScatDisp.shape;
NScat=size(shape,1); Ntype=max(shape);

% rotate T-matrices and solvers
v=0;
for j=1:NScat
    if rot(j) ~= 0
        v=v+1;
        [rotTmat,rotSolver] = RotateTmat(tmat{shape(j)},solver{shape(j)},rot(j));
        shape(j) = Ntype+v;
        tmat{Ntype+v} = rotTmat;
        solver{Ntype+v} = rotSolver;
    end
end

% get expansion order
nmax=tmat{1}.order;


% Solve iteratively multiple particle scattering problem
disp('Solving multiple scattering problem and computing norm')
tic()
% create wave function expansions of plane wave at the centers of
% each of the scatterers
parfor j=1:NScat
    a{j} = regularwavefunctionexpansion(nmax,pos(j),uinc);
end

% apply the T-matrices to the incident field
parfor j=1:NScat

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
nitns = min(100,floor(NScat*(2*nmax+1)/2));

% Solve the linear system using GMRES
[x,~,~,~,~] = gmres(@matrix_product,rhs,nitns,1e-8,1);

% convert coefficients into wavefunction expansions
c = unpack(x);

% Evaluate total field and compute norm
% setup a grid
NEval = 200;
t=linspace(-r,r,NEval); s=linspace(-r,r,NEval);
[x,y]=meshgrid(s,t);
z = x+y*1i;

% get the largest radius of the scatterers
maxrad = tmat{1}.radius;
for j=2:size(tmat,2)
    maxrad=max(maxrad,tmat{j}.radius);
end

mask = abs(z) < r; 

% get the scattered field... this is just the sum of the radiating fields
% from each scatterer
scatfield = c{1}.evaluate(z,mask);
for j=2:NScat
    scatfield = scatfield + c{j}.evaluate(z,mask);
end

% get the total field
totalfield = scatfield + uinc.evaluate(z,mask);
totalfield(isnan(totalfield))=0;

% compute norm of the total field on the circle of radius r
f_squared = abs(totalfield).^2;
integral_val = trapz(s, trapz(t, f_squared, 2), 1);
L2_norm = sqrt(integral_val);

toc()

%-----------------------------------------
% indented function to implement the matrix
% product in GMRES
%-----------------------------------------

    function y = matrix_product(x)
        
        % convert vector of coefficients into wavefunction expansions
        c = unpack(x);
        
        % apply matrix
        for j=1:NScat
            
            % we temporarily set the origin for the T-matrix object to the position
            % of the jth scatterer... this allows the T-matrix to interact with
            % wavefunction expansions with the same origin
            tmat{shape(j)}.setOrigin(pos(j));
            
            % initialize sum
            csum = regularzero(nmax,pos(j),kwave);

            % sum contributions from the other scatterers
            for i=1:NScat
                
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
        vec = zeros(2*nmax+1,NScat);

        % copy the coefficients into the array
        for j=1:NScat
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
        vec = reshape(vec,2*nmax+1,NScat);

        % create radiating wavefunction expansions from the columns of the
        % array
        for j=1:NScat
            a{j} = radiatingwavefunctionexpansion(nmax,pos(j),kwave,vec(:,j));
        end
        
    end

end