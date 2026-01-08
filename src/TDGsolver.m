% Define the TDG solver class, used to solve scattering problem 
% with the DtN-TDG method

classdef TDGsolver < solver
    
    %class properties
    properties
        matrix
        rhs
        matrixPlt
        rhsPlt
        coeffs
        mesh
        param
        scatterer
        uinc
        type
    end
    
    methods
        

        %-----------------------------------------
        % constructor: build the solver
        %-----------------------------------------
        
        function self = TDGsolver(kwave,incidentField,scatterer,param)
 
            % call parent constructor
            self = self@solver(kwave,incidentField)
            
            % inizialize parameters
            self.param = param; self.param.K = kwave;
            self.scatterer=scatterer;
            self.matrix = [];
            self.rhs = [];
            self.coeffs = [];
        end

        %-----------------------------------------
        % setup: build system matrix
        %-----------------------------------------
        
        % This method sets up the solver (system matrix)
        function setup(self)                      
            
            % build mesh
            self.mesh = GenerateMesh(self);
            
            % compute the TDG matrix in the two cases
            if strcmp(self.type,'dir')
                self.matrix = MatrixDtNTDG(self);
            else
                self.matrix = MatrixDtNTDGTrans(self);
            end
            
        end
        
        %--------------------------------------------
        % solve: build rhs and solve linear system
        %--------------------------------------------
        
        % This method solves the scattering problem for every right hand
        % side specified in the self.incidentField cell array.
        function solve(self)

            % check that the setup method has been run
            if isempty(self.matrix)
                error('Must call setup() first') 
            end
            
            % loop through incident directions
            for j=1:length(self.incidentField)
                
                %build rhs in both cases
                if strcmp(self.type,'dir')
                    gd = @(x) -self.incidentField{j}.evaluate(x); % incident wave
                    self.rhs(:,j) = RhsDtNTDG(self,gd);
                else
                    gd = @(x) self.incidentField{j}.evaluate(x); % incident wave
                    gn = @(x) self.incidentField{j}.evaluateGradientVect(x); % incident wave gradient
                    self.rhs(:,j) = RhsDtNTDGTrans(self,gd,gn);
                end
            end

            % solve lynear system
            self.coeffs = self.matrix\self.rhs;
        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------        
        function val = getFarField(self,points,index) 
            
            if nargin<3
                index = 1;
            end
            
            if isempty(self.coeffs)      
                error('Must run solve() first')  
            end
            
            for j=1:length(index)

                % get the far field
                val(:,j) = FarField(self,points,index(j));
            
            end

        end
        
        %===============================================================
        % visualization of the solution near the scatterer
        %===============================================================

        function Visualize(self,uinc,b,center,PlotType)
            self.uinc=uinc;

            % Use the solution obtained with T-matrix to evaluate the 
            % solution inside the radius of the scatterer
            % TDG method with Dirichlet condition on the circular boundary

            disp('Visualing the field near the scatterer')
            tic()
            % distinguish between Dirichlet and transmission case
            if strcmp(self.type,'dir')
                psave=self.mesh.p; u_save = self.uinc; % save original mesh
                % translate the mesh to its position
                self.mesh.p(:,1) = self.mesh.p(:,1)+real(center);
                self.mesh.p(:,2) = self.mesh.p(:,2)+imag(center);
                self.matrixPlt = MatrixTDGDir(self); % compute matrix for plot
                gd = @(x) b.evaluate(x) + self.uinc.evaluate(x); % Dirichlet datum
                self.rhsPlt = RhsTDGDirSingle(self,gd); % compute rhs for plot
                u = self.matrixPlt\self.rhsPlt; % solve linear system
                % Plot either total or scattered field
                if strcmp(PlotType,'tot')
                    SolPlotTot(self,u); % plot total field
                else
                    self.uinc = -self.uinc; % subtract incident wave
                    SolPlotScatt(self,u); % plot scatteres field
                end  
                self.mesh.p=psave; self.uinc=u_save; % restore original mesh
            else % transmission case
                psave=self.mesh.p; u_save = self.uinc; % save original mesh
                % translate the mesh to its position
                self.mesh.p(:,1) = self.mesh.p(:,1)+real(center);
                self.mesh.p(:,2) = self.mesh.p(:,2)+imag(center);
                self.matrixPlt = MatrixTDGDirTrans(self); % compute matrix for plot
                gd = @(x) b.evaluate(x) + self.uinc.evaluate(x); % Dirichlet datum
                self.rhsPlt = RhsTDGDirSingle(self,gd); % compute rhs for plot
                u = self.matrixPlt\self.rhsPlt; % solve linear system
                % Plot either total or scattered field
                if strcmp(PlotType,'tot')
                    SolPlotTotTrans(self,u); % plot total field
                else
                    self.uinc = -self.uinc; % subtract incident wave
                    SolPlotScattTrans(self,u); % plot scattered field
                end 
                self.mesh.p=psave; self.uinc=u_save; % restore original mesh
            end
            toc()
        end

        function VisualizeMult(self,uinc,c,center,numscat)
            self.uinc=uinc;

            % Use the solution obtained with T-matrix to evaluate the 
            % solution inside the radius of the scatterer
            % TDG method with Dirichlet condition on the circular boundary
            % multiple scattering case

            % distinguish between Dirichlet and transmission case
            if strcmp(self.type,'dir')
                psave=self.mesh.p; % save original mesh
                % translate the mesh to its position
                self.mesh.p(:,1) = self.mesh.p(:,1)+real(center);
                self.mesh.p(:,2) = self.mesh.p(:,2)+imag(center);
                self.matrixPlt = MatrixTDGDir(self); % compute matrix for plot
                self.rhsPlt = RhsTDGDirMulti(self,c,numscat); % compute rhs for plot
                u = self.matrixPlt\self.rhsPlt; % solve linear system
                SolPlotTot(self,u); % plot total field
                self.mesh.p=psave; % restore original mesh
            else % transmission case
                psave=self.mesh.p; % save original mesh
                % translate the mesh to its position
                self.mesh.p(:,1) = self.mesh.p(:,1)+real(center);
                self.mesh.p(:,2) = self.mesh.p(:,2)+imag(center);
                self.matrixPlt = MatrixTDGDirTrans(self); % compute matrix for plot
                self.rhsPlt = RhsTDGDirMulti(self,c,numscat); % compute rhs for plot
                u = self.matrixPlt\self.rhsPlt; % solve linear system
                SolPlotTotTrans(self,u); % plot total field
                self.mesh.p=psave; % restore original mesh
            end
        end

    end % end methods

end