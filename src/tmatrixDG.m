classdef tmatrixDG < tmatrix & matlab.mixin.Copyable

    properties
        radius
    end

    methods

        function obj = tmatrixDG(varargin)

            if nargin == 1 && isa(varargin{1}, 'tmatrix')
                base = varargin{1};

                % Extract values first
                order   = base.order;
                kwave   = base.kwave;
                matrix  = base.matrix;
                origin  = base.origin;
                comment = base.comments;

            else
                % Pass-through case
                order   = varargin{1};
                kwave   = varargin{2};
                matrix  = varargin{3};
                origin  = varargin{4};
                comment = varargin{5};
            end

            % call superclass constructor exactly once, unconditionally
            obj@tmatrix(order, kwave, matrix, origin, comment);

        end

        function saveRadius(self, radius)
            self.radius = radius;
        end

    end

end