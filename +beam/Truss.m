classdef Truss < models.beam.StructureElement
    % Truss:
    %
    %
    %
    % @author Daniel Wirtz @date 2011-12-09
    %
    % @new{0,6,dw,2011-12-09} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Constant)
        PlotSplitFactor = 5;
    end
    
    methods
        function this = Truss(model, material, pointsidx)
            this = this@models.beam.StructureElement(model, material, pointsidx);
            this.initialize;
        end
        
        function M = getLocalMassMatrix(this)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines Stabes
            % c kodiert Stoffparameter:
            % c1 = E*A/L      (Federh�rte)
            % c2 = rho*A*L/6
            
            % Transformationsmatrix: nat�rliche Koords -> glob. Koords
            % Sortierung der Variablen:
            % u0, v0, w0, u1, v1, w1
            %  1   2   3   4   5   6
            
            T_block = this.T;
            
            T = zeros(6);
            T([1 2 3], [1 2 3]) = T_block;    % Verschiebung links
            T([4 5 6], [4 5 6]) = T_block;    % Verschiebung rechts
            
            M = zeros(6);
            M([1 4], [1 4]) = this.c(2) * [2 1; 1 2];
            
            M = T * M * T';
        end
        
        function K = getLocalStiffnessMatrix(this)
            % Berechnet lokale Steifigkeitsmatrix eines Stabes
            % c kodiert Stoffparameter:
            % c1 = E*A/L      (Federh�rte)
            % c2 = rho*A*L/6
            
            % Transformationsmatrix: nat�rliche Koords -> glob. Koords
            % Sortierung der Variablen:
            % u0, v0, w0, u1, v1, w1
            %  1   2   3   4   5   6
            T_block = this.T;
            
            T = zeros(6);
            T([1 2 3], [1 2 3]) = T_block;      % Verschiebung links
            T([4 5 6], [4 5 6]) = T_block;    % Verschiebung rechts
            
            K = zeros(6);
            K([1 4], [1 4]) = this.c(1) * [1 -1; -1 1];
            
            K = T * K * T';
        end
        
        function f = getLocalForce(this, gravity)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines Stabes
            % c kodiert Stoffparameter:
            % c1 = E*A/L      (Federh�rte)
            % c2 = rho*A*L/6
            
            % Transformationsmatrix: nat�rliche Koords -> glob. Koords
            % Sortierung der Variablen:
            % u0, v0, w0, u1, v1, w1
            %  1   2   3   4   5   6
            
            T_block = this.T;
            T = zeros(6);
            T([1 2 3], [1 2 3]) = T_block;    % Verschiebung links
            T([4 5 6], [4 5 6]) = T_block;    % Verschiebung rechts
            
            q_lok = this.Material(1) * this.Material(2) * (this.T' * gravity);
            f = zeros(6, 1);
            f([1 4]) = 0.5 * q_lok(1) * this.Length * [1; 1];
            
            f = T * f;
        end
        
        function plot(this, p, u_elem, plot_options)            
            s = 0 : this.Length/(4*this.PlotSplitFactor) : this.Length;
            
            u_p = plot_options.multiplier * (u_elem(2,:)' * s/this.Length + u_elem(1,:)' * (1 - s/this.Length));
            % Parametrisierung des Viertelkreises in lokalen Koords
            y = (this.Length/(2*this.PlotSplitFactor)) * sin( 2*pi*s / (this.Length/this.PlotSplitFactor) );
            x = s;
            z = 0*s;
            % Umrechnung in globale Koords und Verschiebung um globale Koords des
            % lokalen Ursprungs KR(i).pc
            COR = this.T * [x; y; z];
            plot3( p(this.PointsIdx(1), 1) + COR(1,:) + u_p(1,:), ...
                p(this.PointsIdx(1), 2) + COR(2,:) + u_p(2,:), ...
                p(this.PointsIdx(1), 3) + COR(3,:) + u_p(3,:), 'k' );
        end
    end
    
    methods(Access=private)
        function initialize(this)
            m = this.Model;
            % L�ngenberechnung
            dp = (m.Points(this.PointsIdx(2), :) - m.Points(this.PointsIdx(1), :));
            l = norm(dp);
            this.Length = l;
            
            % Lokales Koordinatensystem
            e_x = dp' / l;
            e_y = [-e_x(2) e_x(1) 0]';
            if norm(e_y) == 0
                e_y = [0 e_x(3) -e_x(2)]';
            end
            e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
                e_x(3)*e_y(1) - e_x(1)*e_y(3);
                e_x(1)*e_y(2) - e_x(2)*e_y(1)];
            e_y = e_y / norm(e_y);
            e_z = e_z / norm(e_z);
            this.T = [e_x e_y e_z];
            
            %     Effektive Konstanten
            %     <rho>	<A>     <E>
            %       1    2       3
            material = this.Material;
            this.c(1) = material(3) * material(2) / l;                           % c1 = E*A/L      (Federh�rte)
            this.c(2) = material(1) * material(2) * l / 6;                       % c2 = rho*A*L/6
            
            this.c_theta = material(9);
            this.kappa = material(10);
            this.alphaA = material(11);
        end
    end
    
end