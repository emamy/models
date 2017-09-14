function [KK]= StiffnessMatrix_Vmnh(dim,time_step,dx,D)
% stiffness matrix (V m n h V m n h)
  n=dim;
  KK = zeros(dim,dim);
  
  % introduce factor for compact notation      
  fact = time_step/dx^2 * D;      
      
  for i=2:n-7
      % Wenn Index = V
      if mod(i, 4) == 1
          KK(i,i-4) = fact;
          KK(i,i)   = -2*fact;
          KK(i,i+4) = fact;    
      end
  end
  
%   for j = 1:n
%       if mod(j,4) ~= 1
%           KK(j,j)  = 1;
%       end
%   end
  
  % homogeneous Neumann B.C. - no flow
  % right side
    %KK(n-3,n-3) = -1;
    KK(n-3,n-7) = 1;
  % left side
    %KK(1,1)     = -1;
    KK(1,5)     = 1;  
end
