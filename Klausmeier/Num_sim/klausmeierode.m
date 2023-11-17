function dcombineddt = klausmeierode(t,v,B,nu,d,M,dx,kt,A_vec,x)
% Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
% License: GNU GPL
% Last updated: 17/11/2023


dcombineddt = zeros(2*M,1);
A = interp1(kt,A_vec,t); % find rainfall parameter from input vector of changing values

%% add perturbation to solution at every timestep
% newL = 40;
% pert = sin(2*pi/newL*x);
% v(1:M) = v(1:M) + min(v(1:M))*pert';
% v(M+1:2*M) = v(M+1:2*M) + min(v(M+1:2*M))*pert';


cc=2:(M-1);

% u eq
dcombineddt(1)  = v(1).^2.*v(M+1)   - B*v(1)   +(v(2)-2*v(1)+v(M))/(dx^2);
dcombineddt(cc) = v(cc).^2.*v(M+cc) - B*v(cc)  +(v(cc+1)-2*v(cc)+v(cc-1))/(dx^2);
dcombineddt(M)  = v(M).^2.*v(2*M)   - B*v(M)   +(v(1)-2*v(M)+v(M-1))/(dx^2);


% w eq
dcombineddt(M+1)  = A -  v(M+1) -   v(1).^2.*v(M+1)   + nu*(v(M+2)-v(M+1))/dx      +d*(v(M+2)-2*v(M+1)+v(2*M))/(dx^2);
dcombineddt(M+cc) = A -  v(M+cc)-   v(cc).^2.*v(M+cc) + nu*(v(M+cc+1)-v(M+cc))/dx  +d*(v(M+cc+1)-2*v(M+cc)+v(M+cc-1))/(dx^2);
dcombineddt(2*M)  = A -  v(2*M) - v(M).^2.*v(2*M)     + nu*(v(M+1)-v(2*M))/dx      +d*(v(M+1)-2*v(2*M)+v(2*M-1))/(dx^2);

