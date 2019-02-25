function [y] = RingTotalVariation(x,lam,alp,Center,varargin)

% This code is written based on the code provided in the papers:
% (1)Zhang, H., & Wang, Y. (2013). Edge adaptive directional total variation.
%    The Journal of Engineering, 2013(11), 61-62.

% (2)Bayram, Ilker, and Mustafa E. Kamasak. "A directional total variation."
%   2012 Proceedings of the 20th European Signal Processing Conference (EUSIPCO). IEEE, 2012.

% minimizes wrt z, 0.5*||x - z||_2^2 + lam*RTV_{alpha,theta}(z) (see eq7 in the manuscript)
% Input parameters
% x : noisy image
% lam : lambda parameter in the functional
% alp : length of the major axis of the ellipse
% theta : direction of the ellipse
% varargin : number of iterations (default = 100)
% Output : y : minimizer of the functional above
% Center : The center of the rings

% Modified by: Morteza Salehjahromi
%% Added by me (also x:times converted to .x in lines containing the R)
Lx=size(x,1);
Ly=size(x,2);

Row_c  = Center(1);
Col_c  = Center(2);

y_scan = 1-Row_c : Lx-Row_c;
x_scan = 1-Col_c : Ly-Col_c;

[X,Y] = meshgrid(x_scan,y_scan);

theta=atan2(Y,X);
theta(isnan(theta)==1)=0;

R = exp(1i*theta);%M: cos + i sin

%%
h=[1 -1];
g = h(end:-1:1);

%R = exp(1i*theta);%M: cos + i sin
kappa = alp.^(-2);
kappa = kappa/(8*lam^2);

vx = zeros(size(x)); % vx and vy hold the vector fields described in the manuscript
vy = zeros(size(x));

if isempty(varargin), % number of default iterations
    MAX_ITER = 100;
else
    MAX_ITER = varargin{1};
end

wb = waitbar(0,'Please wait...');
for iter = 1:MAX_ITER,
    waitbar(iter/MAX_ITER,wb)
    % apply 'A' to v (see eq 13 for the definition of A)
    ux = alp*vx; % apply Lambda
    uy = vy;
    u = R.*(ux + 1i*uy); % apply R
    ux = conv2(real(u),g); ux = ux(:,1:end-1); % apply Delta
    uy = conv2(imag(u),g'); uy = uy(1:end-1,:); % M: g=[-1 1]   g'=[-1;1]
    % subtract from x
    u = x - lam*(ux + uy);
    % now apply A'
    ux = conv2(u,h); ux = ux(:,2:end);
    uy = conv2(u,h'); uy = uy(2:end,:);
    u = lam*conj(R).*(ux + 1i*uy);
    vx = vx + kappa*real(u)*alp;
    vy = vy + kappa*imag(u);
    % now apply the threshold

    m = abs(vx + 1i*vy);
    ind = (m > 10^(-10));
    m = max(1,m);
%    m = max(10^(-10),m); % for stability
    vx(ind) = vx(ind)./m(ind);
    vy(ind) = vy(ind)./m(ind);
end
close(wb);

% compute A*v
vx = alp*vx; % apply Lambda
u = R.*(vx + 1i*vy); % apply R
ux = conv2(real(u),g); ux = ux(:,1:end-1); % apply Delta
uy = conv2(imag(u),g');  uy = uy(1:end-1,:);
% subtract from x
y = x - lam*(ux + uy);