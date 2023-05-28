function [C, T] = BSpline_approx(curve)
% B-spline approximation
%
% Input arguments: 
% curve:
%    3-by-m matrix of the original curve data
% 
% Output arguments:
% C:
%    3-by-m matrix of the projected points of the original data on the b-spline
% T: 
%    3-by-m matrix of the normalized tangential vectors of the projected points on the b-spline

M = curve; 
k = 3; % 3-order B-spline

% sample time of data
s = 0:1:size(M,2)-1; 
x = (s-s(1))/(s(end)-s(1)); 

% n is the number of control points
options = optimset('Display', 'iter', 'TolX', 1);
fun = @(n) error(k,n,M,x);
% to get the closest curve
n = fminbnd(fun,k,size(M,2),options); 

% create knot vector
i = k+1:1:n;
t = [zeros(1,k) (i-k)/(floor(n)+1-k) ones(1,k)];

% get control points
D = bspline_approx(k,t,x,M); 

% the projected points on B-spline curve
C = bspline_deboor(k,t,D,x); 

% the knots and control points of the 1st deriv of the curve
[dknots,dctrl] = bspline_deriv(k,t,D); 
% the tangential vectors of the projected points
T = bspline_deboor(k-1,dknots,dctrl,x);
% get normalized tangential vectors
for j = 1:size(T,2)
     T(:,j) = T(:,j)/norm(T(:,j)); 
end

% plot the original curve and its approximated B-spline
% figure(1)
% hold on
% axis equal
% bspline = plot3(C(1,:), C(2,:), C(3,:), 'r', 'LineWidth',1); % B-spline
% ori_curve = plot3(M(1,:), M(2,:), M(3,:), 'kx', 'LineWidth',0.5); % original curve
% legend([bspline,ori_curve],'The approximated b-spline', 'the original curve')
% hold off

end

function E = error(k,n,M,x)
% B-spline approximation averaged error
%
% Input arguments:
% k:
%    3, B-spline order 
% n:
%    number of control points
% M:
%    3-by-m matrix of original curve
% x:
%    B-spline values corresponding to which data points are observed
% 
% Output arguments:
% E
%    averaged pairwise error between the original data and its projected
%    points on the b-spline

    i = k+1:1:n;
    t = [zeros(1,k) (i-k)/(floor(n)+1-k) ones(1,k)];
    D = bspline_approx(k,t,x,M);
    E = bspline_error(k,t,D,M,x);
    E = mean(E);
end

function D = bspline_approx(k,t,x,M)
% B-spline curve control point approximation.
%
% Input arguments:
% k:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% x:
%    B-spline values corresponding to which data points are observed
% M:
%    3-by-m matrix of observed data points, possibly polluted with noise,
%
% Output arguments:
% D:
%    d-by-n matrix of control points

B = bspline_basismatrix(k,t,x);
Q = M * B;
D = Q / (B'*B);
end

function E = bspline_error(k,t,D,M,x)
% B-spline approximation error.
%
% Input arguments:
% k:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% D:
%    control points
% M:
%    points whose error to compute
% x:
%    parameter values w.r.t. which distance from m is minimized
%
% Output arguments:
% E:
%    sum of squares approximation error

Y = bspline_deboor(k,t,D,x);
M = reshape(M, size(M,1), size(M,2), 1);
E = sum((M-Y).^2);
E = sqrt(E); % pairwise distance
end

function [B,x] = bspline_basismatrix(n,t,x)
% B-spline basis function value matrix B(n) for x.
%
% Input arguments:
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% x (optional):
%    an m-dimensional vector of values where the basis function is to be
%    evaluated
%
% Output arguments:
% B:
%    a matrix of m rows and numel(t)-n columns

if nargin > 2
    B = zeros(numel(x),numel(t)-n);
    for j = 0 : numel(t)-n-1
        B(:,j+1) = bspline_basis(j,n,t,x);
    end
else
    [b,x] = bspline_basis(0,n,t);
    B = zeros(numel(x),numel(t)-n);
    B(:,1) = b;
    for j = 1 : numel(t)-n-1
        B(:,j+1) = bspline_basis(j,n,t,x);
    end
end

end

function [y,x] = bspline_basis(j,n,t,x)
% B-spline basis function value B(j,n) at x.
%
% Input arguments:
% j:
%    interval index, 0 =< j < numel(t)-n
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% x (optional):
%    value where the basis function is to be evaluated
%
% Output arguments:
% y:
%    B-spline basis function value, nonzero for a knot span of n

validateattributes(j, {'numeric'}, {'nonnegative','integer','scalar'});
validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
validateattributes(t, {'numeric'}, {'real','vector'});
assert(all( t(2:end)-t(1:end-1) >= 0 ), ...
    'Knot vector values should be nondecreasing.');
if nargin < 4
    x = linspace(t(n), t(end-n+1), 100);  % allocate points uniformly
else
    validateattributes(x, {'numeric'}, {'real','vector'});
end
assert(0 <= j && j < numel(t)-n, ...
    'Invalid interval index j = %d, expected 0 =< j < %d (0 =< j < numel(t)-n).', j, numel(t)-n);

y = bspline_basis_recurrence(j,n,t,x);

function y = bspline_basis_recurrence(j,n,t,x)

y = zeros(size(x));
if n > 1
    b = bspline_basis(j,n-1,t,x);
    dn = x - t(j+1);
    dd = t(j+n) - t(j+1);
    if dd ~= 0  % indeterminate forms 0/0 are deemed to be zero
        y = y + b.*(dn./dd);
    end
    b = bspline_basis(j+1,n-1,t,x);
    dn = t(j+n+1) - x;
    dd = t(j+n+1) - t(j+1+1);
    if dd ~= 0
        y = y + b.*(dn./dd);
    end
elseif t(j+2) < t(end)  % treat last element of knot vector as a special case
    y(t(j+1) <= x & x < t(j+2)) = 1;
else
    y(t(j+1) <= x) = 1;
end

end

end


function [C,U] = bspline_deboor(n,t,P,U)
% Evaluate explicit B-spline at specified locations.
%
% Input arguments:
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% P:
%    control points, typically 2-by-m, 3-by-m or 4-by-m (for weights)
% u (optional):
%    values where the B-spline is to be evaluated, or a positive
%    integer to set the number of points to automatically allocate
%
% Output arguments:
% C:
%    points of the B-spline curve

validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
d = n-1;  % B-spline polynomial degree (1 for linear, 2 for quadratic, etc.)
validateattributes(t, {'numeric'}, {'real','vector'});
assert(all( t(2:end)-t(1:end-1) >= 0 ), 'bspline:deboor:InvalidArgumentValue', ...
    'Knot vector values should be nondecreasing.');
validateattributes(P, {'numeric'}, {'real','2d'});
nctrl = numel(t)-(d+1);
assert(size(P,2) == nctrl, 'bspline:deboor:DimensionMismatch', ...
    'Invalid number of control points, %d given, %d required.', size(P,2), nctrl);
if nargin < 4
    U = linspace(t(d+1), t(end-d), 10*size(P,2));  % allocate points uniformly
elseif isscalar(U) && U > 1
    validateattributes(U, {'numeric'}, {'positive','integer','scalar'});
    U = linspace(t(d+1), t(end-d), U);  % allocate points uniformly
else
    validateattributes(U, {'numeric'}, {'real','vector'});
    assert(all( U >= t(d+1) & U <= t(end-d) ), 'bspline:deboor:InvalidArgumentValue', ...
        'Value outside permitted knot vector value range.');
end

m = size(P,1);  % dimension of control points
t = t(:).';     % knot sequence
U = U(:);
S = sum(bsxfun(@eq, U, t), 2);  % multiplicity of u in t (0 <= s <= d+1)
I = bspline_deboor_interval(U,t);

Pk = zeros(m,d+1,d+1);
a = zeros(d+1,d+1);

C = zeros(size(P,1), numel(U));
for j = 1 : numel(U)
    u = U(j);
    s = S(j);
    ix = I(j);
    Pk(:) = 0;
    a(:) = 0;

    % identify d+1 relevant control points
    Pk(:, (ix-d):(ix-s), 1) = P(:, (ix-d):(ix-s));
    h = d - s;

    if h > 0
        % de Boor recursion formula
        for r = 1 : h
            q = ix-1;
            for i = (q-d+r) : (q-s)
                a(i+1,r+1) = (u-t(i+1)) / (t(i+d-r+1+1)-t(i+1));
                Pk(:,i+1,r+1) = (1-a(i+1,r+1)) * Pk(:,i,r) + a(i+1,r+1) * Pk(:,i+1,r);
            end
        end
        C(:,j) = Pk(:,ix-s,d-s+1);  % extract value from triangular computation scheme
    elseif ix == numel(t)  % last control point is a special case
        C(:,j) = P(:,end);
    else
        C(:,j) = P(:,ix-d);
    end
end

function ix = bspline_deboor_interval(u,t)
% Index of knot in knot sequence not less than the value of u.
% If knot has multiplicity greater than 1, the highest index is returned.

    i = bsxfun(@ge, u, t) & bsxfun(@lt, u, [t(2:end) 2*t(end)]);  % indicator of knot interval in which u is
    [row,col] = find(i);
    [row,ind] = sort(row);  %#ok<ASGLU> % restore original order of data points
    ix = col(ind);
end
end

function [dknots, dctrl] = bspline_deriv(order, knots, ctrl)
% Knots and control points associated with the derivative of B-spline curve.
%
% Input arguments:
% order:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% knots:
%    knot vector
% ctrl:
%    control points, typically 2-by-m, 3-by-m, or 4-by-m (for weights)
%
% Output arguments:
% dctrl:
%    control points of the derivative of the input B-spline curve
% dknots:
%    the new knot vector associated with the derivative B-spline curve

p = order - 1;
tmp = size(ctrl);
n = tmp(2)-1;
dim = tmp(1);

% derivative knots
dknots = knots(2:max(size(knots))-1);

% derivative control points
dctrl = zeros(dim,n);
for i = 1 : n
    dctrl(:,i) = (p / (knots(i+p+1) - knots(i+1))) * (ctrl(:,i+1) - ctrl(:,i));
end

end