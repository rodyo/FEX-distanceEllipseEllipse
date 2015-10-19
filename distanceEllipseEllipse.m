% distanceEllipseEllipse    Computes the distances between two 
%                           arbitrary ellipses (in 3D)
%
% USAGE:
%   [min_dist, fp_min, fs_min] = ...
%       distanceEllipseEllipse(a,b,c,u,v)
%
% The input arguments are 
%
% =============================================================
%  name     description                           size
% =============================================================
%   a       Semi-major axes for both ellipses     (1x2)
%   b       Semi-minor axes for both ellipses     (1x2)
%   c       Locations of ellipse-centers          {(1x3),(1x3)}
%   u       Directions of Ellipse's primary axes  {(1x3),(1x3)}
%           (unit vectors)
%   v       Direction of Ellipse's secondary axis {(1x3),(1x3)}
%           (unit vectors)            
%
% The output arguments are
% =============================================================
%  name     description                           size
% =============================================================
% min_dist  minimum distance between the two      (1x1)  
%           ellipses 
% fp_min    corresponding true anomaly on the     (1x1)
%           "primary" ellipse (-pi <= f <= +pi)         
% fs_min    corresponding true anomaly on the     (1x1)
%           "secondary" ellipse (-pi <= f <= +pi)         
%
% This is a MATLAB-efficient implementation of the algorithm 
% described in:
%
% Ik-Sung Kim: "An algorithm for finding the distance between two 
% ellipses". Commun. Korean Math. Soc. 21 (2006), No.3, pp.559-567
%
% See also distanceEllipsePoint.
function [min_dist, fp_min, fs_min] = distanceEllipseEllipse(a,b,c,u,v)
% Please report bugs and inquiries to: 
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace sàrl
% Licence    : BSD


% ELEMENTARY EXAMPLE:
%{ 
%   Ellipse1   Ellipse2(=circle)
a = [2.0       1.0];
b = [0.5       1.0];
c = {[0,0,0],[-2,2,0]}; % location of centers
u = {[1,0,0], [1,0,0]}; % both oriented in XY-plane
v = {[0,1,0], [0,1,0]}; % to visualize them more easily

% plot the ellipses
u  = cellfun(@(x)x./norm(x), u, 'UniformOutput', false);
v  = cellfun(@(x)x./norm(x), v, 'UniformOutput', false);
f  = linspace(0,2*pi,100);
E1 = a(1)*u{1}'*cos(f) + b(1)*v{1}'*sin(f) + repmat(c{1}',1,numel(f));
E2 = a(2)*u{2}'*cos(f) + b(2)*v{2}'*sin(f) + repmat(c{2}',1,numel(f));
figure(1), clf, hold on
plot(E1(1,:),E1(2,:),'r', E2(1,:),E2(2,:),'b')
axis equal

% run routine
[min_dist, fp_min, fs_min] = ...
    distanceEllipseEllipse(a,b,c,u,v);

% plot the minimum distance returned
p1 = a(1)*u{1}'*cos(fp_min) + b(1)*v{1}'*sin(fp_min) + c{1}';
p2 = a(2)*u{2}'*cos(fs_min) + b(2)*v{2}'*sin(fs_min) + c{2}';
line([p1(1) p2(1)]', [p1(2) p2(2)]', 'color', 'k')
%}

% If you find this work useful, please consider a donation:
% https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6G3S5UYM7HJ3N
    

    % algorithm parameters    
    maxiters = 25;
    tolFun   = 1e-12; % == (1e6)^2, tolerance on distance-squared
    tolX     = 1e-6;  % tolerance on anomaly
    
    % error traps
    assert(isvector(a) && isvector(b) && numel(a)==2 && numel(b)==2,...
        [mfilename ':ab_not_2Dvectors'],...
        'Semi major and minor axes ''%s'' and ''%s'' must be 2-element vectors.',...
        inputname(1), inputname(2));
    assert(iscell(c) && numel(c)==2 && all(cellfun('prodofsize',c)==3),...
        [mfilename ':c_not_cell_or_3Dvector'],...
        ['Coordinates of the ellipse-centers ''%s'' must be given as \n',...
        'two 3-D vectors in a cell-array.'], inputname(3));    
    assert(iscell(u) && numel(u)==2 && all(cellfun('prodofsize',u)==3),...
        [mfilename ':u_not_cell_or_3Dvector'],...
        ['Primary axes ''%s'' of both ellipses must be given as \n',...
        'two 3-D unit-vectors in a cell-array.'], inputname(4));    
    assert(iscell(v) && numel(v)==2 && all(cellfun('prodofsize',v)==3),...
        [mfilename ':v_not_cell_or_3Dvector'],...
        ['Secondary axes ''%s'' of both ellipses must be given as \n',...
        'two 3-D unit-vectors in a cell-array.'], inputname(5))
    assert(all( cellfun(@(x,y) x(:)'*y(:) <= 3*eps(max([x(:); y(:)])), u,v) ),...
        [mfilename ':u_v_not_perpendicular'],...
        'Primary and secondary axes must be perpendicular to each other.');
        
    % make sure everything is correct shape & size
    c{1} = c{1}(:); u{1} = u{1}(:); v{1} = v{1}(:);
    c{2} = c{2}(:); u{2} = u{2}(:); v{2} = v{2}(:);
    
    % make sure [u] and [v] are UNIT-vectors
    u = cellfun(@(x)x./norm(x), u, 'UniformOutput', false);
    v = cellfun(@(x)x./norm(x), v, 'UniformOutput', false);
        
    % initialize some variables to speed up computation
    R{1}  = [u{1}, v{1}, cross(u{1},v{1})]; % rotation matrix to put ellipse in standard form
    R{2}  = [u{2}, v{2}, cross(u{2},v{2})]; % rotation matrix to put ellipse in standard form
    comp0 = [eye(3),[0;0;0]];               % part of a companion matrix for a quartic
    
    % initialize output
    min_dist = inf;
    fp_min   = NaN;
    fs_min   = NaN;
    
    % initial values to try
    % (8 equally distributed random values)    
    f0 = rand + linspace(0,3,8)*pi/2;
    for ii = 1:numel(f0)
    
        % initialize some values
        tinit = f0(ii);     converged  = false; 
        wt    = cell(2,1);  iterations = 0;
        dmin2 = inf;        t_hat      = [0,0];
               
        % initial point on first ellipse        
        wt{2} = R{2}*[a(2)*cos(tinit); b(2)*sin(tinit); 0] + c{2};                      % Cartesian coordinates

        % main loop
        while ~converged
            
            iterations = iterations + 1;
            if (iterations >= maxiters)
                f0(ii) = NaN; break; end
                        
            % save previous t_hat
            t_hatp = t_hat;

            % swap ellipses continuously
            for jj = 1:2

                % find optimal point on the other ellipse
                s = R{jj}.'*(wt{3-jj} - c{jj});% transform current point
                A = a(jj)*s(1);                % The constants A,B and C follow from the
                B = b(jj)*s(2);                %   condition dQ/dt = 0, with Q = Q(s,E,t) the
                C = b(jj)^2 - a(jj)^2;         %   XY-distance between point s and ellipse E

                % we have to find [t_hat], the true anomaly on the ellipse that minimizes
                % the distance between the associated point on the ellipse [E] and the
                % point [s]. The solution depends on the value of [C].
                
                % If C = 0 (ellipse = circle), the solution is easy:
                if (C==0)
                    t_hat(jj) = atan2(B,A);
                    % cos(t_hat), sin(t_hat) (more convenient this way)
                    sinth = sin(t_hat(jj));          
                    costh = cos(t_hat(jj));

                % otherwise, we have to solve a quartic eqution in A,B,C, which is
                % done most quickly by using EIG() on its companion matrix:
                else
                    % solve quartic via associated companion matrix
                    comp  = [-2*A/C, -(A^2+B^2-C^2)/C^2, 2*A/C, (A/C)^2; comp0];
                    Roots = eig(comp);
                    Roots = Roots(imag(Roots)==0);
                    
                    % extract optimal point
                    sint1  = sqrt(max(0,1-Roots.^2));   sint2 = -sint1;
                    sints  = [sint1, sint2];            costs = [Roots, Roots];
                    selld  = (s(1)-a(jj)*costs).^2 + (s(2)-b(jj)*sints).^2;
                    [dummy, tind] = min(selld(:)); %#ok PORT: support MATLAB < R2009a
                    sinth = sints(tind);                costh = costs(tind);
                    t_hat(jj) = atan2(sinth, costh);
                    
                end
                
                % get Cartesian coordinates of the corresponding point 
                % on the jj-th ellipse                
                wt{jj} = R{jj}*[a(jj)*costh; b(jj)*sinth; 0] + c{jj};

            end % for (ellipse swapper)

            % distance-squared between the current optimal points
            diffvec = wt{1} - wt{2};     % difference vector between the two optimal points
            dmin2p  = dmin2;             % store previous value for convergence check
            dmin2   = diffvec.'*diffvec; % distance is magnitude of difference vector
                        
            % check for convergence
            converged = (dmin2 >= dmin2p) || ...    % distance must DECREASE every iteration
               (abs(dmin2-dmin2p) <= tolFun) || ... % diff. between two consecutive distances is smaller than TolFun
                all(abs(t_hatp-t_hat) <= tolX);     % diff. between two consecutive true anomalies is smaller than TolFun
                       
        end % algorithm while-loop
        
        if all(isnan(f0))            
            warning('distanceEllipseEllipse:maxiters_exceeded',...
                ['Maximum number of iterations was exceeded for all initial estimates; ',...
                'results may not be accurate.']);            
        end
        
        % compute the real distance. If this is less than the stored value,
        % replace all corresponding entries
        new_distance = sqrt(dmin2);        
        if (new_distance < min_dist)
            min_dist = new_distance;
            fp_min = t_hat(1);
            fs_min = t_hat(2);
        end
            
    end % for all initial values
            
end % function (Kim's method)
