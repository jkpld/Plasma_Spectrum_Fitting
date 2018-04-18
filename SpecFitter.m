classdef SpecFitter < handle
    properties
        bg
        x
        y
        fitBG
        extraCostTerms
        
        BG_knts0 = [ ...
          182       0.0059;
          200       0.0148;
          270       0.1668;
          320       0.3965;
          409       0.6882;
          448       0.6365;
          500            1;
          570       0.7562;
          650       0.3793;
          740       0.1091;
          820       0.0174;
          875       0.0032];
      
        BG_ly = [0, 0, 0.075, 0.255, 0.535, 0.525, 0.99, 0.735, 0.355, 0.085, 0, 0].';
        BG_uy = [0.035, 0.035, 0.195, 0.485, 0.715, 0.665, 1.201, 0.785, 0.405, 0.135, 0.045, 0.025].';
%         BG_ux = [192, 210, 280, 330, 419, 453, 515, 580, 660, 750, 830, 876].';
%         BG_lx = [181, 190, 260, 310, 399, 443, 500, 560, 640, 730, 810, 865].';
        
%         BG_ly = [0, 0, 0.075, 0.385, 0.485, 0.525, 0.99, 0.735, 0.355, 0.085, 0, 0].';
%         BG_uy = [0.035, 0.035, 0.365, 0.435, 0.715, 0.665, 1.001, 0.785, 0.405, 0.135, 0.045, 0.025].';
        BG_ux = [192, 210, 280, 330, 419, 473, 515, 580, 660, 750, 830, 876].';
        BG_lx = [181, 190, 260, 310, 399, 443, 490, 560, 640, 730, 810, 865].';
        
        BG_fitter
        BG_knts
        
        robust_spec_fit = false;
    end
    properties (Hidden)
        iteration = 0;
        myeps
        y_i
    end
    methods
        function obj = SpecFitter(x,y,fitBG)
            if nargin == 0
                obj.x = 180:900;
                obj.BG_knts = [1;obj.BG_knts0(:)];
                return;
            end
            obj.x = x(:)';
            obj.y = y(:);
            
            obj.myeps = eps(max(abs(y)));

            options = optimoptions('fmincon',...
                'display','none',...
                'OptimalityTolerance',1e-3,...
                'StepTolerance',1e-3);
            
            obj.BG_fitter = @(p0,y) fmincon(@(p) BG_cost(obj,p,y,true), p0, [],[],[],[],[0;obj.BG_lx;obj.BG_ly], [inf; obj.BG_ux; obj.BG_uy],[],options);
            
            if fitBG
%                 min_y = mean(y(x<220 | x>800));
%                 tmp = max(ordfilt2(y(:), 15, ones(91,1),'symmetric'), 0);
%                 tmp = imerode(max(min(y(:), ordfilt2(y(:), ceil(51*0.3), ones(51,1),'symmetric')), 0),ones(11,1));
%                 obj.myeps = 2.22044604925031e-16;
                % Smooth spectrum
                y = imfilter(y,fspecial('average',[5,1]));
%                 [mad(y,1), std(y), std(y)/mad(y,1)]
%                 apply rolling-ball filter
                xR=200;

                yR = 3/sqrt(std(y)/mad(y,1));
                th = linspace(0,pi,xR*2);
                se = offsetstrel(yR*sin(th)'); 
                y = imdilate(imerode(y,se),se);

%                 tmp = y(:);%-min_y;
%                 figure(1)
%                 clf(1)
%                 line(x,tmp)
                tmp = y(:);
                [~,idx] = min(abs(x-500));
                obj.BG_knts = obj.BG_fitter([y(idx); obj.BG_knts0(:)], tmp);
                obj.bg = evaluateBG(obj);
%                 obj.BG_knts
%                 line(x,obj.bg,'color','r')

            else
                obj.bg = zeros(size(obj.y));
            end
            obj.y_i = obj.y - obj.bg;
            obj.fitBG = fitBG;

            
        end
        
        function bg = evaluateBG(obj,knts)
            if nargin < 2
                p = obj.BG_knts;
            else
                p = knts;
            end
                
            A = p(1);
            p = reshape(p(2:end),[],2);
            bg = A*pchip(p(:,1),p(:,2),obj.x).';
        end
        
        function cost = BG_cost(obj,p,y,robust)
            A = p(1);
            p = reshape(p(2:end),[],2);
            yh = A*pchip(p(:,1),p(:,2),obj.x).';
            err = y - yh;
            if robust
                w = iBisquare(err, obj.myeps);
            else
                w = ones(size(err));
            end
            w(obj.x>490 & obj.x<540) = w(obj.x>490 & obj.x<540)*5.0;
            cost = err.' * (w.*err)+ 100*p(3:4,2).'*p(3:4,2) + 2*p(5,2)^2;% + 100*A.^2;
        end
        
        function [cost, grad] = fitSpec(obj, p)
            % model parameters
            n = numel(p)/3;
            p = reshape(p,n,3);

            % gaussian argument
            arg = (obj.x - p(:,2))./p(:,3); % n x X
            arg_e = exp(-arg.^2/2); % n x X

            % individual gaussians
            yh =  arg_e .* p(:,1); % n x X

            % error
%             if obj.fitBG
%                 Y = obj.y - obj.bg; % X x 1
%             else
%                 Y = obj.y;
%             end
            err = obj.y_i - sum(yh)';
            
            if obj.robust_spec_fit
                w = iBisquare_spec(err, obj.myeps);
%                 w(~robustInd) = 1;
            else
                w = ones(size(err));
            end
            cost = (err'*(w.*err)); % 1 x 1
            err = w.*err;
            % Gradiant
            tmp2 = arg.*arg_e.*p(:,1)./p(:,3);

            dA = arg_e * err; % n x 1
            dmu = tmp2 * err; % n x 1
            dsig = (arg.*tmp2) * err; % n x 1

            grad = -2* ([dA; dmu; dsig]);
            
            [ca,ga] = obj.extraCostTerms(p(:));
            
%             size(cost)
%             size(ga)
%             size(ca)
%             size(ga)
            cost = cost + ca;
            grad = grad + ga;
            
%             if obj.fitBG
%                 if obj.iteration > 1
%                     tmp = max(ordfilt2(BG, 25, ones(91,1),'symmetric'), 0);
%                     obj.BG_knts = obj.BG_fitter(obj.BG_knts, tmp);
%                     obj.bg = evaluateBG(obj);
%                 end
%                 obj.iteration = obj.iteration + 1;
%             end
        end
    end
end

function w = iBisquare(r,myeps)

% Only use non-NaN residuals to compute median
valid = ~isnan( r );
% And bound the median away from zero
s = max( 1e8 * myeps, median( abs( r(valid) ) ) );
% r = r/(20*s);
r = r/(20*s);

% Covert the residuals to wrights
w = zeros( size( r ) );
idx = abs( r ) < 1;
w(idx) = (1 - r(idx).^2).^2;
w(valid & r<0) = w(valid & r<0)*1.1;
% Everything with NaN residual should have zero weight
w(~valid) = 0;
end

function w = iBisquare_spec(r,myeps)

% Only use non-NaN residuals to compute median
valid = ~isnan( r );
% And bound the median away from zero
s = max( 1e8 * myeps, median( abs( r(valid) ) ) );
% r = r/(20*s);
r = r/(20*s);

% Covert the residuals to wrights
w = zeros( size( r ) );
idx = abs( r ) < 1;
w(idx) = (1 - r(idx).^2).^2;
% w(valid & r>0) = w(valid & r>0)*1.1;
% Everything with NaN residual should have zero weight
w(~valid) = 0;
end