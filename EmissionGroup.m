classdef EmissionGroup
    properties
        % Name of emission group
        name(1,1) string 
        
        % Emission lines of group
        lines(1,:) double
        
        % Function for non-linear contraints
        nonLinFitConstraint = []
        
        % Function applied to initial amplitudes before fitting
        a0_PostProcess = @(x) x;
        
        % Allowed movement of line centers whild fitting
        dc = 0.3;
        
        % Initial line width (Gaussian sigma)
        s0 = 2;
        
        % Allowed change in line width while fitting
        ds = 1;
        
%         % Initial amplitudes. If empty, then they will be automatically
%         % estimated.
%         a0 = [];
        
        % removeSmallPeaks : logical value. If true, then small peaks will
        % not attempt to be fitted.
        removeSmallPeaks = true;
    end
    
    properties %(Hidden)
        % Additional cost : function that returns an additional cost to be
        % added to the total fitting cost. This allows for additional
        % constraints to be applied by minimization..
        additionalCost = [];
    end
    methods 
        function obj = EmissionGroup(name,lines)
            if nargin ~= 0
                obj.name = name;
                obj.lines = unique(lines);
            end
        end
        
        function obj = addLine(obj,lines)
            obj.lines = unique([obj.lines, lines]);
        end
        
        function [p0,p0lb,p0ub] = getInitialParameters(obj,x,y,threshold)
            c0 = obj.lines(:);
            dc = obj.dc(:).*ones(size(c0));
            
            c0lb = c0 - dc;
            c0ub = c0 + dc;
            
            s0 = obj.s0(:).*ones(size(c0));
            ds = obj.ds(:).*ones(size(c0));
            
            s0lb = s0 - ds;
            s0ub = s0 + ds;
            
            if nargin < 4
                threshold = 0;
            end

            a0 = max(interp1(x,medfilt1(y,9,'truncate'),c0),0);
            
            if ~isempty(obj.nonLinFitConstraint)
                try
                    [~,~,~,~,corr] = obj.nonLinFitConstraint([a0;c0;s0],true(size(c0)));
                    idx = corr>0;
                    a0(idx) = a0(idx) - corr(idx);
                catch
                end
            end
            
            if obj.removeSmallPeaks
                toSmall = a0 < threshold;
            else
                toSmall = false(size(a0));
            end
            a0(~toSmall) = obj.a0_PostProcess(a0(~toSmall));
            
            a0lb = zeros(size(c0));
            a0ub = max(a0*5,0.0001);
            
            % combine
            p0 = [a0(:),c0(:),s0(:)];
            p0ub = [a0ub(:),c0ub(:),s0ub(:)];
            p0lb = [a0lb(:),c0lb(:),s0lb(:)];

            
            % remove small peaks
            p0(toSmall,:) = [];
            p0ub(toSmall,:) = [];
            p0lb(toSmall,:) = [];
            
        end
        
        function str = show(obj)
            if ~isempty(obj.lines)
                s = string(obj.lines).join(", ");
            else
                s = '';
            end
            s = sprintf('%s :\n %s\n',obj.name, s);
            if nargout > 0
                str = s;
            else
                fprintf('%s',s);
            end
        end
        
        function plot(obj,ax,col)
            if nargin < 2
                ax = gca;
            end
            if nargin < 3
                col = 'k';
            end
            
            y = [ax.YLim(:); nan] .* ones(size(obj.lines));
            x = obj.lines .* [1;1;nan];
            
            line(x(:),y(:),'Color',col,'Parent',ax,'ButtonDownFcn',@displayInfo,'UserData',obj.name,'Tag','SpecLine')
            line(obj.lines,ax.YLim(2)*ones(size(obj.lines)),'Marker','o','Color',col,'LineStyle','none','Parent',ax,'LineWidth',2,'ButtonDownFcn',@displayInfo,'UserData',obj.name,'Tag','SpecLine')
        end
        
        function a = initFitAmplitude(obj,x,y,threshold)
            if nargin < 3
                threshold = 0;
            end
            valid = ~isnan(y);
            a = interp1(x(valid),y(valid),obj.lines);
            toRemove = a < threshold;
            a(~toRemove) = obj.a0_PostProcess(a(~toRemove));
            a(toRemove) = nan;
        end
    end
end

function displayInfo(obj,evt)


xt = evt.IntersectionPoint(1);
[~,idx] = min(abs(obj.XData-xt));
x = obj.XData(idx);

y = obj.Parent.YLim(2)+0.02*diff(obj.Parent.YLim);
s = strrep(obj.UserData,'_',' ');
% s = strrep(s,'p','+');
s = s + " (" + string(x) + ")";
try delete(obj.Parent.UserData), catch ME, rethrow(ME); end
obj.Parent.UserData = text(x,y,s,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',obj.Color,'FontWeight','bold','Parent',obj.Parent,'Tag','SpecLineName');

end