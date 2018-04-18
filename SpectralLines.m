classdef SpectralLines
    properties
        emissionGroups(1,:) EmissionGroup
    end
    properties (Dependent)
        names(1,:) string
    end
    
    properties (Hidden)
        robust_spec_fit = false;
    end
    methods
        function obj = SpectralLines()
        end

        function names = get.names(obj)
            names = [obj.emissionGroups.name];
        end

        function obj = addGroup(obj,name,lines)
            obj.emissionGroups(end+1) = EmissionGroup(name,lines);
        end

        function obj = addLine(obj,name,lines)
            idx = obj.names == name;
            if any(idx)
                obj.emissionGroups(idx) = obj.emissionGroups(idx).addLine(lines);
            else
                obj.emissionGroups(end+1) = EmissionGroup(name,lines);
            end
        end

        function show(obj,name)
            if nargin < 2
                for i = 1:numel(obj.emissionGroups)
                    obj.emissionGroups(i).show()
                end
            else
                obj.emissionGroups(obj.names==name).show()
            end
        end



        function plot(obj,ax,name)
            if nargin < 2 || isempty(ax)
                ax = gca;
            end
            if nargin < 3
                name = obj.names;
            end
            % Remove any old SpecLines plotted
            toDelete = false(size(ax.Children));
            for i = 1:numel(ax.Children)
                try
                    if startsWith(ax.Children(i).Tag,'SpecLine')
                        toDelete(i) = true;
                    end
                catch
                end
            end
            delete(ax.Children(toDelete))

%             cmap = brewermap(ceil(numel(obj.emissionGroups)/64)*64,'*Spectral');
            cmap = flipud(myColorMap(ceil(numel(obj.emissionGroups)/64)*64,'v1'));
            cols = cmap(round(linspace(1,size(cmap,1),numel(obj.emissionGroups))),:);

            % sort the groups by mean wavelength for determining color
            groups = obj.emissionGroups;
            groupMean = arrayfun(@(x) median(x.lines),groups);
            [~,idx] = sort(groupMean);
            groups = groups(idx);

            for i = 1:numel(groups)
                if ismember(groups(i).name,name)
                    groups(i).plot(ax,cols(i,:))
                end
            end

            for i = 1:numel(groups)
                if ismember(groups(i).name,name)
                    text(ax.XLim(2)+0.02*diff(ax.XLim),ax.YLim(2) - (i-1)*diff(ax.YLim)/numel(groups), strrep(groups(i).name,'_',' '),'color',cols(i,:),'VerticalAlignment','top','FontWeight','bold','Tag','SpecLineLegend')
                end
            end
        end

        function [p,specGroupIdx,fitter,noiseLvl] = fit(obj, x, y, groupsToInclude, fitBG, noiseLvl,p0_user)

            if nargin < 4 || isempty(groupsToInclude)
                toInclude = 1:numel(obj.emissionGroups);
            else
                if islogical(groupsToInclude)
                    toInclude = find(groupsToInclude);
                elseif isstring(groupsToInclude)
                    toInclude = find(ismember(obj.names,groupsToInclude));
                else
                    error('fit:badInput','groupsToInclude should be a logical array or a string array with the the groups to include.')
                end
            end

            notValid = isnan(y);
            x(notValid) = [];
            y(notValid) = [];
            x = x(:);
            y = y(:);

            if nargin < 5
                fitBG = false;
            else
                fitBG = logical(fitBG(1));
            end


            if nargin < 6 || isempty(noiseLvl)
%                 bg = medfilt1(y,131,'truncate');
%                 bg = ordfilt2(y,21,ones(131,1),'symmetric');
%                 ytmp = y-bg;
%                 noiseLvl = 2*mad(ytmp(ytmp<mad(ytmp,1)),1);
                noiseLvl = 2*mad(y(x>800 | x<220),1);
            end
            
%             noiseLvl
            y = y./noiseLvl;
            % Extract groups toInclude
            groups = obj.emissionGroups(toInclude);

            % Initialize SpecFitter
            fitter = SpecFitter(x,y,fitBG);
            fitter.robust_spec_fit = obj.robust_spec_fit;
                

%             figure
%             line(x,y)
%             line(x,fitter.bg,'color','r','linewidth',2)
%             axis tight
%
%             size(fitter.bg)
%             size(y)


            % ------- Get initial Gaussian parameters

            p0user =  nargin == 7 && ~isempty(p0_user);

            threshold = 1;%~fitBG;%1;%noiseLvl;

            p0 = cell(numel(groups),1);
            p0lb = cell(numel(groups),1);
            p0ub = cell(numel(groups),1);
            groupIdx = cell(numel(groups),1);
            for i = 1:numel(groups)
                if p0user
                    groups(i).removeSmallPeaks = false;
                end
                [p0{i},p0lb{i},p0ub{i}] = groups(i).getInitialParameters(x,fitter.y-fitter.bg,threshold);
%                 if groups(i).name == "H"
%                     p0{i}
%                     p0lb{i}
%                     p0ub{i}
%                 end
                groupIdx{i} = toInclude(i)*ones(size(p0{i},1),1);
            end
% error('some err')
            p0 = cat(1,p0{:});
            p0lb = cat(1,p0lb{:});
            p0ub = cat(1,p0ub{:});
            specGroupIdx = cat(1,groupIdx{:});

            
            if p0user
                
                p0_user(:,1) = p0_user(:,1)/noiseLvl;
                p0(:,1:size(p0_user,2)) = p0_user;
                
                p0ub(:,1) = max(p0_user(:,1)*2,0);
                p0lb(:,1) = -0.0001;

            end
           

%             [p0,p0lb,p0ub]
%             error('some err')
%             [c,~,g,~,corr] = obj.emissionGroups(1).nonLinFitConstraint(p0,specGroupIdx==1)

%             [c,~,g] = obj.nonLinConstraints(p0,specGroupIdx)
%             p0 = reshape(p0,[],3);
%             [p0, specGroupIdx, obj.names(specGroupIdx)']
%             error('some err')

            if isempty(p0)
                error('fit:noLinesToFit','There are no lines to fit!')
            end

            % Non-linear constraints
            nonlcon = @(p) nonLinConstraints(obj,p,specGroupIdx);

            % Fit options
            options = optimoptions('fmincon',...
                'Algorithm','interior-point',...
                'HonorBounds',true,...
                'OptimalityTolerance',1e-4,...
                'ConstraintTolerance',1e-4,...
                'StepTolerance',1e-15,...
                'SpecifyObjectiveGradient',true,...
                'SpecifyConstraintGradient',true,...
                'display','none',...
                'MaxFunctionEvaluations',8000,...
                'MaxIterations',2000,...
                'CheckGradients',false,...
                'ScaleProblem','none');%'obj-and-constr');

            fitter.extraCostTerms = @(p) additionalCost(obj,p,specGroupIdx);
%             p0(specGroupIdx==9,:)
            fun = @(p) fitter.fitSpec(p);
            p = fmincon(fun,p0,[],[],[],[],p0lb,p0ub,nonlcon,options);


%             [c,~,g] = obj.nonLinConstraints(p,specGroupIdx)
%             size(p)
%             p = reshape(p,[],3);



            % Add the noise level back in.
            p(:,1) = p(:,1)*noiseLvl;
            fitter.y = fitter.y*noiseLvl;
            fitter.bg = fitter.bg*noiseLvl;
            
            try 
                fitter.BG_knts(1) = fitter.BG_knts(1)*noiseLvl;
            catch
            end
        end

        function [c, grad] = additionalCost(obj,p,groupIdx)
            groups = unique(groupIdx);
            c = 0;
            grad = zeros(size(p));

            for i = 1:numel(groups)
                if ~isempty(obj.emissionGroups(groups(i)).additionalCost)
                    [ci, gi] = obj.emissionGroups(groups(i)).additionalCost(p,groupIdx==groups(i));
                    c = c + ci;
                    grad = grad + gi;
                end
            end
        end

        function [c,ceq,g,geq] = nonLinConstraints(obj,p,groupIdx)
            % Non-linear constraints only apply the the gaussian amplitudes

            groups = unique(groupIdx);
            c_idx = 0;
            N = numel(groups);
            c = cell(N,1);
            ceq = cell(N,1);
            g = cell(N,1);
            geq = cell(N,1);

            for i = 1:numel(groups)
                if ~isempty(obj.emissionGroups(groups(i)).nonLinFitConstraint)
                    [c{i},ceq{i},g{i},geq{i}] = obj.emissionGroups(groups(i)).nonLinFitConstraint(p,groupIdx==groups(i));


                    if ~isempty(g{i})
                        g{i}(:,2) = g{i}(:,2) + c_idx;
                    end
                    if ~isempty(geq{i})
                        geq{i}(:,2) = geq{i}(:,2) + c_idx;
                    end

                    c_idx = c_idx + numel(c{i});
                end
            end
            nEmpt = ~cellfun(@isempty, c);

            c = cat(1,c{nEmpt});
            ceq = cat(1,ceq{nEmpt});
            g = cat(1,g{nEmpt});
            geq = cat(1,geq{nEmpt});

            if ~isempty(g)
                g = sparse(g(:,1),g(:,2),g(:,3),numel(p),numel(c));
            end
            if ~isempty(geq)
                geq = sparse(geq(:,1),geq(:,2),geq(:,3),numel(p),numel(c));
            end
        end

        function plotfit(obj,x,p,idx,bg,ax)
            if nargin < 6
                ax = gca;
            end
            if nargin < 5 || isempty(bg)
                bg = zeros(size(x));
            end

            % Remove any old SpecLines plotted
            toDelete = false(size(ax.Children));
            for i = 1:numel(ax.Children)
                try
                    if startsWith(ax.Children(i).Tag,'FitSpecLine')
                        toDelete(i) = true;
                    end
                catch
                end
            end
            delete(ax.Children(toDelete))

            groups = obj.emissionGroups;
            cmap = flipud(myColorMap(ceil(numel(groups)/64)*64,'v1'));
            cols = cmap(round(linspace(1,size(cmap,1),numel(groups))),:);

            % sort the groups by mean wavelength for determining color
            groupMean = arrayfun(@(x) median(x.lines),groups);
            [~,sidx] = sort(groupMean);
            [~,sidx] = sort(sidx);
            cols = cols(sidx,:);

            gs = zeros(size(x));
            for i = size(p):-1:1
                g = p(i,1)*exp(-(x - p(i,2)).^2/(2*p(i,3)^2));
                gs = gs + g;
                valid = x-p(i,2) < 4.3*p(i,3) & x-p(i,2) > -4.3*p(i,3); %g>p(i,1)/100;
                lh(i) = line(x(valid),g(valid)+bg(valid),'color',cols(idx(i),:),'linestyle','-','linewidth',0.5,'Parent',ax,'Tag','FitSpecLineIndividual');
            end
            line(x,gs+bg,'color','r','linestyle','-','Parent',ax,'Tag','FitSpecLineSum');
%             uistack(lh,'top');
            if ~all(bg==0)
                line(x,bg,'color','k','linewidth',1,'Parent',ax,'Tag','FitSpecLineBG');
            end
        end

        function create_latex_table_body(obj,option)
            fprintf('\\begin{table}\n')
            fprintf('\\caption{\\label{}} \n')
            fprintf('\\begin{indented}\n')
            fprintf('\\lineup\n')

            if nargin == 1 || strcmp(option,'down')

                fprintf('\\item[]\\begin{tabular}{@{}*{%d}{l}}\n',numel(obj.emissionGroups))
                fprintf('\\br\n')

                fprintf('%s\n', obj.names.join(' & ') + "\cr")
                fprintf('\\mr\n')

                for i = 1:numel(obj.emissionGroups)
                    l = round(obj.emissionGroups(i).lines,1);
                    lines(i,1:numel(l)) = string(l);
                end
                lines(ismissing(lines)) = " ";
                lines = lines';
                lines = lines.join(" & ");
                for i = 1:numel(lines)
                    fprintf('%s\\cr\n', lines(i));
                end
            elseif strcmp(option,'over')
                fprintf('\\item[]\\begin{tabular}{@{}p{3cm}p{\\textwidth - 3cm}}\n')
                fprintf('\\br\n')
                fprintf('Specie & Lines\\cr\n')
                fprintf('\\mr\n')

                for i = 1:numel(obj.emissionGroups)
                    fprintf('%s & %s \\cr\n', obj.names(i), string(round(obj.emissionGroups(i).lines,1)).join(", "))
                end
            end

            fprintf('\\br\n')
            fprintf('\\end{tabular}\n')
            fprintf('\\end{indented}\n')
            fprintf('\\end{table}\n')

        end
    end

    methods (Static)
        function y = evaluateModel(x,p)
            y = zeros(size(x));
            for i = size(p):-1:1
                g = p(i,1)*exp(-(x - p(i,2)).^2/(2*p(i,3)^2));
                y = y + g;
            end
        end
    end
end
