% Load the example SpectralLines object
%  - This object has emission groups for the following species
%      H, He , O, O2+, NO, OH, N2_1PS, N2_2PS, N2_1NS, OH2
%    as well as two groups with lines that were not identified
%      Unknown_N2_related, Unknown_O2_related.
%  - Note the emission lines in these groups are simply the lines that were
%    observed in the emission spectrum of our plasma source.
load('specLines.mat')
show(specLines) % Display the different EmissionGroups and their emission lines

% Load two different example spectrum.
%  - The first spectrum is from a He plasma jet surrounded with N2 and in
%    contact with a liquid surface.
%  - The second spectrum is from a He plasma jet surrounded with O2 and in
%    contact with a liquid surface.
load('exampleSpectra.mat')

% Fit the spectra with N2, and do not include the O2 related emission groups
% while fitting
groups_to_use = ~contains(specLines.names, "O2");
[p, group_idx] = specLines.fit(x1,y1, groups_to_use);

% Plot the spectrum and the fit
%  - This will result in the figure on the left above.
%
figure('color','w')
line(x1,y1,'Color',0.5*[1 1 1],'LineWidth',2)
specLines.plotfit(x1, p, group_idx)
xlim([200,870])
xlabel('Wavelength (nm)')
ylabel('Intensity')
ax2 = axes('Parent',gcf,'Position',[0.45,0.45,0.5,0.5],'Box','on','XGrid','on','YGrid','on','XLim',[565,800]);
line(x1,y1,'Color',0.5*[1 1 1],'LineWidth',2,'Parent',ax2)
specLines.plotfit(x1, p, group_idx,[],ax2)
setTheme(gcf,'light')

% Fit the spectra with O2, and do not include several of the N2 related
% emission groups while fitting
% - Note: The O2 spectra observed have a continuous background emission.
%   This background must be taken into account to correctly fit the spectrum.
%   The fit function already has an empirical profile for this continuous
%   background. Note that the background is specifically for the spectrum
%   observed in our experiments using O2 surrounding the He jet, and is not
%   necessarily general.
groups_to_use = ~contains(specLines.names, ["N2_1PS", "Unknown_N2", "N2+"]);
fit_O2_background = true;
[p, group_idx, fitter] = specLines.fit(x2,y2, groups_to_use, fit_O2_background);

% Plot the spectrum and the fit
%  - This will result in the figure on the right above.
figure
line(x2,y2,'Color',0.5*[1 1 1],'LineWidth',2)
background = fitter.bg; % extract the fit background.
specLines.plotfit(x2, p, group_idx, background)

% Adjust the scaling
ax = gca;
thrsh = 0.2e4;
for i = 1:numel(ax.Children)
    idx = ax.Children(i).YData > thrsh;
    ax.Children(i).YData(idx) = (ax.Children(i).YData(idx) - thrsh).^(0.7) + thrsh;
end

ax.XLim = [200,870];
ax.YLim(1) = -20;
ax.YTick = [0,1000,2000,18000^0.7 + 2000];
ax.YTickLabel = string([0, 1, 2, 20]);
ax.YRuler.SecondaryLabel.String = '\times10^3';
ax.YRuler.SecondaryLabel.Visible = 'on';
line(ax.XLim,thrsh*[1,1],'LineStyle','--','Color',0.4*[1 1 1])
xlabel('Wavelength (nm)')
ylabel('Intensity')
setTheme(gcf,'light')