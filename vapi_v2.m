%% VAPI METHOD - TOSCANA-EMILIA-MARCHE
% Evaluation of rainfall heights of assigned return period T, procedure:
% zone A = 1; zone B = 2; zone C = 3; zone D = 4; zone E = 5; zone F = 6; zone G = 7.     

clear all
close all
warning('off','all')

%% Loading data and variables 

location = '\bf Zena, Monte Ceresa (Emilia-Romagna, Italy)';
load prpb
load time_int
load tr

% Zone parameter
r      = 0.89;    % -,  for Toscana-Emilia-Marche (Rep GNDCI-Line 1). 
m_hg   = 65;      % mm, average daily rainfall (isometric map).              % modify
m_hh   = 25;      % mm, average hourly rainfall (isometric map).             % modify
z      = 3;       % -,  type of zone                                         % modify

% Precipitation

% Senio (Valsenio)          % m_hg = 70
% rain3  = 36;      % mm                                                       % modify
% rain12 = 112.8;   % mm                                                       % modify
% rain24 = 205.8;   % mm                                                       % modify
% rain36 = 242.8;   % mm                                                       % modify
% rain48 = ??;      % mm                                                       % modify

% Santerno (Castel del Rio) % m_hg = 80
% rain3  = 30;      % mm                                                       % modify
% rain12 = 102.6;   % mm                                                       % modify
% rain24 = 169.8;   % mm                                                       % modify
% rain36 = 225.6;   % mm                                                       % modify
% rain48 = ??;      % mm                                                       % modify

% Zena (Monte Ceresa) % m_hg = 65
rain1  = 10.8;    % mm                                                       % modify
rain2  = 21.2;    % mm                                                       % modify
rain3  = 30.8;    % mm                                                       % modify
rain4  = 39.8;    % mm                                                       % modify
rain5  = 47.2;    % mm                                                       % modify
rain6  = 51.8;    % mm                                                       % modify
rain12 = 78.2;    % mm                                                       % modify
rain24 = 147.4;   % mm                                                       % modify
rain36 = 182.2;   % mm                                                       % modify
% rain48 = ??;      % mm                                                       % modify

% Lamone (San Cassiano)     % m_hg = 75
% rain3  = 41;      % mm                                                       % modify
% rain12 = 132.4;   % mm                                                       % modify
% rain24 = 207.2;   % mm                                                       % modify
% rain36 = 258;     % mm                                                       % modify
% rain48 = ??;      % mm                                                       % modify

% Marzeno (Trebbio)     % m_hg = 75
% rain3  = 38;      % mm                                                       % modify
% rain12 = 138.4;   % mm                                                       % modify
% rain24 = 213.2;   % mm                                                       % modify
% rain36 = 250.4;   % mm                                                       % modify
% rain48 = ??;      % mm                                                       % modify

%% Programming

for i = 1:7:length(prpb)
    zone(:,:,i) = prpb(i:(i+6),:);        % R³ matrix with 7 layer of zone
end
zone = zone(:,:,1:7:length(prpb));
n = (log(m_hg)-log(m_hh)-log(r))/log(24);
for i = 1:length(time_int)
    m_hd(i) = m_hh*(time_int(i)^n);
end
for i = 1:length(zone)                                                               
    for j = 1:length(tr)
        for k = z
            T = @(x) (1/(1-exp(-zone(i,3,k)*exp(-x*zone(i,4,k))-zone(i,1,k)*zone(i,3,k)^(1/zone(i,2,k))*exp((-x*zone(i,4,k))/zone(i,2,k)))))-tr(j);
            kt(i,j) = fzero(T,0);
        end
    end
end
kt = kt(1:length(zone), 1:length(tr))';                                              
ht = kt.*m_hd;
a = (sum(log(time_int)'.^2)*sum(log(ht)')-sum(log(time_int)')*sum(log(time_int)'.*log(ht)'))/(length(time_int)*sum(log(time_int)'.^2)-sum(log(time_int)')^2);
b = (length(time_int)*sum(log(time_int)'.*log(ht)')-sum(log(time_int)')*sum(log(ht)'))/(length(time_int)*sum(log(time_int)'.^2)-sum(log(time_int)')^2);
at = exp(a);
for i = 1:length(at)
    for j = 1:length(time_int)
        ht_d(i,j) = at(i)*(time_int(j)^b(i));
        it_d(i,j) = at(i)*(time_int(j)^(b(i)-1));  
    end
end

%% Plots

% Date settings
lw = 1.15;
minx = min(time_int)-1;  
maxx = max(time_int);
miny = 0;
maxy = round(max(max(ht_d)),-2)+50;

% Figure settings
figure(1)
set(1,'Units','centimeters','Position',[1,1,35,19]);
set(1,'defaultaxesfontname', 'CMU Serif');
set(1,'defaulttextfontname', 'CMU Serif');
set(gca,'FontSize',20,'Box','on');
plot(time_int,ht_d,'LineWidth',4*lw), hold on
plot(1,  rain1, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(2,  rain2, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(3,  rain3, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(4,  rain4, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(5,  rain5, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(6,  rain6, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(12, rain12,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(24, rain24,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(36, rain36,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
% plot(48, rain48,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
vxt = minx:3:maxx;
vyt = miny:50:maxy;
set(gca,'xtick',vxt,'ytick',vyt,'XMinorTick','on','YMinorTick','on');
axis([minx maxx miny maxy]);
ylabel('rainfall depth, h (mm)');
xlabel('rainfall duration, d (h)');
nw = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.15;
text(nw(1)-2,maxy-10,location,'Color','k','FontSize',20,'VerticalAlignment','Top','HorizontalAlignment','Left')
legend('T = 5 y','T = 10 y','T = 20 y','T = 50 y','T = 100 y','T = 200 y','T = 500 y','T = 1000 y','Location','bestoutside','NumColumns',1);
legend('boxon')
grid on
set(gca,'PlotBoxAspectRatio',[2.5 1.7 2],'FontSize',20,'Ticklength',[0.015 0.00],'LineWidth',1.5,'Layer','top');

orient(figure(1),'landscape')
saveas(1,'figone','svg')
system('inkscape --export-filename=FIG1_vapi_ink.pdf figone.svg')
system('rm figone.svg')
system('pdfcrop FIG1_vapi_ink.pdf FIG_vapi_depth.pdf')
system('rm FIG1_vapi_ink.pdf')

miny = 0;
% maxy = round(max(max(it_d)),-1)+5;
maxy = 60;

figure(2)
set(2,'Units','centimeters','Position',[1,1,35,19]);
set(2,'defaultaxesfontname', 'CMU Serif');
set(2,'defaulttextfontname', 'CMU Serif');
set(gca,'FontSize',20,'Box','on');
plot(time_int,it_d,'LineWidth',4*lw), hold on
plot(1,  rain1/1, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(2,  rain2/2, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(3,  rain3/3, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(4,  rain4/4, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(5,  rain5/5, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(6,  rain6/6, 'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(12, rain12/12,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(24, rain24/24,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
plot(36, rain36/36,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
% plot(48, rain48/48,'LineWidth',lw,'Marker','o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10), hold on
vxt = minx:3:maxx;
vyt = miny:5:maxy;
set(gca,'xtick',vxt,'ytick',vyt,'XMinorTick','on','YMinorTick','on');
axis([minx maxx miny maxy]);
ylabel('rainfall intensity, i (mm h^{-1})');
xlabel('rainfall duration, d (h)');
nw = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.15;
% xx = (minx+maxx)/2;
text(nw(1)-2,maxy-1.8,location,'Color','k','FontSize',20,'VerticalAlignment','Top','HorizontalAlignment','Left')
legend('T = 5 y','T = 10 y','T = 20 y','T = 50 y','T = 100 y','T = 200 y','T = 500 y','T = 1000 y','Location','bestoutside','NumColumns',1);
legend('boxon')
grid on
set(gca,'PlotBoxAspectRatio',[2.5 1.7 2],'FontSize',20,'Ticklength',[0.015 0.00],'LineWidth',1.5,'Layer','top');

orient(figure(2),'landscape')
saveas(2,'figtwo','svg')
system('inkscape --export-filename=FIG2_vapi_ink.pdf figtwo.svg')
system('rm figtwo.svg')
system('pdfcrop FIG2_vapi_ink.pdf FIG_vapi_intensity.pdf')
system('rm FIG2_vapi_ink.pdf')