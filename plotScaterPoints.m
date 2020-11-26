function [rmse,pcc] = plotScaterPoints(x,y,x_label, y_label)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rmse = sqrt(mean((y-x).^2));
pcc = corr(y,x); 
% h1 = figure;
% set(h1,'Position',[300 300 1200 320]);
% ax_min = min(round(min(x)),round(min(y)));
% ax_max = max(ceil(max(x)), ceil(max(y)));
ax_min = -250;
ax_max = 250;
% scatter(x,y,sz,'filled');
dscatter(x,y);
hold on;
x1 = [ax_min ax_max];
y1 = [ax_min ax_max];
line(x1,y1,'Color','r','LineWidth',1);
x2 = [ax_min ax_max];
y2 = [ax_max ax_max];
line(x2,y2,'Color','k');%,'LineWidth',1);
x3 = [ax_max ax_max];
y3 = [ax_min ax_max];
line(x3,y3,'Color','k');%,'LineWidth',1);
xlim([ax_min ax_max]);
ylim([ax_min ax_max]); 
xlabel(x_label);
ylabel(y_label);
TXT_1=['RMSD  = ';'PCC   = '];
TXT_2=strjust(num2str([rmse;pcc],'%.2f'),'right');
% text(0.05*(ax_max-ax_min),0.80*(ax_max-ax_min),[TXT_1,TXT_2],'FontName', 'Courier','FontSize',10,'FontWeight','bold','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1],'Color',[0 0 0]);
text((ax_max-250),(ax_min+50),[TXT_1,TXT_2],'FontName', 'Courier','FontSize',10,'FontWeight','bold','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1],'Color',[0 0 0]);

end

