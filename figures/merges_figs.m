% Load saved figures
c=hgload('newanno.fig');
k=hgload('newbc.fig');
% Prepare subplots
figure
h(1)=subplot(1,2,1);
h(2)=subplot(1,2,2);
% Paste figures on the subplots
copyobj(allchild(get(c,'CurrentAxes')),h(1));
copyobj(allchild(get(k,'CurrentAxes')),h(2));