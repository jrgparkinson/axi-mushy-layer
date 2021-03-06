function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 11-Apr-2015 14:04:26

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1);
set(plot1(1),'Color',[0 0 1]);
set(plot1(2),'Color',[1 0 0]);

% Create xlabel
xlabel('a');

% Create ylabel
ylabel('q');

% Create title
title('q.\nabla \theta(z=H/2) at r=b (blue) and r=a (red),Da =5e-05, Rm =60, R =0.25, H =0.25b =0.035');

