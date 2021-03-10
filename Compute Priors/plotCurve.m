function plotCurve(X,fig,linew,msize,viewvec)

% Plot a 2D or 3D curve. For 3D plots, show a cubic spline interpolation 
% between points for a smoother plot, and also colorize according to 
% increasing sample point index.
%
% Inputs: 
%   X - An n x p matrix representing an n-tuple in R^p.
%   fig - The figure number
%   linew - The linewidth
%   msize - The marker size. If msize = 0 then do not overlay markers on
%       top of the line.
%   viewvec - [Alt Az] view vector for 3D plot.

[n,p]=size(X);
figure(fig); 

if p==2
    
    if msize>0
        plot(X(:,1),X(:,2),'ko-','LineWidth',linew,'MarkerSize',msize);
    else
        plot(X(:,1),X(:,2),'k','LineWidth',linew);
    end
    axis equal off tight;
    
elseif p==3
    
    N=10*n;
    t1=linspace(0,1,n);
    t2=linspace(0,1,N);
    x=spline(t1,X(:,1)',t2);
    y=spline(t1,X(:,2)',t2);
    z=spline(t1,X(:,3)',t2);
    surface([x;x],[y;y],[z;z],[t2;t2],'facecol','no','edgecol','interp','linew',linew);
    if msize>0
        hold on;
        plot3(X(:,1),X(:,2),X(:,3),'ko','MarkerSize',msize);
        hold off;
    end
    axis equal tight off;
    view(viewvec);
    
end
