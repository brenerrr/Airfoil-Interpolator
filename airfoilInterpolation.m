close all
clear all
clc

% ******************************************************************************
% ******************************************************************************
% ******************************************************************************

%% User inputs  
filename = 'sd7003.dat'; 
isClockWise = false;         % Points distribution direction 
startsFromTE =  true;        % true if coordinates start from trailing edge               
nOut = 441;                  % Final number of points
percentageSuctionSide = 0.7; % Percentage of points in the suction side surface;
dsMinLE = 0.0010;            % Minimum distance between points at leading edge
dsMinTE = 0.00001;           % Minimum distance between points at trailing edge
xStartCoarse = 0.99;         % Where to start distribution coarsening 
xEndCoarse = 0.00001;        % Where to end distribution coarsening
roundFlag = true;            % true if one desires to round trailing edge
rTE = 0.0008;                % Trailing edge radius (in chord's length) - only required when rounding 
tip = 1;                     % Higher values stretch trailing edge making it become a tip. Default = 1
openTE = false;              % true for airfoils without closed trailing edges

% Increase these coefficients to concentrate more points in that region (default=1)
alphaTE = 2;
alphaLE = 2;
alphaPressureSide = 1;
alphaSuctionSide = 1;

% ******************************************************************************
% ******************************************************************************
% ******************************************************************************

%% Pre processing
nCoarse = round((1-percentageSuctionSide) * nOut);
coordCoarse = load(filename); 

% Coordinates need to be in clock wise orientation
if isClockWise == false
    coordCoarse(:,1) = coordCoarse(end:-1:1,1);
    coordCoarse(:,2) = coordCoarse(end:-1:1,2);
end

% Coordinates need to start from trailing edge
if (startsFromTE==false)    
    
    % Rearrange points
    [~,i] = max(coordCoarse(:,1));    
    if (openTE==true)  
        i = i + 1;
        coordCoarse = [coordCoarse(i:end,:); coordCoarse(2:i-1,:)];
    else
        coordCoarse = [coordCoarse(i:end,:); coordCoarse(2:i,:)];
    end        
end

% Close open trailing edges
if (openTE)
    % Find x location of TE point
    x = coordCoarse(:,1);
    y = coordCoarse(:,2);    
    tanAlpha = (y(1)-y(2))/(x(1)-x(2));
    xTE = (tanAlpha*x(1)-y(1))/tanAlpha;
    yTE = 0;    
    
    % Increase size of coordinates dimensions
    coordCoarse = [[0,0] ; coordCoarse; [0,0] ];
    coordCoarse(:,1) = [ xTE; x; xTE ];
    coordCoarse(:,2) = [ yTE; y; yTE ];   
        
end

% Make sure first and last points are the same 
if abs(coordCoarse(end,2)-coordCoarse(1,2))>1e-9
    coordCoarse = [coordCoarse; coordCoarse(1,:)];
end

% Make sure xstart and xend coarse are in the right format 
if xStartCoarse < xEndCoarse
  temp = xEndCoarse;
  xEndCoarse = xStartCoarse;
  xStartCoarse = temp;
end


% Normalize points so chord is unitary 
[~,iLE] = min(coordCoarse(:,1));
coordCoarse(:,1) = coordCoarse(:,1) - coordCoarse(iLE,1);
coordCoarse = coordCoarse / coordCoarse(1,1);

% Rounding
if (roundFlag == true)
  
  nInt = 20001; % Number of points when interpolating airfoil surface 
  
  % Add more points to surface to facilitate rounding later 
  [xOut, yOut, ~, ~] = interpolateSurface ( ...
                                coordCoarse(:,1), coordCoarse(:,2), nInt, 0.000025, ...
                                0.000025, round(nInt*0.5), 0.99, xEndCoarse, alphaTE,...
                                alphaLE, alphaPressureSide, alphaSuctionSide );
                        
  % Do it again to make sure everything is smooth                              
  [xOut, yOut, ~, ~] = interpolateSurface ( ...
                              xOut, yOut, nInt, 0.000025, ...
                              0.000025, round(nInt*0.5), 0.99, xEndCoarse, alphaTE,...
                              alphaLE, alphaPressureSide, alphaSuctionSide );
                  
                                      
  % Round trailing edge 
  [xOut,yOut] = roundTE(xOut, yOut, rTE, tip);
%  [xOut,yOut] = roundTE(coordCoarse(:,1), coordCoarse(:,2), rTE, tip);
end

% Normalize 
[~,iLE] = min(xOut(:,1));
xOut = xOut - xOut(iLE);
xOut = xOut / xOut(1);
yOut = yOut / xOut(1);

% Distribute points through surface 
[xOut, yOut, ~, ~] = interpolateSurface ( ...
                                xOut, yOut, nOut, dsMinTE, ...
                                dsMinLE, nCoarse, xStartCoarse, xEndCoarse, alphaTE,...
                                alphaLE, alphaPressureSide, alphaSuctionSide );

% Do it again to make sure everything is smooth
[xOut, yOut, sOut, dsOut] = interpolateSurface ( ...
                                xOut, yOut, nOut, dsMinTE, ...
                                dsMinLE, nCoarse, xStartCoarse, xEndCoarse, alphaTE,...
                                alphaLE, alphaPressureSide, alphaSuctionSide );                                

% Normalize 
[~,iLE] = min(xOut(:,1));
xOut = xOut - xOut(iLE);
xOut = xOut / xOut(1);
yOut = yOut / xOut(1);

%% Plot
% Airfoil 
figure(1);
hold on 
plot(coordCoarse(:,1),coordCoarse(:,2),'bo-','MarkerSize',10)   
plot(xOut,yOut,'ro--','MarkerFaceColor','r','MarkerSize',2)    
axis([0, 1, -0.5, 0.5])   
pbaspect([1 1 1])
grid on 
legend('Original data', 'Interpolated data') 
%set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'Color','white')
set(gcf, 'Position',  [100, 100, 800, 800])
%export_fig airfoil.png 

% X distribution
figure(2)
plot(sOut,xOut,'-o')
hold on
xlabel('s','FontSize',20)
ylabel('x','FontSize',20)
%set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'Color','white')
set(gcf, 'Position',  [100, 100, 800, 800])
%export_fig x_distribution.png

% Y distribution
figure(3)
plot(sOut,yOut,'-o')
hold on
xlabel('s')
ylabel('y')
%set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'Color','white')
set(gcf, 'Position',  [100, 100, 800, 800])
%export_fig y_distribution.png

figure(4)
plot(sOut(2:end),dsOut,'-o') 
hold on 
xlabel('s')
ylabel('ds')
%set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'Color','white')
set(gcf, 'Position',  [100, 100, 800, 800])
%export_fig ds.png

figure(5)
plot(sOut,'-o') 
hold on 
xlabel('Index','fontsize',18)
ylabel('s')
%set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'Color','white')
set(gcf, 'Position',  [100, 100, 800, 800])
%export_fig s.png


% Make output array start from leading edge and follow clock wise
% distribution
[~, iTemp] = min(xOut);
xOut = [ xOut(iTemp:end); xOut(2:iTemp)];
yOut = [ yOut(iTemp:end); yOut(2:iTemp)];
xOut = fliplr(xOut);
yOut = fliplr(yOut);

% Create output array
nOut = length(xOut);
output = zeros(nOut,2);
output(:,1)=xOut;
output(:,2)=yOut;

% Output chord length (unitary)
c = max(xOut)-min(xOut);

save airfoilRefined.dat nOut -ascii
save airfoilRefined.dat output -ascii -double -append

fprintf(['\n Output file: airfoilRefined.dat \n ', ...
         num2str(nOut), ' points in new geometry \n' ...
         ' Output s chord length: ', num2str(c), '\n']);

         
