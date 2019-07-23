% First and last point must be the trailing edge tip with clock wise point distribution 
function [xOut,yOut] = roundTE(x,y,r, beta)

  % Extract points at pressure and suction side  (excluding airfoil tip)
  [_,iLE] = min(x);
  suctionSideX = x(iLE:end-1); 
  suctionSideY = y(iLE:end-1);
  pressureSideX = x(1:iLE-1);
  pressureSideY = y(1:iLE-1);
  
  % Trailing edge tip 
  Px = pressureSideX(1);
  Py = pressureSideY(1);
  
  % Work only with points close to trailinge edge  
  pressureSideCutY = pressureSideY( pressureSideX>0.95 );
  pressureSideCutX = pressureSideX( pressureSideX>0.95 );
  suctionSideCutY = suctionSideY ( suctionSideX>0.95 );
  suctionSideCutX = suctionSideX ( suctionSideX>0.95 );   

  % Find initial values that are at distant by at least 2r 
  % Start from trailing edge and sweep airfoil until d is approximattely 2r 
  d = 0;
  xx = Px;
  dx = 0.00001;
  while (d < 2*r)
    xx = xx - dx;
        
    % B is the point at the pressure side 
    Bx = xx;
    By = spline(pressureSideCutX, pressureSideCutY, xx);

    % A is the point at the suction 
    Ax = xx;
    Ay = spline(suctionSideCutX, suctionSideCutY, xx);
 
    d = sqrt((Bx-Ax)^2 + (By-Ay)^2);
        
  end

  % Throw away points to the right of xx on pressure and suction sides 
  pressureSideCutY = [ By; pressureSideCutY( pressureSideCutX < xx)];
  pressureSideCutX = [ xx; pressureSideCutX( pressureSideCutX < xx)];  
  suctionSideCutY = [ suctionSideCutY( suctionSideCutX < xx); Ay];  
  suctionSideCutX = [ suctionSideCutX( suctionSideCutX < xx); xx];
 
  % Before continuing, bring trailing edge tip to the left
  % M is the point between first pressure and last suction sides points
  Mx = ( pressureSideCutX(1) + suctionSideCutX(end) ) / 2;
  My = ( pressureSideCutY(1) + suctionSideCutY(end) ) / 2;
    
  d = sqrt((Mx-Px)^2 + (My-Py)^2);
  PNewX = beta*r / d * ( Px-Mx ) + Mx;
  PNewY = beta*r / d * ( Py-My ) + My;
   
  pressureSideCutX = [PNewX; pressureSideCutX];
  pressureSideCutY = [PNewY; pressureSideCutY];
 
  % Calculate s from pressure side
  s = zeros(length(pressureSideCutX)+length(suctionSideCutX),1);
  counter = 2;
  for i = length(pressureSideCutX):-1:2
    ds = sqrt( (pressureSideCutX(i)-pressureSideCutX(i-1))^2 + (pressureSideCutY(i)-pressureSideCutY(i-1))^2 );
    s(counter) = s(counter-1) + ds;
    counter = counter + 1;
  end 
  
  sArc0 = s(counter-1);
  ds = sqrt( (pressureSideCutX(1)-suctionSideCutX(end))^2 + (pressureSideCutY(i)-suctionSideCutY(end))^2 );
  s(counter) = s(counter-1) + ds;
  sArcN = s(counter);
  counter = counter + 1;
  
  % Calculate s from suction side
  for i = length(suctionSideCutX):-1:2
    ds = sqrt( (suctionSideCutX(i)-suctionSideCutX(i-1))^2 + (suctionSideCutY(i)-suctionSideCutY(i-1))^2 );
    s(counter) = s(counter-1) + ds;
    counter = counter + 1;
  end 
  
  % Interpolate points from suction and pressure sides
  sInterp = linspace(s(1), s(end), 3000);   
  xx = [pressureSideCutX(end:-1:1); suctionSideCutX(end:-1:1) ];
  yy = [pressureSideCutY(end:-1:1); suctionSideCutY(end:-1:1) ];
  xInterp = pchip(s, xx, sInterp); 
  yInterp = pchip(s, yy, sInterp);   
  
  % Append interpolated data to geometry 
  [_,iTE] = max(xInterp); 
  
  xOut = [ xInterp(iTE:-1:1)'; ...
    pressureSideX( pressureSideX <= 0.95 ); ...
    suctionSideX( suctionSideX <= 0.95); ...
    xInterp(end:-1:iTE)' ];
    
   yOut = [ yInterp(iTE:-1:1)'; ...
    pressureSideY( pressureSideX <= 0.95 ); ...
    suctionSideY( suctionSideX <= 0.95); ...
    yInterp(end:-1:iTE)' ];
  
 

end
