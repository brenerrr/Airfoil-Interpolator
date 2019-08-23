% Points distribution must be clockwise
function [xOut,yOut, sOut, dsOut] = interpolateSurface ( ...
                                x, y, nTotal, dsMinTE, dsMinLE, nCoarse, ...
                                xStartCoarse, xEndCoarse, alphaTE, alphaLE,...
                                alphaPressure, alphaSuction )

  nRefined = nTotal - nCoarse;
  
  %% Create a distribution that goes from xStartCoarse to xEndCoarse with 
  %% minimum dsMinTE, sMinLE and nTotal points with nCoarse points in the coarse
  %% region
  dsCoarseSum = GetSizeOfArc(x, y, xStartCoarse, xEndCoarse);
    
  % Find ds of coarse part of airfoil using newthon's method to change dsMax  
  i = 0;
  dsMax = 0.1;
  while true 
    i = i + 1;
    if mod(nCoarse-1,2) ~= 0 
      ds1 = tanhDistribution(dsMinTE, dsMax, (nCoarse-1)/2+1, alphaTE, alphaPressure);  
    else
      ds1 = tanhDistribution(dsMinTE, dsMax, (nCoarse-1)/2, alphaTE, alphaPressure);
    end
    ds2 = tanhDistribution(dsMinLE, dsMax, (nCoarse-1)/2, alphaLE, alphaPressure);     
    ds = [ds1, ds2(end:-1:1)];

    % In case this is the first iteration, the next dsMax needs to be guessed 
    if i == 1      
      previousDsMax = dsMax; 
      previousF = sum(ds) - dsCoarseSum;   
      dsMax = dsMax*0.999; 
      if mod((nCoarse-1),2) ~= 0 
        ds1 = tanhDistribution(dsMinTE, dsMax, (nCoarse-1)/2+1, alphaTE, alphaPressure);  
      else
        ds1 = tanhDistribution(dsMinTE, dsMax, (nCoarse-1)/2, alphaTE, alphaPressure);
      end
      ds2 = tanhDistribution(dsMinLE, dsMax, (nCoarse-1)/2, alphaLE, alphaPressure);     
      ds = [ds1, ds2(end:-1:1)];         
    end

    
    % Perform step with newthon's method 
    f = sum(ds) - dsCoarseSum;  
    df = (f - previousF) / (dsMax - previousDsMax);    
    previousDsMax = dsMax;
    dsMax = dsMax - f/df;
    previousF = f; 

    if abs(dsMax - previousDsMax) < 1e-6 
      break
    end
    
    if isnan(dsMax)
      error ('dsMax = NaN')
    end

  end
  
  

  if dsMax < 0 
    error ('dsMax < 0 - Try different dsMinLE or dsMinTE')
  end
  dsCoarse = ds;
  
  
  % Do the same but know for suction side of airfoil
  % But first, find total arc of airfoil surface 
  dsTotal = 0;
  for i=1:length(x)-1
    dsTotal = dsTotal + sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 ) ;
  end
  dsOutSum = dsTotal - dsCoarseSum;
  dsMax = 0.1;
  while true 
    i = i + 1;
        
    % Here number of poins is not nRefined-1 on purpose to make final size be nPoints-1
    if mod(nRefined,2) ~= 0 
      ds1 = tanhDistribution(dsMinLE, dsMax, (nRefined)/2+1, alphaLE, alphaSuction);    
    else
      ds1 = tanhDistribution(dsMinLE, dsMax, (nRefined)/2, alphaLE, alphaSuction);   
    end    
    ds2 = tanhDistribution(dsMinTE, dsMax, (nRefined)/2, alphaTE, alphaSuction);     
    ds = [ds1, ds2(end:-1:1)];
    
    % In case this is the first iteration, the next dsMax needs to be guessed 
    if i == 1      
      previousDsMax = dsMax; 
      previousF = sum(ds) - dsOutSum;   
      dsMax = dsMax*0.999; 
      if mod((nRefined),2) ~= 0 
        ds1 = tanhDistribution(dsMinLE, dsMax, (nRefined)/2+1, alphaLE, alphaSuction);    
      else
        ds1 = tanhDistribution(dsMinLE, dsMax, (nRefined)/2, alphaLE, alphaSuction);   
      end    
      ds2 = tanhDistribution(dsMinTE, dsMax, (nRefined)/2, alphaTE, alphaSuction);     
      ds = [ds1, ds2(end:-1:1)];         
    end
    
    % Newthon's method 
    f = sum(ds) - dsOutSum;
    df = (f - previousF) / (dsMax - previousDsMax);
    previousDsMax = dsMax;
    dsMax = dsMax - f/df;
    previousF = f; 

    if abs(dsMax - previousDsMax) < 1e-6 
      break
    end 
    
    if isnan(dsMax)
      error ('dsMax = NaN')
    endif

  end
  
  if dsMax < 0 
    error ('dsMax < 0 - Try different dsMinLE or dsMinTE')
  end
  
  % Calculate refined s from ds 
  dsOut = ds; 
  ds = [dsCoarse, dsOut];  
  sOut = zeros(nTotal,1);
  for i = 2:nTotal
    sOut(i) = sum(ds(1:i-1));
  end

  % Calculate coarse s from inputs x and y   
  s = zeros(length(x),1);
  ds = zeros(length(s)-1,1);
  for i = 2:length(x)
    ds(i-1) = sqrt( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 );
    s(i) = s(i-1) + ds(i-1);
  end 
  
  % Calculate coordinates from refined s
  xOut = pchip(s, x, sOut);
  yOut = pchip(s, y, sOut);
  
  dsOut=zeros(length(xOut)-1,1);
  for i = 1:length(xOut)-1
    dsOut(i) = sqrt( (xOut(i) - xOut(i+1))^2 + (yOut(i) - yOut(i+1))^2 );
  end       
  
%  min(dummy)
%    figure(1)
%  plot(x,y,'o-') 
%  axis equal
%  figure(2)
%  plot(ds,'o-') 
%  figure(3)
%  plot(s,'o-')
% 
%  figure(10)
%  plot(sOut(2:end),dummy,'o-') 
%%  figure(11)
%%%  plot(xOut, yOut,'o-')
%%  axis equal
%  hold on 
%  plot(s,x,'ro-')
%  plot(sOut,xOut)
%%    figure(12)
%  plot(s,x,'o-')
%  

 
end

% Bigger alpha means more points close to minValue, while bigger betas means 
% more points closer to maxValue
function [y] = tanhDistribution(minValue, maxValue, nPoints, alpha, beta)
  x = linspace(-3*alpha, 3*beta,nPoints);
  y = (tanh(x) + 1) / 2 * (maxValue - minValue) + minValue ;

end

% Get size of arc starting at xStart to xEnd (arc is admitted to be in clock wise)
function [arcSize] = GetSizeOfArc(x, y, xStart, xEnd)

  [~,iLE] = min(x);
  n = length(x);
  arcSize = 0;
  shouldStartAdding = false;
  for i = 2:n
  
    % Loop until xStart is surpassed and then start storing arc segments sum 
    if x(i) < xStart && shouldStartAdding == false   
    
      % Find yStart through weighted average
      yStart = (xStart - x(i-1)) / (x(i)-x(i-1)) * (y(i)-y(i-1)) + y(i-1);
      arcSize = sqrt( (x(i)-xStart)^2 + (y(i)-yStart)^2 );
      shouldStartAdding = true;

    % Loop until xEnd or leading edge then get out of loop
    elseif x(i) < xEnd || i == iLE
      yEnd = (xEnd - x(i-1)) / (x(i)-x(i-1)) * (y(i)-y(i-1)) + y(i-1);
      arcSize = arcSize + sqrt( (xEnd-x(i-1))^2 + (yEnd-y(i-1))^2 );
      break
    
    % Summing arcs until xEnd is reached
    elseif shouldStartAdding
      arcSize = arcSize + sqrt( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 );
    end
    
  end

end