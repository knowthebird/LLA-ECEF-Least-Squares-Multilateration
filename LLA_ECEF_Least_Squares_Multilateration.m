 format long g;

%if station locations LLD format, WGS84, alt in meters
station = DataPoints(:,1:3);
numToSub = 0;
%subtract for LLD to LLA
%see http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/intptW.html
station(:,3) = station(:,3) - numToSub;

range = DataPoints(:,4);
ECEFstation = lla2ecef(station, 'WGS84');

%Station locations ECEF format, WGS84, range in meters
ECEFguess = MinEst(ECEFstation,range,300,.1,0);

LLA = ecef2lla(ECEFguess, 'WGS84')
ECEFguess(:,3) = ECEFguess(:,3) + numToSub;

function minxHat = MinEst(RefLocations, Ranges, Maxiter, TolX, DispOpt)
   
    
    if  size(RefLocations,1) == size(Ranges,1)
        if size(RefLocations,2) ~= 3
            disp('Locations should be [x 3] matrix, for x,y,z values.');
            return;
        elseif size(Ranges,2) ~= 1
            disp('Ranges should be [x 1] matrix, for range values.');
            return;
        else
            clear dimLoc dimRang;
            NumOfStations = size(Ranges,1);
        end
    else
        disp('Number of locations and ranges is not the same.');
        return;
    end
    H = zeros(NumOfStations,4);
    H(:,1) = 1;
    H(:,2) = -2*RefLocations(:,1);
    H(:,3) = -2*RefLocations(:,2);
    H(:,4) = -2*RefLocations(:,3);
    
    b = zeros(NumOfStations,1);
    b(:,1) = Ranges(:).^2 - RefLocations(:,1).^2 - RefLocations(:,2).^2 - RefLocations(:,3).^2;
    
    xHat = (transpose(H)*H)\(transpose(H)*b); 
    Est = [xHat(2,1) xHat(3,1) xHat(4,1)];

    if DispOpt
        options = optimset('PlotFcns',@optimplotfval,'MaxIter',Maxiter,'TolX',TolX);
    else
        options = optimset('MaxIter',Maxiter,'TolX',TolX);
    end
    
    lb = [min(RefLocations(:,1)),min(RefLocations(:,2)),min(RefLocations(:,3))];
    ub = [max(RefLocations(:,1)),max(RefLocations(:,2)),max(RefLocations(:,3))];

    A = []; b = []; Aeq = []; beq = [];
    
    %minxHat = fminsearch(@three_var,Est,options);
    minxHat = fmincon(@three_var,Est,A,b,Aeq,beq,lb,ub);
    function e = three_var(Est)
        e=0;
        for i=1:size(Ranges,1);
            e = e + (sqrt((RefLocations(i,1)-Est(1)).^2 + (RefLocations(i,2)-Est(2)).^2 + (RefLocations(i,3)-Est(3)).^2) - Ranges(i))^2;
        end
    end
end
