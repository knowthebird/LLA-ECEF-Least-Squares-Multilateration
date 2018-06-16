clc; clearvars; format compact; format long g;

n = 10; % generate n ground stations, minimum is 4
%n = input('Number of Ground Stations?')

MaxD = 1000; % max distance for ground station from origin
%MaxD = input('Maximum distance of Ground Station from Origin?')

UserLoc = randi(MaxD, 1, 3); % Generate a random location for Rover
GSLoc = randi(MaxD, n, 3); % Generate Ground Station Locations
% determin distance from ground stations to rover
GSLocSpheres = [GSLoc sqrt((UserLoc(1)-GSLoc(:,1)).^2 + (UserLoc(2)-GSLoc(:,2)).^2 + (UserLoc(3)-GSLoc(:,3)).^2)];

%enter an error multiplier as a percent
NoiseMultiplier = 0;
NoiseMultiplier = (NoiseMultiplier / 100) * MaxD;
GSLocSpheres(:,4) =  GSLocSpheres(:,4) + randn(n,1)*NoiseMultiplier;


% to display spheres with radius and intersections uncomment code....

% figure
% [x y z] = sphere;
% hold on
% for c = 1:n
%     s(c) = surf(x*GSLocSpheres(c,4)+GSLocSpheres(c,1),y*GSLocSpheres(c,4)+GSLocSpheres(c,2),z*GSLocSpheres(c,4)+GSLocSpheres(c,3));
% end
% daspect([1 1 1])
% scatter3(UserLoc(1),UserLoc(2),UserLoc(3), 'x')
% view(30,10)
% grid on
% hold off
% shg

intersections = zeros(n-3,3);


%figure
for i= 1:n-3
    DistR = GSLocSpheres(i,4);
    Distr = GSLocSpheres(i+1,4);
    DistC = sqrt((GSLocSpheres(i+1,1) - GSLocSpheres(i,1))^2+(GSLocSpheres(i+1,2) - GSLocSpheres(i,2))^2+(GSLocSpheres(i+1,3) - GSLocSpheres(i,3))^2);
    x = (DistC^2 - Distr^2 + DistR^2) / (2 * DistC);
    CircleCenter1 = [((x/DistC)*(GSLocSpheres(i+1,1) - GSLocSpheres(i,1)))+GSLocSpheres(i,1),((x/DistC)*(GSLocSpheres(i+1,2) - GSLocSpheres(i,2)))+GSLocSpheres(i,2),((x/DistC)*(GSLocSpheres(i+1,3) - GSLocSpheres(i,3)))+GSLocSpheres(i,3)];
    CircleRadius1 = (1/(2*DistC))*sqrt(4*DistC^2*DistR^2-(DistC^2-Distr^2+DistR^2)^2);
    CircleNormalVector1 = [GSLocSpheres(i+1,1) - GSLocSpheres(i,1),GSLocSpheres(i+1,2) - GSLocSpheres(i,2),GSLocSpheres(i+1,3) - GSLocSpheres(i,3)];
    
    
    theta=0:0.01:2*pi;
    v=null(CircleNormalVector1);
    points1=repmat(CircleCenter1',1,size(theta,2))+CircleRadius1*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    %plot3(points1(1,:),points1(2,:),points1(3,:),'r-');
    %hold on

    DistR = GSLocSpheres(i,4);
    Distr = GSLocSpheres(i+2,4);
    DistC = sqrt((GSLocSpheres(i+2,1) - GSLocSpheres(i,1))^2+(GSLocSpheres(i+2,2) - GSLocSpheres(i,2))^2+(GSLocSpheres(i+2,3) - GSLocSpheres(i,3))^2);
    x = (DistC^2 - Distr^2 + DistR^2) / (2 * DistC);
    CircleCenter2 = [((x/DistC)*(GSLocSpheres(i+2,1) - GSLocSpheres(i,1)))+GSLocSpheres(i,1),((x/DistC)*(GSLocSpheres(i+2,2) - GSLocSpheres(i,2)))+GSLocSpheres(i,2),((x/DistC)*(GSLocSpheres(i+2,3) - GSLocSpheres(i,3)))+GSLocSpheres(i,3)];
    CircleRadius2 = (1/(2*DistC))*sqrt(4*DistC^2*DistR^2-(DistC^2-Distr^2+DistR^2)^2);
    CircleNormalVector2 = [GSLocSpheres(i+2,1) - GSLocSpheres(i,1),GSLocSpheres(i+2,2) - GSLocSpheres(i,2),GSLocSpheres(i+2,3) - GSLocSpheres(i,3)];

    v=null(CircleNormalVector2);
    points2=repmat(CircleCenter2',1,size(theta,2))+CircleRadius2*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    %plot3(points2(1,:),points2(2,:),points2(3,:),'r-');

    DistR = GSLocSpheres(i,4);
    Distr = GSLocSpheres(i+3,4);
    DistC = sqrt((GSLocSpheres(i+3,1) - GSLocSpheres(i,1))^2+(GSLocSpheres(i+3,2) - GSLocSpheres(i,2))^2+(GSLocSpheres(i+3,3) - GSLocSpheres(i,3))^2);
    x = (DistC^2 - Distr^2 + DistR^2) / (2 * DistC);
    CircleCenter2 = [((x/DistC)*(GSLocSpheres(i+3,1) - GSLocSpheres(i,1)))+GSLocSpheres(i,1),((x/DistC)*(GSLocSpheres(i+3,2) - GSLocSpheres(i,2)))+GSLocSpheres(i,2),((x/DistC)*(GSLocSpheres(i+3,3) - GSLocSpheres(i,3)))+GSLocSpheres(i,3)];
    CircleRadius2 = (1/(2*DistC))*sqrt(4*DistC^2*DistR^2-(DistC^2-Distr^2+DistR^2)^2);
    CircleNormalVector2 = [GSLocSpheres(i+3,1) - GSLocSpheres(i,1),GSLocSpheres(i+3,2) - GSLocSpheres(i,2),GSLocSpheres(i+3,3) - GSLocSpheres(i,3)];

    v=null(CircleNormalVector2);
    points3=repmat(CircleCenter2',1,size(theta,2))+CircleRadius2*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    %plot3(points3(1,:),points3(2,:),points3(3,:),'r-');

    
    
    
    %scatter3(UserLoc(1),UserLoc(2),UserLoc(3), 'x')
%     grid on
%     hold off
    %find minimum distance between points to determine intersections
    distance = 1;
    for c = 1:size(points1,2)
        for d = 1:size(points2,2)
            if (d > 1) | (c > 1)
                lastdistance = distance;
                distance = norm(sqrt( (points1(1,c)-points2(1,d))^2 + (points1(2,c)-points2(2,d))^2 + (points1(3,c)-points2(3,d))^2));
                if distance < MinDist1(1,3);
                    MinDist1 = [c, d, distance];
                end

            else
                distance = norm(sqrt( (points1(1,c)-points2(1,d))^2 + (points1(2,c)-points2(2,d))^2 + (points1(3,c)-points2(3,d))^2));
                MinDist1 = [c, d, distance];
            end
        end
    end

    intersection1 = [ (points1(1,MinDist1(1,1)) + points2(1,MinDist1(1,2)))/2, (points1(2,MinDist1(1,1)) + points2(2,MinDist1(1,2)))/2, (points1(3,MinDist1(1,1)) + points2(3,MinDist1(1,2)))/2];

    for c = 1:size(points1,2)
        for d = 1:size(points2,2)
            if (d > 1) | (c > 1)
                lastdistance = distance;
                distance = norm(sqrt( (points1(1,c)-points2(1,d))^2 + (points1(2,c)-points2(2,d))^2 + (points1(3,c)-points2(3,d))^2));
                if distance < MinDist2(1,3);
                    MinDist2Temp = [c, d, distance];
                    intersection2 = [ (points1(1,MinDist2Temp(1,1)) + points2(1,MinDist2Temp(1,2)))/2, (points1(2,MinDist2Temp(1,1)) + points2(2,MinDist2Temp(1,2)))/2, (points1(3,MinDist2Temp(1,1)) + points2(3,MinDist2Temp(1,2)))/2];

                    if intersection1 ~= intersection2
                        if (c ~= size(points1,2))&& (c~= 1) && (d~= 1)&&(d ~= size(points2,2))
                            highdistance1 = norm(sqrt( (points1(1,c+1)-points2(1,d))^2 + (points1(2,c+1)-points2(2,d))^2 + (points1(3,c+1)-points2(3,d))^2));
                            lowdistance1 = norm(sqrt( (points1(1,c-1)-points2(1,d))^2 + (points1(2,c-1)-points2(2,d))^2 + (points1(3,c-1)-points2(3,d))^2));
                            highdistance2 = norm(sqrt( (points1(1,c)-points2(1,d+1))^2 + (points1(2,c)-points2(2,d+1))^2 + (points1(3,c)-points2(3,d+1))^2));
                            lowdistance2 = norm(sqrt( (points1(1,c)-points2(1,d-1))^2 + (points1(2,c)-points2(2,d-1))^2 + (points1(3,c)-points2(3,d-1))^2));
                            if (distance < highdistance1) && (distance < lowdistance1) && (distance < highdistance2) && (distance < lowdistance2)
                                MinDist2 = [c, d, distance];
                            end

                        elseif (c == size(points1,2))
                            highdistance1 = norm(sqrt( (points1(1,1)-points2(1,d))^2 + (points1(2,1)-points2(2,d))^2 + (points1(3,1)-points2(3,d))^2));
                            lowdistance1 = norm(sqrt( (points1(1,c-1)-points2(1,d))^2 + (points1(2,c-1)-points2(2,d))^2 + (points1(3,c-1)-points2(3,d))^2));
                            if (d == 1)
                                highdistance2 = norm(sqrt( (points1(1,c)-points2(1,d+1))^2 + (points1(2,c)-points2(2,d+1))^2 + (points1(3,c)-points2(3,d+1))^2));
                                lowdistance2 = norm(sqrt( (points1(1,c)-points2(1,size(points2,2)))^2 + (points1(2,c)-points2(2,size(points2,2)))^2 + (points1(3,c)-points2(3,size(points2,2)))^2));
                            elseif (d == size(points2,2))
                                highdistance2 = norm(sqrt( (points1(1,c)-points2(1,1))^2 + (points1(2,c)-points2(2,1))^2 + (points1(3,c)-points2(3,1))^2));
                                lowdistance2 = norm(sqrt( (points1(1,c)-points2(1,d-1))^2 + (points1(2,c)-points2(2,d-1))^2 + (points1(3,c)-points2(3,d-1))^2));
                            else
                                highdistance2 = norm(sqrt( (points1(1,c)-points2(1,d+1))^2 + (points1(2,c)-points2(2,d+1))^2 + (points1(3,c)-points2(3,d+1))^2));
                                lowdistance2 = norm(sqrt( (points1(1,c)-points2(1,d-1))^2 + (points1(2,c)-points2(2,d-1))^2 + (points1(3,c)-points2(3,d-1))^2));
                            end
                            
                            if (distance < highdistance1) && (distance < lowdistance1) && (distance < highdistance2) && (distance < lowdistance2)
                                MinDist2 = [c, d, distance];
                            end
                        else
                            highdistance1 = norm(sqrt( (points1(1,c+1)-points2(1,d))^2 + (points1(2,c+1)-points2(2,d))^2 + (points1(3,c+1)-points2(3,d))^2)); 
                            if c == 1
                                lowdistance1 = norm(sqrt( (points1(1,size(points1,2))-points2(1,d))^2 + (points1(2,size(points1,2))-points2(2,d))^2 + (points1(3,size(points1,2))-points2(3,d))^2));
                            else
                                lowdistance1 = norm(sqrt( (points1(1,c-1)-points2(1,d))^2 + (points1(2,c-1)-points2(2,d))^2 + (points1(3,c-1)-points2(3,d))^2));
                            end
                            if (d == 1)
                                highdistance2 = norm(sqrt( (points1(1,c)-points2(1,d+1))^2 + (points1(2,c)-points2(2,d+1))^2 + (points1(3,c)-points2(3,d+1))^2));
                                lowdistance2 = norm(sqrt( (points1(1,c)-points2(1,size(points2,2)))^2 + (points1(2,c)-points2(2,size(points2,2)))^2 + (points1(3,c)-points2(3,size(points2,2)))^2));
                            elseif (d == size(points2,2))
                                highdistance2 = norm(sqrt( (points1(1,c)-points2(1,1))^2 + (points1(2,c)-points2(2,1))^2 + (points1(3,c)-points2(3,1))^2));
                                lowdistance2 = norm(sqrt( (points1(1,c)-points2(1,d-1))^2 + (points1(2,c)-points2(2,d-1))^2 + (points1(3,c)-points2(3,d-1))^2));
                            else
                                highdistance2 = norm(sqrt( (points1(1,c)-points2(1,d+1))^2 + (points1(2,c)-points2(2,d+1))^2 + (points1(3,c)-points2(3,d+1))^2));
                                lowdistance2 = norm(sqrt( (points1(1,c)-points2(1,d-1))^2 + (points1(2,c)-points2(2,d-1))^2 + (points1(3,c)-points2(3,d-1))^2));
                            end
                            
                            if (distance < highdistance1) && (distance < lowdistance1) && (distance < highdistance2) && (distance < lowdistance2)
                                MinDist2 = [c, d, distance];
                            end
                        end
                    end

                end

            else
                distance = norm(sqrt( (points1(1,c)-points2(1,d))^2 + (points1(2,c)-points2(2,d))^2 + (points1(3,c)-points2(3,d))^2));
                MinDist2 = [c, d, distance];
            end
        end
    end

    intersection2 = [ (points1(1,MinDist2(1,1)) + points2(1,MinDist2(1,2)))/2, (points1(2,MinDist2(1,1)) + points2(2,MinDist2(1,2)))/2, (points1(3,MinDist2(1,1)) + points2(3,MinDist2(1,2)))/2];

    %find which intersection is closest to third ring
    for c = 1:size(points3,2)
        if c > 1
            Intsect1DistThirdRing = norm(sqrt( (points3(1,c)-intersection1(1,1))^2 + (points3(2,c)-intersection1(1,2))^2 + (points3(3,c)-intersection1(1,3))^2));
            if Intsect1DistThirdRing < LowestOneDist
                LowestOneDist = Intsect1DistThirdRing;
            end
        else
            LowestOneDist = norm(sqrt( (points3(1,c)-intersection1(1,1))^2 + (points3(2,c)-intersection1(1,2))^2 + (points3(3,c)-intersection1(1,3))^2));
        end
    end


    %find which intersection is closest to third ring
    for c = 1:size(points3,2)
        if c > 1
            Intsect2DistThirdRing = norm(sqrt( (points3(1,c)-intersection2(1,1))^2 + (points3(2,c)-intersection2(1,2))^2 + (points3(3,c)-intersection2(1,3))^2));
            if Intsect2DistThirdRing < LowestTwoDist
                LowestTwoDist = Intsect2DistThirdRing;
            end
        else
            LowestTwoDist = norm(sqrt( (points3(1,c)-intersection2(1,1))^2 + (points3(2,c)-intersection2(1,2))^2 + (points3(3,c)-intersection2(1,3))^2));
        end
    end


    if LowestOneDist < LowestTwoDist
        intersections(i,:) = intersection1;
    else
        intersections(i,:) = intersection2;
    end
    
    
end



intersections = round(abs(intersections))
%MeanIntersection = mean(intersections);
[idx,C] = kmeans(intersections,1);
KMeanIntersection = C
TrueIntersection = UserLoc
PercentError = 100*(abs(KMeanIntersection - TrueIntersection)/TrueIntersection)




% %plot spheres from locations and radius
% figure
% scatter3(GSLocSpheres(:,1),GSLocSpheres(:,2),GSLocSpheres(:,3))
% hold on
% scatter3(UserLoc(1),UserLoc(2),UserLoc(3), 'x')
% scatter3(((x/DistC)*(GSLocSpheres(2,1) - GSLocSpheres(i,1)))+GSLocSpheres(i,1),((x/DistC)*(GSLocSpheres(2,2) - GSLocSpheres(i,2)))+GSLocSpheres(i,2),((x/DistC)*(GSLocSpheres(2,3) - GSLocSpheres(i,3)))+GSLocSpheres(i,3),'x')
% % r=radius;
% % teta=-pi:0.01:pi;
% % x=r*cos(teta);
% % y=r*sin(teta);
% % plot3(x,y,zeros(1,numel(y)))
% hold off


