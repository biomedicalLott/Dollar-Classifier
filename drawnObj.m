% Written by Robert Lott for Rubine Classification Project
classdef  drawnObj < handle
    %drawnObj Will hold and process values directly related to
    %classification of data samples
    %   Ground truth is assembled from a simple assignment of "groundTruth" which
    %   acts as an ID for the ground truth, it should be between 0 and 9 at
    %   present
    
    properties
        %For Template
        groundTruth     %Classification ID related to what it represents
        ID              % Sample ID of the sample, optional arg
        oldX            %oldX position of points
        oldY            %oldY position of points
        X               %interpolated oldX points
        Y               %interpolated oldY points
        center          %center of adjusted image
        t               %time of when each point ws placed
        theta           %angle of points
        euDist          %euclidian distance along points
        distance        %distance between each point
        cumulativeDist  %cumulative distance between points
        totalDistance   %Total distance traveled
        covariance      %covariance of the position vectors
        covMean         %the mean of a set of observations
    end
    
    methods
        
        function drawnObj = drawnObj(Load,points,groundTruth,ID)
%             and if i did, it'll check to see if i wanted to just make an
%             empty one
            if(exist('Load','var') && Load == 1)
                %                    drawnObj.loadObj(savedObj)
                return
            else
                Load = 0;
            end
            %             Checks to make sure i didn't call an empty shape
            if ~exist('points','var')
                %           If there aren't any points, you did a bad thing
                error('You forgot some points.')
                return
            end
            if exist('groundTruth','var')
                if ~exist('ID','var')
                    %                     If it doesn't exist, just generate a random ID
                    ID = rand()*50204918;
                end
                drawnObj.storeAsTemplate(groundTruth,ID);
                %                 disp(['Storing class ' , num2str(groundTruth)])
            end
            %             step 0 store and remove bad points
            drawnObj.processPoints(points);
            
            %             step 1 - resample points path
            drawnObj.getPointDistance();
            drawnObj.resamplePath(64);
            drawnObj.findCenter();
            %             Step 2 - Rotate to 0 degrees
            drawnObj.rotateToZero();
            %             Step 3 - Scale Points
            scaleSize = [64,64];
            drawnObj.scaleToSquare(scaleSize);
            drawnObj.translateToOrigin();
            %             Step 4 - Match Points, done outside of this.
            
        end
        
        %%
 
% % % % %       START STEP 2
        function rotateToZero(drawnObj)
%               corrects rotation from center point to be 0
            theta = drawnObj.findAngleFromCenter();
            points = drawnObj.rotateBy(-theta);
            %            drawnObj.comparePlots([drawnObj.X,drawnObj.Y],points,['rotation'])
            drawnObj.X = points(:,1);
            drawnObj.Y = points(:,2);
            drawnObj.theta = drawnObj.findAngleFromCenter();
        end
          function theta = findAngleFromCenter(drawnObj)
%               finds rotation from center point 
            theta = atan2((drawnObj.center(2) - drawnObj.Y(1)),...
                (drawnObj.center(1)-drawnObj.X(1))) ;
          end      
        function points = rotateBy(drawnObj, theta)
%             rotates an image 
            c = drawnObj.center;
            x = drawnObj.X;
            y = drawnObj.Y;
            adjustedX = (x - c(1));
            adjustedY = (y - c(2));
            x2 = adjustedX*cos(theta) - adjustedY*sin(theta) + c(1);
            y2 = adjustedX*sin(theta) + adjustedY*cos(theta) + c(2);
            siz = size(x2);
            points = zeros(siz(1),2,siz(2));
            points(:,1,:) = x2;
            points(:,2,:) = y2;
        end       
% % % % %       END STEP 2
% % % % %       START STEP 3

        function scaleToSquare(drawnObj, scaleSize)
% %             Scales the image to a bounding box
            Box = drawnObj.boundingBox();
            width = abs(Box(2,1) - Box(1,1));
            height = abs(Box(2,2) - Box(1,2));
            qx =  drawnObj.X.* (scaleSize(1)./width);
            qy = drawnObj.Y.* (scaleSize(2)./height);
            %                drawnObj.comparePlots([drawnObj.X,drawnObj.Y],[qx,qy],['scale'])
            drawnObj.X = qx;
            drawnObj.Y = qy;
        end
        function Box =  boundingBox(drawnObj)
%             creates a bounding box of the drawn object
            Box = [...
                min(drawnObj.X), min(drawnObj.Y)
                max(drawnObj.X), max(drawnObj.Y)
                ];
        end
        function translateToOrigin(drawnObj)
%             creates a translation to origin
            drawnObj.findCenter();
            qx = drawnObj.X - drawnObj.center(1);
            qy = drawnObj.Y- drawnObj.center(2);
            %             drawnObj.comparePlots([drawnObj.X,drawnObj.Y],[qx,qy],['translation'])
            drawnObj.X = qx;
            drawnObj.Y = qy;
            
        end
% % % % %       END STEP 3

        function comparePlots(drawnObj,oldPoints,newPoints, description)
%             for debugging, compares changes visually
            if(~exist('description','var'))
                description = 'Comparing points';
            end
            figure
            plot(newPoints(:,1),newPoints(:,2),'ro')
            hold on
            plot(oldPoints(:,1),oldPoints(:,2),'kx')
            hold off
            title(description)
        end

        
% % % % %       START STEP 1
                function getPointDistance(drawnObj)
%                     grabs the distance, cumulative distance, and total
%                     distances 
            points = [drawnObj.oldX, drawnObj.oldY];
            d1 = points(2:end,:) - points(1:end-1,:);
            d2 = sqrt(d1(:,1).^2+d1(:,2).^2);
            drawnObj.distance = [0;d2];
            drawnObj.cumulativeDist = cumsum(d2);
            drawnObj.totalDistance = drawnObj.cumulativeDist(end);
        end
        function resamplePath(drawnObj, n)
            %           n points to split the path into
            pathLength = (drawnObj.cumulativeDist);
            I = pathLength(end)/ (n-1);
            D = 0;
            points =[drawnObj.oldX,drawnObj.oldY];
            points1 = points;
            newPoints = [drawnObj.oldX(1),drawnObj.oldY(1)];
            %          for i = 1:length(drawnObj.oldX)
            index = 2;
            %             while length(points(:,1)) >= index
            for index = 2:length(points(:,1))
                d = sqrt(sum((points (index,:) - points (index-1,:)).^2));
                if D + d >= I
                    qx = points(index-1,1) +...
                        ((I - D)/d * (points(index,1) - points(index-1,1)));
                    qy = points(index-1,2) +...
                        ((I- D)/d * (points(index,2) - points(index-1,2)));
                    q = [qx,qy];
                    newPoints = [newPoints;q];
                    points(index,:) = q;
                    D = 0;
                    continue
                end
                D = D+d;
                %                 index = index+1;
            end
            drawnObj.X = newPoints(:,1);
            drawnObj.Y = newPoints(:,2);
         end
 function findCenter(drawnObj)
            meanX = mean(drawnObj.X );
            meanY = mean(drawnObj.Y);
            drawnObj.center = [meanX, meanY];
        end
% % % % %       END STEP 1

% % % % %         START STEP 0
        function storeAsTemplate(drawnObj,groundTruth,ID)
            %             Why is there an ID? just in case.
            drawnObj.ID = ID;
            drawnObj.groundTruth = groundTruth;
        end
        function storeValues(drawnObj,points)
            drawnObj.oldX = points(:,1);
            drawnObj.oldY = points(:,2);
            drawnObj.t = points(:,3);
        end
       % *STEP 0 - STORE AND REMOVE BAD POINTS *
        function processPoints(drawnObj,points)
            drawnObj.calculateDistanceBetweenPoints(points);
            points = drawnObj.removeClosePoints(points);
            drawnObj.storeValues(points);
        end
        function removeFeatures(drawnObj,featureToRemove)
            drawnObj(:,:).features(featureToRemove) = [];
        end
        function calculateDistanceBetweenPoints(drawnObj,points)
            drawnObj.euDist = sqrt(...
                sum((points(3:end,1:2)...
                - points(1:end-2,1:2)).^2,2));
        end
        function [points] = removeClosePoints(drawnObj,points)
            % remove points that are too close to one another that
            % could throw off the processing
            %               badPoints = find(drawnObj.euDist(2:end-1) < 0.01)
            len=length(drawnObj.euDist);
            badCounter = 1;
            badPoints = 1;
            for i=2:len-2
                if drawnObj.euDist(i) < 0.01
                    drawnObj.euDist(i+1)=drawnObj.euDist(i+1)+drawnObj.euDist(i);
                    badPoints(badCounter )=i;
                    badCounter =badCounter +1;
                end
            end
            if isempty(badPoints)
                return
            end
            if(length(badPoints) > length(drawnObj.euDist(:,1))-3)
                badPoints = badPoints(badPoints < length(drawnObj.euDist(:,1))-3);
            end
            drawnObj.euDist(badPoints,:)=[];
            points(badPoints+2,:) = [];
        end
% % % % %         END STEP 0 

% % % % %  START STEP 4 
        function [bestTemplate, Score] = Recognize(drawnObj, templateObjects)
            tempLen = length(templateObjects);
            score = ones(tempLen,1)*9999999;
            tic
            drawnObj.findCenter();
            drawnObj.findAngleFromCenter();
%             pre-declaring the range shifting
            rangeMax = 5;
            rangeInc = 1;
            angleShift = ([-rangeMax:rangeInc:rangeMax] .* pi/180);
%             getting the rotated series of points
            rotatedPoints = drawnObj.generateAngleShifts(angleShift);
%             compare the drawn object to the group visually (plot)
%             drawnObj.compareToGroup(rotatedPoints,drawnObj.center, angleShift);
%             now determine the best class fit!
            for i = 1:length(templateObjects)
                templateObj = templateObjects(i);
                templateObj.findCenter();
                templateObj.findAngleFromCenter();
                score(i) = drawnObj.DistanceAtBestAngle(templateObj,rotatedPoints);
            end
            
            toc
            dSize = length(drawnObj.X);
%             score = 1 - score./(0.5*sqrt(dSize^2+ dSize^2));
            
            %             Get the best fit to the line
            [max1,class1] = max(score);
            %             then get the second best by storing the distance and removing
            %             the previous minimum
            score2 = score;
            score2(class1) = [];
            [max2,class2] = max(score2);
            class2 = find(score == max2);
            %             Get classification from the respective objects
            class1 = templateObjects(class1).groundTruth;
            class2 = templateObjects(class2).groundTruth;
            %             see if they're close or not, if they aren't just pass the
            %             best.
            if max1 * 1.1 < max2
                bestTemplate = {class1;class2};
                Score = [max1;max2];
                return
            end
            bestTemplate = class1;
            Score = max1;
            
        end
        function rotatedPoints = generateAngleShifts(drawnObj, angleShift)
            % drawnObj.theta +
            rotLen = length(drawnObj.X);
            angleLen = length(angleShift);
            rotatedPoints = zeros(rotLen,2,angleLen);
            rotatedPoints = drawnObj.rotateBy(-angleShift);
            
        end
        function dist = DistanceAtBestAngle(drawnObj,templateObj, rotatedPoints)
            dist = PathDistance(rotatedPoints,templateObj);
        end
        function dist = PathDistance(rotatedPoints,templateObj)
%             finds the distances between template points and drawn points
            Tx = templateObj.X;
            Ty = templateObj.Y;
            rotSiz = size(rotatedPoints);
            tLen = length(Tx);
            x = rotatedPoints(:,1,:);
            y = rotatedPoints(:,2,:);
            if tLen > rotSiz(1)
                Tx = Tx(1:rotSiz(1));
                Ty = Ty(1:rotSiz(1));
                distSize = rotSiz(1);
            else
                x = x(1:tLen,:,:);
                y = y(1:tLen,:,:);
                distSize= tLen;
            end
            points = zeros([length(x),rotSiz(2), rotSiz(3)]);
            distance = drawnObj.euclidean(x,y,Tx,Ty);
            %             distance = drawnObj.mahalanobis(points, x,y,Tx, Ty, templateObj.center);
            score = 1 - distance./(0.5*sqrt(distSize^2+ distSize^2));
            
            dist = max(max(score));
        end
    end
    methods(Static)
        function loadedObj = loadobj(saved)
%             loads object from save file
            if isstruct(saved)
                newObj = drawnObj(1);
                newObj.groundTruth = saved.groundTruth;
                newObj.ID = saved.ID;
                newObj.X = saved.X;
                newObj.Y = saved.Y;
                newObj.center = saved.center;
                newObj.theta = saved.theta;
                loadedObj = newObj;
            else
                loadedObj = saved;
            end
        end
        function mDistance = mahalanobis(points,x,y,Tx,Ty, center)
            %             finds mahalanobis distance of point from
            %             distribution within templates
            points(:,1,:) = x - Tx;
            points(:,2,:) = y - Ty;
            coVar = (Tx*Ty');
            meanT = center;
            rotSiz = size(points);
            mDistance = zeros(rotSiz(3));
            for i = 1:rotSiz(3)
                XX = (points(:,:,i) - meanT);
                mDistance(i) = sum(sqrt(((XX' * inv(coVar) * XX))),'all');
            end
        end
        function eDist = euclidean(x1,y1,x2,y2)
%             finds euclidean distance and sums it
            x = x1 - x2;
            y = y1 - y2;
            eDist = sum(sqrt(x.^2 + y.^2));
            
        end
        function compareToGroup(points, center, titles)
%             Creates a set of plots to compare against the whole group 
            figure(1);
            clf;
            n = round(sqrt(size(points,3)))+1;
            
            for i = 1:size(points,3)
                subplot(n,n,i+1)
                %                  template = templateObjects(i)
                plot(points(:,1,i),points(:,2,i))
                hold on
                plot(center(1),center(2),'gs')
                title(titles(i))
            end
            hold off
            
        end
        function compareToTemplates(drawnObj,templateObjects)
%             compares the templates to what you drew
            figure(2)
            clf
            objectCount = length(templateObjects) + 1;
            objectSize = round(sqrt(objectCount))+1
            subplot(objectSize,objectSize,1)
            plot(drawnObj.X,drawnObj.Y,'rh')
            hold on
            title(drawnObj.theta);
            for i = 1:length(templateObjects)
                subplot(objectSize,objectSize,i+1)
                template = templateObjects(i);
                plot(template.X,template.Y)
                title(template.theta);
            end
            hold off
            
        end
    end
end

