%% STARTING

function main
clear all, close all

%importing labels from txt
load vesselLabels.txt;


% ----------------------- CONST ------------------------- %
RegionBuffer = [];
stepRoi = 4;

%baseBkg = 13; % Initial Frame: 0 %
baseNum = 13;

% To use txt values use nVesselLabels = nFrames + 1 %
% nVesselLabels start in 1 and nFrames starts in 0  %
nTotalFrames = 1433; % Total: 1433
nInitialFrame = 0000;  % Initial Boat: 12

thr_global = 180; % 180

thr_diff = 25;    % 18 or 60 fails detecting the boat sometimes

minArea = 50;  % 50 or 100
maxArea = 1000; % 1000? maybe

%1cm=58pixes(units)
distanceBetweenVessels = 100; % 58? or 80? 

bufferArr = [];
reglist = [];
rectangleAux = [];
arrAllIndsRectangleAux = [];
indsTemp = [];
mainFigure = figure(1);
vesselTrail = [];
vesselTrailSREShift = [];
vesselTrailNow = [];
%matrix  
% 0 1 0  _
% 1 1 1   |
% 0 1 0   v
se = strel('disk',3);
vector=[];
maxBufferNum = 7;
numFrameIterations = 0;
numFrameIterationsAux = 0;

%Creating Output_labelling.txt but we don't put anything inside
fileID = fopen('Output_labelling.txt','wb');
fprintf(fileID,'%6s %2s %6s %10s %9s\n','Frame Number','X','Y','Width','Height');


% ---------------------- END Const ---------------------- %

% ------------------------------------------------------- %

% ---------------------- ROI ---------------------------- %
% Remove object intersection
% Faz as caixinhas

% ------------------------------------------------------- %
% --------------------- BUFFER -------------------------- %
% ------------------------------------------------------- %

%	Temporal validation is equal to 8 Spacial validation

bufferStructNames = ['a'; ...
    'b'; ...
    'c'; ...
    'd'; ...
    'e'; ...
    'f'; ...
    'g'];

bufferStruct = struct('a', {}, ...
    'b', {}, ...
    'c', {}, ...
    'd', {}, ...
    'e', {}, ...
    'f', {}, ...
    'g', {});



% ------------------------------------------------------- %
% ------------------- END BUFFER ------------------------ %
% ------------------------------------------------------- %

for f = nInitialFrame : stepRoi : nTotalFrames
    %if f > 100
       %break; 
    %end
    array_inds = [];
    labelDraw=[];
    
    imgfrNew = imread(sprintf('../Frames/frame%.4d.jpg', ...
        baseNum + f));
    
    % ------------------------------------------------------ %
    % Buffer Shift Lines
    % ------------------------------------------------------ %
    
    % if the buffer is not full, use only the length of the buffer
	% so we can put more frames inside the buffer

    if numFrameIterations < 8
        numFrameIterationsAux = maxBufferNum;
        
        
    else
        %if iterations > 7 the buffer is full and it's possible to shift
        numFrameIterationsAux = numFrameIterations;
        
        bufferStruct(1).g = bufferStruct(1).f;
        bufferStruct(1).f = bufferStruct(1).e;
        bufferStruct(1).e = bufferStruct(1).d;
        bufferStruct(1).d = bufferStruct(1).c;
        bufferStruct(1).c = bufferStruct(1).b;
        bufferStruct(1).b = bufferStruct(1).a;
        
        
        
    end
    
    % ------------------------------------------------------ %
    % END Buffer Shift Lines
    % ------------------------------------------------------ %
    hold on    
    imgdif = (abs(double(imgfrNew(:,:,1)))>thr_global) | ...
        (abs(double(imgfrNew(:,:,2))-double(imgfrNew(:,:,1)))>thr_diff) | ...
        (abs(double(imgfrNew(:,:,3))-double(imgfrNew(:,:,1)))>thr_diff);
    
    
    bw = imclose(imgdif,se);
    bw = imclose(imgdif,se);
    bw = imclose(imgdif,se);
    bw = imclose(imgdif,se);
    str = sprintf('Frame: %d',f);
    title(str);
    
	% HOW TO SEE THE MOVIE IN FRAMES? IN BALCK&WHITE OR REAL VIDEO?
	
    % ----------------------------------------------------------- %
    imshow(bw);  %%Mete Background preto, bom para ver threshold
    % ----------------------------------------------------------- %
    
    % ----------------------------------------------------------- %
    %imshow(imgfrNew);  %%Faz o mesmo video mas com cor
    % ----------------------------------------------------------- %
    
    [lb num]=bwlabel(bw);
    regionProps = regionprops(lb,'area','FilledImage','Centroid');
    
    inds = [];
    for k = 1 : length(regionProps)
        if find([regionProps(k).Area] < maxArea & [regionProps(k).Area] > minArea)
            inds = [inds k];
        end
    end
    
    regnum = length(inds);
    
    %if number of regions in one image accepted > 0 
    if regnum
        
        % ------------------------------------------------------ %
        % Spatial Validation Algorithm
        % ------------------------------------------------------ %
        
        for k=1:regnum
            for m=1:regnum
                if k ~= m
                    %run all inds to search for inds that are too close
                    vesselAX = regionProps(inds(k)).Centroid(1,1);
                    vesselAY = regionProps(inds(k)).Centroid(1,2);
                    vesselBX = regionProps(inds(m)).Centroid(1,1);
                    vesselBY = regionProps(inds(m)).Centroid(1,2);
                    
                    distBetweenVessels = [vesselAX, vesselAY; ...
                        vesselBX, vesselBY];
                    pdistBetweenVessels = pdist(distBetweenVessels, 'euclidean');
                    
                    if pdistBetweenVessels < distanceBetweenVessels

                        arrayDetection = ismember([k m],array_inds);

                        %k is not find on array_inds
                        if arrayDetection(1,1) == 0
                            array_inds = [array_inds k];
                        end
                        
                        %m is not find on array_inds
                        if arrayDetection(1,2) == 0
                            array_inds = [array_inds m];
                        end
                    end
                end
            end
        end

        
        % ------------------------------------------------------ %
        % END Spatial Validation Algorithm
        % ------------------------------------------------------ %
        
      
        allInds = inds;
 
        allInds(array_inds) = [];
        
        % NOW allInds have only vessels aproved by spacial validation algoritm
        
        %%USE if spacial validation is not used, to check if is true
        %array_inds = inds;
        
        % ----------------------------------------------------------- %
        %%%Temporal Buffer
        % ----------------------------------------------------------- %
        %%change buffer lines from 1-6 to 2-7
        %%frame 7 will be overwriten
        
        
        
        % ------------------------------------------------------ %
        % Temporal Validation Algorithm
        % ------------------------------------------------------ %
        
        % Converting from lin col to rectangle format
        
        regnumAllInds = length(allInds);
        
        %if existes regions to filter with temporal algorithm
        if regnumAllInds
            arrAllIndsRectangleAux = [];
            for j = 1 : regnumAllInds
                [lin, col] = find(lb == allInds(j));
                upLPoint = min([lin col]);
                dWindow  = max([lin col]) - upLPoint + 1;
                
                rectangleAux = [fliplr(upLPoint) fliplr(dWindow)];
                
                % add the rectangleAux to arrAllIndsRectangleAux
                arrAllIndsRectangleAux = [arrAllIndsRectangleAux; rectangleAux];

            end
            % add the arrAllIndsRectangleAux to buffer first line

            
            %if numFrameIterations < 7 buffer is not full
            if numFrameIterations < 7
                numFrameIterationsAux = maxBufferNum;
                
                
                %%%%%%%%%%%%To REcode
                if numFrameIterations == 0
                    bufferStruct(1).g = arrAllIndsRectangleAux;
                end
                if numFrameIterations == 1
                    bufferStruct(1).f = arrAllIndsRectangleAux;
                end
                if numFrameIterations == 2
                    bufferStruct(1).e = arrAllIndsRectangleAux;
                end
                if numFrameIterations == 3
                    bufferStruct(1).d = arrAllIndsRectangleAux;
                end
                if numFrameIterations == 4
                    bufferStruct(1).c = arrAllIndsRectangleAux;
                end
                if numFrameIterations == 5
                    bufferStruct(1).b = arrAllIndsRectangleAux;
                end
            else
                
                bufferStruct(1).a = arrAllIndsRectangleAux;
                
            end
                       
            %%%%%%%%%%%%EndTo REcode
            
            % ----------------------------------------------------------- %
            %   filtering vessels to know what to print
            % ----------------------------------------------------------- %
            % bufferStruct(1).a has the possible vessels to print
            
            %%calculates the number of occurencies of vessels in Layers A on
            %%other Layers
            
            [colA,n] = size(bufferStruct(1).a);
            vesselOcurrencies = zeros(1,colA);
            
            
            vesselOcurrencies = vesselOcurrencies + foundOnBufferLayer(bufferStruct(1).a,bufferStruct(1).b);
            
            vesselOcurrencies = vesselOcurrencies + foundOnBufferLayer(bufferStruct(1).a,bufferStruct(1).c);
            
            vesselOcurrencies = vesselOcurrencies + foundOnBufferLayer(bufferStruct(1).a,bufferStruct(1).d);
            
            vesselOcurrencies = vesselOcurrencies + foundOnBufferLayer(bufferStruct(1).a,bufferStruct(1).e);
           
            vesselOcurrencies = vesselOcurrencies + foundOnBufferLayer(bufferStruct(1).a,bufferStruct(1).f);
            
            vesselOcurrencies = vesselOcurrencies + foundOnBufferLayer(bufferStruct(1).a,bufferStruct(1).g);
            
            %%NOW vesselOcurrencies has counter with numbers of ocurrencies in
            %%of the vessels in first layer with the rest of the buffer
            
            indsTemp = [];
            for r=1:colA
                  % -----------------3 out of 5-------------------- %
                if vesselOcurrencies(1,r) > 2
                    indsTemp = [indsTemp r];
                end
            end
            
            %BufferStruct is incremented because buffer will be incremented
            numFrameIterations = numFrameIterations + 1;         
            
            % ------------------------------------------------------ %
            % END Temporal Validation Algorithm
            % ------------------------------------------------------ %
            
            % ----------------------------------------------------------- %
            %doing boxes on approved inds
            % ----------------------------------------------------------- %
            
            %regnumAllInds = length(allInds); % change variables
            
            %number of yellow boxes to print
            [nIndsTemp,m] = size(indsTemp);

            
            for j=1: 1: nIndsTemp % change variables
                vesselTrailNow = bufferStruct(1).a(indsTemp(1,j),:);
            end
        end
    end
    
    [linLabel colLabel] = find (vesselLabels(:,1) == f+14);
    
    if colLabel == 1
        labelDraw = [labelDraw vesselLabels(linLabel,2:5)];
    end
    isnotempty = 0;
    if ~isempty(labelDraw)
		%se o width for negativo
		if labelDraw(3) < 0
			labelDraw(1)=labelDraw(1) + labelDraw(3);
			labelDraw(3) = abs(labelDraw(3));
		end
		%se o height for negativo
		if labelDraw(4) < 0
			labelDraw(2)=labelDraw(2) + labelDraw(4);
			labelDraw(4) = abs(labelDraw(4));
		end
		rectangle('Position', labelDraw,'EdgeColor',[0 1 0],'linewidth',2);
        if ~isempty(vesselTrailNow)
            isnotempty = 1;
            
            vesselTrail = [vesselTrail; f+1 labelDraw];
            vesselTrailSREShift = [vesselTrailSREShift; f * 1.10 + 1 labelDraw];
            
			fprintf(fileID,'%6d %9d %6d %9d %9d\n',f,vesselTrailNow(1),vesselTrailNow(2),vesselTrailNow(3),vesselTrailNow(4));
           
            vector=[vector bboxOverlapRatio(labelDraw, vesselTrailNow)];            
            
            
            labelDraw1R = labelDraw;
            labelDraw1R(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw1R(3) = labelDraw(3) * 0.80;
            labelDraw1R(4) = labelDraw(4) * 0.80;
            
            labelDraw1U = labelDraw;
            labelDraw1U(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw1U(3) = labelDraw(3) * 0.80;
            labelDraw1U(4) = labelDraw(4) * 0.80;
            
            labelDraw1D = labelDraw;
            labelDraw1D(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw1D(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw1D(3) = labelDraw(3) * 0.80;
            labelDraw1D(4) = labelDraw(4) * 0.80;
            
            labelDraw2R = labelDraw;
            labelDraw2R(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw2R(3) = labelDraw(3) * 0.90;
            labelDraw2R(4) = labelDraw(4) * 0.90;
            
            labelDraw2U = labelDraw;
            labelDraw2U(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw2U(3) = labelDraw(3) * 0.90;
            labelDraw2U(4) = labelDraw(4) * 0.90;
            
            labelDraw2D = labelDraw;
            labelDraw2D(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw2D(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw2D(3) = labelDraw(3) * 0.90;
            labelDraw2D(4) = labelDraw(4) * 0.90;
            
            labelDraw3R = labelDraw;
            labelDraw3R(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw3R(3) = labelDraw(3) * 1.10;
            labelDraw3R(4) = labelDraw(4) * 1.10;
            
            labelDraw3U = labelDraw;
            labelDraw3U(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw3U(3) = labelDraw(3) * 1.10;
            labelDraw3U(4) = labelDraw(4) * 1.10;
            
            labelDraw3D = labelDraw;
            labelDraw3D(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw3D(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw3D(3) = labelDraw(3) * 1.10;
            labelDraw3D(4) = labelDraw(4) * 1.10;
            
            labelDraw4R = labelDraw;
            labelDraw4R(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw4R(3) = labelDraw(3) * 1.20;
            labelDraw4R(4) = labelDraw(4) * 1.20;
            
            labelDraw4U = labelDraw;
            labelDraw4U(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw4U(3) = labelDraw(3) * 1.20;
            labelDraw4U(4) = labelDraw(4) * 1.20;
            
            labelDraw4D = labelDraw;
            labelDraw4D(1) = labelDraw(1) + labelDraw(1) * 0.10;
            labelDraw4D(2) = labelDraw(2) + labelDraw(2) * 0.10;
            labelDraw4D(3) = labelDraw(3) * 1.20;
            labelDraw4D(4) = labelDraw(4) * 1.20;
            
            rectangle('Position',vesselTrailNow,'EdgeColor',[1 1 0],...
                    'linewidth',2);
                
            
            
        else
            vector=[vector 0];
        end 
    else
        vector=[vector 0];
    end
        
    drawnow
   
end

mFigure = figure('Name','IoU')
title('Graphic')
xlabel('Frames') % x-axis label
ylabel('Ratio') % y-axis label
plot(vector);
grid on
grid minor
xlim([0 400]); % x-axis limits
ylim([-0.4 1]); % y-axis limits

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw1R);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw1U);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw1D);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw2R);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw2U);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw2D);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw3R);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw3U);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw3D);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw4R);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw4U);
grid on
grid minor

mFigureSRE = figure('Name','SRE: Success Plot')
title('Graphic')
xlabel('Distance')    % x-axis label
ylabel('Frames')      % y-axis label
plot(labelDraw4D);
grid on
grid minor

fclose(fileID);
end


% ------------------------------------------------------ %
% Function foundOnBufferLayer( layerA, layerN )
% ------------------------------------------------------ %

%input layer a and other layer and as output vector with number of
%proximities found
function bufferCount = foundOnBufferLayer( layerA, layerN )
%%%Layer A is bufferStruct(1).a
%%%if bufferStruct(1).a is a matrix
[p,q] = size(layerA);
[r,s] = size(layerN);
colLayerA = p;
colLayerN = r;
bufferCount = zeros(1,p);
distanceBetweenVessels = 80;
%LayerA is a Matrix

%k index in buffer.a
for k=1:colLayerA
    isFirstLayer = false;
    for m=1:colLayerN
        
        vesselA = layerA(k,:);
        vesselN = layerN(m,:);
        %search if A(k) and N(m) vessels are close
        vesselAX = vesselA(1)+vesselA(3)/2;
        vesselAY = vesselA(2)+vesselA(4)/2;
        vesselNX = vesselN(1)+vesselN(3)/2;
        vesselNY = vesselN(2)+vesselN(4)/2;
        distBetweenVessels = [vesselAX, vesselAY; ...
            vesselNX, vesselNY];
        pdistBetweenVessels = pdist(distBetweenVessels, 'euclidean');
        
        if pdistBetweenVessels < distanceBetweenVessels
            isFirstLayer = true;
            if isFirstLayer == true
                bufferCount(k) = bufferCount(k) + 1;
            end
        end
    end
end

end

