%   Spectral Imaging Toolbox 1.0 by Miles Aron
%   Copyright (C) 2017  Miles Aron
%   Please cite the corresponding publication if you publish work inspired by or using this code
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   See included user manual for detailed instructions.

close all
format long

% add dependencies to path so user doesn't have to
addpath('sc/')
addpath('PATCH_3Darray')
addpath('Demo')
addpath(genpath('Demo'))
addpath('bfmatlab')

inputGUI() % Launch GUI to get inputs
uiwait % Wait until inputs are submitted

dirParts = strsplit(inputDir, filesep); % Parts of input directory path
dirName = dirParts{end}; % Name of input directory
dirData = dir([inputDir filesep '*']); % Input directory contents - must contain only images with ome-compatible meta-data
fileList = {dirData.name}'; % Names of files in input directory

% Remove files that are clearly not images
for i=1:length(fileList)
    if dirData(i).isdir
        fileList{i}=[];
    end
     if length(fileList{i}) < 5
         fileList{i}=[];
     end
end
fileList = fileList(~cellfun('isempty',fileList)); 

% Name and make directories for saving
parentDir = [outputDir filesep dirName '-results']; % parent directory for output
mkdir(parentDir);

% Pre-segmentation results output directories
resultsDir = [parentDir filesep 'Results - no segmentation']; % directory for pre-segmentation results
GPDir = [resultsDir filesep 'GP']; % directory for pre-segmentation GP data
histogramDir = [resultsDir filesep 'GP Histograms']; % directory for pre-segmentation GP histograms
spectraDir = [resultsDir filesep 'Spectra']; % directory for pre-segmentation spectra
mkdir(resultsDir);
mkdir(GPDir);
mkdir(histogramDir);
mkdir(spectraDir);

if ismember(1,options) % Segmented membrane results output directories
    resultsMembranesDir = [parentDir filesep 'Results - segmented membranes']; % directory for segmented membrane results
    GPmembranesDir = [resultsMembranesDir filesep 'GP']; % directory for segmented membrane GP data
    histogramMembranesDir = [resultsMembranesDir filesep 'GP Histograms']; % directory for segmented membrane GP histograms
    spectraMembranesDir = [resultsMembranesDir filesep 'Spectra']; % directory for segmented membrane spectra
    mkdir(resultsMembranesDir);
    mkdir(GPmembranesDir);
    mkdir(histogramMembranesDir);
    mkdir(spectraMembranesDir);
end
if ismember(1,options) || ismember(2,options) % Segmented object results output directories
    resultsObjectsDir = [parentDir filesep 'Results - segmented objects']; % directory for segmented object results
    GPobjectsDir = [resultsObjectsDir filesep 'GP']; % directory for segmented object GP data
    histogramObjectsDir = [resultsObjectsDir filesep 'GP Histograms']; % directory for segmented object GP histograms
    spectraObjectsDir = [resultsObjectsDir filesep 'Spectra']; % directory for segmented object spectra
    regionPropsObjectsDir = [resultsObjectsDir filesep 'Object properties']; % directory for segmented object region properties
    mkdir(resultsObjectsDir);
    mkdir(GPobjectsDir);
    mkdir(histogramObjectsDir);
    mkdir(spectraObjectsDir);
    mkdir(regionPropsObjectsDir);
end
if ismember(3,options) % Segmented spherical object results output directories
    resultsCirclesDir = [parentDir filesep 'Results - segmented spherical objects']; % directory for segmented spherical object results
    GPcirclesDir = [resultsCirclesDir filesep 'GP']; % directory for segmented spherical object GP data
    histogramCirclesDir = [resultsCirclesDir filesep 'GP Histograms']; % directory for segmented spherical object GP histograms
    spectraCirclesDir = [resultsCirclesDir filesep 'Spectra']; % directory for segmented spherical object spectra
    sizeDistributionCirclesDir = [resultsCirclesDir filesep 'Size distribution']; % directory for segmented spherical object size distribution
    circleDetectionCirclesDir = [resultsCirclesDir filesep 'Circle detection']; % directory for segmented spherical object detection
    mkdir(resultsCirclesDir);
    mkdir(GPcirclesDir);
    mkdir(histogramCirclesDir);
    mkdir(spectraCirclesDir);
    mkdir(sizeDistributionCirclesDir);
    mkdir(circleDetectionCirclesDir);
end

% Initialize variables
nImg = numel(fileList); % number of images
spectraWeighted = cell(nImg,2); % mean (col 1) and std (col 2) of intensity  for pixels with signal (after applying thresholding mask) at each wavelength
[membraneSpectra, diameters, meanGP, stdGP, medGP, pxCountGP, meanObjectGP,...
    stdObjectGP, medObjectGP, pxCountObjectGP, meanMembraneGP, stdMembraneGP, ...
    medMembraneGP, pxCountMembraneGP, membraneCentroids, membraneDiameters, ...
    wavelengthsAll, GPfull] = deal(cell(nImg,1));
sizeRangeGlobal = 0;
imgObjNameInd = 1;
imgMembraneNameInd = 1;
imgObjNames = cell(1,1);
imgMembraneNames = cell(1,1);

if(~isnan(GPref) && ~isempty(refIm)) % load reference image specified in input GUI
    if (isnan(GPref) || isempty(refIm))
        GPcorrectionFactorWarning = warndlg('GP correction factor will not be used because either reference GP or reference image was not correctly specified. Please try again if GP correction factor is desired.');
        waitfor(GPcorrectionFactorWarning)
    else
        data=bfopen(refIm); % opens reference image
        imgData = data{1, 1}; % image data
        metadata = data{1, 2}; % raw metadata
        omeMeta = data{1, 4}; % standardized metadata
        nLambda = omeMeta.getChannelCount(0); % number of C slices
        
        % Get wavelengths from metadata
        wavelengthsRef = zeros(nLambda,1);

        for i=1:nLambda-1 % because TD is the last channel and we will ignore it
            chKey = append('Global Name #',string(i+1)); % because Name #1 is (erroneously) called TD
            chVal = metadata.get(chKey); % a string including " nm"
            chWave = strsplit(chVal, " "); % strip the nm if present
            wavelengthsRef(i)=str2double(chWave{1});
            %wavelengthsRef(i)=str2double(omeMeta.getChannelName(0,(i-1)));
        end
        
        % Set wavelengths to values chosen in GUI. If not exact, take closest.
        % In event of two options that are closest, take lower wavelength.
        [~, wavelengthsInd(1,1)] = min(abs(wavelengthsRef-lambdas(1)));
        [~, wavelengthsInd(2,1)] = min(abs(wavelengthsRef-lambdas(2)));
        
        rawRef = imgData(:,1);
        I_low_ref=im2double(rawRef{wavelengthsInd(1)});
        I_high_ref=im2double(rawRef{wavelengthsInd(2)});
        G = (nanmean(I_low_ref(:))*(1-GPref))/(nanmean(I_high_ref(:))*(1+GPref));
    end
end

% Loop through images
i = 1; 
while i <= nImg          
    % Clear memory and close figures
    clear croppedObjectsGP
    close all
    
    % [~,imageID,~] = fileparts(fileList{i}); % for saving
    imageID = strsplit(fileList{i}, '.'); 
    imageID = imageID{1}; % remove file extension for saving
    disp(['processing image: ' imageID]);

    if startsWith(imageID, '._')
    % Skip this file if it starts with dot underline.
        continue; % Jump to the bottom of the loop.
    end

    try
        data=bfopen([inputDir filesep fileList{i}]); % opens image files sequentially from directory
    catch
        disp(['file ' imageID ' is not a bio-formats image. The Spectral Imaging Toolbox requires bio-format-compatible images.']);
        i = i+1;
        imageID = strsplit(fileList{i}, '.'); 
        imageID = imageID{1}; % remove file extension for saving
        disp(['processing image: ' imageID]);
    end
    
    imgData = data{1, 1}; % image data
    metadata = data{1, 2}; % raw metadata
    omeMeta = data{1, 4}; % standardized metadata
    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM); % in µm
    voxelSizeY = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM); % in µm
    stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
    nLambda = omeMeta.getChannelCount(0); % number of C slices
    
    % Get wavelengths from metadata

    wavelengths = zeros(nLambda,1);
    for cc=1:nLambda-1 % because TD is the last channel and we will ignore it
        chKey = append('Global Name #',string(cc+1)); % because Name #1 is (erroneously) called TD
        chVal = metadata.get(chKey); % a string including " nm"
        chWave = strsplit(chVal, " "); % strip the nm
        wavelengths(cc)=str2double(chWave{1});
        %wavelengths(cc)=str2double(omeMeta.getChannelName(0,(cc-1)));
    end
    wavelengthsAll{i} = wavelengths;
    
    % Set wavelengths to values chosen in GUI. If not exact, take closest.
    % In event of two options that are closest, take lower wavelength.
    [~, wavelengthsInd(1,1)] = min(abs(wavelengths-lambdas(1)));
    [~, wavelengthsInd(2,1)] = min(abs(wavelengths-lambdas(2)));

    if stackSizeZ>1 % Z-stack
        voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM); % in µm
        % Split stack every n_wavelengths
        rawData = cell(stackSizeZ,1);
        for ii=1:stackSizeZ
            ind1 = (ii-1)*nLambda+1;
            ind2 = ii*nLambda;
            rawData{ii} = imgData(ind1:ind2,1); % raw image data
        end
    else % Open spectral image stack
        rawData = imgData(:,1); % raw image data
    end
    
    GPslices = cell(nImg, stackSizeZ);
    for k=1:stackSizeZ % process z-stack
        if stackSizeZ>1
            raw = rawData{k}; % Slice of Z-stack
            imageTag = [imageID '-slice' num2str(k)]; % Include slice number in image tag for saving
        else
            raw = rawData; % Single slice spectral image stack
            imageTag = imageID;
        end
        
        % Initialize variables
        summedImages = zeros(size(raw{1,1}));
        
        % Sum images in stack for selecting appropriate thresholding level
        for c=1:nLambda % for each wavelength
            raw{c}= im2double(raw{c}); % convert to double
            raw{c}(isnan(raw{c}))=0; % removal of NaN values from wavelength image for future processing
            summedImages=summedImages+raw{c}; % summed intensity of all wavelength images
        end
        
        % Generate thresholding mask
        I = im2bw(summedImages, graythresh(summedImages)); % Determine threshold level by Otsu's method with the summed image
        
        % Apply thresholding and generate spectrum
        thresholded = raw; % create thresholded image object array
        cSliceWeighted = zeros(nLambda,1); 
        
        for c=1:nLambda % for each wavelength
            % Apply median filtering if enabled
            if preProcess
                thresholded{c} = medfilt2(thresholded{c});
            end
            
            % Apply thresholding
            thresholded{c}(I==0) = NaN; 
            cSliceWeighted(c) = nansum(thresholded{c}(:))/numel(I(I==1)); % mean of c-slice in c-stack image
        end
        
        % Plot spectrum and save
        figure,plot(wavelengths, cSliceWeighted, '-m.', 'MarkerSize', 30, 'LineWidth', 3, 'Color', 'b', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
        spectraSettings();
        save([spectraDir filesep 'spectraWeighted'], 'cSliceWeighted') % save spectral data
        savefig([spectraDir filesep imageTag '-spectra.fig']) % Save spectra figure
        I_low = thresholded{wavelengthsInd(1)};
        I_high = thresholded{wavelengthsInd(2)};
        
        if(~isnan(GPref) && ~isempty(refIm)) % calculate GP with correction factor
            GP = (I_low-G.*I_high)./(I_low+G.*I_high);
        else % calculate GP without correction factor
            GP=(I_low-I_high)./(I_low+I_high); 
        end
        
        GP(I_high==0)=NaN; % Removal of 1's from 0 values in 490nm
        GP(I_low==0)=NaN; % Removal of -1's from 0 values in 440nm
        
        % Generate and save pseudocolored GP maps
        figure()
        sc(GP, 'jet','k') % GP map
        GPmapSettings();
        save([GPDir filesep imageTag '-GP'], 'GP') % save whole image GP data
        savefig([GPDir filesep imageTag '-GP.fig']) % save whole image GP map
        GPslices{k}=GP;
        
        meanGP{i} = nanmean(GP(:));
        stdGP{i} = nanstd(GP(:));
        medGP{i} = nanmedian(GP(:));
        pxCountGP{i} = numel(GP(:));

        % Generate and save GP histogram figure and fit
        histogramData = histogramSettings(GP, -1:0.01:1, twoTerm);
        xlim([-1 1])
        a1 = xlabel('GP', 'FontSize', 24);
        a2 = ylabel('N_{pixels}', 'FontSize', 24);
        set([a1 a2], 'interpreter', 'tex');
        save([histogramDir filesep imageTag '-histogramData'], 'histogramData') % Save histogram data and fit parameters
        savefig(gcf,[histogramDir filesep imageTag '-histogram.fig']) % Save histogram figure
        
        if ismember(3,options) % Segmenting spherical objects
            if sizeRangeGlobal==0
                tryMe = GP;
                tryMe(~isnan(tryMe)) = 1;
                tryMe(isnan(tryMe)) = 0;
                findCirclesGUI(tryMe); % launches GUI that allows user to interactively choose size range for imfindcircles and finds the circles
                uiwait
            end
            
            % Find circles
            [centers, radii, metric] = imfindcircles(tryMe,[sizeRange(1) sizeRange(2)],'ObjectPolarity', 'bright', 'Sensitivity', 0.95);
            
            decentCircles = find(metric>0.065);
            centers = centers(decentCircles, :);
            radii = radii(decentCircles);
            
            % Label circles on image
            figure()
            sc(GP,'jet','k');
            hold on
            viscircles(centers, radii,'EdgeColor','b');
            movegui(gcf,'center')
            [objMembraneGP, objProps, meanObjGP, stdObjGP, medObjGP, pxCountObjGP] = deal(cell(length(radii),1));
            objMembraneThresh = cell(length(radii),length(wavelengths));
            nObj = length(radii);
            for jj = 1:nObj
                % Draw circle ROI and create mask
                t = 0:1/ceil(2*radii(jj)):2*pi;
                xi = radii(jj)*cos(t)+centers(jj,1);
                yi = radii(jj)*sin(t)+centers(jj,2);
                circleMask = poly2mask(xi, yi, size(GP,1), size(GP,2));
                
                % Perform image dilation
                se = strel('diamond',4); % structuring element for dilation
                circleMask = imdilate(circleMask,se);
                
                % Get stats for saving later on
                objProps{jj} = regionprops(circleMask, 'all');
                
                % Segment circles
                extents = regionprops(circleMask, 'BoundingBox');
                circle = GP;
                circle(~circleMask) = nan;
                objMembraneGP{jj} = imcrop(circle, extents(1).BoundingBox);
                meanObjGP{jj} = nanmean(objMembraneGP{jj}(:));
                stdObjGP{jj} = nanstd(objMembraneGP{jj}(:));
                medObjGP{jj} = nanmedian(objMembraneGP{jj}(:));
                pxCountObjGP{jj} = numel(objMembraneGP{jj}(:));

                % Apply segmentation to pre-processed images
                for c=1:length(wavelengths) % for each wavelength in c-stack
                    circle = thresholded{c};
                    circle(~circleMask) = nan;
                    objMembraneThresh{jj}{c} = imcrop(circle, extents(1).BoundingBox);
                end
            end
            savefig(gcf,[circleDetectionCirclesDir filesep imageTag '-circles.fig'])
            diameters{i,1} = radii*2;
            meanObjectGP{i} = meanObjGP;
            stdObjectGP{i} = stdObjGP;
            medObjectGP{i} = medObjGP;
            pxCountObjectGP{i} = pxCountObjGP;
        elseif ismember(1,options) || ismember(2,options) % Segmenting objects
            % Segment images into individual objects
            [objMembraneGP, objMembraneThresh, croppedObjectsGP,...
                croppedObjectsThresh, croppedObjectsBoundingBoxes, objProps, membraneSpectrum] = deal(cell(1,1));
            objectSeparatorGUI(thresholded, GP); % Crops objects and provides GUI for separating connected objects
            uiwait % Wait until segmentation into single object images is complete 
            nObj = numel(croppedObjectsGP);
            [meanObjGP, stdObjGP, medObjGP, pxCountObjGP, objCentroids] = deal(cell(nObj,1));
            for ii=1:nObj
                
                objN = num2str(ii);
                objTag = ['-object' objN];
                
                % Object tag for saving
                imageName = fileList{i};
                fullObjTag = [imageName ' object' objN];
                imgObjNames{imgObjNameInd,1} = fullObjTag;
                imgObjNameInd = imgObjNameInd+1;
                    
                % Generate and save pseudocolored cropped object GP maps
                objectGP = croppedObjectsGP{ii};
                close all
                figure()
                sc(objectGP, 'jet','k') % cropped object GP map
                GPmapSettings();
                savefig(gcf,[GPobjectsDir filesep imageTag objTag ' - GP.fig']) % cropped object GP map
                save([GPobjectsDir filesep imageTag objTag '-GP'], 'objectGP') % cropped object GP data
                
                meanObjGP{ii} = nanmean(objectGP(:));
                stdObjGP{ii} = nanstd(objectGP(:));
                medObjGP{ii} = nanmedian(objectGP(:));
                pxCountObjGP{ii} = numel(objectGP(:));
                
                % Store membrane histogram and spectrum figure and data
                cSliceObject = zeros(nLambda,1);
                objectThresh = croppedObjectsThresh{ii};
                for c=1:nLambda % for each wavelength
                    cSliceObject(c) = nansum(objectThresh{c}(:))/numel(~isnan(objectThresh{c}(:))); % mean of c-slice in c-stack
                end
                figure,plot(wavelengths, cSliceObject(:), '-m.', 'MarkerSize', 30, 'LineWidth', 3, 'Color', 'b', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
                spectraSettings();
                save([spectraObjectsDir filesep imageTag objTag '-spectrum'], 'cSliceObject') % Save object spectral data
                savefig([spectraObjectsDir filesep imageTag objTag '-spectrum.fig']) % Save object spectra figure
                histogramData = histogramSettings(objectGP, -1:0.01:1, twoTerm);
                xlim([-1 1])
                a1 = xlabel('GP', 'FontSize', 24);
                a2 = ylabel('N_{pixels}', 'FontSize', 24);
                set([a1 a2], 'interpreter', 'tex');
                save([histogramObjectsDir filesep imageTag objTag '-histogramData'], 'histogramData') % Save histogram data and fit parameters
                savefig(gcf,[histogramObjectsDir filesep imageTag objTag '-histogram.fig']) % Save histogram figure
                
                % Segment membranes
                if ismember(1,options) % 1 if performing membrane segmentation
                    % Make mask from cropped object GP map
                    mask = croppedObjectsGP{ii};
                    mask(~isnan(croppedObjectsGP{ii})) = 1;
                    mask(isnan(croppedObjectsGP{ii})) = 0;
                    
                    % Remove residual disconnected areas leftover from cropping
                    stats = sort(cell2mat(struct2cell(regionprops(mask,'Area'))));
                    if length(stats)>1
                        blob = min(stats(end));
                        mask = bwareaopen(mask, blob-1);
                    end
                    
                    % Fill holes
                    mask = imfill(mask, 'holes');
                    
                    % Get stats for saving later on
                    objPropTemp = regionprops(mask, 'EquivDiameter');
                    objProps{ii} = objPropTemp.EquivDiameter*voxelSizeX.doubleValue(); % diameter
                    
                    % Centroids to aid comparing before and after images
                    objCentroidTemp = regionprops(mask, 'Centroid');
                    objCentroidTemp = objCentroidTemp.Centroid; % centroid
                    bbTemp = croppedObjectsBoundingBoxes{ii};
                    
                    % Segment membranes
                    [~, threshold] = edge(mask, 'sobel');
                    fudgeFactor = .5;
                    mask = edge(mask,'sobel', threshold * fudgeFactor);
                    se90 = strel('line', 5, 90); % Valerio 
                    se0 = strel('line', 5, 0);
                    mask = imdilate(mask, [se90 se0]);
                    
                    % Clean-up as necessary
                    stats = sort(cell2mat(struct2cell(regionprops(mask,'Area'))));
                    if length(stats)>1
                        blob = min(stats(end));
                        mask = bwareaopen(mask, blob-1);
                    end
                    
                    % Apply segmentation to GP maps
                    objMembraneGP{ii}=croppedObjectsGP{ii};
                    objMembraneGP{ii}(~mask)=nan;
                    
                    % Apply segmentation to pre-processed images
                    for c=1:nLambda % for each wavelength in c-stack
                        objMembraneThresh{ii}{c} = croppedObjectsThresh{ii}{c};
                        objMembraneThresh{ii}{c}(~mask) = nan;
                    end
                end
            end
            meanObjectGP{i} = meanObjGP;
            stdObjectGP{i} = stdObjGP;
            medObjectGP{i} = medObjGP;
            pxCountObjectGP{i} = pxCountObjGP;
            membraneDiameters{i} = objProps';
            membraneCentroids{i} = objCentroids;
        end
        if ismember(1,options) || ismember(3,options) % Segmenting membranes or spherical objects
            try
                % Membrane segmentation verification
                membraneSegmentationGUI(objMembraneGP, objMembraneThresh, objProps);
                uiwait % Wait until segmentation verification is complete
                
                % Update indeces to keep labelling consistent
                objIndeces = 1:nObj;
                objIndeces(imgsRemoved) = [];
                nMembranes = numel(objMembraneGP);
                [meanMemGP, stdMemGP, medMemGP, pxCountMemGP, objNames] = deal(cell(nMembranes,1));
                for ii=1:nMembranes
                    close all
                    
                    % Object tag for saving
                    objN = num2str(objIndeces(ii));
                    imageName = fileList{i};
                    fullObjTag = [imageName ' object' objN];
                    imgMembraneNames{imgMembraneNameInd,1} = fullObjTag;
                    imgMembraneNameInd = imgMembraneNameInd+1;
                    objTag = ['-obj' objN];
                    
                    % Initialize temporary variables for saving
                    membraneGP = objMembraneGP{ii};
                    membraneThresh = objMembraneThresh{ii};
                    
                    % Generate and save pseudocoloured membrane GP maps
                    figure()
                    sc(membraneGP, 'jet','k') % membrane GP map
                    GPmapSettings();
                      
                    meanMemGP{ii} = nanmean(membraneGP(:));
                    stdMemGP{ii} = nanstd(membraneGP(:));
                    medMemGP{ii} = nanmedian(membraneGP(:));
                    pxCountMemGP{ii} = numel(membraneGP(:));
                
                    if ismember(1,options) % membranes
                        save([GPmembranesDir filesep imageTag objTag '-GP'], 'membraneGP') % membrane GP data
                        savefig(gcf,[GPmembranesDir filesep imageTag objTag '-GP.fig']) % membrane GP map
                    else % spherical objects
                        save([GPcirclesDir filesep imageTag objTag '-GP'], 'membraneGP') % spherical object GP data
                        savefig(gcf,[GPcirclesDir filesep imageTag objTag '-GP.fig']) % spherical object GP map
                    end
                    
                    % Store membrane histogram and spectrum figure and data
                    cSliceMembrane = zeros(nLambda,1);
                    for c=1:nLambda % for each wavelength
                        cSliceMembrane(c) = nansum(membraneThresh{c}(:))/numel(~isnan(membraneThresh{c}(:))); % mean of c-slice in c-stack
                    end
                    figure,plot(wavelengths, cSliceMembrane(:), '-m.', 'MarkerSize', 30, 'LineWidth', 3, 'Color', 'b', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
                    spectraSettings();
                    if ismember(1,options) % membranes
                        save([spectraMembranesDir filesep imageTag objTag '-spectrum'], 'cSliceMembrane') % save membrane spectral data
                        savefig([spectraMembranesDir filesep imageTag objTag '-spectrum.fig']) % Save membrane spectra figure
                    else % spherical objects
                        save([spectraCirclesDir filesep imageTag objTag '-spectrum'], 'cSliceMembrane') % save spherical object spectral data
                        savefig([spectraCirclesDir filesep imageTag objTag '-spectrum.fig']) % Save spherical object spectra figure
                    end
                    histogramData = histogramSettings(membraneGP, -1:0.01:1, twoTerm);
                    xlim([-1 1])
                    a1 = xlabel('GP', 'FontSize', 24);
                    a2 = ylabel('N_{pixels}', 'FontSize', 24);
                    set([a1 a2], 'interpreter', 'tex');
                    if ismember(1,options) % membranes
                        save([histogramMembranesDir filesep imageTag objTag '-histogramData'], 'histogramData') % Save histogram data and fit parameters
                        savefig(gcf,[histogramMembranesDir filesep imageTag objTag '-histogram.fig']) % Save histogram figure
                    else % spherical objects
                        save([histogramCirclesDir filesep imageTag objTag '-histogramData'], 'histogramData') % Save histogram data and fit parameters
                        savefig(gcf,[histogramCirclesDir filesep imageTag objTag '-histogram.fig']) % Save histogram figure
                    end
                    close all
                end
                meanMembraneGP{i} = meanMemGP;
                stdMembraneGP{i} = stdMemGP;
                medMembraneGP{i} = medMemGP;
                pxCountMembraneGP{i} = pxCountMemGP;
            catch err
                disp('Membrane segmentation failed, likely due to poor image quality. If z-stack this could be an out-of-focus slice.')
                warning(err.identifier, err.message);
            end
        end
    end
    GPfull{i} = GPslices{1};
    if stackSizeZ > 1 % z-stack 3D reconstruction
        try
            disp('3D reconstruction underway...')
            close all
            
            % Plot each image as a layer, indicated by color with each pixel as a separate patch
            figure;
            set(gca, 'XLim', [0 512]);
            set(gca, 'YLim', [0 512]);
            hold on;
            GP3D = cat(3,GPslices{:}); % to 3D array
            gridX = 1:size(GP3D,1);
            gridY = 1:size(GP3D,2);
            gridZ = 1:(size(GP3D,3));
            gridX = gridX*voxelSizeX.doubleValue();
            gridY = gridY*voxelSizeY.doubleValue();
            gridZ = gridZ*voxelSizeZ.doubleValue();
            colorbar('vertical');
            clim = [-1 1];
            hpat = PATCH_3Darray(GP3D,gridX,gridY,gridZ,clim,'col');
            view(3)
            colormap(jet(512));
            xlabel('\mum')
            ylabel('\mum')
            zlabel('\mum')
            set(get(gca,'xlabel'),'rotation',25);
            set(get(gca,'ylabel'),'rotation',-47.5);
            set(gca, 'FontSize',24);
            disp('3D reconstruction complete...')
            savefig(gcf,[GPDir filesep imageTag '-3Dreconstruction.fig']) % save 3D reconstruction figure
        catch err
            disp('3D reconstruction failed...')
            warning(err.identifier, err.message);
        end
    end
    i=i+1;
end

% Save wavelengths in pre-segmentation results spectra directory
save([spectraDir filesep 'wavelengths'], 'wavelengthsAll') % Save wavelength labels

% Save mean, std, median, and pixel count at all levels analyzed
save([GPDir filesep 'mean GP'], 'meanGP') % mean GP for pre-segmentation images
save([GPDir filesep 'std GP'], 'stdGP') % std GP for pre-segmentation images
save([GPDir filesep 'median GP'], 'medGP') % median GP for pre-segmentation images
save([GPDir filesep 'pixel count GP'], 'pxCountGP') % number of pixels in GP calculations for pre-segmentation images

if ismember(1,options) || ismember(2,options) % Segmenting objects    
    save([GPobjectsDir filesep 'mean GP'], 'meanObjectGP') % mean GP for segmented objects
    save([GPobjectsDir filesep 'std GP'], 'stdObjectGP') % std GP for segmented objects
    save([GPobjectsDir filesep 'median GP'], 'medObjectGP') % median GP for segmented objects
    save([GPobjectsDir filesep 'pixel count GP'], 'pxCountObjectGP') % number of pixels in GP calculations for segmented objects
    
    % For output to excel 
    a = cell(1,1);
    colNames = {'Filename - Object ID', 'Mean GP', 'Standard deviation GP', 'Median GP', 'Number of Pixels', 'Equivalent diameter (microns)'}';
    a(1,1:6) = colNames;
    a(2:numel(imgObjNames)+1,1) = imgObjNames;
    meanObjectGPCol = cell2mat(cellfun(@cell2mat, meanObjectGP, 'UniformOutput', false));
    a(2:length(meanObjectGPCol)+1,2) = num2cell(meanObjectGPCol);
    stdObjectGPCol = cell2mat(cellfun(@cell2mat, stdObjectGP, 'UniformOutput', false));
    a(2:length(stdObjectGPCol)+1,3) = num2cell(stdObjectGPCol);
    medObjectGPCol = cell2mat(cellfun(@cell2mat, medObjectGP, 'UniformOutput', false));
    a(2:length(medObjectGPCol)+1,4) = num2cell(medObjectGPCol);
    pxCountObjectGPCol = cell2mat(cellfun(@cell2mat, pxCountObjectGP, 'UniformOutput', false));
    a(2:length(pxCountObjectGPCol)+1,5) = num2cell(pxCountObjectGPCol);
    membraneDiametersCol = cell2mat(cellfun(@cell2mat, membraneDiameters, 'UniformOutput', false));
    a(2:length(membraneDiametersCol)+1,6) = num2cell(membraneDiametersCol);
    xlswrite([parentDir filesep 'Results_summary_objects.xlsx'],a)
end

if ismember(1,options) % Segmenting membranes 
    save([GPmembranesDir filesep 'mean GP'], 'meanMembraneGP') % mean GP for segmented membranes
    save([GPmembranesDir filesep 'std GP'], 'stdMembraneGP') % std GP for segmented membranes
    save([GPmembranesDir filesep 'median GP'], 'medMembraneGP') % median GP for segmented membranes
    save([GPmembranesDir filesep 'pixel count GP'], 'pxCountMembraneGP') % mean GP for segmented membranes
    save([GPobjectsDir filesep 'diameters'], 'membraneDiameters') % Object properties
    save([GPobjectsDir filesep 'imgObjNames'], 'imgObjNames') % Object names
    
    % For output to excel
    a = cell(1,1);
    colNames = {'Filename - Object ID', 'Mean GP', 'Standard deviation GP', 'Median GP', 'Number of Pixels', 'Equivalent diameter (microns)'}';
    a(1,1:6) = colNames;
    a(2:numel(imgMembraneNames)+1,1) = imgMembraneNames;
    meanMembraneGPCol = cell2mat(cellfun(@cell2mat, meanMembraneGP, 'UniformOutput', false));
    a(2:length(meanMembraneGPCol)+1,2) = num2cell(meanMembraneGPCol);
    stdMembraneGPCol = cell2mat(cellfun(@cell2mat, stdMembraneGP, 'UniformOutput', false));
    a(2:length(stdMembraneGPCol)+1,3) = num2cell(stdMembraneGPCol);
    medMembraneGPCol = cell2mat(cellfun(@cell2mat, medMembraneGP, 'UniformOutput', false));
    a(2:length(medMembraneGPCol)+1,4) = num2cell(medMembraneGPCol);
    pxCountMembraneGPCol = cell2mat(cellfun(@cell2mat, pxCountMembraneGP, 'UniformOutput', false));
    a(2:length(pxCountMembraneGPCol)+1,5) = num2cell(pxCountMembraneGPCol);
    membraneDiametersCol = cell2mat(cellfun(@cell2mat, membraneDiameters, 'UniformOutput', false));
    a(2:length(membraneDiametersCol)+1,6) = num2cell(membraneDiametersCol);
   
    xlswrite([parentDir filesep 'Results_summary_membranes.xlsx'],a)
end

if ismember(3,options) % Segmenting spherical objects    
    diameters = cell2mat(diameters);
    sizeHistogramData = histogramSettings(diameters, 0:1:ceil(max(diameters)), twoTerm);
    xlim([0 ceil(max(diameters))])
    a1 = xlabel('Diameter (pixels)', 'FontSize', 24);
    a2 = ylabel('Microbubble number', 'FontSize', 24);
    save([sizeDistributionCirclesDir filesep imageTag '-diametersInPixels'], 'diameters') % histogram data
    save([sizeDistributionCirclesDir filesep imageTag '-histogramData'], 'sizeHistogramData') % histogram data
    savefig(gcf,[sizeDistributionCirclesDir filesep imageTag '-sizeHistogram.fig']) % save histogram figure
    save([GPcirclesDir filesep 'mean GP'], 'meanObjectGP') % mean GP for spherical objects
    save([GPcirclesDir filesep 'std GP'], 'stdObjectGP') % std GP for spherical objects
    save([GPcirclesDir filesep 'median GP'], 'medObjectGP') % median GP for spherical objects
    save([GPcirclesDir filesep 'pixel count GP'], 'pxCountObjectGP') % mean GP for spherical objects
    
    % For output to excel
    a = cell(1,1);
    colNames = {'Filename - Object ID', 'Mean GP', 'Standard deviation GP', 'Median GP', 'Number of Pixels', 'Diameter (microns)'}';
    a(1,1:6) = colNames;
    a(2:numel(imgMembraneNames)+1,1) = imgMembraneNames;
    meanObjectGPCol = cell2mat(cellfun(@cell2mat, meanObjectGP, 'UniformOutput', false));
    a(2:length(meanObjectGPCol)+1,2) = num2cell(meanObjectGPCol);
    stdObjectGPCol = cell2mat(cellfun(@cell2mat, stdObjectGP, 'UniformOutput', false));
    a(2:length(stdObjectGPCol)+1,3) = num2cell(stdObjectGPCol);
    medObjectGPCol = cell2mat(cellfun(@cell2mat, medObjectGP, 'UniformOutput', false));
    a(2:length(medObjectGPCol)+1,4) = num2cell(medObjectGPCol);
    pxCountObjectGPCol = cell2mat(cellfun(@cell2mat, pxCountObjectGP, 'UniformOutput', false));
    a(2:length(pxCountObjectGPCol)+1,5) = num2cell(pxCountObjectGPCol);
    a(2:length(diameters(:))+1,6) = num2cell(diameters(:));
    xlswrite([parentDir filesep 'Results_summary_spherical_objects.xlsx'],a)
end