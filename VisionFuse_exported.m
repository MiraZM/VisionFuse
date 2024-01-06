classdef VisionFuse_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ImageRegistratorLabel           matlab.ui.control.Label
        TabGroup                        matlab.ui.container.TabGroup
        IMG1Tab                         matlab.ui.container.Tab
        ClearButton                     matlab.ui.control.Button
        ArtifactRemovalPanel            matlab.ui.container.Panel
        ArtifactRemovalButton           matlab.ui.control.Button
        ImageResamplingPanel            matlab.ui.container.Panel
        AngleSlider                     matlab.ui.control.Slider
        AngleSliderLabel                matlab.ui.control.Label
        ResampleSlider                  matlab.ui.control.Slider
        ResampleSliderLabel             matlab.ui.control.Label
        IntensityNormalizationPanel     matlab.ui.container.Panel
        ContrastStrechingDropDown       matlab.ui.control.DropDown
        ContrastStrechingDropDownLabel  matlab.ui.control.Label
        HistogramMatchingButton         matlab.ui.control.Button
        NoiseReductionPanel             matlab.ui.container.Panel
        NonLocSlider                    matlab.ui.control.Slider
        NonLocSliderLabel               matlab.ui.control.Label
        MedianSlider                    matlab.ui.control.Slider
        MedianSliderLabel               matlab.ui.control.Label
        GaussianSlider                  matlab.ui.control.Slider
        GaussianSliderLabel             matlab.ui.control.Label
        Import1Button                   matlab.ui.control.Button
        IMG2Tab                         matlab.ui.container.Tab
        ClearButton_2                   matlab.ui.control.Button
        ArtifactRemovalPanel_2          matlab.ui.container.Panel
        ArtifactRemovalButton_2         matlab.ui.control.Button
        ImageResamplingPanel_2          matlab.ui.container.Panel
        AngleSlider_2                   matlab.ui.control.Slider
        AngleSlider_2Label              matlab.ui.control.Label
        ResampleSlider_2                matlab.ui.control.Slider
        ResampleSlider_2Label           matlab.ui.control.Label
        IntensityNormalizationPanel_2   matlab.ui.container.Panel
        ContrastStrechingDropDown_2     matlab.ui.control.DropDown
        ContrastStrechingDropDown_2Label  matlab.ui.control.Label
        HistogramMatchingButton_2       matlab.ui.control.Button
        NoiseReductionPanel_2           matlab.ui.container.Panel
        NonLocSlider_2                  matlab.ui.control.Slider
        NonLocSlider_2Label             matlab.ui.control.Label
        MedianSlider_2                  matlab.ui.control.Slider
        MedianSlider_2Label             matlab.ui.control.Label
        GaussianSlider_2                matlab.ui.control.Slider
        GaussianSlider_2Label           matlab.ui.control.Label
        Import2Button_2                 matlab.ui.control.Button
        ExtrensicTab                    matlab.ui.container.Tab
        Resize                          matlab.ui.control.Button
        ApplyExtButton                  matlab.ui.control.Button
        TranslationPanel                matlab.ui.container.Panel
        YaxisSlider                     matlab.ui.control.Slider
        AngelLabel_2                    matlab.ui.control.Label
        XaxisSlider                     matlab.ui.control.Slider
        AngelLabel                      matlab.ui.control.Label
        RotationPanel                   matlab.ui.container.Panel
        BoundingboxDropDown             matlab.ui.control.DropDown
        BoundingboxDropDownLabel        matlab.ui.control.Label
        MethodDropDown                  matlab.ui.control.DropDown
        MethodDropDownLabel             matlab.ui.control.Label
        AngelSlider                     matlab.ui.control.Slider
        AngelSliderLabel                matlab.ui.control.Label
        FeatureBasedMethodsTab          matlab.ui.container.Tab
        DropDown                        matlab.ui.control.DropDown
        DropDownLabel                   matlab.ui.control.Label
        Slider                          matlab.ui.control.Slider
        SliderLabel                     matlab.ui.control.Label
        DetectORBFeaturesButton         matlab.ui.control.Button
        SIFTButton                      matlab.ui.control.Button
        SURFButton                      matlab.ui.control.Button
        IntensityBasedMethodsTab        matlab.ui.container.Tab
        MMIButton                       matlab.ui.control.Button
        MIPanel                         matlab.ui.container.Panel
        EpsilonSliderLabel_3            matlab.ui.control.Label
        EpsilonSliderLabel_2            matlab.ui.control.Label
        EpsilonSlider                   matlab.ui.control.Slider
        EpsilonSliderLabel              matlab.ui.control.Label
        MaxIterationsSlider             matlab.ui.control.Slider
        MaxIterationsSliderLabel        matlab.ui.control.Label
        GrowthFactorSlider              matlab.ui.control.Slider
        InitialRadiusSlider             matlab.ui.control.Slider
        MIButton                        matlab.ui.control.Button
        Switch                          matlab.ui.control.Switch
        PhasecorrelationLabel           matlab.ui.control.Label
        NCCButton                       matlab.ui.control.Button
        UIAxes4                         matlab.ui.control.UIAxes
        UIAxes3                         matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        I;
        I2;
        Image1;
        Image2;
        matchedPoints1;
        matchedPoints2;


    end % Description
    
    

    % Callbacks that handle component events
    methods (Access = private)

        % Callback function
        function SURFButtonPushed(app, event)

            sizeImg1 = size(app.Image1);
            sizeImg2 = size(app.Image2);

            % Choose the smaller dimension as the target size
            targetSize = min([sizeImg1; sizeImg2]);

            % Resize both images to the target size
            img1 = im2single(imresize(rgb2gray(app.Image1), targetSize(1:2)));
            img2 = im2single(imresize(rgb2gray(app.Image2), targetSize(1:2)));

            % Detect SURF features
            pointsA = detectSURFFeatures(img1);
            pointsB = detectSURFFeatures(img2);


            p1 = pointsA.selectStrongest(10);

            % Display the original image in the specified axes
            imshow(img1, 'Parent', app.UIAxes);

            % Hold on to the current axes to plot over the image
            hold(app.UIAxes, 'on');

            % Plot the feature points on the image
            plot(app.UIAxes, p1.Location(:,1), p1.Location(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 10);

            % Release the hold off the axes
            hold(app.UIAxes, 'off');

            p2 = pointsB.selectStrongest(10);

            % Display the original image in the specified axes
            imshow(img2, 'Parent', app.UIAxes2);

            % Hold on to the current axes to plot over the image
            hold(app.UIAxes2, 'on');

            % Plot the feature points on the image
            plot(app.UIAxes2, p1.Location(:,1), p1.Location(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 10);

            % Release the hold off the axes
            hold(app.UIAxes2, 'off');


            % Extract features
            [featuresA, pointsA] = extractFeatures(img1, pointsA);
            [featuresB, pointsB] = extractFeatures(img2, pointsB);


            % Match features
            indexPairs = matchFeatures(featuresA, featuresB, 'Method', 'Threshold');
            matchedPointsA = pointsA(indexPairs(:, 1), :);
            matchedPointsB = pointsB(indexPairs(:, 2), :);

            % Show matched features in a montage
            axes(app.UIAxes5);
            showMatchedFeatures(img1, img2, matchedPointsA, matchedPointsB, 'montage', 'Parent', app.UIAxes5);



            % Match features
            indexPairs = matchFeatures(featuresA, featuresB, 'Unique', true);
            matchedPointsA = pointsA(indexPairs(:, 1), :);
            matchedPointsB = pointsB(indexPairs(:, 2), :);

            % Estimate geometric transformation
            [tform, ~, ~] = estimateGeometricTransform(matchedPointsB, matchedPointsA, 'affine');

            % Apply the transformation to align the moving image
            alignedImg2 = imwarp(img2, tform, 'OutputView', imref2d(size(img1)));

            % Display the fixed image
            imshow(img1, 'Parent', app.UIAxes5);
            hold(app.UIAxes5, 'on');

            % Overlay the aligned moving image
            imshow(alignedImg2, 'Parent', app.UIAxes5);

            % Optionally, you can adjust the properties of the second imshow for better visibility
            % For example, setting the 'AlphaData' for transparency
            imshow(alignedImg2, 'Parent', app.UIAxes5, 'AlphaData', 0.5);

            hold(app.UIAxes5, 'off');
        end

        % Callback function
        function SIFTButtonPushed(app, event)
            vlfeatPath = "C:\Users\miraz\OneDrive\Desktop\4th year 1st sem\CV\FINALZ\vlfeat-0.9.21\toolbox";
            run(fullfile(vlfeatPath, 'vl_setup'));


            sizeImg1 = size(app.Image1);
            sizeImg2 = size(app.Image2);

            % Choose the smaller dimension as the target size
            targetSize = min([sizeImg1; sizeImg2]);

            % Resize both images to the target size
            img1 = im2single(imresize(rgb2gray(app.Image1), targetSize(1:2)));
            img2 = im2single(imresize(rgb2gray(app.Image2), targetSize(1:2)));

            % Compute SIFT features
            [frames1, descriptors1] = vl_sift(img1);
            [frames2, descriptors2] = vl_sift(img2);

            % Match features
            [matches, scores] = vl_ubcmatch(descriptors1, descriptors2);

            % Extract matched points
            app.matchedPoints1 = frames1(1:2, matches(1, :))';
            app.matchedPoints2 = frames2(1:2, matches(2, :))';
            % Display the image in UIAxes2
            imshow(app.Image1, 'Parent', app.UIAxes);

            % Hold on to the current axes, so you can plot over the image
            hold(app.UIAxes, 'on');

            % Plot matched feature points on the image
            plot(app.UIAxes, app.matchedPoints1(:,1), app.matchedPoints1(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 5);

            % Release hold off the axes
            hold(app.UIAxes, 'off');


            % Display the image in UIAxes2
            imshow(app.Image2, 'Parent', app.UIAxes2);

            % Hold on to the current axes, so you can plot over the image
            hold(app.UIAxes2, 'on');

            % Plot matched feature points on the image
            plot(app.UIAxes2, app.matchedPoints2(:,1), app.matchedPoints2(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 5);

            % Release hold off the axes
            hold(app.UIAxes2, 'off');

            showMatchedFeatures(app.Image1, app.Image2, app.matchedPoints1, app.matchedPoints2, 'montage', 'Parent', app.UIAxes3);


            % Get the number of features to display from the slider
            numFeatures = round(app.Slider.Value);
            % Estimate the transformation matrix
            [tform, inlierIdx] = estimateGeometricTransform(app.matchedPoints2, app.matchedPoints1, 'similarity');

            % Check if the transformation is close to an identity matrix
            if norm(tform.T - eye(3)) < 0.01
                % The transformation is negligible. Display original images without applying transformation.
                imshow(app.Image1, 'Parent', app.UIAxes4);
            else
                % Apply the transformation to the moving image
                movingReg = imwarp(app.Image2, tform, 'OutputView', imref2d(size(app.Image1)));

                % Display the overlay in UIAxes_4
                imshowpair(app.Image1, movingReg, 'blend', 'Parent', app.UIAxes4);
    end  
        end

        % Callback function
        function SwitchValueChanged(app, event)
            value = app.Switch.Value;
            if value == "On"
                fixed = app.ImageFile;

                moving = app.image2;
                %Estimate the registration required to bring these two images into alignment.
                % imregcorr returns a simtform2d object that defines the transformation.
                tformEstimate = imregcorr(moving,fixed)
                %Apply the estimated geometric transform to the misaligned image.
                Rfixed = imref2d(size(fixed));
                %Specify the OutputView name-value argument to ensure the registered image is the same size as the reference image.
                movingReg = imwarp(moving,tformEstimate,OutputView=Rfixed);
                imshowpair(fixed,movingReg,"montage");
                %View the aligned image overlaid on the original image, using imshowpair.
                % In this view, imshowpair uses color to highlight areas of misalignment.
                imshowpair(fixed,movingReg,"falsecolor");
                % imregcorr can handle rotation and scale distortions well, but not shear distortion.
            end
        end

        % Button pushed function: HistogramMatchingButton
        function HistogramMatchingButtonPushed(app, event)
                    if isempty(app.Image1)
                    % Display an alert message
                    uialert(app.UIFigure, 'No image has been loaded.', 'Image Required');
                    return;
                    end
                
                    % Open a file dialog for the user to select a reference image
                    [file, path] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp;*.tif;*.tiff', 'Image Files'}, 'Select a Reference Image');
                
                    % Check if the user selected a file
                    if isequal(file, 0)
                        return; % User canceled the operation
                    end
                
                    % Load the selected reference image
                    referenceImagePath = fullfile(path, file);
                    referenceImage = imread(referenceImagePath);
                
                    % Display the reference image in the appropriate axes (UIAxes3)
                    imshow(referenceImage, 'Parent', app.UIAxes3);
                
                    % Perform histogram matching between Image1 and the reference image
                    matchedImage = imhistmatch(app.Image1, referenceImage);
                
                    % Display the processed image in the appropriate axes (UIAxes2)
                    imshow(matchedImage, 'Parent', app.UIAxes);
        end

        % Value changed function: GaussianSlider
        function GaussianSliderValueChanged(app, event)
           % Get the current value of the Gaussian blur slider
            sliderValue = app.GaussianSlider.Value;

            % Apply Gaussian blur to Image1
            blurredImage = imgaussfilt(app.Image1, sliderValue);

            % Display the blurred image in the appropriate axes
            imshow(blurredImage, 'Parent', app.UIAxes);
        end

        % Value changed function: MedianSlider
        function MedianSliderValueChanged(app, event)
        sliderValue = round(app.MedianSlider.Value);

        % Check if Image1 is empty
        if isempty(app.Image1)
            % Handle the case where Image1 is not loaded
            % You can display a message or take appropriate action
            disp('Image1 is not loaded.');
            return;
        end
    
        % Apply the median filter to Image1
        if size(app.Image1, 3) == 3 % Check if it's a color image
            % Separate the color channels
            redChannel = app.Image1(:,:,1);
            greenChannel = app.Image1(:,:,2);
            blueChannel = app.Image1(:,:,3);
            
            % Apply the median filter to each channel
            filteredRedChannel = medfilt2(redChannel, [sliderValue sliderValue]);
            filteredGreenChannel = medfilt2(greenChannel, [sliderValue sliderValue]);
            filteredBlueChannel = medfilt2(blueChannel, [sliderValue sliderValue]);
            
            % Combine the filtered channels into the final color image
            filteredImage = cat(3, filteredRedChannel, filteredGreenChannel, filteredBlueChannel);
        else % Grayscale image
            % Apply the median filter directly to the grayscale image
            filteredImage = medfilt2(app.Image1, [sliderValue sliderValue]);
        end
    
        % Display the filtered image in the appropriate axes
        imshow(filteredImage, 'Parent', app.UIAxes);

        end

        % Value changed function: NonLocSlider
        function NonLocSliderValueChanged(app, event)
           if isempty(app.Image1)
                % Handle the case where Image1 is not loaded
                % You can display a message or take appropriate action
                disp('Image1 is not loaded.');
                return;
            end

            % Convert Image1 to LAB color space
            noisyLAB = rgb2lab(app.Image1);

            % Estimate noise level from a region of the image
            roi = [210, 24, 52, 41]; 
            patch = imcrop(noisyLAB, roi);
            patchSq = patch.^2;
            edist = sqrt(sum(patchSq, 3));
            patchSigma = sqrt(var(edist(:)));

            % Degree of smoothing is based on the slider value
            DoS = app.NonLocSlider.Value * patchSigma;

            % Apply non-local means filter
            denoisedLAB = imnlmfilt(noisyLAB, 'DegreeOfSmoothing', DoS);

            % Convert the denoised LAB image back to RGB
            denoisedRGB = lab2rgb(denoisedLAB, 'OutputType', 'uint8');

            % Display the denoised RGB image in the appropriate axes
            imshow(denoisedRGB, 'Parent', app.UIAxes);

        end

        % Value changed function: ContrastStrechingDropDown
        function ContrastStrechingDropDownValueChanged(app, event)
          selectedOption = app.ContrastStrechingDropDown.Value;

        % Check the selected option and set contrast stretching parameters accordingly
        switch selectedOption
            case 'Low Contrast'
                minIntensity = 0.2;
                maxIntensity = 0.8;
            case 'Medium Contrast'
                minIntensity = 0.1;
                maxIntensity = 0.9;
            case 'High Contrast'
                minIntensity = 0;
                maxIntensity = 1;
            otherwise
                % Handle any other selected option or errors
                disp('Invalid contrast stretching option selected.');
                return; % You can choose to return or take other appropriate action
        end
    
        % Assuming you have an image loaded in 'app.Image1', perform contrast stretching
        if isempty(app.Image1)
            % Handle the case where Image1 is not loaded
            % You can display a message or take appropriate action
            disp('Image1 is not loaded.');
            return;
        end
    
        % Perform contrast stretching on 'app.Image1'
        stretchedImage = imadjust(app.Image1, [minIntensity, maxIntensity], []);
    
        % Display the stretched image in the appropriate axes
        imshow(stretchedImage, 'Parent', app.UIAxes);

        end

        % Value changed function: ResampleSlider
        function ResampleSliderValueChanged(app, event)
     
            scalingFactor = app.ResampleSlider.Value;

            % Check if Image1 is empty
            if isempty(app.Image1)
                % Handle the case where Image1 is not loaded
                % You can display a message or take appropriate action
                disp('Image1 is not loaded.');
                return;
            end

            % Resize (resample) Image1 based on the scaling factor
            % You can specify the interpolation method if needed (e.g., 'bilinear')
            resizedImage = imresize(app.Image1, scalingFactor, 'Method', 'bilinear');

            % Display the resized image in the appropriate axes
            imshow(resizedImage, 'Parent', app.UIAxes);
        end

        % Button pushed function: ArtifactRemovalButton
        function ArtifactRemovalButtonPushed(app, event)
              if isempty(app.Image1)
                % Display an alert message
                uialert(app.UIFigure, 'No image has been loaded.', 'Image Required');
                return;
             end
        
            % Apply an artifact removal technique (e.g., Gaussian filtering)
            sigma = 2; % Standard deviation for Gaussian filter
            outputImage = imgaussfilt(app.Image1, sigma);
        
            % Display the processed image in the specified axes (w1)
            imshow(outputImage, 'Parent', app.UIAxes);
        end

        % Button pushed function: DetectORBFeaturesButton
        function DetectORBFeaturesButtonPushed2(app, event)
             if isempty(app.Image1)
                % Handle the case where Image1 is not loaded
                % You can display a message or take appropriate action
                disp('Image1 is not loaded.');
                return;
            end

            % Convert Image1 to grayscale if it's not already
            if size(app.Image1, 3) == 3
                grayImage = rgb2gray(app.Image1);
            else
                grayImage = app.Image1;
            end

            % Detect ORB features in the grayscale image
            orbPoints = detectORBFeatures(grayImage);

            % Visualize the detected ORB features on the image
            % You can customize the visualization as needed
            annotatedImage = insertMarker(app.Image1, orbPoints.Location);

            % Display the annotated image in the appropriate axes
            imshow(annotatedImage, 'Parent', app.UIAxes3);
        end

        % Button pushed function: ApplyExtButton
        function ApplyExtButtonPushed(app, event)
            
            importedImage= app.Image1;
            rotation_angle = app.AngelSlider.Value;
            rotation_method = app.MethodDropDown.Value;
            rotation_BBox = app.BoundingboxDropDown.Value;
            translation_xasix = app.XaxisSlider.Value;
            translation_yasix = app.YaxisSlider.Value;

            % Check if any value is empty and provide default values if needed
            if isempty(rotation_angle)
                rotation_angle = 0;
            end

             if isempty(rotation_method)
                 rotation_method = 'bilinear';
             end

            if isempty(rotation_BBox)
                rotation_BBox = 'loose';
            end

            if isempty(translation_xasix)
                translation_xasix = 0;
            end

            if isempty(translation_yasix)
                translation_yasix = 0;
            end


            Rotval = imrotate(importedImage,rotation_angle,rotation_method,rotation_BBox);


            tranval = imtranslate(Rotval,[translation_xasix,translation_yasix ]);
            app.Image1 = tranval;
            imshow(app.Image1,'Parent',app.UIAxes);
        end

        % Value changed function: MethodDropDown
        function MethodDropDownValueChanged(app, event)
            value = app.MethodDropDown.Value;
            
        end

        % Button pushed function: Resize
        function ResizeButtonPushed(app, event)
          scalingFactor = app.ResampleSlider.Value;

            % Check if Image1 is empty
            if isempty(app.Image1)
                % Handle the case where Image1 is not loaded
                % You can display a message or take appropriate action
                disp('Image1 is not loaded.');
                return;
            end

            % Resize (resample) Image1 based on the scaling factor
            % You can specify the interpolation method if needed (e.g., 'bilinear')
            resizedImage = imresize(app.Image1, scalingFactor, 'Method', 'bilinear');

            % Display the resized image in the appropriate axes
            imshow(resizedImage, 'Parent', app.UIAxes2);
        end

        % Button pushed function: Import1Button
        function Import1ButtonPushed(app, event)
                startingFolder = 'C:\Program Files\MATLAB';
                %addpath(genpath('path_to_nifti_toolbox'));

                if ~exist(startingFolder, 'dir')
                    startingFolder = pwd;
                end
                
                % Get the name of the first file that the user wants to use.
                defaultFileName = fullfile(startingFolder, '*.*');
                [baseFileName1, folder1] = uigetfile(defaultFileName, 'Select the First Image');
                
                if baseFileName1 == 0
                    % User clicked the Cancel button for the first image.
                    return;
                end
                
                fullFileName1 = fullfile(folder1, baseFileName1);
                %if isNifti(fullFileName1)  % You might need to define this function to check NIfTI files
            
                %app.Image1 = niftiread(fullFileName1)
                
                % Select a slice to display
                % slice = app.Image1(:, :, some_slice_number); % Replace 'some_slice_number' with a valid slice number
                % imshow(slice, [], 'parent', app.UIAxes);
                % Check if the file is DICOM
                if isdicom(fullFileName1)
                    app.Image1 = dicomread(fullFileName1);
                else
                    app.Image1 = imread(fullFileName1);
                end
                
                imshow(app.Image1,[], 'parent', app.UIAxes);  % Display the image
    
                 
                % % Open file dialog to select an image including DICOM files
                % [file, path] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.dcm', ...
                %     'Image Files (*.png, *.jpg, *.jpeg, *.tif, *.dcm)'; '*.*', 'All Files (*.*)'}, ...
                %     'Select an Image');
                % 
                % % Check if a file was selected
                % if isequal(file, 0)
                %     disp('File selection canceled');
                %     return;
                % end
                % 
                % % Full path to the selected file
                % fullPath = fullfile(path, file);
                % 
                % % Check the file extension
                % [~, ~, ext] = fileparts(file);
                % 
                % % Read the selected image
                % if strcmpi(ext, '.dcm') % If DICOM file
                %     image = dicomread(fullPath);
                %     % Optional: Convert to uint8 if needed for display
                %     image = uint8(255 * mat2gray(image));
                % else % For other image formats
                %     image = imread(fullPath);
                % end
                % 
                % % Display the image in the appropriate UIAxes
                % imshow(image, 'Parent', app.UIAxes);
                % 
                % % Store the image for use in cross-correlation
                % % You need to ensure that the app has properties image1 and image2 defined
                % if isempty(app.image1)
                %     app.image1 = image;
                % elseif isempty(app.image2)
                %     app.image2 = image;
                % end
        end

        % Callback function
        function Import2ButtonPushed(app, event)
      startingFolder = 'C:\Program Files\MATLAB';

                if ~exist(startingFolder, 'dir')
                    startingFolder = pwd;
                end
                
                % Get the name of the first file that the user wants to use.
                defaultFileName = fullfile(startingFolder, '*.*');
                [baseFileName1, folder1] = uigetfile(defaultFileName, 'Select the First Image');
                
                if baseFileName1 == 0
                    % User clicked the Cancel button for the first image.
                    return;
                end
                
                fullFileName1 = fullfile(folder1, baseFileName1);
                
                % Check if the file is DICOM
                if isdicom(fullFileName1)
                    app.Image2 = dicomread(fullFileName1);
                else
                    app.Image2 = imread(fullFileName1);
                end
                
                imshow(app.Image2,[], 'parent', app.UIAxes2);  % Display the image

        end

        % Value changed function: AngleSlider
        function AngleSliderValueChanged(app, event)
             value = app.AngleSlider.Value;
            

            % Assuming the original image is stored in app.currentImageOriginal
            imOriginal = app.Image1;
        
            % Check if the image is loaded
            if isempty(imOriginal)
                uialert(app.UIAxes, 'Load an image first', 'Error');
                return;
            end
        
            % Create a random affine transformation based on the slider value
            tform = randomAffine2d('Rotation',[-value value]);
            outputView = affineOutputView(size(imOriginal), tform);
        
            % Apply the transformation
            app.Image2 = imwarp(imOriginal, tform, 'OutputView', outputView);
        
            % Display the augmented image
            imshow(app.Image2, 'Parent', app.UIAxes); 
        end

        % Button pushed function: MIButton
        function MIButtonPushed(app, event)
            % Read DICOM images 
            fixed = app.Image1;
            moving = app.Image2;
            
            % Check if images are loaded and non-empty
            if isempty(fixed) || isempty(moving)
                error('One or both images are empty.');
            end
            
            % Convert to grayscale if they are RGB images
            if size(fixed, 3) == 3
                fixed = rgb2gray(fixed);
            end
            if size(moving, 3) == 3
                moving = rgb2gray(moving);
            end
            
            % Ensure images are of type uint8 for consistency
            fixed = im2uint8(fixed);
            moving = im2uint8(moving);
            
            % Resize images if they are too small
            minSize = 5; % Minimum size required for imregister
            if any(size(fixed) < minSize) || any(size(moving) < minSize)
                newSize = max([size(fixed); size(moving); minSize minSize]);
                fixed = imresize(fixed, newSize);
                moving = imresize(moving, newSize);
            end

            % Get parameters from sliders (example: app.InitialRadiusSlider.Value)
            optimizer = registration.optimizer.OnePlusOneEvolutionary;
            optimizer.InitialRadius = app.InitialRadiusSlider.Value;
            optimizer.Epsilon = app.EpsilonSlider.Value;
            optimizer.GrowthFactor = app.GrowthFactorSlider.Value;
            optimizer.MaximumIterations = app.MaxIterationsSlider.Value;
        
            metric = registration.metric.MattesMutualInformation;
           
            % Perform registration
            movingRegistered = imregister(moving, fixed, 'affine', optimizer, metric);
        
            % Display registered images
            imshowpair(fixed, movingRegistered, 'Scaling', 'joint', 'Parent', app.UIAxes3);
        end

        % Value changed function: InitialRadiusSlider
        function InitialRadiusSliderValueChanged(app, event)
            value = app.InitialRadiusSlider.Value;
            
        end

        % Value changed function: GrowthFactorSlider
        function GrowthFactorSliderValueChanged(app, event)
            value = app.GrowthFactorSlider.Value;
            
        end

        % Value changed function: EpsilonSlider
        function EpsilonSliderValueChanged(app, event)
            value = app.EpsilonSlider.Value;
            
        end

        % Value changed function: MaxIterationsSlider
        function MaxIterationsSliderValueChanged(app, event)
            value = app.MaxIterationsSlider.Value;
            
        end

        % Value changed function: Switch
        function SwitchValueChanged2(app, event)
            value = app.Switch.Value;
            fixed = app.Image1;
            moving = app.Image2;
            %Estimate the registration required to bring these two images into alignment.
            % imregcorr returns a simtform2d object that defines the transformation.
            tformEstimate = imregcorr(moving,fixed)
            %Apply the estimated geometric transform to the misaligned image.
            Rfixed = imref2d(size(fixed));
            %Specify the OutputView name-value argument to ensure the registered image is the same size as the reference image.
            movingReg = imwarp(moving,tformEstimate,OutputView=Rfixed);
            %View the aligned image overlaid on the original image, using imshowpair. 
            % In this view, imshowpair uses color to highlight areas of misalignment.
            imshowpair(fixed,movingReg,"falsecolor",'Parent',app.UIAxes3);
            % imregcorr can handle rotation and scale distortions well, but not shear distortion
        end

        % Button pushed function: NCCButton
        function NCCButtonPushed2(app, event)
         
            % Check if both images are loaded
            if isempty(app.Image1) || isempty(app.Image2)
                disp('Both images need to be loaded first.');
                return;
            end

            % Convert images to grayscale if they are RGB
            if size(app.Image1, 3) == 3
                image1_gray = rgb2gray(app.Image1);
            else
                image1_gray = app.Image1;
            end

            if size(app.Image2, 3) == 3
                image2_gray = rgb2gray(app.Image2);
            else
                image2_gray = app.Image2;
            end

            % Check and display sizes for debugging
            disp('Size of image1_gray:');
            disp(size(image1_gray));
            disp('Size of image2_gray:');
            disp(size(image2_gray));

            % Proceed with NCC if sizes are correct
            if all(size(image1_gray) > 0) && all(size(image2_gray) > 0)
                % Perform Normalized Cross-Correlation
                nccResult = normxcorr2(image2_gray, image1_gray);

                % Find the peak location representing the best match
                [ypeak, xpeak] = find(nccResult == max(nccResult(:)));

                % Display the result
                imshow(nccResult, 'Parent', app.UIAxes4); % Assuming UIAxes4 is the axes for displaying NCC result

                % You can also display the peak location or take appropriate action
                % For example, displaying a marker at the peak location
                hold(app.UIAxes4, 'on');
                plot(app.UIAxes4, xpeak, ypeak, 'ro', 'MarkerSize', 10);
                hold(app.UIAxes4, 'off');

            else
                error('One of the images is empty or not correctly loaded.');
            end
        end

        % Value changed function: XaxisSlider
        function XaxisSliderValueChanged(app, event)
            value = app.XaxisSlider.Value;
            
        end

        % Value changed function: YaxisSlider
        function YaxisSliderValueChanged(app, event)
            value = app.YaxisSlider.Value;
            
        end

        % Value changed function: GaussianSlider_2
        function GaussianSlider_2ValueChanged(app, event)
            sliderValue = app.GaussianSlider_2.Value;
           
            % Apply Gaussian blur to Image1
            blurredImage = imgaussfilt(app.Image2, sliderValue);

            % Display the blurred image in the appropriate axes
            imshow(blurredImage, 'Parent', app.UIAxes2);
        end

        % Value changed function: MedianSlider_2
        function MedianSlider_2ValueChanged(app, event)
           
        sliderValue = round(app.MedianSlider_2.Value);

        % Check if Image1 is empty
        if isempty(app.Image2)
            % Handle the case where Image1 is not loaded
            % You can display a message or take appropriate action
            disp('Image1 is not loaded.');
            return;
        end
    
        % Apply the median filter to Image1
        if size(app.Image1, 3) == 3 % Check if it's a color image
            % Separate the color channels
            redChannel = app.Image2(:,:,1);
            greenChannel = app.Image2(:,:,2);
            blueChannel = app.Image2(:,:,3);
            
            % Apply the median filter to each channel
            filteredRedChannel = medfilt2(redChannel, [sliderValue sliderValue]);
            filteredGreenChannel = medfilt2(greenChannel, [sliderValue sliderValue]);
            filteredBlueChannel = medfilt2(blueChannel, [sliderValue sliderValue]);
            
            % Combine the filtered channels into the final color image
            filteredImage = cat(3, filteredRedChannel, filteredGreenChannel, filteredBlueChannel);
        else % Grayscale image
            % Apply the median filter directly to the grayscale image
            filteredImage = medfilt2(app.Image2, [sliderValue sliderValue]);
        end
    
        % Display the filtered image in the appropriate axes
        imshow(filteredImage, 'Parent', app.UIAxes2);
        end

        % Value changed function: NonLocSlider_2
        function NonLocSlider_2ValueChanged(app, event)
            value = app.NonLocSlider_2.Value;
               if isempty(app.Image2)
                % Handle the case where Image1 is not loaded
                % You can display a message or take appropriate action
                disp('Image1 is not loaded.');
                return;
            end

            % Convert Image1 to LAB color space
            noisyLAB = rgb2lab(app.Image2);

            % Estimate noise level from a region of the image
            roi = [210, 24, 52, 41]; 
            patch = imcrop(noisyLAB, roi);
            patchSq = patch.^2;
            edist = sqrt(sum(patchSq, 3));
            patchSigma = sqrt(var(edist(:)));

            % Degree of smoothing is based on the slider value
            DoS = app.NonLocSlider_2.Value * patchSigma;

            % Apply non-local means filter
            denoisedLAB = imnlmfilt(noisyLAB, 'DegreeOfSmoothing', DoS);

            % Convert the denoised LAB image back to RGB
            denoisedRGB = lab2rgb(denoisedLAB, 'OutputType', 'uint8');

            % Display the denoised RGB image in the appropriate axes
            imshow(denoisedRGB, 'Parent', app.UIAxes2);
   

        end

        % Value changed function: ResampleSlider_2
        function ResampleSlider_2ValueChanged(app, event)
           
            scalingFactor = app.ResampleSlider_2.Value;

            % Check if Image1 is empty
            if isempty(app.Image2)
                % Handle the case where Image1 is not loaded
                % You can display a message or take appropriate action
                disp('Image1 is not loaded.');
                return;
            end

            % Resize (resample) Image1 based on the scaling factor
            % You can specify the interpolation method if needed (e.g., 'bilinear')
            resizedImage = imresize(app.Image2, scalingFactor, 'Method', 'bilinear');

            % Display the resized image in the appropriate axes
            imshow(resizedImage, 'Parent', app.UIAxes2);
        end

        % Value changed function: AngleSlider_2
        function AngleSlider_2ValueChanged(app, event)
            value = app.AngleSlider_2.Value;
            
          
            % Assuming the original image is stored in app.currentImageOriginal
            imOriginal = app.Image2;
        
            % Check if the image is loaded
            if isempty(imOriginal)
                uialert(app.UIAxes2, 'Load an image first', 'Error');
                return;
            end
        
            % Create a random affine transformation based on the slider value
            tform = randomAffine2d('Rotation',[-value value]);
            outputView = affineOutputView(size(imOriginal), tform);
        
            % Apply the transformation
            app.Image2 = imwarp(imOriginal, tform, 'OutputView', outputView);
        
            % Display the augmented image
            imshow(app.Image2, 'Parent', app.UIAxes2); 
        end

        % Button pushed function: HistogramMatchingButton_2
        function HistogramMatchingButton_2Pushed(app, event)
             if isempty(app.Image2)
                    % Display an alert message
                    uialert(app.UIFigure, 'No image has been loaded.', 'Image Required');
                    return;
             end
                
                    % Open a file dialog for the user to select a reference image
                    [file, path] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp;*.tif;*.tiff', 'Image Files'}, 'Select a Reference Image');
                
                    % Check if the user selected a file
                    if isequal(file, 0)
                        return; % User canceled the operation
                    end
                
                    % Load the selected reference image
                    referenceImagePath = fullfile(path, file);
                    referenceImage = imread(referenceImagePath);
                
                    % Display the reference image in the appropriate axes (UIAxes3)
                    imshow(referenceImage, 'Parent', app.UIAxes3);
                
                    % Perform histogram matching between Image1 and the reference image
                    matchedImage = imhistmatch(app.Image2, referenceImage);
                
                    % Display the processed image in the appropriate axes (UIAxes2)
                    imshow(matchedImage, 'Parent', app.UIAxes22);
        end

        % Value changed function: ContrastStrechingDropDown_2
        function ContrastStrechingDropDown_2ValueChanged(app, event)
            
        selectedOption = app.ContrastStrechingDropDown_2.Value;

        % Check the selected option and set contrast stretching parameters accordingly
        switch selectedOption
            case 'Low Contrast'
                minIntensity = 0.2;
                maxIntensity = 0.8;
            case 'Medium Contrast'
                minIntensity = 0.1;
                maxIntensity = 0.9;
            case 'High Contrast'
                minIntensity = 0;
                maxIntensity = 1;
            otherwise
                % Handle any other selected option or errors
                disp('Invalid contrast stretching option selected.');
                return; % You can choose to return or take other appropriate action
        end
    
        % Assuming you have an image loaded in 'app.Image1', perform contrast stretching
        if isempty(app.Image2)
            % Handle the case where Image1 is not loaded
            % You can display a message or take appropriate action
            disp('Image2 is not loaded.');
            return;
        end
    
        % Perform contrast stretching on 'app.Image1'
        stretchedImage = imadjust(app.Image2, [minIntensity, maxIntensity], []);
    
        % Display the stretched image in the appropriate axes
        imshow(stretchedImage, 'Parent', app.UIAxes2);
        end

        % Button pushed function: ArtifactRemovalButton_2
        function ArtifactRemovalButton_2Pushed(app, event)
            if isempty(app.Image2)
                % Display an alert message
                uialert(app.UIFigure, 'No image has been loaded.', 'Image Required');
                return;
             end
        
            % Apply an artifact removal technique (e.g., Gaussian filtering)
            sigma = 2; % Standard deviation for Gaussian filter
            outputImage = imgaussfilt(app.Image2, sigma);
        
            % Display the processed image in the specified axes (w1)
            imshow(outputImage, 'Parent', app.UIAxes2);
        end

        % Button pushed function: Import2Button_2
        function Import2Button_2Pushed(app, event)
            startingFolder = 'C:\Program Files\MATLAB';

                if ~exist(startingFolder, 'dir')
                    startingFolder = pwd;
                end
                
                % Get the name of the first file that the user wants to use.
                defaultFileName = fullfile(startingFolder, '*.*');
                [baseFileName1, folder1] = uigetfile(defaultFileName, 'Select the First Image');
                
                if baseFileName1 == 0
                    % User clicked the Cancel button for the first image.
                    return;
                end
                
                fullFileName1 = fullfile(folder1, baseFileName1);
                
                % Check if the file is DICOM
                if isdicom(fullFileName1)
                    app.Image2 = dicomread(fullFileName1);
                else
                    app.Image2 = imread(fullFileName1);
                end
                
                imshow(app.Image2,[], 'parent', app.UIAxes2);  % Display the image
        end

        % Button pushed function: MMIButton
        function MMIButtonPushed(app, event)
            fixed = app.Image1;
            moving = app.Image2;
            
            % Determine a common size
            minSize = 5; % Minimum size for each dimension
            
            % Determine the target size
            targetSize = max([size(fixed); size(moving); repmat(minSize, 1, ndims(fixed))], [], 1);
            
            % Resize the fixed image
            fixed = imresize(fixed, targetSize(1:2));
            
            % Resize the moving image
            moving = imresize(moving, targetSize(1:2));
            
            % Convert images to appropriate type if needed
            fixed = im2double(fixed);
            moving = im2double(moving);
            
            % Debugging: Check the size of images
            disp(['Size of fixed image: ', mat2str(size(fixed))]);
            disp(['Size of moving image: ', mat2str(size(moving))]);
            
            % Continue with the registration process
            [optimizer, metric] = imregconfig("multimodal");
            registeredDefault = imregister(moving, fixed, "affine", optimizer, metric);
            imshowpair(registeredDefault, fixed, "falsecolor", 'Parent', app.UIAxes3);
        end

        % Button pushed function: SURFButton
        function SURFButtonPushed2(app, event)
sizeImg1 = size(app.Image1);
sizeImg2 = size(app.Image2);

% Choose the smaller dimension as the target size
targetSize = min([sizeImg1(1:2); sizeImg2(1:2)]);

% Check if Image1 is RGB and convert to grayscale if it is
if size(app.Image1, 3) == 3
    img1 = im2single(imresize(rgb2gray(app.Image1), targetSize));
else
    img1 = im2single(imresize(app.Image1, targetSize));
end

% Check if Image2 is RGB and convert to grayscale if it is
if size(app.Image2, 3) == 3
    img2 = im2single(imresize(rgb2gray(app.Image2), targetSize));
else
    img2 = im2single(imresize(app.Image2, targetSize));
end
 original = img1;
 distorted = img2;
% 
% % Detect SURF features
 pointsA = detectSURFFeatures(original);
 pointsB = detectSURFFeatures(distorted);
% 
% 
 p1 = pointsA.selectStrongest(10);
% 
% % Display the original image in the specified axes
 imshow(app.Image1, 'Parent', app.UIAxes);
% 
% % Hold on to the current axes to plot over the image
 hold(app.UIAxes, 'on');
% 
% % Plot the feature points on the image
 plot(app.UIAxes, p1.Location(:,1), p1.Location(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 10);
% 
% % Release the hold off the axes
 hold(app.UIAxes, 'off');
% 
 p2 = pointsB.selectStrongest(10);
% 
% % Display the original image in the specified axes
 imshow(app.Image2, 'Parent', app.UIAxes2);
 % Hold on to the current axes to plot over the image
 hold(app.UIAxes2, 'on');
% Plot the feature points on the image
 plot(app.UIAxes2, p1.Location(:,1), p1.Location(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 10);
 % Release the hold off the axes
 hold(app.UIAxes2, 'off');


ptsOriginal  = detectSURFFeatures(original);
ptsDistorted = detectSURFFeatures(distorted);

[featuresOriginal,validPtsOriginal] = extractFeatures(original,ptsOriginal);
[featuresDistorted,validPtsDistorted] = extractFeatures(distorted,ptsDistorted);


% Detect SURF features
ptsOriginal = detectSURFFeatures(original);


% Plot the detected features
plot(validPtsOriginal, 'showOrientation', true);
title('Detected SURF Features');


index_pairs = matchFeatures(featuresOriginal,featuresDistorted);
matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));

showMatchedFeatures(original, distorted, matchedPtsOriginal, matchedPtsDistorted, 'montage', 'Parent', app.UIAxes3);

    transformationType = app.DropDown.Value;



[tform,inlierIdx] = estimateGeometricTransform2D(matchedPtsDistorted,matchedPtsOriginal,transformationType);
inlierPtsDistorted = matchedPtsDistorted(inlierIdx,:);
inlierPtsOriginal  = matchedPtsOriginal(inlierIdx,:);
% axes(app.UIAxes4); 
% showMatchedFeatures(original, distorted, inlierPtsDistorted, inlierPtsOriginal, 'montage', 'Parent', app.UIAxes3);


 % Check if the transformation is close to an identity matrix
    if norm(tform.T - eye(3)) < 0.01
        % The transformation is negligible. Display original images without applying transformation.
        imshow(original, 'Parent', app.UIAxes4);
    else
        % Apply the transformation to the moving image
        movingReg = imwarp(distorted, tform, 'OutputView', imref2d(size(app.Image1)));

        % Display the overlay in UIAxes_4
        imshowpair(original, movingReg, 'blend', 'Parent', app.UIAxes4);
    end  

        end

        % Button pushed function: SIFTButton
        function SIFTButtonPushed2(app, event)
vlfeatPath = "C:\Users\miraz\OneDrive\Desktop\4th year 1st sem\CV\FINALZ\vlfeat-0.9.21\toolbox"; 
run(fullfile(vlfeatPath, 'vl_setup'));


sizeImg1 = size(app.Image1);
sizeImg2 = size(app.Image2);

% Choose the smaller dimension as the target size
targetSize = min([sizeImg1; sizeImg2]);

% Resize both images to the target size
img1 = im2single(imresize(rgb2gray(app.Image1), targetSize(1:2)));
img2 = im2single(imresize(rgb2gray(app.Image2), targetSize(1:2)));

    % Compute SIFT features
    [frames1, descriptors1] = vl_sift(img1);
    [frames2, descriptors2] = vl_sift(img2);

    % Match features
    [matches, scores] = vl_ubcmatch(descriptors1, descriptors2);

    % Extract matched points
    app.matchedPoints1 = frames1(1:2, matches(1, :))';
    app.matchedPoints2 = frames2(1:2, matches(2, :))';
% Display the image in UIAxes2
imshow(app.Image1, 'Parent', app.UIAxes);

% Hold on to the current axes, so you can plot over the image
hold(app.UIAxes, 'on');

% Plot matched feature points on the image
plot(app.UIAxes, app.matchedPoints1(:,1), app.matchedPoints1(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 5);

% Release hold off the axes
hold(app.UIAxes, 'off');

    
% Display the image in UIAxes2
imshow(app.Image2, 'Parent', app.UIAxes2);

% Hold on to the current axes, so you can plot over the image
hold(app.UIAxes2, 'on');

% Plot matched feature points on the image
plot(app.UIAxes2, app.matchedPoints2(:,1), app.matchedPoints2(:,2), 'y*', 'LineWidth', 1, 'MarkerSize', 5);

% Release hold off the axes
hold(app.UIAxes2, 'off');

     showMatchedFeatures(app.Image1, app.Image2, app.matchedPoints1, app.matchedPoints2, 'montage', 'Parent', app.UIAxes3);

transformationType = app.DropDown.Value;
    % Get the number of features to display from the slider
    numFeatures = round(app.Slider.Value);
% Estimate the transformation matrix
[tform, inlierIdx] = estimateGeometricTransform(app.matchedPoints2, app.matchedPoints1, transformationType);

 % Check if the transformation is close to an identity matrix
    if norm(tform.T - eye(3)) < 0.01
        % The transformation is negligible. Display original images without applying transformation.
        imshow(app.Image1, 'Parent', app.UIAxes4);
    else
        % Apply the transformation to the moving image
        movingReg = imwarp(app.Image2, tform, 'OutputView', imref2d(size(app.Image1)));

        % Display the overlay in UIAxes_4
        imshowpair(app.Image1, movingReg, 'blend', 'Parent', app.UIAxes4);
    end  

        end

        % Value changed function: DropDown
        function DropDownValueChanged(app, event)
            value = app.DropDown.Value;
            
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
 
    % Clear the first axes
    cla(app.UIAxes, 'reset');
    
    % Clear the second axes
    cla(app.UIAxes2, 'reset');

    % Clear the third axes
    cla(app.UIAxes3, 'reset');
    cla(app.UIAxes4, 'reset');

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.902 0.8471 0.8941];
            app.UIFigure.Position = [100 100 1141 726];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'First Image ')
            app.UIAxes.FontName = 'Kristen ITC';
            app.UIAxes.FontWeight = 'bold';
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
            app.UIAxes.GridColor = [1 1 1];
            app.UIAxes.Position = [22 328 310 258];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'Second Image')
            app.UIAxes2.FontName = 'Kristen ITC';
            app.UIAxes2.FontWeight = 'bold';
            app.UIAxes2.XTick = [];
            app.UIAxes2.XTickLabel = '';
            app.UIAxes2.YTick = [];
            app.UIAxes2.YTickLabel = '';
            app.UIAxes2.GridColor = [1 1 1];
            app.UIAxes2.Position = [342 332 309 256];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.UIFigure);
            app.UIAxes3.FontName = 'Kristen ITC';
            app.UIAxes3.FontWeight = 'bold';
            app.UIAxes3.XTick = [];
            app.UIAxes3.YTick = [];
            app.UIAxes3.GridColor = [1 1 1];
            app.UIAxes3.Position = [22 47 308 259];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.UIFigure);
            title(app.UIAxes4, 'Output Image')
            app.UIAxes4.FontName = 'Kristen ITC';
            app.UIAxes4.FontWeight = 'bold';
            app.UIAxes4.XTick = [];
            app.UIAxes4.YTick = [];
            app.UIAxes4.GridColor = [1 1 1];
            app.UIAxes4.Position = [342 40 309 260];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [666 17 456 613];

            % Create IMG1Tab
            app.IMG1Tab = uitab(app.TabGroup);
            app.IMG1Tab.Title = 'IMG1';
            app.IMG1Tab.BackgroundColor = [0.7804 0.7333 0.7804];

            % Create Import1Button
            app.Import1Button = uibutton(app.IMG1Tab, 'push');
            app.Import1Button.ButtonPushedFcn = createCallbackFcn(app, @Import1ButtonPushed, true);
            app.Import1Button.BackgroundColor = [0.8706 0.851 0.8667];
            app.Import1Button.FontName = 'Kristen ITC';
            app.Import1Button.Position = [31 492 178 55];
            app.Import1Button.Text = 'Import1';

            % Create NoiseReductionPanel
            app.NoiseReductionPanel = uipanel(app.IMG1Tab);
            app.NoiseReductionPanel.TitlePosition = 'centertop';
            app.NoiseReductionPanel.Title = 'Noise Reduction';
            app.NoiseReductionPanel.BackgroundColor = [0.8706 0.851 0.8667];
            app.NoiseReductionPanel.FontName = 'Kristen ITC';
            app.NoiseReductionPanel.FontWeight = 'bold';
            app.NoiseReductionPanel.Position = [18 244 211 194];

            % Create GaussianSliderLabel
            app.GaussianSliderLabel = uilabel(app.NoiseReductionPanel);
            app.GaussianSliderLabel.HorizontalAlignment = 'right';
            app.GaussianSliderLabel.FontName = 'Kristen ITC';
            app.GaussianSliderLabel.Position = [4 138 57 22];
            app.GaussianSliderLabel.Text = 'Gaussian';

            % Create GaussianSlider
            app.GaussianSlider = uislider(app.NoiseReductionPanel);
            app.GaussianSlider.Limits = [0 5];
            app.GaussianSlider.ValueChangedFcn = createCallbackFcn(app, @GaussianSliderValueChanged, true);
            app.GaussianSlider.MinorTicks = [0.5 1.5 2.5 3.5 4.5];
            app.GaussianSlider.Position = [81 157 119 3];
            app.GaussianSlider.Value = 1;

            % Create MedianSliderLabel
            app.MedianSliderLabel = uilabel(app.NoiseReductionPanel);
            app.MedianSliderLabel.FontName = 'Kristen ITC';
            app.MedianSliderLabel.Position = [12 82 49 22];
            app.MedianSliderLabel.Text = 'Median';

            % Create MedianSlider
            app.MedianSlider = uislider(app.NoiseReductionPanel);
            app.MedianSlider.Limits = [3 15];
            app.MedianSlider.ValueChangedFcn = createCallbackFcn(app, @MedianSliderValueChanged, true);
            app.MedianSlider.MinorTicks = [];
            app.MedianSlider.Position = [76 101 119 3];
            app.MedianSlider.Value = 3;

            % Create NonLocSliderLabel
            app.NonLocSliderLabel = uilabel(app.NoiseReductionPanel);
            app.NonLocSliderLabel.FontName = 'Kristen ITC';
            app.NonLocSliderLabel.Position = [11 22 57 22];
            app.NonLocSliderLabel.Text = 'Non Loc';

            % Create NonLocSlider
            app.NonLocSlider = uislider(app.NoiseReductionPanel);
            app.NonLocSlider.Limits = [1 7];
            app.NonLocSlider.MajorTicks = [1 2 3 4 5 6 7];
            app.NonLocSlider.ValueChangedFcn = createCallbackFcn(app, @NonLocSliderValueChanged, true);
            app.NonLocSlider.MinorTicks = [];
            app.NonLocSlider.Position = [75 41 120 3];
            app.NonLocSlider.Value = 5;

            % Create IntensityNormalizationPanel
            app.IntensityNormalizationPanel = uipanel(app.IMG1Tab);
            app.IntensityNormalizationPanel.TitlePosition = 'centertop';
            app.IntensityNormalizationPanel.Title = 'Intensity Normalization';
            app.IntensityNormalizationPanel.BackgroundColor = [0.8706 0.851 0.8667];
            app.IntensityNormalizationPanel.FontName = 'Kristen ITC';
            app.IntensityNormalizationPanel.FontWeight = 'bold';
            app.IntensityNormalizationPanel.Position = [18 19 211 178];

            % Create HistogramMatchingButton
            app.HistogramMatchingButton = uibutton(app.IntensityNormalizationPanel, 'push');
            app.HistogramMatchingButton.ButtonPushedFcn = createCallbackFcn(app, @HistogramMatchingButtonPushed, true);
            app.HistogramMatchingButton.FontName = 'Kristen ITC';
            app.HistogramMatchingButton.FontWeight = 'bold';
            app.HistogramMatchingButton.Position = [33 101 150 38];
            app.HistogramMatchingButton.Text = 'Histogram Matching';

            % Create ContrastStrechingDropDownLabel
            app.ContrastStrechingDropDownLabel = uilabel(app.IntensityNormalizationPanel);
            app.ContrastStrechingDropDownLabel.FontName = 'Kristen ITC';
            app.ContrastStrechingDropDownLabel.FontWeight = 'bold';
            app.ContrastStrechingDropDownLabel.Position = [47 69 122 22];
            app.ContrastStrechingDropDownLabel.Text = 'Contrast  Streching';

            % Create ContrastStrechingDropDown
            app.ContrastStrechingDropDown = uidropdown(app.IntensityNormalizationPanel);
            app.ContrastStrechingDropDown.Items = {'Low Contrast', 'Medium Contrast', 'High Contrast', ''};
            app.ContrastStrechingDropDown.ValueChangedFcn = createCallbackFcn(app, @ContrastStrechingDropDownValueChanged, true);
            app.ContrastStrechingDropDown.FontName = 'Kristen ITC';
            app.ContrastStrechingDropDown.Position = [33 20 150 38];
            app.ContrastStrechingDropDown.Value = 'Low Contrast';

            % Create ImageResamplingPanel
            app.ImageResamplingPanel = uipanel(app.IMG1Tab);
            app.ImageResamplingPanel.TitlePosition = 'centertop';
            app.ImageResamplingPanel.Title = 'Image Resampling';
            app.ImageResamplingPanel.BackgroundColor = [0.8706 0.851 0.8667];
            app.ImageResamplingPanel.FontName = 'Kristen ITC';
            app.ImageResamplingPanel.FontWeight = 'bold';
            app.ImageResamplingPanel.Position = [239 245 203 193];

            % Create ResampleSliderLabel
            app.ResampleSliderLabel = uilabel(app.ImageResamplingPanel);
            app.ResampleSliderLabel.FontName = 'Kristen ITC';
            app.ResampleSliderLabel.Position = [68 136 59 22];
            app.ResampleSliderLabel.Text = 'Resample';

            % Create ResampleSlider
            app.ResampleSlider = uislider(app.ImageResamplingPanel);
            app.ResampleSlider.Limits = [0 2];
            app.ResampleSlider.MajorTicks = [0.5 1 1.5 2];
            app.ResampleSlider.ValueChangedFcn = createCallbackFcn(app, @ResampleSliderValueChanged, true);
            app.ResampleSlider.MinorTicks = [];
            app.ResampleSlider.FontName = 'Kristen ITC';
            app.ResampleSlider.Position = [35 124 128 3];

            % Create AngleSliderLabel
            app.AngleSliderLabel = uilabel(app.ImageResamplingPanel);
            app.AngleSliderLabel.HorizontalAlignment = 'right';
            app.AngleSliderLabel.FontName = 'Kristen ITC';
            app.AngleSliderLabel.Position = [76 67 39 22];
            app.AngleSliderLabel.Text = 'Angle';

            % Create AngleSlider
            app.AngleSlider = uislider(app.ImageResamplingPanel);
            app.AngleSlider.ValueChangedFcn = createCallbackFcn(app, @AngleSliderValueChanged, true);
            app.AngleSlider.Position = [35 53 127 3];

            % Create ArtifactRemovalPanel
            app.ArtifactRemovalPanel = uipanel(app.IMG1Tab);
            app.ArtifactRemovalPanel.TitlePosition = 'centertop';
            app.ArtifactRemovalPanel.Title = 'Artifact Removal';
            app.ArtifactRemovalPanel.BackgroundColor = [0.8706 0.851 0.8667];
            app.ArtifactRemovalPanel.FontName = 'Kristen ITC';
            app.ArtifactRemovalPanel.FontWeight = 'bold';
            app.ArtifactRemovalPanel.Position = [239 20 203 174];

            % Create ArtifactRemovalButton
            app.ArtifactRemovalButton = uibutton(app.ArtifactRemovalPanel, 'push');
            app.ArtifactRemovalButton.ButtonPushedFcn = createCallbackFcn(app, @ArtifactRemovalButtonPushed, true);
            app.ArtifactRemovalButton.FontName = 'Kristen ITC';
            app.ArtifactRemovalButton.FontWeight = 'bold';
            app.ArtifactRemovalButton.Position = [35 100 132 37];
            app.ArtifactRemovalButton.Text = 'Artifact Removal';

            % Create ClearButton
            app.ClearButton = uibutton(app.IMG1Tab, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.BackgroundColor = [0.902 0.902 0.902];
            app.ClearButton.FontName = 'Kristen ITC';
            app.ClearButton.Position = [252 492 168 55];
            app.ClearButton.Text = 'Clear';

            % Create IMG2Tab
            app.IMG2Tab = uitab(app.TabGroup);
            app.IMG2Tab.Title = 'IMG2';
            app.IMG2Tab.BackgroundColor = [0.7804 0.7294 0.7804];

            % Create Import2Button_2
            app.Import2Button_2 = uibutton(app.IMG2Tab, 'push');
            app.Import2Button_2.ButtonPushedFcn = createCallbackFcn(app, @Import2Button_2Pushed, true);
            app.Import2Button_2.BackgroundColor = [0.8706 0.851 0.8667];
            app.Import2Button_2.FontName = 'Kristen ITC';
            app.Import2Button_2.Position = [36 492 178 55];
            app.Import2Button_2.Text = 'Import2';

            % Create NoiseReductionPanel_2
            app.NoiseReductionPanel_2 = uipanel(app.IMG2Tab);
            app.NoiseReductionPanel_2.TitlePosition = 'centertop';
            app.NoiseReductionPanel_2.Title = 'Noise Reduction';
            app.NoiseReductionPanel_2.BackgroundColor = [0.8706 0.851 0.8667];
            app.NoiseReductionPanel_2.FontName = 'Kristen ITC';
            app.NoiseReductionPanel_2.FontWeight = 'bold';
            app.NoiseReductionPanel_2.Position = [18 257 211 194];

            % Create GaussianSlider_2Label
            app.GaussianSlider_2Label = uilabel(app.NoiseReductionPanel_2);
            app.GaussianSlider_2Label.HorizontalAlignment = 'right';
            app.GaussianSlider_2Label.FontName = 'Kristen ITC';
            app.GaussianSlider_2Label.Position = [4 138 57 22];
            app.GaussianSlider_2Label.Text = 'Gaussian';

            % Create GaussianSlider_2
            app.GaussianSlider_2 = uislider(app.NoiseReductionPanel_2);
            app.GaussianSlider_2.Limits = [0 5];
            app.GaussianSlider_2.ValueChangedFcn = createCallbackFcn(app, @GaussianSlider_2ValueChanged, true);
            app.GaussianSlider_2.MinorTicks = [0.5 1.5 2.5 3.5 4.5];
            app.GaussianSlider_2.Position = [81 157 119 3];
            app.GaussianSlider_2.Value = 1;

            % Create MedianSlider_2Label
            app.MedianSlider_2Label = uilabel(app.NoiseReductionPanel_2);
            app.MedianSlider_2Label.FontName = 'Kristen ITC';
            app.MedianSlider_2Label.Position = [12 82 49 22];
            app.MedianSlider_2Label.Text = 'Median';

            % Create MedianSlider_2
            app.MedianSlider_2 = uislider(app.NoiseReductionPanel_2);
            app.MedianSlider_2.Limits = [3 15];
            app.MedianSlider_2.ValueChangedFcn = createCallbackFcn(app, @MedianSlider_2ValueChanged, true);
            app.MedianSlider_2.MinorTicks = [];
            app.MedianSlider_2.Position = [76 101 119 3];
            app.MedianSlider_2.Value = 3;

            % Create NonLocSlider_2Label
            app.NonLocSlider_2Label = uilabel(app.NoiseReductionPanel_2);
            app.NonLocSlider_2Label.FontName = 'Kristen ITC';
            app.NonLocSlider_2Label.Position = [11 22 57 22];
            app.NonLocSlider_2Label.Text = 'Non Loc';

            % Create NonLocSlider_2
            app.NonLocSlider_2 = uislider(app.NoiseReductionPanel_2);
            app.NonLocSlider_2.Limits = [1 7];
            app.NonLocSlider_2.MajorTicks = [1 2 3 4 5 6 7];
            app.NonLocSlider_2.ValueChangedFcn = createCallbackFcn(app, @NonLocSlider_2ValueChanged, true);
            app.NonLocSlider_2.MinorTicks = [];
            app.NonLocSlider_2.Position = [75 41 120 3];
            app.NonLocSlider_2.Value = 5;

            % Create IntensityNormalizationPanel_2
            app.IntensityNormalizationPanel_2 = uipanel(app.IMG2Tab);
            app.IntensityNormalizationPanel_2.TitlePosition = 'centertop';
            app.IntensityNormalizationPanel_2.Title = 'Intensity Normalization';
            app.IntensityNormalizationPanel_2.BackgroundColor = [0.8706 0.851 0.8667];
            app.IntensityNormalizationPanel_2.FontName = 'Kristen ITC';
            app.IntensityNormalizationPanel_2.FontWeight = 'bold';
            app.IntensityNormalizationPanel_2.Position = [18 32 211 178];

            % Create HistogramMatchingButton_2
            app.HistogramMatchingButton_2 = uibutton(app.IntensityNormalizationPanel_2, 'push');
            app.HistogramMatchingButton_2.ButtonPushedFcn = createCallbackFcn(app, @HistogramMatchingButton_2Pushed, true);
            app.HistogramMatchingButton_2.FontName = 'Kristen ITC';
            app.HistogramMatchingButton_2.FontWeight = 'bold';
            app.HistogramMatchingButton_2.Position = [32 96 150 38];
            app.HistogramMatchingButton_2.Text = 'Histogram Matching';

            % Create ContrastStrechingDropDown_2Label
            app.ContrastStrechingDropDown_2Label = uilabel(app.IntensityNormalizationPanel_2);
            app.ContrastStrechingDropDown_2Label.FontName = 'Kristen ITC';
            app.ContrastStrechingDropDown_2Label.FontWeight = 'bold';
            app.ContrastStrechingDropDown_2Label.Position = [47 61 122 22];
            app.ContrastStrechingDropDown_2Label.Text = 'Contrast  Streching';

            % Create ContrastStrechingDropDown_2
            app.ContrastStrechingDropDown_2 = uidropdown(app.IntensityNormalizationPanel_2);
            app.ContrastStrechingDropDown_2.Items = {'Low Contrast', 'Medium Contrast', 'High Contrast', ''};
            app.ContrastStrechingDropDown_2.ValueChangedFcn = createCallbackFcn(app, @ContrastStrechingDropDown_2ValueChanged, true);
            app.ContrastStrechingDropDown_2.FontName = 'Kristen ITC';
            app.ContrastStrechingDropDown_2.Position = [33 12 150 38];
            app.ContrastStrechingDropDown_2.Value = 'Low Contrast';

            % Create ImageResamplingPanel_2
            app.ImageResamplingPanel_2 = uipanel(app.IMG2Tab);
            app.ImageResamplingPanel_2.TitlePosition = 'centertop';
            app.ImageResamplingPanel_2.Title = 'Image Resampling';
            app.ImageResamplingPanel_2.BackgroundColor = [0.8706 0.851 0.8667];
            app.ImageResamplingPanel_2.FontName = 'Kristen ITC';
            app.ImageResamplingPanel_2.FontWeight = 'bold';
            app.ImageResamplingPanel_2.Position = [239 258 203 193];

            % Create ResampleSlider_2Label
            app.ResampleSlider_2Label = uilabel(app.ImageResamplingPanel_2);
            app.ResampleSlider_2Label.FontName = 'Kristen ITC';
            app.ResampleSlider_2Label.Position = [68 136 59 22];
            app.ResampleSlider_2Label.Text = 'Resample';

            % Create ResampleSlider_2
            app.ResampleSlider_2 = uislider(app.ImageResamplingPanel_2);
            app.ResampleSlider_2.Limits = [0 2];
            app.ResampleSlider_2.MajorTicks = [0.5 1 1.5 2];
            app.ResampleSlider_2.ValueChangedFcn = createCallbackFcn(app, @ResampleSlider_2ValueChanged, true);
            app.ResampleSlider_2.MinorTicks = [];
            app.ResampleSlider_2.FontName = 'Kristen ITC';
            app.ResampleSlider_2.Position = [35 124 128 3];

            % Create AngleSlider_2Label
            app.AngleSlider_2Label = uilabel(app.ImageResamplingPanel_2);
            app.AngleSlider_2Label.HorizontalAlignment = 'right';
            app.AngleSlider_2Label.FontName = 'Kristen ITC';
            app.AngleSlider_2Label.Position = [76 67 39 22];
            app.AngleSlider_2Label.Text = 'Angle';

            % Create AngleSlider_2
            app.AngleSlider_2 = uislider(app.ImageResamplingPanel_2);
            app.AngleSlider_2.ValueChangedFcn = createCallbackFcn(app, @AngleSlider_2ValueChanged, true);
            app.AngleSlider_2.Position = [35 53 127 3];

            % Create ArtifactRemovalPanel_2
            app.ArtifactRemovalPanel_2 = uipanel(app.IMG2Tab);
            app.ArtifactRemovalPanel_2.TitlePosition = 'centertop';
            app.ArtifactRemovalPanel_2.Title = 'Artifact Removal';
            app.ArtifactRemovalPanel_2.BackgroundColor = [0.8706 0.851 0.8667];
            app.ArtifactRemovalPanel_2.FontName = 'Kristen ITC';
            app.ArtifactRemovalPanel_2.FontWeight = 'bold';
            app.ArtifactRemovalPanel_2.Position = [239 33 203 174];

            % Create ArtifactRemovalButton_2
            app.ArtifactRemovalButton_2 = uibutton(app.ArtifactRemovalPanel_2, 'push');
            app.ArtifactRemovalButton_2.ButtonPushedFcn = createCallbackFcn(app, @ArtifactRemovalButton_2Pushed, true);
            app.ArtifactRemovalButton_2.FontName = 'Kristen ITC';
            app.ArtifactRemovalButton_2.FontWeight = 'bold';
            app.ArtifactRemovalButton_2.Position = [35 100 132 37];
            app.ArtifactRemovalButton_2.Text = 'Artifact Removal';

            % Create ClearButton_2
            app.ClearButton_2 = uibutton(app.IMG2Tab, 'push');
            app.ClearButton_2.BackgroundColor = [0.902 0.902 0.902];
            app.ClearButton_2.FontName = 'Kristen ITC';
            app.ClearButton_2.Position = [254 492 168 55];
            app.ClearButton_2.Text = 'Clear';

            % Create ExtrensicTab
            app.ExtrensicTab = uitab(app.TabGroup);
            app.ExtrensicTab.Title = 'Extrensic';
            app.ExtrensicTab.BackgroundColor = [0.7804 0.7294 0.7804];

            % Create RotationPanel
            app.RotationPanel = uipanel(app.ExtrensicTab);
            app.RotationPanel.Title = 'Rotation';
            app.RotationPanel.FontName = 'Kristen ITC';
            app.RotationPanel.Position = [76 399 268 171];

            % Create AngelSliderLabel
            app.AngelSliderLabel = uilabel(app.RotationPanel);
            app.AngelSliderLabel.HorizontalAlignment = 'right';
            app.AngelSliderLabel.FontName = 'Kristen ITC';
            app.AngelSliderLabel.Position = [5 112 39 22];
            app.AngelSliderLabel.Text = 'Angel';

            % Create AngelSlider
            app.AngelSlider = uislider(app.RotationPanel);
            app.AngelSlider.Limits = [0 180];
            app.AngelSlider.MajorTicks = [0 45 90 135 180];
            app.AngelSlider.Position = [59 131 182 3];

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.RotationPanel);
            app.MethodDropDownLabel.FontName = 'Kristen ITC';
            app.MethodDropDownLabel.Position = [38 62 53 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create MethodDropDown
            app.MethodDropDown = uidropdown(app.RotationPanel);
            app.MethodDropDown.Items = {'', 'Nearest', 'Bilinear', 'Bicubic'};
            app.MethodDropDown.ValueChangedFcn = createCallbackFcn(app, @MethodDropDownValueChanged, true);
            app.MethodDropDown.FontName = 'Kristen ITC';
            app.MethodDropDown.Position = [98 56 114 33];
            app.MethodDropDown.Value = '';

            % Create BoundingboxDropDownLabel
            app.BoundingboxDropDownLabel = uilabel(app.RotationPanel);
            app.BoundingboxDropDownLabel.HorizontalAlignment = 'center';
            app.BoundingboxDropDownLabel.FontName = 'Kristen ITC';
            app.BoundingboxDropDownLabel.Position = [26 10 82 22];
            app.BoundingboxDropDownLabel.Text = 'Boundingbox';

            % Create BoundingboxDropDown
            app.BoundingboxDropDown = uidropdown(app.RotationPanel);
            app.BoundingboxDropDown.Items = {'Loose', 'Crop'};
            app.BoundingboxDropDown.FontName = 'Kristen ITC';
            app.BoundingboxDropDown.Position = [115 4 126 33];
            app.BoundingboxDropDown.Value = 'Loose';

            % Create TranslationPanel
            app.TranslationPanel = uipanel(app.ExtrensicTab);
            app.TranslationPanel.Title = 'Translation';
            app.TranslationPanel.FontName = 'Kristen ITC';
            app.TranslationPanel.Position = [77 205 268 175];

            % Create AngelLabel
            app.AngelLabel = uilabel(app.TranslationPanel);
            app.AngelLabel.HorizontalAlignment = 'right';
            app.AngelLabel.FontName = 'Kristen ITC';
            app.AngelLabel.Position = [7 119 40 22];
            app.AngelLabel.Text = 'x_axis';

            % Create XaxisSlider
            app.XaxisSlider = uislider(app.TranslationPanel);
            app.XaxisSlider.ValueChangedFcn = createCallbackFcn(app, @XaxisSliderValueChanged, true);
            app.XaxisSlider.Position = [68 128 150 3];

            % Create AngelLabel_2
            app.AngelLabel_2 = uilabel(app.TranslationPanel);
            app.AngelLabel_2.HorizontalAlignment = 'right';
            app.AngelLabel_2.FontName = 'Kristen ITC';
            app.AngelLabel_2.Position = [11 59 39 22];
            app.AngelLabel_2.Text = 'y_axis';

            % Create YaxisSlider
            app.YaxisSlider = uislider(app.TranslationPanel);
            app.YaxisSlider.ValueChangedFcn = createCallbackFcn(app, @YaxisSliderValueChanged, true);
            app.YaxisSlider.Position = [71 68 150 3];

            % Create ApplyExtButton
            app.ApplyExtButton = uibutton(app.ExtrensicTab, 'push');
            app.ApplyExtButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyExtButtonPushed, true);
            app.ApplyExtButton.IconAlignment = 'center';
            app.ApplyExtButton.BackgroundColor = [0.8 0.8 0.8];
            app.ApplyExtButton.FontName = 'Kristen ITC';
            app.ApplyExtButton.FontSize = 16;
            app.ApplyExtButton.FontWeight = 'bold';
            app.ApplyExtButton.FontColor = [1 1 1];
            app.ApplyExtButton.Position = [81 153 269 29];
            app.ApplyExtButton.Text = 'ApplyExt';

            % Create Resize
            app.Resize = uibutton(app.ExtrensicTab, 'push');
            app.Resize.ButtonPushedFcn = createCallbackFcn(app, @ResizeButtonPushed, true);
            app.Resize.IconAlignment = 'center';
            app.Resize.BackgroundColor = [0.8 0.8 0.8];
            app.Resize.FontName = 'Kristen ITC';
            app.Resize.FontSize = 16;
            app.Resize.FontWeight = 'bold';
            app.Resize.FontColor = [1 1 1];
            app.Resize.Position = [81 101 269 29];
            app.Resize.Text = 'Resize';

            % Create FeatureBasedMethodsTab
            app.FeatureBasedMethodsTab = uitab(app.TabGroup);
            app.FeatureBasedMethodsTab.Title = 'Feature Based Methods';
            app.FeatureBasedMethodsTab.BackgroundColor = [0.7804 0.7333 0.7804];

            % Create SURFButton
            app.SURFButton = uibutton(app.FeatureBasedMethodsTab, 'push');
            app.SURFButton.ButtonPushedFcn = createCallbackFcn(app, @SURFButtonPushed2, true);
            app.SURFButton.FontName = 'Kristen ITC';
            app.SURFButton.Position = [89 473 257 74];
            app.SURFButton.Text = 'SURF';

            % Create SIFTButton
            app.SIFTButton = uibutton(app.FeatureBasedMethodsTab, 'push');
            app.SIFTButton.ButtonPushedFcn = createCallbackFcn(app, @SIFTButtonPushed2, true);
            app.SIFTButton.FontName = 'Kristen ITC';
            app.SIFTButton.Position = [89 358 257 74];
            app.SIFTButton.Text = 'SIFT';

            % Create DetectORBFeaturesButton
            app.DetectORBFeaturesButton = uibutton(app.FeatureBasedMethodsTab, 'push');
            app.DetectORBFeaturesButton.ButtonPushedFcn = createCallbackFcn(app, @DetectORBFeaturesButtonPushed2, true);
            app.DetectORBFeaturesButton.FontName = 'Kristen ITC';
            app.DetectORBFeaturesButton.Position = [89 218 257 74];
            app.DetectORBFeaturesButton.Text = 'DetectORBFeatures';

            % Create SliderLabel
            app.SliderLabel = uilabel(app.FeatureBasedMethodsTab);
            app.SliderLabel.HorizontalAlignment = 'right';
            app.SliderLabel.Position = [37 175 36 22];
            app.SliderLabel.Text = 'Slider';

            % Create Slider
            app.Slider = uislider(app.FeatureBasedMethodsTab);
            app.Slider.Position = [94 184 184 3];

            % Create DropDownLabel
            app.DropDownLabel = uilabel(app.FeatureBasedMethodsTab);
            app.DropDownLabel.HorizontalAlignment = 'right';
            app.DropDownLabel.Position = [51 100 65 22];
            app.DropDownLabel.Text = 'Drop Down';

            % Create DropDown
            app.DropDown = uidropdown(app.FeatureBasedMethodsTab);
            app.DropDown.Items = {'affine', 'similarity', 'rigid', 'projective'};
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
            app.DropDown.Position = [131 93 99 35];
            app.DropDown.Value = 'affine';

            % Create IntensityBasedMethodsTab
            app.IntensityBasedMethodsTab = uitab(app.TabGroup);
            app.IntensityBasedMethodsTab.Title = 'Intensity Based Methods';
            app.IntensityBasedMethodsTab.BackgroundColor = [0.7804 0.7333 0.7804];

            % Create NCCButton
            app.NCCButton = uibutton(app.IntensityBasedMethodsTab, 'push');
            app.NCCButton.ButtonPushedFcn = createCallbackFcn(app, @NCCButtonPushed2, true);
            app.NCCButton.BackgroundColor = [0.8392 0.8196 0.8392];
            app.NCCButton.FontName = 'Kristen ITC';
            app.NCCButton.Position = [83 19 150 36];
            app.NCCButton.Text = 'NCC';

            % Create PhasecorrelationLabel
            app.PhasecorrelationLabel = uilabel(app.IntensityBasedMethodsTab);
            app.PhasecorrelationLabel.FontName = 'Kristen ITC';
            app.PhasecorrelationLabel.FontSize = 18;
            app.PhasecorrelationLabel.FontWeight = 'bold';
            app.PhasecorrelationLabel.Position = [93 81 166 24];
            app.PhasecorrelationLabel.Text = 'Phase correlation:';

            % Create Switch
            app.Switch = uiswitch(app.IntensityBasedMethodsTab, 'slider');
            app.Switch.ValueChangedFcn = createCallbackFcn(app, @SwitchValueChanged2, true);
            app.Switch.Position = [303 81 50 22];

            % Create MIPanel
            app.MIPanel = uipanel(app.IntensityBasedMethodsTab);
            app.MIPanel.Title = 'MI';
            app.MIPanel.BackgroundColor = [0.8353 0.8157 0.8392];
            app.MIPanel.Position = [109 128 233 442];

            % Create MIButton
            app.MIButton = uibutton(app.MIPanel, 'push');
            app.MIButton.ButtonPushedFcn = createCallbackFcn(app, @MIButtonPushed, true);
            app.MIButton.FontName = 'Kristen ITC';
            app.MIButton.Position = [33 25 177 30];
            app.MIButton.Text = 'MI';

            % Create InitialRadiusSlider
            app.InitialRadiusSlider = uislider(app.MIPanel);
            app.InitialRadiusSlider.Limits = [0.0001 0.01];
            app.InitialRadiusSlider.ValueChangedFcn = createCallbackFcn(app, @InitialRadiusSliderValueChanged, true);
            app.InitialRadiusSlider.Position = [34 375 162 3];
            app.InitialRadiusSlider.Value = 0.0001;

            % Create GrowthFactorSlider
            app.GrowthFactorSlider = uislider(app.MIPanel);
            app.GrowthFactorSlider.Limits = [1.0001 1.05];
            app.GrowthFactorSlider.ValueChangedFcn = createCallbackFcn(app, @GrowthFactorSliderValueChanged, true);
            app.GrowthFactorSlider.Position = [31 218 162 3];
            app.GrowthFactorSlider.Value = 1.0001;

            % Create MaxIterationsSliderLabel
            app.MaxIterationsSliderLabel = uilabel(app.MIPanel);
            app.MaxIterationsSliderLabel.HorizontalAlignment = 'right';
            app.MaxIterationsSliderLabel.FontName = 'Kristen ITC';
            app.MaxIterationsSliderLabel.Position = [58 140 93 22];
            app.MaxIterationsSliderLabel.Text = 'Max Iterations';

            % Create MaxIterationsSlider
            app.MaxIterationsSlider = uislider(app.MIPanel);
            app.MaxIterationsSlider.Limits = [300 5000];
            app.MaxIterationsSlider.ValueChangedFcn = createCallbackFcn(app, @MaxIterationsSliderValueChanged, true);
            app.MaxIterationsSlider.FontName = 'Kristen ITC';
            app.MaxIterationsSlider.Position = [30 122 162 3];
            app.MaxIterationsSlider.Value = 300;

            % Create EpsilonSliderLabel
            app.EpsilonSliderLabel = uilabel(app.MIPanel);
            app.EpsilonSliderLabel.HorizontalAlignment = 'right';
            app.EpsilonSliderLabel.FontName = 'Kristen ITC';
            app.EpsilonSliderLabel.Position = [86 302 47 22];
            app.EpsilonSliderLabel.Text = 'Epsilon';

            % Create EpsilonSlider
            app.EpsilonSlider = uislider(app.MIPanel);
            app.EpsilonSlider.Limits = [1e-06 0.001];
            app.EpsilonSlider.ValueChangedFcn = createCallbackFcn(app, @EpsilonSliderValueChanged, true);
            app.EpsilonSlider.FontName = 'Kristen ITC';
            app.EpsilonSlider.Position = [35 293 162 3];
            app.EpsilonSlider.Value = 1e-06;

            % Create EpsilonSliderLabel_2
            app.EpsilonSliderLabel_2 = uilabel(app.MIPanel);
            app.EpsilonSliderLabel_2.HorizontalAlignment = 'right';
            app.EpsilonSliderLabel_2.Position = [76 387 70 22];
            app.EpsilonSliderLabel_2.Text = 'InitialRadius';

            % Create EpsilonSliderLabel_3
            app.EpsilonSliderLabel_3 = uilabel(app.MIPanel);
            app.EpsilonSliderLabel_3.HorizontalAlignment = 'right';
            app.EpsilonSliderLabel_3.Position = [75 230 78 22];
            app.EpsilonSliderLabel_3.Text = 'GrowthFactor';

            % Create MMIButton
            app.MMIButton = uibutton(app.IntensityBasedMethodsTab, 'push');
            app.MMIButton.ButtonPushedFcn = createCallbackFcn(app, @MMIButtonPushed, true);
            app.MMIButton.BackgroundColor = [0.8392 0.8196 0.8392];
            app.MMIButton.FontName = 'Kristen ITC';
            app.MMIButton.Position = [251 20 150 36];
            app.MMIButton.Text = 'MMI';

            % Create ImageRegistratorLabel
            app.ImageRegistratorLabel = uilabel(app.UIFigure);
            app.ImageRegistratorLabel.HorizontalAlignment = 'center';
            app.ImageRegistratorLabel.FontName = 'Kristen ITC';
            app.ImageRegistratorLabel.FontSize = 36;
            app.ImageRegistratorLabel.FontWeight = 'bold';
            app.ImageRegistratorLabel.Position = [281 651 581 50];
            app.ImageRegistratorLabel.Text = 'Image Registrator ';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = VisionFuse_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end