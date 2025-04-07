% Define parameters
timePoints = 260; % Number of time points
zLayers = 47;     % Number of z layers
inputFolder = 'H:\data\PLOS data\Experiment\Experiment\Figure 2A\Data S1\DataS1_RawImage_Membrane\'; % Path to the folder with the images
outputFolder = fullfile(inputFolder, 'image'); % Save output in a new "image" folder within the input folder

% Create the output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Loop through each time point
 for t = 1:timePoints
    % Initialize an array to hold the images for each Z layer
    imgStack = [];

    % Loop through each z layer
    for z = 1:zLayers
        % Construct the filename
        filename = sprintf('%sMembrane_t%03d-p%02d.tif', inputFolder, t, z);
        % Load the image
        img = imread(filename);

        % Store in the stack
        imgStack(:,:,z) = img;
    end

    % Perform maximum intensity projection
    imgProjection = max(imgStack, [], 3);

    % Create an RGB image where white parts become red and black remains black
    imgRGB = zeros([size(imgProjection), 3], 'uint8'); % Initialize RGB image
    imgRGB(:,:,1) = imgProjection; % Red channel set to projection (red for white parts)

    % Rotate the image 180 degrees
    imgRGB = imrotate(imgRGB, 180);
    
    % Save the projected image as PNG in the new "image" folder
    outputFilename = sprintf('%s\\Projection_t%03d.png', outputFolder, t);
    imwrite(imgRGB, outputFilename);
    end