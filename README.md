# VisionFuse MATLAB Application

#Description
VisionFuse is a MATLAB application designed for image processing and registration. It offers a range of functionalities including noise reduction, intensity normalization, image resampling, artifact removal, feature-based methods (like SURF and SIFT), and intensity-based methods (like NCC and MMI).

#Features

▹Noise Reduction: Gaussian, Median, Non-local means filtering.
▹Intensity Normalization: Histogram Matching, Contrast Stretching.
▹Image Resampling: Resizing and Rotation.
▹Artifact Removal: Techniques to remove or reduce artifacts in images.
▹Feature-Based Methods: SURF, SIFT, and ORB feature detection.
▹Intensity-Based Methods: Normalized Cross-Correlation (NCC), Mutual Information (MI), Multimodal Image registration (MMI).

#Prerequisites
To run the VisionFuse application, the following MATLAB toolboxes are required:

▹Image Processing Toolbox
▹Computer Vision Toolbox
Note: The code includes a dependency on VLFeat (for SIFT features), which needs to be installed separately.

#Installation
1-Clone or download the VisionFuse repository to your local machine.
2-Ensure that the required MATLAB toolboxes are installed.
3-Run the VisionFuse.mlapp file in MATLAB.

#Usage
▹Import images using the provided buttons and view the images on axis1/axis2.
▹Select the desired processing methods and adjust the parameters as needed.
▹View the processed images and results in the designated axes (output on axis3 or axis4).
