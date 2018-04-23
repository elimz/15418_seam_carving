# 15418_seam_carving
15418 Team Project: Seam carving in Parallel. Implementing content-aware images resizing in real-time

Procedures:
- Converting image formats form png/jpeg to ppm in command line "convert trial.png -compress none dest.ppm". The ppm file contains image info such as width and height, as well as a char array of RGB values at each pixel.
- seam.c function takes in the ppm file and works on image processing.

Images contains photos from National Geography, and are used to test the algorithm;

