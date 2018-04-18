
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

int build_matrix(pixel_t** image_pixel_array, int *rows, int *cols);
void compute_E(pixel_t** imagePixelArray, double** E, int num_rows, int num_cols);
void compute_M(double** E, double** M, int num_rows, int num_cols);
void find_seam(double** E, int* seam_path, int num_rows, int num_cols);
void color_seam(pixel_t** imagePixelArray, int* seam_path, int num_rows, int num_cols);

