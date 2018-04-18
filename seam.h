
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

pixel_t** build_matrix(int *rows, int *cols, int *max_px);
void find_energy_map(pixel_t** image_pixel_array, double** energy_array, int num_rows, int num_cols);
void output_image(pixel_t** image_pixel_array, char* output_file, int num_cols, int num_rows, int max_px_val);
void compute_E(pixel_t** imagePixelArray, double** E, int num_rows, int num_cols);
void compute_M(double** E, double** M, int num_rows, int num_cols);
void find_seam(double** E, int* seam_path, int num_rows, int num_cols);
void color_seam(pixel_t*** imagePixelArray, double** M, int* seam_path, int num_rows, int num_cols, int max_px_val);
void remove_seam(pixel_t*** image_pixel_array, int* seam_path, int* rows, int* cols);
