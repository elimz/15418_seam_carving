
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <omp.h>
#include "cycletimer.h"

// debug 
#include <assert.h>


#ifndef OMP
#define OMP 0
#endif

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

pixel_t** build_matrix(int *rows, int *cols, int *max_px, char* file);
void find_energy_map(pixel_t** image_pixel_array, double** energy_array, int num_rows, int num_cols);
void output_image(pixel_t** image_pixel_array, char* output_file, int num_cols, int num_rows, int max_px_val);
void compute_E(pixel_t** imagePixelArray, double* E, int num_rows, int num_cols, int NTHREAD);
// void compute_M(double** E, double** M, int num_rows, int num_cols);
double find_seam(double* E, int* seam_path, int num_rows, int num_cols);
void color_seam(pixel_t*** imagePixelArray, int* seam_paths, int num_rows, int num_cols, int max_px_val, char* seam_file);
void remove_seam(pixel_t*** image_pixel_array, int* seam_paths, int* rows, int* cols);

// helper functions: 
int start_pos_partition (int N, int P, int i);
// debug / visualization functions
void intermediary_img(double ** matrix, char* output_file,  \
    int num_rows, int num_cols, int max_px_val, int min_px_val);

// // outer wrapper
// int main_support(int nthread);

