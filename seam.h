
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

int build_matrix(int *rows, int *cols, int *max_px, pixel_t ***mx);
void find_energy_map(pixel_t** image_pixel_array, double** energy_array, int num_rows, int num_cols);
void output_image(pixel_t** image_pixel_array, char* output_file, int num_cols, int num_rows, int max_px_val);