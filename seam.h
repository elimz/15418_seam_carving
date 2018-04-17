
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

int build_matrix(void);
void find_energy_map(pixel_t** image_pixel_array, double** energy_array, int num_rows, int num_cols);
