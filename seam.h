
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

int build_matrix(void);
void compute_E(pixel_t** imagePixelArray, double** E, int num_rows, int num_cols);
void find_seam(double** energy_array, int num_rows, int num_cols)
