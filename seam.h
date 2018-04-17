
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

int build_matrix(void);
void find_energy_map(pixel_t** imagePixelArray, double** energy_array, num_rows, num_cols);
