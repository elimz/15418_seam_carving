
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct pixel {
    int R;
    int G; 
    int B; 
} pixel_t;

int build_matrix(void);
int find_energy_map(pixel_t** imagePixelArray);
