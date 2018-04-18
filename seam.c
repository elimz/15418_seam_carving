
// seam carving application in parallel
// steps: 
// 1. command line, convert jpg/png format into ppm; 
//  "convert source.png -compress none dest.ppm" // require imageMagick

#include <math.h>

#ifndef SEAM_H
#include "seam.h"
#define SEAM_H
#endif

#define MAX_ENERGY 550.0;

int main(){
    int num_rows, num_cols, max_px_val;
    pixel_t** image_pixel_array;
    // TODO: pass in file paths into all the functions 
    // char *file_path = "./images/tower.ppm";
    // pass in pointers for image array pointer, width, height
    build_matrix(&num_rows, &num_cols, &max_px_val, &image_pixel_array);

    // alloc mem for energy array
    double** energy_array = malloc(sizeof(double*) * num_rows);
    int row;
    for (row = 0; row < num_rows; row++) {
        energy_array[row] = malloc(sizeof(double) * num_cols);
    }
    // compute energy map and store in engery_array
    find_energy_map(image_pixel_array, energy_array, num_rows, num_cols);

    // output image;        TODO: change file name to match file input name
    char* out_file_name = "tower_out.ppm";
    output_image(image_pixel_array, out_file_name, num_rows, num_cols, max_px_val);
    printf("Finished - Image Processing Finished! \n");
    return 0;
}


// this function parses the ppm file and put all the pixels RGB values in a 
//  2D matrix. returns the pointer to this 2D array; 
// int build_matrix(char *file_path){
int build_matrix(int *rows, int *cols, int *max_px, pixel_t ***mx){

    FILE *ppm_file = fopen("./images/tower.ppm", "r");
    if (ppm_file == NULL) {
        printf("Error: failed to open file!\n");
        fclose(ppm_file);
        exit(1);
    }

    int bufsize = 1500;
    char buf[bufsize]; 
    int max_px_val;
    int num_rows, num_cols; 

    // first 2 lines include image mode and dimension; 
    fgets(buf, bufsize, ppm_file);
    // reads in image dimension
    if (fgets(buf, bufsize, ppm_file) != NULL){
        // take in image height and width
        sscanf(buf, "%d %d\n", &num_cols, &num_rows);
    }

    // pass this global value to other functions;
    *rows = num_rows; 
    *cols = num_cols;

    // reads in max pixel value;
    if (fgets(buf, bufsize, ppm_file) != NULL){
        // take in image height and num_cols
        sscanf(buf, "%d\n", &max_px_val);
        // printf("max_px_val = %d\n", max_px_val);
    }
    *max_px = max_px_val;


    printf("Reading in image with size (%d, %d), max pixel value = %d\n", \
        num_cols, num_rows, max_px_val);


    // now parse all the rgb info into a matrix;
    // image is represented by a list of structs, each struct has R, G, B field 
    int image_size = num_cols * num_rows; 
    
    // --------------------------------TODO: allocating a 1D array; --------------------------------
    // pixel_t *matrix = (pixel_t *)malloc(sizeof(pixel_t) * num_cols * num_rows); 
    // --------------------------------TODO: allocating a 1D array; --------------------------------

    // allocate a 2D array ;
    pixel_t **matrix = (pixel_t **) malloc(sizeof(pixel_t *) * num_rows);
    if (matrix == NULL){
        printf("ERROR: Malloc failed- matrix\n");
        exit(1);
    }
    int i;
    for (i = 0; i < num_rows; i ++){
        matrix[i] = (pixel_t *)malloc(sizeof(pixel_t) * num_cols);
        if (matrix[i] == NULL){
            printf("ERROR: Malloc failed - matrix[i]\n");
            exit(1);
        }
    }

    int curr_row, curr_col; 
    int pixel_count = 0;

    while (fgets(buf, bufsize ,ppm_file) != NULL){
        if(pixel_count == image_size){
            printf("ERROR: getting more pixels than image_size. pixel_count = %d, image_size =%d\n", \
                pixel_count, image_size);
        }
        char *buf_ptr = buf;
        int R, G, B;
        int offset; 

        // each line contains 8 groups of RGB pixel values; reading them 3 at a 
        //  time, and store them to RGB fields;
        while (sscanf(buf_ptr, "%d %d %d%n",  \
                          &R, &G, &B, &offset) == 3){

            // calculate curr_row, curr_col;
            curr_row = pixel_count / num_cols; 
            curr_col = pixel_count % num_cols;

            (matrix[curr_row][curr_col]).R = R;
            (matrix[curr_row][curr_col]).G = G;
            (matrix[curr_row][curr_col]).B = B;

            buf_ptr = buf_ptr + offset; 
            pixel_count ++;
        }
    }
    printf("Finished - convert input image into matrix\n");

    // allowing other funcitons to access matrix: 
    *mx = matrix;
    fclose(ppm_file);
    return 0;
}   

double pixel_difference(pixel_t pU, pixel_t pD, pixel_t pL, pixel_t pR) {
    // find partial derivatives for x
    int dxR = abs(pR.R - pL.R);
    int dxG = abs(pR.G - pL.G);
    int dxB = abs(pR.B - pL.B);
    int dx = (dxR + dxG + dxB) / 2;

    // find partial derivatives for y
    int dyR = abs(pD.R - pU.R);
    int dyG = abs(pD.G - pU.G);
    int dyB = abs(pD.B - pU.B);
    int dy = (dyR + dyG + dyB) / 2;

    // return magnitude
    return sqrt((dx * dx) + (dy * dy));

}

void output_image(pixel_t** matrix, char* output_file,  int num_rows, int num_cols, int max_px_val){

    FILE *fp = fopen(output_file, "wb");
    if (fp == NULL) {
        printf("Error: failed to open file! file - %s\n", output_file);
        fclose(fp);
        exit(1);
    }

    // write file header first; 
    fprintf(fp, "P3\n%d %d\n%d\n", num_cols, num_rows, max_px_val);

    // write the rest of data;
    int row, col;
    pixel_t curr_pixel;
    for (row = 0; row < num_rows; row ++){
        // write each row as a row in the ppm file as well;
        for (col = 0; col < num_cols; col ++) {
            // unpack matrix values and write to file; 
            curr_pixel = matrix[row][col];

            if (col > 0){
                fprintf(fp, " ");   // appending space to later items
            }

            fprintf(fp, "%d %d %d", curr_pixel.R, curr_pixel.G, curr_pixel.B);
            // last item appends a new line at the end;
            if (col == num_cols - 1){
                fprintf(fp, "\n");
            }
        }
    }
    printf("Finished - writing image to output.\n");
}

// takes in an image, and return an energy map calculated by gradient magnitude
void find_energy_map(pixel_t** image_pixel_array, double** energy_array, int num_rows, int num_cols) {
    int i;
    for (i = 0; i < num_rows; i++) {
        int j;
        for (j = 0; j < num_cols; j++) {
            // don't want to remove the edges
            if (i == 0 || j == 0 ||
                i == num_rows - 1 || j == num_cols - 1) {
                energy_array[i][j] = MAX_ENERGY;
            } else {
                energy_array[i][j] = pixel_difference(image_pixel_array[i-1][j],
                                                     image_pixel_array[i+1][j],
                                                     image_pixel_array[i][j-1],
                                                     image_pixel_array[i][j+1]);
            }
        }
    }
}

