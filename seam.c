
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
    int num_rows;
    int num_cols;
    pixel_t*** image_pixel_array;
    // char *file_path = "./images/tower.ppm";
    // pass in pointers for image array pointer, width, height
    build_matrix();

    // alloc mem for energy array
    double** energy_array = malloc(sizeof(double*) * num_rows);
    int row;
    // for (row = 0; row < NUMBER OF ROWS; row++) {
    for (row = 0; row < num_rows; row++) {
        energy_array[row] = malloc(sizeof(double) * num_cols);
    }

    // compute energy map and store in engery_array
    find_energy_map(image_pixel_array, energy_array, num_rows, num_cols);



    return 0;
}


// this function parses the ppm file and put all the pixels RGB values in a 
//  2D matrix. returns the pointer to this 2D array; 
// int build_matrix(char *file_path){
int build_matrix(void){

    FILE *ppm_file = fopen("./images/tower.ppm", "r");
    if (ppm_file == NULL) {
        printf("Error: failed to open file!\n");
        fclose(ppm_file);
        return 1;
    }

    printf("entering function build_matrix\n");
    int bufsize = 100;
    char buf[bufsize];
    int width, height; 
    int max_px_val;

    // first 2 lines include image mode and dimension; 
    fgets(buf, bufsize, ppm_file);
    // reads in image dimension
    if (fgets(buf, bufsize, ppm_file) != NULL){
        // take in image height and width
        sscanf(buf, "%d %d\n", &width, &height);
        printf("image width = %d, height = %d\n", width, height);
    }
    // reads in max pixel value;
    if (fgets(buf, bufsize, ppm_file) != NULL){
        // take in image height and width
        sscanf(buf, "%d\n", &max_px_val);
        printf("max_px_val = %d\n", max_px_val);
    }


    // now parse all the rgb info into a matrix;
    // image is represented by a list of structs, each struct has R, G, B field 
    int image_size = width * height; 

    // allocating a list of structs; 
    pixel_t *matrix = (pixel_t *)malloc(image_size * sizeof(pixel_t));
    // parse ppm image to fill in the matrix; 
    char pixel_buf[1500];   // TODO: only need to fit in24 numbers separated by space; 
    char *buf_ptr = pixel_buf;

    // PARSE pixel info from pixel_buf and fill the matrix
    int i; 

    for (i = 0; i < image_size; i ++){
        // go to a new struct, read 3 numbers; 
        pixel_t curr_pixel = matrix[i];
        int offset; 

        // every 8 groups of RGB values, need to read in another line; 
        if (i % 8 == 0){
            printf("   entering -- \n");
            fgets(pixel_buf, 1500 , ppm_file);
            printf("pixel_buf = %s\n", pixel_buf);
            char *buf_ptr = pixel_buf;
        }

        sscanf(buf_ptr, "%d %d %d%n",  \
                          &curr_pixel.R, &curr_pixel.G, &curr_pixel.B, &offset); 
        int R = curr_pixel.R;
        int G = curr_pixel.G;
        int B = curr_pixel.B;
        printf("offset = %d, (%d, %d, %d)\n", offset, R, G, B); // test
        buf_ptr = buf_ptr + offset; // move pointer to next number;

        if (i == 10)
            break;
    }


    fclose(ppm_file);
    // return matrix;
    return 0;
}   

// // this function parses the ppm file and put all the pixels RGB values in a 
// //  2D matrix
// int build_matrix(FIL){
//     return 0;
// }   

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

