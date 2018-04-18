
// seam carving application in parallel
// steps: 
// 1. command line, convert jpg/png format into ppm; 
//  "convert source.png -compress none dest.ppm" // require imageMagick

#include <math.h>

#ifndef SEAM_H
#include "seam.h"
#define SEAM_H
#endif

#define MAX_ENERGY 550.0
#define NUM_SEAMS_TO_REMOVE 1

int main(){

    int num_rows, num_cols, max_px_val;
    pixel_t** image_pixel_array;
    // TODO: pass in file paths into all the functions 
    // char *file_path = "./images/tower.ppm";
    // pass in pointers for image array pointer, width, height
    build_matrix(&num_rows, &num_cols, &max_px_val, &image_pixel_array);

    // alloc for E array
    int row; 
    double** E = malloc(sizeof(double*) * num_rows);
    for (row = 0; row < num_rows; row++) {
        E[row] = malloc(sizeof(double) * num_cols);
    }

    // alloc mem for M array
    double** M = malloc(sizeof(double*) * num_rows);
    for (row = 0; row < num_rows; row++) {
        M[row] = malloc(sizeof(double) * num_cols);
    }

    // alloc mem for seam_path array
    int* seam_path = malloc(sizeof(int) * num_rows);

    // remove NUM_SEAMS_TO_REMOVE number of lowest cost seams
    int seam_num;
    for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
        // compute energy map and store in engery_array
        compute_E(image_pixel_array, E, num_rows, num_cols);

        //b TODO This function causes seg fault
        // // compute M from E
        // compute_M(E, M, num_rows, num_cols);
        // printf("here? 53  \n");

        // // find seam to remove
        // find_seam(M, seam_path, num_rows, num_cols);

        // // color the seam and output the image
        // color_seam(image_pixel_array, seam_path, num_rows, num_cols);

        // // remove the seam from the image, also sets new values for num_rows and num_cols
        // remove_seam(image_pixel_array, seam_path, &num_rows, &num_cols);
    }

    // output image;        TODO: change file name to match file input name
    char* out_file_name = "tower_out.ppm";
    output_image(image_pixel_array, out_file_name, num_rows, num_cols, max_px_val);
    
    // free data structures 
    for (row = 0; row < num_rows; row++) {
        free(E[row]);
    }

    for (row = 0; row < num_rows; row++) {
        free(M[row]);
    }

    for (row = 0; row < num_rows; row++) {
        free(image_pixel_array[row]);
    }

    free(E);
    free(M);
    free(seam_path);
    free(image_pixel_array);
    printf("Finished - Image Processing Finished! \n");
    return 0;
}


// this function parses the ppm file and put all the pixels RGB values in a 
//  2D matrix. returns the pointer to this 2D array; 
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

    // allocate a 2D array for the pixel matrix 
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

    // allowing other funcitons to access matrix: 
    *mx = matrix;
    fclose(ppm_file);
    printf("Finished - convert input image into matrix\n");
    return 0;
}   

double pixel_difference(pixel_t pU, pixel_t pD, pixel_t pL, pixel_t pR) {
    // find pixel difference for x
    int dxR = (pR.R - pL.R);
    int dxG = (pR.G - pL.G);
    int dxB = (pR.B - pL.B);
    int deltx2 = (dxR*dxR + dxG*dxG + dxB*dxB);

    // find pixel difference for y
    int dyR = (pD.R - pU.R);
    int dyG = (pD.G - pU.G);
    int dyB = (pD.B - pU.B);
    int delty2 = (dyR*dyR + dyG*dyG + dyB*dyB);

    // return magnitude for dual-gradient
    return sqrt(deltx2 + delty2);

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
void compute_E(pixel_t** image_pixel_array, double** E, int num_rows, int num_cols) {
    int i;
    for (i = 0; i < num_rows; i++) {
        int j;
        for (j = 0; j < num_cols; j++) {
            // don't want to remove the edges
            if (i == 0 || j == 0 ||
                i == num_rows - 1 || j == num_cols - 1) {
                E[i][j] = MAX_ENERGY;
            } else {
                E[i][j] = pixel_difference(image_pixel_array[i-1][j],
                                           image_pixel_array[i+1][j],
                                           image_pixel_array[i][j-1],
                                           image_pixel_array[i][j+1]);
            }
        }
    }
}

// M(i, j) = E(i, j) + min(M(i - 1, j - 1), M(i - 1, j), M(i - 1, j + 1))
void compute_M(double** E, double** M, int num_rows, int num_cols) {
    int i;
    for (i = 0; i < num_rows; i++) {
        int j;
        for (j = 0; j < num_cols; j++) {
            // printf("- computeM - row, col = (%d, %d); max (%d, %d)\n", i, j, num_rows, num_cols);
            int middle = E[i][j];
            
            int left;
            int right;
            if (j - 1 < 0) {
                left = MAX_ENERGY;
            } else {
                left = E[i][j-1];
            }

            if (j + 1 >= num_rows) {
                right = MAX_ENERGY;
            } else {
                right = E[i][j+1];
            }

            // TODO: added lines to avoid seg fault
            if (i < num_cols - 1){
                M[i+1][j] = M[i+1][j] + fmin(left, fmin(middle, right));
            }
        }
    }
}

// finds the seam from M
void find_seam(double** M, int* seam_path, int num_rows, int num_cols) {
    // find the min seam cost col in the last row
    double* last_row = M[num_rows - 1];

    int j;
    int min_cost_col = 0;
    int min_cost = last_row[0];
    for (j = 1; j < num_cols; j++) {
        int current_cost = last_row[j];
        if (current_cost < min_cost) {
            min_cost = current_cost;
            min_cost_col = j;
        } 
    }

    // go up from the bottom and find the small cost path
    int i;
    seam_path[num_rows - 1] = min_cost_col;
    for (i = num_rows - 2; i >= 0; i++) {
        int prev_col = seam_path[i + 1];
        int middle = M[i][prev_col];

        int left;
        int right;
        if (j - 1 < 0) {
            left = MAX_ENERGY;
        } else {
            left = M[i][prev_col - 1];
        }
        if (j + 1 >= num_rows) {
            right = MAX_ENERGY;
        } else {
            right = M[i][prev_col + 1];
        }
        int current_min_cost = fmin(left, fmin(middle, right));


        if (current_min_cost == left) {
            seam_path[i] = j - 1;
        } else if (current_min_cost == middle) {
            seam_path[i] = j;
        } else {
            seam_path[i] = j + 1;
        }
    }
}

// colors the seam pixels red and outputs the image
void color_seam(pixel_t** image_pixel_array, int* seam_path, int num_rows, int num_cols) {
    int i;
    for (i = 0; i < num_rows; i++) {
        image_pixel_array[i][seam_path[i]].R = 255;
        image_pixel_array[i][seam_path[i]].G = 0;
        image_pixel_array[i][seam_path[i]].B = 0;
    }

    // call function to output image
}

// remove the seam from the image pixel array and get new dimensions
void remove_seam(pixel_t** image_pixel_array, int* seam_path, int* rows, int* cols) {
    int num_rows = *rows;
    int num_cols = *cols;

    // need to shift every pixel after the seam over to the "left"
    int i;
    for (i = 0; i < num_rows; i++) {
        int j;
        int seam_col = seam_path[i];
        // seam col guaranteed to be at least 0
        for (j = seam_col + 1; j < num_cols; j++) {
            image_pixel_array[i][j - 1] = image_pixel_array[i][j];
        }
    }

    *cols = num_cols - 1;
}
