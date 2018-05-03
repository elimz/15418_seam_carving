
// seam carving application in parallel
// steps: 
// 1. command line, convert jpg/png format into ppm; 
//  "convert source.png -compress none dest.ppm" // require imageMagick


// TODO next: 
// - measure timing; 
// - use valgrind to measaure memory allocated, before and after changing to a 1D array; 
// - handle image input and output naming inside main function
// Stretch goal:
// - different data structure to store RGB info. 
//      instead of storing all of (RGB) in one huge matrix, separate into 3 matrices; 
//      this might benefit parallelizing later 

#include "seam.h"

#define MAX_ENERGY 9999999.0
#define NUM_SEAMS_TO_REMOVE 70

#ifndef OMP
#define OMP 1
#endif

// debug - batch write for compute_E; 
#ifndef BATCH 
#define BATCH 0
#endif


char* input_file = "../images/tower.ppm";
char* output_file = "output.ppm";
char* seam_file = "output_seam.ppm";

// TODO: 
int nthread = 8;

int main(){
    
    #if TIMING
    double t0 = currentSeconds();
    #endif

    int num_rows, num_cols, original_cols, max_px_val;
    // TODO: pass in file paths into all the functions 
    // char *file_path = "./images/tower.ppm";
    // pass in pointers for image array pointer, width, height
    pixel_t** image_pixel_array = build_matrix(&num_rows, &num_cols, &max_px_val, input_file);
    original_cols = num_cols;
    // alloc for E array
    int row; 
    double** E = malloc(sizeof(double*) * num_rows);
    for (row = 0; row < num_rows; row++) {
        E[row] = malloc(sizeof(double) * num_cols);
    }
    // double* E = malloc(sizeof(double) * (num_rows * num_cols));

    // alloc mem for seam_path array
    // int* seam_path = malloc(sizeof(int) * num_rows);
    int seam_num;
    int** seam_paths = malloc(sizeof(int*) * NUM_SEAMS_TO_REMOVE);
    for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
        seam_paths[seam_num] = malloc(sizeof(int) * num_rows);
    }

    double start = currentSeconds();

    // compute energy map and store in engery_array
    compute_E(image_pixel_array, E, num_rows, num_cols);
    printf("Finished computing E\n");

    // set up arrays to sort
    double* first_row = E[0];
    double first_row_pixel[num_cols];
    int first_row_indices[num_cols];
    memcpy(first_row_pixel, first_row, num_cols * sizeof(double));
    int idx;
    for (idx = 0; idx < num_cols; idx++) {
        first_row_indices[idx] = idx;
    }

    // sort to find the indices of least energy in the first row
    quickSort_double(first_row_pixel, first_row_indices, 0, num_cols - 1);
    int ordered_indicies[NUM_SEAMS_TO_REMOVE];
    memcpy(ordered_indicies, first_row_indices, NUM_SEAMS_TO_REMOVE * sizeof(int));
    quickSort_int(ordered_indicies, 0, NUM_SEAMS_TO_REMOVE - 1);

    // remove NUM_SEAMS_TO_REMOVE number of lowest cost seams
    for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
        printf("%d, ", ordered_indicies[seam_num]);
        seam_paths[seam_num][0] = ordered_indicies[seam_num];
        find_seam(E, seam_paths[seam_num], num_rows, num_cols);
        // printf("Finished finding seam\n");
    }

    // color the seam and output the image
    color_seam(&image_pixel_array, seam_paths, num_rows, num_cols, max_px_val, seam_file);
    printf("Finished coloring seam\n");
    
    // remove the seam from the image, also sets new values for num_rows and num_cols
    remove_seam(&image_pixel_array, seam_paths, &num_rows, &num_cols);
    printf("Finished removing seam\n");


    double delta = currentSeconds() - start;
    printf("%d seams of a %dx%d image removed in %.3f seconds\n", NUM_SEAMS_TO_REMOVE, original_cols, num_rows, delta);

    output_image(image_pixel_array, output_file, num_rows, num_cols, max_px_val);
    
    // // free data structures 
    // for (row = 0; row < num_rows; row++) {
    //     free(E[row]);
    // }

    for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
        free(seam_paths[seam_num]);
    }

    for (row = 0; row < num_rows; row++) {
        free(image_pixel_array[row]);
    }

    free(E);
    free(seam_paths);
    free(image_pixel_array);
    printf("Finished - Image Processing Finished! \n");
    return 0;
}  

double pixel_difference(pixel_t pU, pixel_t pD, pixel_t pL, pixel_t pR) {
    // find partial derivative for x
    int dxR = abs(pR.R - pL.R);
    int dxG = abs(pR.G - pL.G);
    int dxB = abs(pR.B - pL.B);
    double deltx2 = (dxR + dxG + dxB) / 2.0;
    // double deltx2 = (dxR*dxR + dxG*dxG + dxB*dxB);

    // find partial derivative for y
    int dyR = abs(pD.R - pU.R);
    int dyG = abs(pD.G - pU.G);
    int dyB = abs(pD.B - pU.B);
    double delty2 = (dyR + dyG + dyB) / 2.0;
    // double delty2 = (dyR*dyR + dyG*dyG + dyB*dyB);

    // return magnitude for gradient magnitude
    // return sqrt(deltx2 + delty2);
    return deltx2 + delty2;
}

// takes in an image, and return an energy map calculated by gradient magnitude

void compute_E(pixel_t** image_pixel_array, double** E, int num_rows, int num_cols) {
    // debug: 
    // int max_pix_val = 255; 
    // int min_pix_val = 0;

    // int i, j;
    int i ;
    for (i = 0; i < num_rows; i++) {
        int j;
        for (j = 0; j < num_cols; j++) {
            // don't want to remove the edges
            if (i == num_rows - 1 || j == num_cols - 1) {
                E[i][j] = MAX_ENERGY;
            } else {
                E[i][j] = pixel_difference(image_pixel_array[i][j],
                                           image_pixel_array[i + 1][j],
                                           image_pixel_array[i][j],
                                           image_pixel_array[i][j + 1]);
            }
            // debug 
            // if (E[i][j] > max_pix_val){
            //     max_pix_val = E[i][j];    
            // } if (E[i][j] < min_pix_val){
            //     min_pix_val = E[i][j];
            // } 
        }
    }
    // // debug: print out gradient file;
    // char* output_file1 = "gradient.ppm";
    // intermediary_img(E, output_file1,  num_rows, num_cols, max_pix_val, min_pix_val);
}

// finds the seam from E
void find_seam(double** E, int* seam_path, int num_rows, int num_cols) {

    // find min seams with least energy cost pixel in first row
    // traverse downwards and greedily pick the least energy cost neighbor

    int i;
    for (i = 1; i < num_rows; i++) {
        int prev_col = seam_path[i - 1];
        double middle = E[i][prev_col];
        double left;
        double right;
        if (prev_col - 1 < 0) {
            left = INT_MAX;
        } else {
            left = E[i][prev_col - 1];
        }
        if (prev_col + 1 >= num_cols) {
            right = INT_MAX;
        } else {
            right = E[i][prev_col + 1];
        }

        double current_min_cost = fmin(middle, fmin(left, right));

        if (current_min_cost == INT_MAX) {
            printf("left: %f, middle: %f, right: %f, prev_col: %d\n", left, middle, right, prev_col);
        }

        if (current_min_cost == middle) {
            seam_path[i] = prev_col;
            E[i][prev_col] = INT_MAX;
        } else if (current_min_cost == left) {
            seam_path[i] = prev_col - 1;
            E[i][prev_col - 1] = INT_MAX;
        } else {
            seam_path[i] = prev_col + 1;
            E[i][prev_col + 1] = INT_MAX;
        }
    }

    // for (i = 0; i < num_rows; i++) {
    //     if (seam_path[i] < 0) {
    //         printf("seam index: %d with value %d\n", i, seam_path[i]);
    //     }
    // }
    printf("\n\n");
}

// colors the seam pixels red and outputs the image
void color_seam(pixel_t*** image_pixel_array, int** seam_paths, int num_rows, int num_cols, int max_px_val, char* seam_file) {
    int seam_num;
    for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
        int* current_seam_path = seam_paths[seam_num];

        int i;
        for (i = 0; i < num_rows; i++) {
            (*image_pixel_array)[i][current_seam_path[i]].R = 255;
            (*image_pixel_array)[i][current_seam_path[i]].G = 0;
            (*image_pixel_array)[i][current_seam_path[i]].B = 0;
        }
    }

    // call function to output image
    output_image(*image_pixel_array, seam_file, num_rows, num_cols, max_px_val);
}

// remove the seam from the image pixel array and get new dimensions
void remove_seam(pixel_t*** image_pixel_array_pt, int** seam_paths, int* rows, int* cols) {
    int num_rows = *rows;
    int num_cols = *cols;

    // need to shift every pixel after the seam over to the "left"
    int seams_in_row[NUM_SEAMS_TO_REMOVE];
    int row;
    for (row = 0; row < num_rows; row++) { 
        int seam_num;
        for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
            seams_in_row[seam_num] = seam_paths[seam_num][row];
        }
        quickSort_int(seams_in_row, 0, NUM_SEAMS_TO_REMOVE - 1);
        int prev = 0;
        for (seam_num = 1; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
            // printf("%d, ", seams_in_row[seam_num]);
            prev = seams_in_row[seam_num - 1];
            if (prev == seams_in_row[seam_num]) {
                printf("AAAAAA\n");
            }
        }
        // printf("\n\n");

        int num_to_remove = 1;
        int next_seam_idx = 1;

        int j;
        for (j = seams_in_row[0] + 1; j < num_cols; j++) {
            (*image_pixel_array_pt)[row][j - num_to_remove] = (*image_pixel_array_pt)[row][j];
             if (j == seams_in_row[next_seam_idx]) {
                num_to_remove++;
                next_seam_idx++;
            }
        }
    }


    *cols = num_cols - NUM_SEAMS_TO_REMOVE;


            // if (seam_col < 0 || seam_col >= num_cols ) {
        //     printf("LESS THAN 0 \n");
        // }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////


// this function parses the ppm file and put all the pixels RGB values in a 
//  2D matrix. returns the pointer to this 2D array; 
pixel_t** build_matrix(int *rows, int *cols, int *max_px, char* file){
    FILE *ppm_file = fopen(file, "r");
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
    // *mx = matrix;
    fclose(ppm_file);
    printf("Finished - convert input image into matrix\n");
    return matrix;
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
    // printf("Finished - writing image to output.\n");
}

void intermediary_img(double** matrix, char* output_file,
    int num_rows, int num_cols, int max_px_val, int min_px_val){

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
    double curr_value;
    for (row = 0; row < num_rows; row ++){
        // write each row as a row in the ppm file as well;
        for (col = 0; col < num_cols; col ++) {
            // // unpack matrix values and write to file; 
            curr_value = matrix[row][col];

            // scale up so it's bright enough to see; 
            int rand_offset = 1000;
            // convert to a scale of 0 to 255
            curr_value = (curr_value + rand_offset)* 255  / (max_px_val - min_px_val) ;


            // cast to int in order to display
            int int_curr_value = (int)curr_value;
            if (int_curr_value <  0 ){
                int_curr_value = 0;
            }


            if (col > 0){
                fprintf(fp, " ");   // appending space to later items
            }

            fprintf(fp, "%d %d %d", int_curr_value, int_curr_value, int_curr_value);
            // fprintf(fp, "%d", int_curr_value);
            // last item appends a new line at the end;
            if (col == num_cols - 1){
                fprintf(fp, "\n");
            }
            // printf("%lf\n", curr_value);
        }
    }
    printf("Finished - writing intermediary image to output.\n");
}

///////////////////////////////////////////////////////////////////////////////////////////

