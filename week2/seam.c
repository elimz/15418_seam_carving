
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

#ifndef TIMING 
#define TIMING 1
#endif

#define MAX_ENERGY 9999999
#define NUM_SEAMS_TO_REMOVE 300

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
    // double** E = malloc(sizeof(double*) * num_rows);
    // for (row = 0; row < num_rows; row++) {
    //     E[row] = malloc(sizeof(double) * num_cols);
    // }
    double* E = malloc(sizeof(double) * (num_rows * num_cols));

    // alloc mem for M array
    double** M = malloc(sizeof(double*) * num_rows);
    for (row = 0; row < num_rows; row++) {
        M[row] = malloc(sizeof(double) * num_cols);
    }

    // alloc mem for seam_path array
    int* seam_path = malloc(sizeof(int) * num_rows);

    // remove NUM_SEAMS_TO_REMOVE number of lowest cost seams
    #if TIMING
    double c_time_malloc = currentSeconds(); 
    double t_malloc = c_time_malloc - t0;
    #endif

    int seam_num;

    // timer array for different sections 
    #if TIMING
    double timing[4] = {0.0, 0.0, 0.0, 0.0 };
    #define T_COMP_E        0
    #define T_COMP_M        1
    #define T_FIND_SEAM     2
    #define T_REMOVE_SEAM   3
    double start = currentSeconds();        // after malloc
    #endif

    
    #if OMP
        printf("++++++++ TEST, OMP is ON\n");
        omp_set_num_threads(nthread);
    #endif

    #if BATCH
        printf("++++++++ TEST, BATCH is ON\n");
    #endif

    double t_start_E, t_start_M, t_start_f_seam, t_start_rm_seam; 
    double t_loop_start;

    for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {
        t_loop_start = currentSeconds();
        // compute energy map and store in engery_array
        compute_E(image_pixel_array, E, num_rows, num_cols);
        // printf("Finished computing E\n");
        #if TIMING
        t_start_E = currentSeconds();
        double t_compute_E = t_start_E - t_loop_start;
        timing[T_COMP_E] += t_compute_E;
        #endif

        // compute M from E
        compute_M(E, M, num_rows, num_cols);
        // printf("Finished computing M\n");
        #if TIMING
        t_start_M = currentSeconds();
        double t_compute_M = t_start_M - t_start_E;
        timing[T_COMP_M] += t_compute_M;
        #endif
        
        // find seam to remove
        find_seam(M, seam_path, num_rows, num_cols);
        // printf("Finished finding seam\n");
        #if TIMING
        t_start_f_seam = currentSeconds();
        double t_find_seam = t_start_f_seam - t_start_M;
        timing[T_FIND_SEAM] += t_find_seam;
        #endif
        
        // // color the seam and output the image
        // color_seam(&image_pixel_array, M, seam_path, num_rows, num_cols, max_px_val, seam_file);
        // printf("Finished coloring seam\n");

        // remove the seam from the image, also sets new values for num_rows and num_cols
        remove_seam(&image_pixel_array, seam_path, &num_rows, &num_cols);
        // printf("Finished removing seam\n");
        #if TIMING
        t_start_rm_seam = currentSeconds();
        double t_remove_seam = t_start_rm_seam - t_start_f_seam;
        timing[T_REMOVE_SEAM] += t_remove_seam;
        #endif
    }

    #if TIMING
    double delta = currentSeconds() - start;
    printf("------  TIMING ------ %d seams of a %dx%d image removed in %.3f seconds\n", NUM_SEAMS_TO_REMOVE, original_cols, num_rows, delta);
    #endif

    #if OMP
        printf("    OMP IS ON\n");
    #endif

    #if TIMING
    printf("------ TIMING SPLITS ------\n");
    printf("    T_MALLOC = %f\n", t_malloc);
    printf("    T_COMP_E = %f\n    T_COMP_M = %f\n    T_FIND_SEAM = %f\n    T_REMOVE_SEAM = %f\n", \
        timing[T_COMP_E], timing[T_COMP_M], timing[T_FIND_SEAM], timing[T_REMOVE_SEAM]);
    #endif
    // output image;        TODO: change file name to match file input name
    // char* out_file_name = "tower_out.ppm";
    output_image(image_pixel_array, output_file, num_rows, num_cols, max_px_val);
    
    // // free data structures 
    // for (row = 0; row < num_rows; row++) {
    //     free(E[row]);
    // }

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

double pixel_difference(pixel_t pU, pixel_t pD, pixel_t pL, pixel_t pR) {
    // find partial derivative for x
    int dxR = abs(pR.R - pL.R);
    int dxG = abs(pR.G - pL.G);
    int dxB = abs(pR.B - pL.B);
    double deltx2 = (dxR + dxG + dxB) / 2.0;

    // find partial derivative for y
    int dyR = abs(pD.R - pU.R);
    int dyG = abs(pD.G - pU.G);
    int dyB = abs(pD.B - pU.B);
    double delty2 = (dyR + dyG + dyB) / 2.0;

    // return magnitude for gradient magnitude
    return sqrt(deltx2 + delty2);
}

// takes in an image, and return an energy map calculated by gradient magnitude
void compute_E(pixel_t** image_pixel_array, double* E, int num_rows, int num_cols) {
    // TODO: flatten the array; need to work on indexing; 

    // int i, j;
    int i ;
    double current_val;
    double temp_array[8];       // store 8 numbers
    int temp_counter = 0;       // 8 number counter; 
    int start_idx;              // starting index of the 8 numbers
    for (i = 0; i < num_rows; i++) {
        int j;

        #if OMP 
            #pragma omp parallel for
        #endif
        // TODO: need to separate this into better chunks 

            for (j = 0; j < num_cols; j++) {
                // don't want to remove the edges
                if (i == num_rows - 1 || j == num_cols - 1) {
                    temp_array[temp_counter] = MAX_ENERGY;
                } else {
                    temp_array[temp_counter] = pixel_difference(image_pixel_array[i][j],
                                           image_pixel_array[i + 1][j],
                                           image_pixel_array[i][j],
                                           image_pixel_array[i][j + 1]);
                }
                
                #if BATCH
                    if (temp_counter == 0){
                        start_idx = i * num_cols + j; 
                    }
                    // store thigns at location rep by current counter;
                    // temp_array[temp_counter] = current_val;
                    temp_counter ++; 
                    if (temp_counter == 8){
                        // flush to memory; 
                        int offset; 
                        for (offset = 0; offset < 8; offset ++){
                            E[offset + start_idx] = temp_array[offset];
                        }
                        // printf(" ===== flushing - %ld bytes total\n", sizeof(double) * 8);
                        // reset counter ;
                        temp_counter = 0;
                    }
                #else 
                    E[i * num_cols + j] = current_val; 
                #endif    
            } 
    }
    // // debug: print out gradient file;
    // char* output_file = "gradient.ppm";
    // intermediary_img(E, output_file,  num_rows, num_cols, max_pix_val, min_pix_val);
}

// M(i, j) = E(i, j) + min(M(i - 1, j - 1), M(i - 1, j), M(i - 1, j + 1))
void compute_M(double* E, double** M, int num_rows, int num_cols) {
    int i;
    memcpy(M[0], E, sizeof(double) * num_cols);// copy the first row;

    for (i = 0; i < num_rows; i++) {
        int j;

        #if OMP 
        #pragma omp parallel for
        #endif

        for (j = 0; j < num_cols; j++) {
            double middle = M[i][j];
            
            double left;
            double right;
            if (j - 1 < 0) {
                left = MAX_ENERGY;
            } else {
                left = M[i][j - 1];
            }

            if (j + 1 >= num_cols) {
                right = MAX_ENERGY;
            } else {
                right = M[i][j + 1];
            }

            if (i < num_rows - 1){
                // M[i + 1][j] = E[i + 1][j] + fmin(middle, fmin(left, right));
                int new_row  = i + 1;
                M[i + 1][j] = E[new_row * num_cols + j] + fmin(middle, fmin(left, right));
            }
        }
    }
}

// finds the seam from M
void find_seam(double** M, int* seam_path, int num_rows, int num_cols) {
    // find the min seam cost col in the last row
    double* last_row = M[num_rows - 1];
    double MAX_ENERGY_SUM = MAX_ENERGY*num_rows;

    int j;
    int min_cost_col = 0;
    int min_cost = last_row[0];
    #if OMP 
        #pragma omp parallel for
    #endif
    for (j = 1; j < num_cols; j++) {
        int current_cost = last_row[j];
        if (current_cost < min_cost) {
            min_cost = current_cost;
            min_cost_col = j;
        } 
    }
    M[num_rows - 1][min_cost_col] = MAX_ENERGY_SUM;

    // go up from the bottom and find the small cost path
    int i;
    seam_path[num_rows - 1] = min_cost_col;
    #if OMP 
        #pragma omp parallel for
    #endif
    for (i = num_rows - 2; i >= 0; i--) {
        int prev_col = seam_path[i + 1];
        double middle = M[i][prev_col];

        double left;
        double right;
        if (prev_col - 1 < 0) {
            left = MAX_ENERGY_SUM;
        } else {
            left = M[i][prev_col - 1];
        }
        if (prev_col + 1 >= num_cols) {
            right = MAX_ENERGY_SUM;
        } else {
            right = M[i][prev_col + 1];
        }
        double current_min_cost = fmin(left, fmin(middle, right));

        // printf("left: %f, middle: %f, right: %f, prev_col: %d\n", left, middle, right, prev_col);

        if (current_min_cost == middle) {
            seam_path[i] = prev_col;
            M[i][prev_col] = MAX_ENERGY_SUM;
        } else if (current_min_cost == left) {
            seam_path[i] = prev_col - 1;
            M[i][prev_col - 1] = MAX_ENERGY_SUM;
        } else {
            seam_path[i] = prev_col + 1;
            M[i][prev_col + 1] = MAX_ENERGY_SUM;
        }
    }

    // for (i = 0; i < num_rows; i++) {
    //     printf("seam index: %d with value %d\n", i, seam_path[i]);
    // }
}

// colors the seam pixels red and outputs the image
void color_seam(pixel_t*** image_pixel_array, double** M, int* seam_path, int num_rows, int num_cols, int max_px_val, char* seam_file) {
    int i;
    for (i = 0; i < num_rows; i++) {
        (*image_pixel_array)[i][seam_path[i]].R = 255;
        (*image_pixel_array)[i][seam_path[i]].G = 0;
        (*image_pixel_array)[i][seam_path[i]].B = 0;
    }

    // int x;
    // for (x = 0; x < 300; x++) {
    //     find_seam(M, seam_path, num_rows, num_cols);
    //     for (i = 0; i < num_rows; i++) {
    //         (*image_pixel_array)[i][seam_path[i]].R = 255;
    //         (*image_pixel_array)[i][seam_path[i]].G = 0;
    //         (*image_pixel_array)[i][seam_path[i]].B = 0;
    //     }
    // }

    // call function to output image
    output_image(*image_pixel_array, seam_file, num_rows, num_cols, max_px_val);
}

// remove the seam from the image pixel array and get new dimensions
void remove_seam(pixel_t*** image_pixel_array_pt, int* seam_path, int* rows, int* cols) {
    int num_rows = *rows;
    int num_cols = *cols;

    // need to shift every pixel after the seam over to the "left"
    int i;
    #if OMP 
        #pragma omp parallel for
    #endif
    for (i = 0; i < num_rows; i++) {
        int j;
        int seam_col = seam_path[i];
        #if OMP 
            #pragma omp parallel for
        #endif
        // seam col guaranteed to be at least 0
        for (j = seam_col + 1; j < num_cols; j++) {
            (*image_pixel_array_pt)[i][j - 1] = (*image_pixel_array_pt)[i][j];
        }
    }

    *cols = num_cols - 1;
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

void intermediary_img(double ** matrix, char* output_file,  \
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
