
// seam.c on seq_1D branch, inside week3_para folder 
// with the most up-to-date algorithm from branch week 3
// now adding timing, and omp parallel


// (before changing to 1D array, batch write and 1D block assignment)
// seq baseline:            3.370 seconds
// seq new compute_E (1D, batch write, 1D block assignemtn)
//                          3.263 seconds
// 

#include "seam.h"

#define MAX_ENERGY 9999999.0
#define NUM_SEAMS_TO_REMOVE 300
// #define NUM_SEAMS_TO_TRY 20

// #ifndef OMP  
// #define OMP 1
// #endif

#ifndef TIMING 
#define TIMING 1
#endif

#ifndef BLOCK_ASSIGN 
#define BLOCK_ASSIGN 1
#endif

char* input_file = "../images/tower.ppm";
char* output_file = "output.ppm";
char* seam_file = "output_seam.ppm";

// number of threads
int nthread = 8;

int main(){

    int num_rows, num_cols, original_cols, max_px_val;
    time_t t;

    // pass in pointers for image array pointer, width, height
    pixel_t** image_pixel_array = build_matrix(&num_rows, &num_cols, &max_px_val, input_file);
    original_cols = num_cols;

    // alloc for E array - 1D
    int row; 
    double* E = malloc(sizeof(double) * (num_rows * num_cols));

    // num seams to try in each iteration
    // limit max num seams to try to 200
    int NUM_SEAMS_TO_TRY = num_cols / 3;
    if (NUM_SEAMS_TO_TRY > 200) {
        NUM_SEAMS_TO_TRY = 200;
    }

    // array of random indicies
    int first_row_indices[NUM_SEAMS_TO_TRY];

    // alloc mem for seam_path array
    int seam_num;
    int** seam_paths = malloc(sizeof(int*) * NUM_SEAMS_TO_TRY);
    for (seam_num = 0; seam_num < NUM_SEAMS_TO_TRY; seam_num++) {
        seam_paths[seam_num] = malloc(sizeof(int) * num_rows);
    }

    srand((unsigned) time(&t));

    // array of timing data for different sections
    #if TIMING
    double start = currentSeconds();
    double timing[2] = {0.0, 0.0};
    #define T_COMP_E            0       // compute energy map
    #define T_FIND_SEAM         1       // find seam
    #endif

    #if OMP 
        omp_set_num_threads(nthread);
    #endif


    // remove NUM_SEAMS_TO_REMOVE number of lowest cost seams
    for (seam_num = 0; seam_num < NUM_SEAMS_TO_REMOVE; seam_num++) {

        // compute energy map and store in engery_array
        #if TIMING
        double t_start_E = currentSeconds();
        #endif

        compute_E(image_pixel_array, E, num_rows, num_cols);
        
        #if TIMING
        double t_end_E = currentSeconds();
        timing[T_COMP_E] += (t_end_E - t_start_E);
        #endif

        // randomly select NUM_SEAMS_TO_TRY number of seams, find their energy, 
        //  and remove the seam with lowest energy
        int idx;
        for (idx = 0; idx < NUM_SEAMS_TO_TRY; idx++) {
            first_row_indices[idx] = rand() % (num_cols - 1);
        }

        // New: each thread keeps a copy of the smallest cost path, and use omp 
        // reduce to find the smallest of all
        // int start_pos_partition (int N, int P, int i)
        int trial_num;
        double smallest_cost_path = INT_MAX * 1.0;
        int smallest_cost_index = -1;

        // #if !OMP 
        //     // each thread keeps a copy of the smallest cost path, and use omp
        //     // reduce to find the smallest of all
        //     // #pragma omp parallel {
        //         int my_tid = omp_get_thread_num();
        //         int start_num = start_pos_partition(NUM_SEAMS_TO_TRY, nthread, my_tid);
        //         int end_num = (my_tid == nthread - 1)? NUM_SEAMS_TO_TRY : start_pos_partition(NUM_SEAMS_TO_TRY, nthread, my_tid + 1);
        //         printf("-------- tid = %d, start = %d, end = %d\n");

        //         for (trial_num = start_num; trial_num < end_num; trial_num ++){
        //             seam_paths[trial_num][0] = first_row_indices[trial_num];
        //             double path_cost = find_seam(E, seam_paths[trial_num], num_rows, num_cols);
        //             if (path_cost < smallest_cost_path) {
        //                 smallest_cost_index = trial_num;
        //             }
        //         }
        //     // }
        // #else
            // non-OMP version, 1 thread calculates all of the seams;
            for (trial_num = 0; trial_num < NUM_SEAMS_TO_TRY; trial_num++) {
                seam_paths[trial_num][0] = first_row_indices[trial_num];
                double path_cost = find_seam(E, seam_paths[trial_num], num_rows, num_cols);
                if (path_cost < smallest_cost_path) {
                    smallest_cost_index = trial_num;
                }
            }
        // #endif


        #if TIMING
        double t_end_seam = currentSeconds();
        timing[T_FIND_SEAM] += (t_end_seam - t_end_E);
        #endif

        // color the seam and output the image
        // color_seam(&image_pixel_array, seam_paths[smallest_cost_index], num_rows, num_cols, max_px_val, seam_file);
        // printf("Finished coloring seam\n");
        
        // remove the seam from the image, also sets new values for num_rows and num_cols
        remove_seam(&image_pixel_array, seam_paths[smallest_cost_index], &num_rows, &num_cols);
        // printf("Finished removing seam\n");
    }


    double delta = currentSeconds() - start;
    printf("%d seams of a %dx%d image removed in %.3f seconds\n", NUM_SEAMS_TO_REMOVE, original_cols, num_rows, delta);

    #if TIMING
    printf("------ TIMING SPLITS ------\n");
    printf("    T_COMP_E = %f\n    T_FIND_SEAM = %f\n ", timing[T_COMP_E], timing[T_FIND_SEAM]);
    #endif

    output_image(image_pixel_array, output_file, num_rows, num_cols, max_px_val);
    
    // free data structures 
    // for (row = 0; row < num_rows; row++) {
    //     free(E[row]);
    // }

    for (seam_num = 0; seam_num < NUM_SEAMS_TO_TRY; seam_num++) {
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

// code from 15418 Spring 2018 assignment 4 write up
int start_pos_partition (int N, int P, int i){
    int base = N / P; 
    int extra = N % P; 
    if (i < extra)
        return i * (base + 1);
    else
        return i * base + extra;
}

// // takes in an image, and return an energy map calculated by gradient magnitude
// void compute_E(pixel_t** image_pixel_array, double** E, int num_rows, int num_cols) {
//     // debug: 
//     // int max_pix_val = 255; 
//     // int min_pix_val = 0;

//     int i;
//     for (i = 0; i < num_rows; i++) {
//         int j;
//         for (j = 0; j < num_cols; j++) {
//             // don't want to remove the edges
//             if (i == num_rows - 1 || j == num_cols - 1) {
//                 E[i][j] = MAX_ENERGY;
//             } else {
//                 E[i][j] = pixel_difference(image_pixel_array[i][j],
//                                            image_pixel_array[i + 1][j],
//                                            image_pixel_array[i][j],
//                                            image_pixel_array[i][j + 1]);
//             }
//             // debug 
//             // if (E[i][j] > max_pix_val){
//             //     max_pix_val = E[i][j];    
//             // } if (E[i][j] < min_pix_val){
//             //     min_pix_val = E[i][j];
//             // } 
//         }
//     }

//     // // debug: print out gradient file;
//     // char* output_file1 = "gradient.ppm";
//     // intermediary_img(E, output_file1,  num_rows, num_cols, max_pix_val, min_pix_val);
// }

// takes in an image, and return an energy map calculated by gradient magnitude
void compute_E(pixel_t** image_pixel_array, double* E, int num_rows, int num_cols) {
    // find my thread id, and work region
    #if OMP
        int i, j, my_tid;
        for (i = 0; i < num_rows; i++) {
            // New partition: block assignment, break down the task by total num_cols / nthread;
            #pragma omp parallel num_threads(nthread) shared (i) private (j, my_tid) 
            {
                my_tid = omp_get_thread_num();
                double current_val;
                double temp_array[8];       // store 8 numbers
                int temp_counter = 0;       // 8 number counter; 
                int store_start_idx;              // starting index of the 8 numbers
                
                int my_start = start_pos_partition(num_cols, nthread, my_tid);
                int my_end = (my_tid == nthread - 1)? num_cols: start_pos_partition(num_cols, nthread, my_tid + 1);
                // printf(" compute_E - my tid = %d, start = %d, end = %d\n",  my_tid, my_start, my_end);

                assert((my_start < my_end) && (my_start >= 0));
                for (j = my_start; j < my_end; j++) {
                    // don't want to remove the edge
                    assert(temp_counter >= 0 && temp_counter <= 7);
                    if (i == num_rows - 1 || j == num_cols - 1) {
                        temp_array[temp_counter] = MAX_ENERGY;
                    } else {
                        temp_array[temp_counter] = pixel_difference(image_pixel_array[i][j],
                                               image_pixel_array[i + 1][j],
                                               image_pixel_array[i][j],
                                               image_pixel_array[i][j + 1]);
                    }

                    if (temp_counter == 0){
                        store_start_idx = i * num_cols + j; 
                    }
                    // store thigns at location rep by current counter;
                    temp_counter ++; 
                    if (temp_counter == 8){
                        // flush to memory; 
                        int offset; 
                        for (offset = 0; offset < 8; offset ++){
                            E[offset + store_start_idx] = temp_array[offset];
                        }
                        // reset counter ;
                        temp_counter = 0;
                    }
                }
            }
        }
    #else 
        int i, j; 
        for (i = 0; i < num_rows; i++) {
            for (j = 0; j < num_cols; j++) {
                // don't want to remove the edges
                if (i == num_rows - 1 || j == num_cols - 1) {
                    // temp_array[temp_counter] = MAX_ENERGY;
                    E[i * num_cols + j] = MAX_ENERGY;
                } else {
                    E[i * num_cols + j] = pixel_difference(image_pixel_array[i][j],
                                           image_pixel_array[i + 1][j],
                                           image_pixel_array[i][j],
                                           image_pixel_array[i][j + 1]);
                }
            }
        }
    #endif
}




// finds the seam from E
double find_seam(double* E, int* seam_path, int num_rows, int num_cols) {

    // find min seams with least energy cost pixel in first row
    // traverse downwards and greedily pick the least energy cost neighbor

    int i;
    // double sum = E[0][seam_path[0]];
    double sum = E[seam_path[0]];
    for (i = 1; i < num_rows; i++) {
        int prev_col = seam_path[i - 1];
        // double middle = E[i][prev_col];
        double middle = E[i * num_cols + prev_col];

        double left;
        double right;
        if (prev_col - 1 < 0) {
            left = INT_MAX;
        } else {
            // left = E[i][prev_col - 1];
            left = E[i * num_cols + (prev_col - 1)];
        }
        if (prev_col + 1 >= num_cols) {
            right = INT_MAX;
        } else {
            // right = E[i][prev_col + 1];
            right = E[i * num_cols + (prev_col + 1)];
        }
        // record which of the 3 neighbors should be part of the seam
        double current_min_cost = fmin(middle, fmin(left, right));
        if (current_min_cost == middle) {
            seam_path[i] = prev_col;
        } else if (current_min_cost == left) {
            seam_path[i] = prev_col - 1;
        } else {
            seam_path[i] = prev_col + 1;
        }

        // sum += E[i][seam_path[i]];
        sum += E[i * num_cols + seam_path[i]];
    }

    return sum;

    // for (i = 0; i < num_rows; i++) {
    //     if (seam_path[i] < 0) {
    //         printf("seam index: %d with value %d\n", i, seam_path[i]);
    //     }
    // }
    // printf("\n\n");
}

// colors the seam pixels red and outputs the image
void color_seam(pixel_t*** image_pixel_array, int* seam_path, int num_rows, int num_cols, int max_px_val, char* seam_file) {
    int i;
    for (i = 0; i < num_rows; i++) {
        (*image_pixel_array)[i][seam_path[i]].R = 255;
        (*image_pixel_array)[i][seam_path[i]].G = 0;
        (*image_pixel_array)[i][seam_path[i]].B = 0;
    }

    // call function to output image
    output_image(*image_pixel_array, seam_file, num_rows, num_cols, max_px_val);
}

// remove the seam from the image pixel array and get new dimensions
void remove_seam(pixel_t*** image_pixel_array_pt, int* seam_path, int* rows, int* cols) {
    int num_rows = *rows;
    int num_cols = *cols;

    // need to shift every pixel after the seam over to the "left"
    int i;
    for (i = 0; i < num_rows; i++) {
        int j;
        int seam_col = seam_path[i];
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
    printf("Finished - writing image to output.\n");
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

