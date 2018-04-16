
// seam carving application in parallel
// steps: 
// 1. command line, convert jpg/png format into ppm; 
//  "convert source.png -compress none dest.ppm" // require imageMagick

#ifndef SEAM_H
#include "seam.h"
#define SEAM_H
#endif

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
    // fgets(pixel_buf, 1500 , ppm_file);
   
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

        // if (!strcmp(buf_ptr, "\n")){
        //     printf("here--- egde \n");
        //     buf_ptr ++;
        // }
        if (i == 10)
            break;
    }


    fclose(ppm_file);
    // return matrix;
    return 0;
}   


int main(){
    // char *file_path = "./images/tower.ppm";
    build_matrix();
    
    return 0;
}

// // this function parses the ppm file and put all the pixels RGB values in a 
// //  2D matrix
// int build_matrix(FIL){
//     return 0;
// }   


// // takes in an image, and return an energy map calculated by gradient magnitude
// int find_energy_map(){

// }
