

#include <stdio.h>
#include <stdlib.h>

int main(){
    FILE *fp = fopen("output.txt","wb");
    int num_cols = 1428;
    int num_rows = 968;
    int max_px_val = 255;

    // write file header first; 
    fprintf(fp, "P3\n%d %d\n%d\n", num_cols, num_rows, max_px_val);
    // fprintf(fp, "writing - name:%s, score = %d\n", "bar", 90);
    fclose(fp);
    printf("done\n");
    return 1;
}   