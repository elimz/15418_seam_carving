#include "sort.h"

void swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

void swap_int(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

 
int partition_double (double pixel_arr[], int index_arr[], int low, int high)
{
    int pivot = pixel_arr[high];    // pivot
    int i = (low - 1);  // Index of smaller element
    
    int j;
    for (j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (pixel_arr[j] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&pixel_arr[i], &pixel_arr[j]);
            swap_int(&index_arr[i], &index_arr[j]);
        }
    }
    swap(&pixel_arr[i + 1], &pixel_arr[high]);
    swap_int(&index_arr[i + 1], &index_arr[high]);
    return (i + 1);
}
 
void quickSort_double (double pixel_arr[], int index_arr[], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition_double(pixel_arr, index_arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort_double(pixel_arr, index_arr, low, pi - 1);
        quickSort_double(pixel_arr, index_arr, pi + 1, high);
    }
}


int partition_int (int index_arr[], int low, int high)
{
    int pivot = index_arr[high];    // pivot
    int i = (low - 1);  // Index of smaller element
    
    int j;
    for (j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (index_arr[j] <= pivot)
        {
            i++;    // increment index of smaller element
            swap_int(&index_arr[i], &index_arr[j]);
        }
    }
    swap_int(&index_arr[i + 1], &index_arr[high]);
    return (i + 1);
}
 
void quickSort_int (int index_arr[], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition_int(index_arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort_int(index_arr, low, pi - 1);
        quickSort_int(index_arr, pi + 1, high);
    }
}