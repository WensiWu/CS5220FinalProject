#include<stdio.h>


void write_array_double(const char* fname,int n,  double* a);

void write_array_int(const char* fname,int n,  int* a);



void write_array_double(const char* fname,int n,  double* a)
{
    FILE* fp = fopen(fname, "w+");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", fname);
        exit(-1);
    }
    for (int i = 0; i < n; ++i) {
            fprintf(fp, "%g ", a[i]);
    }
    fclose(fp);
}

void write_array_int(const char* fname, int n,  int* a)
{
    FILE* fp = fopen(fname, "w+");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", fname);
        exit(-1);
    }
    for (int i = 0; i < n; ++i) {
            fprintf(fp, "%d ", a[i]);
    }
    fclose(fp);
}


