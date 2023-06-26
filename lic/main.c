#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


void ScanfVector(int xres, int yres, float*  pVectr, int* pVectrFlag) {
    FILE *fp;
    if ((fp = fopen("./fielddataglobal.txt", "r")) == NULL) {
        printf("error in reading file !\n");
        exit(1);
    }
    float f1, f2, f3, f4;
    int index = 0;
    while (fscanf(fp, "%f %f %f %f", &f1, &f2, &f3, &f4) != EOF) {
        pVectr[index] = f3;
        pVectrFlag[index] = ((int)(f3-9999)==0)?1:0;
        index++;
        pVectr[index] = f4;
        pVectrFlag[index] = ((int)(f4-9999)==0)?1:0;
        index++;
    }
    fclose(fp);
}

long int GetFileSize(FILE* fp){
    long int size;
    fseek(fp,0,SEEK_END);
    size = ftell(fp);
    rewind(fp);
    return size;
}

void ReadVector(int xres, int yres, float*  pVectr, int* pVectrFlag){
    FILE *fp;
    if( (fp = fopen("./fielddataglobal.txt","r")) == NULL ){
        printf("open file failed!\n");
        exit(1);
    }
    long int fileSize = GetFileSize(fp);

    char* buffer = malloc(fileSize+1);
    buffer[fileSize] = 0;       //添加结束符
    int readSize = fread(buffer,sizeof(char),fileSize,fp);
    fclose(fp);
    float f[4];
    char*p = strtok(buffer,"\t");
    int x = 0,y = 0;
    int index = 0;
    float f1, f2, f3, f4;
    while(p){
        f1 = atof(p); p = strtok(NULL,"\t");
        f2 = atof(p); p = strtok(NULL,"\t");
        f3 = atof(p); p = strtok(NULL,"\n");
        f4 = atof(p); p = strtok(NULL,"\t");
        pVectr[index] = f3;
        pVectrFlag[index] = ((int)(f3-9999)==0)?1:0;
        index++;
        pVectr[index] = f4;
        pVectrFlag[index] = ((int)(f4-9999)==0)?1:0;
        index++;
    }
    free(fp);
}


int main(int argc, char** argv) {
    clock_t end,begin;
    float left = -180, right = 179.75, low =-80, high = 79.75;
    float res=0.25;
    int n_xres = (right-left)/res+1, n_yres = (high-low)/res+1;
    float* pVectr = (float*) malloc(sizeof(float)*n_xres*n_yres*2);
    int* pVectrFlag = (int*) calloc(sizeof(int),n_xres*n_yres*2);

    begin = clock();
    ScanfVector(n_xres,n_yres,pVectr,pVectrFlag);
    end = clock();
    printf("ScanfVector Time= %ld\n",end-begin);

    begin = clock();
    ReadVector(n_xres,n_yres,pVectr,pVectrFlag);
    end = clock();
    printf("ReadVector  Time= %ld\n",end-begin);



    return 0;
}
