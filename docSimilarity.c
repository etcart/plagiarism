#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include <time.h>
#define RIVSIZE 5000
#include "RIVtools.h"

#define THRESHOLD .90f
void directoryToL2s(char *rootString, sparseRIV** fileRIVs, int *fileCount);
	
int main(int argc, char *argv[]){
	int fileCount = 0;
	setKeyData();
	sparseRIV *fileRIVs = (sparseRIV*) malloc(1*sizeof(sparseRIV));
	char rootString[2000];
	if(argc <2){ 
		printf("give me a directory");
		return 1;
	}
	strcpy(rootString, argv[1]);
	strcat(rootString, "/");

	directoryToL2s(rootString, &fileRIVs, &fileCount);
	printf("fileCount: %d\n", fileCount);
	getMagnitudes(fileRIVs, fileCount); 

	for(int i=0; i<fileCount; i++){
		if(fileRIVs[i].boolean){
			cosineCompare(fileRIVs[i], fileRIVs+i+1, fileCount-(i+1), THRESHOLD);
		}
		
	}
	
	printf("number of pairs found: %d", RIVKey.thing);
	free(fileRIVs);
	
return 0;
}
void directoryToL2s(char *rootString, sparseRIV** fileRIVs, int *fileCount){
	int length;
	char pathString[2000];
	DIR *directory;
    struct dirent *files = 0;
	
	if(!(directory = opendir(rootString))){
		printf("location not found, %s\n", rootString);
		return;
	}
	
	while((files=readdir(directory))){
		if(*(files->d_name) == '.') continue;
	
		if(files->d_type == DT_DIR){
			strcpy(pathString, rootString);
			
			strcat(pathString, files->d_name);
			strcat(pathString, "/");
			directoryToL2s(pathString, fileRIVs, fileCount);
		}
			
		length = strlen(files->d_name);
		if(strcmp((files->d_name+length-2),".s")){
			continue;
		}
		strcpy(pathString, rootString);
		strcat(pathString, files->d_name);
		FILE *input = fopen(pathString, "r");
		if(!input){
			printf("file %s doesn't seem to exist, breaking out of loop", pathString);
			return;
		}else{
			
			(*fileRIVs) = (sparseRIV*)realloc((*fileRIVs), ((*fileCount)+1)*sizeof(sparseRIV));
			(*fileRIVs)[(*fileCount)] = fileToL2(input);
			strcpy((*fileRIVs)[(*fileCount)].name, pathString);
			
			fclose(input);
			(*fileCount)++;
			
		}	
	}
}
