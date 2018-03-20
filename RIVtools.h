#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef RIVSIZE
#define RIVSIZE 25000
#endif

#ifndef NONZEROS
#define NONZEROS 2
#endif

#ifndef CACHESIZE
#define CACHESIZE 20
#endif

/* the sparseRIV is a RIV form optimized for RIVs that will be mostly 0s
 * as this is often an ideal case, it is adviseable as the default 
 * unless we are doing long term RIV aggregation.
 * specifically, a sparseRIV contains a pair of arrays, 
 * containing locations and values, where pairs are found in like array 
 * indices.
 */
typedef struct{
	char name[100];
	int *values;
	int *locations;
	size_t count;
	unsigned int frequency;
	float magnitude;
	int boolean;
}sparseRIV;
/* the denseRIV is a RIV form optimized for overwhelmingly non-0 vectors
 * this is rarely the case, but its primary use is for performing vector
 * math, as comparisons and arithmetic between vectors are ideally 
 * performed between sparse and dense (hetero-arithmetic)
 */
typedef struct{
	char name[100];
	int* values;
	int* frequency;
	float magnitude;
}denseRIV;
/*
 * RIVKey, holds globally important data that should not be changed partway through
* first function call in the program should always be: 
* setKeyData(RIVsize, nonZeros, blocksize)
* this will set these variables, check for incompatible choices, and open up 
* memory blocks which the program will use
*/
struct RIVData{
	size_t RIVsize;
	int nonZeros;
	int *h_tempBlock;
	int tempSize;
	int thing;
	denseRIV* RIVCache;
	int cacheSize;
}static RIVKey;

/*setKeyData should be the first function called in any usage of this library
 * it sets global variables that practically all functions will reference,
 * it checks that your base parameters are valid, and allocates memory for
 * the functions to use, so that we can move fast with rare allocations.
 * it should only be called once, as changing the RIVsize and nonZeros, variables 
 * mid operation will cause a breakdown of the underlying logic. 
 * do so at your own peril
 */
void setKeyData();
/*consolidateD2S takes a denseRIV value-set input, and returns a sparse RIV with
 * all 0s removed. it does nto automatically carry metadata, which must be assigned
 * to a denseRIV after the fact.  often denseRIVs are only temporary, and don't
 * need to carry metadata
 */

sparseRIV consolidateD2S(int *denseInput);  //fixthis
/* mapS2D expands a sparseRIV out to denseRIV values, filling array locations
 * based on location-value pairs 
 */
int* mapS2D(int* destination, sparseRIV input);
/* makeSparseLocations produces a series of locations out of a series of seeds
 * this produces an "implicit" RIV which can be used with the mapI2D function
 * to create a denseRIV.  seeds should be produced using the makeSeeds function
 */
int* makeSparseLocations(int *seeds, size_t seedCount);
void makeSeed(char* word, int *seeds, size_t seedCount);
int* mapI2D(int *locations, size_t seedCount);

sparseRIV fileToL2(FILE *input);
/* fileToL2Clean operates the same as fileToL2 butkeeps only words 
 * containing lowercase letters and the '_' symbol
 * this is important if you will be lexPush-ing those words later
 */

void cosineCompare(sparseRIV baseRIV, sparseRIV *multipliers, size_t multiplierCount, float threshold);
void getMagnitudes(sparseRIV *inputs, size_t RIVCount);

 sparseRIV fileToL2(FILE *data){
	unsigned int blockSize;
	
	char word[100] = {0};
	int *seeds = RIVKey.h_tempBlock;
	int seedCount = 0;
	while(fgets(word, 100, data)){
		
		if(feof(data)){
			break;
		}
		if(!(*word)){
			break;
		}	
		blockSize = ((seedCount+1)* RIVKey.nonZeros);
		if(blockSize>RIVKey.tempSize){
			RIVKey.h_tempBlock = (int*) realloc(RIVKey.h_tempBlock, blockSize*sizeof(int));
			seeds = RIVKey.h_tempBlock;
			RIVKey.tempSize+=RIVKey.nonZeros;
		}
		
		makeSeed(word, seeds, seedCount);
		seedCount++;
		
	}
	
	int *locations = makeSparseLocations(seeds, seedCount);
	
	int *L2dense;
	L2dense = mapI2D(locations, seedCount);
		
	sparseRIV output = consolidateD2S(L2dense);	
	free(L2dense);
	output.frequency = seedCount/RIVKey.nonZeros;
	output.boolean = 1;
	return output;
}

void cosineCompare(sparseRIV baseRIV, sparseRIV *multipliers, size_t multiplierCount, float threshold){
	
	int *baseDenseRIV = RIVKey.h_tempBlock;
	mapS2D(baseDenseRIV, baseRIV);
	float cosSim;
	sparseRIV *multipliersStop = multipliers+multiplierCount;
	float minsize = baseRIV.magnitude * .85;
	float maxsize = baseRIV.magnitude * 1.15;
	int dot = 0;
	int *values;
	int *locations;
	int *locations_Stop;
	while(multipliers<multipliersStop){
		if(((*multipliers).boolean) 
		&& (((*multipliers).magnitude < maxsize) 
		&& ((*multipliers).magnitude > minsize))){
			dot = 0;
			values = (*multipliers).values;
			locations = (*multipliers).locations;
			locations_Stop = locations+(*multipliers).count;
			
			while(locations<locations_Stop){
				
				dot += (*values)*(*(baseDenseRIV+(*locations)));
				locations++;
				values++;
			}
			cosSim= dot/((baseRIV.magnitude)*((*multipliers).magnitude));
			if(cosSim>=threshold){
				printf("%s\t%s\n%f\n", (*multipliers).name, baseRIV.name, cosSim);
				(*multipliers).boolean = 0;
				RIVKey.thing ++;
			}
		}
		multipliers++;
		
	}
}

void getMagnitudes(sparseRIV *inputs, size_t RIVCount){
	for(int i=0; i<RIVCount; i++){
		unsigned int temp = 0;
		int *values = inputs[i].values;
		int *values_stop = values+inputs[i].count;
		while(values<values_stop){	
			temp += (*values)*(*values);
			values++;
			
		}
		float magnitude = sqrt(temp);
		inputs[i].magnitude = magnitude;
	}
}

int* mapS2D(int* destination, sparseRIV input){
	memset(destination, 0, RIVKey.RIVsize*sizeof(int));
	
	int *locations_slider = input.locations;
	int *values_slider = input.values;
	int *locations_stop = locations_slider+input.count;
	
	while(locations_slider<locations_stop){
		destination[*locations_slider] = *values_slider;
		locations_slider++;
		values_slider++;
	}
	
	return destination;
}

int* mapI2D(int *locations, size_t valueCount){
	int *destination = (int*)calloc(RIVKey.RIVsize,sizeof(int));
	int *locations_slider = locations;
	int *locations_stop = locations_slider+(valueCount*RIVKey.nonZeros);
	while(locations_slider<locations_stop){
	
		destination[*locations_slider] +=1;
		locations_slider++;
		destination[*locations_slider] -= 1;
		locations_slider++;
	}
	
	
	return destination;
}
sparseRIV consolidateD2S(int *denseInput){
	sparseRIV output;
	output.count = 0;
	int* locations = RIVKey.h_tempBlock;
	int* values = RIVKey.h_tempBlock+RIVKey.RIVsize;
	int* locations_slider = locations;
	int* values_slider = values;
	for(int i=0; i<RIVKey.RIVsize; i++){
		if(denseInput[i]){
			*(locations_slider++) = i;
			*(values_slider++) = denseInput[i];
			output.count++;
		}
	}
	output.locations = (int*) malloc(output.count*sizeof(int));
	if(!output.locations){
		printf("broke at locations malloc");
	}
	memcpy(output.locations, locations, output.count*sizeof(int));

	output.values = (int*) malloc(output.count*sizeof(int));
	if(!output.locations){
		printf("broke at values malloc");
	}
	memcpy(output.values, values, output.count*sizeof(int));

	return output;
}


void setKeyData(){
	RIVKey.RIVsize = RIVSIZE;
	RIVKey.nonZeros = NONZEROS;
	if(RIVKey.nonZeros%2){
		printf("your NONZEROS value must be an even number");
		RIVKey.nonZeros++;
		printf(", changed to %d", RIVKey.nonZeros);
	}

	RIVKey.h_tempBlock = (int*)malloc(3*RIVKey.RIVsize*sizeof(int));
	RIVKey.tempSize = 3*RIVKey.RIVsize;
	RIVKey.thing = 0;
	RIVKey.cacheSize = CACHESIZE;
	RIVKey.RIVCache = (denseRIV*)calloc(RIVKey.cacheSize,sizeof(denseRIV));
}


void makeSeed(char* word,  int *seeds, size_t seedCount){

	int i=0;
	int seed = 0;
	while(*word){
		seed += (*(word))<<(i*5);
		word++;
		i++;
	}
	*(seeds + (seedCount*RIVKey.nonZeros)) = seed;
	return;
}
int* makeSparseLocations(int* seeds, size_t seedCount){

	int *locations = RIVKey.h_tempBlock;
	int *locations_slider = locations;
	int *seeds_stop = seeds+seedCount*RIVKey.nonZeros;

	while(seeds<seeds_stop){
		srand(*seeds);
		for(int i=0; i<RIVKey.nonZeros; i++){
			*locations_slider = rand() % RIVKey.RIVsize;
			locations_slider++;
		}
		seeds+=RIVKey.nonZeros;
	}


	return locations;
}


