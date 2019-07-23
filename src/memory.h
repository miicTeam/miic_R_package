#ifndef MEMORY_H
#define MEMORY_H

typedef struct {
	int maxlevel;
	int** sample;
	int** sortedSample;
	int** Opt_sortedSample;
	int* orderSample;
	int* sampleKey;
	int* Nxyuiz;
	int* Nyuiz;
	int* Nuiz;
	int* Nz;
	int* Ny;
	int* Nxui;
	int* Nx;
	int** Nxuiz;
	int* bridge;
	double* Pxyuiz;

	//continuous data
	int* samplesToEvaluate; 
	int* samplesToEvaluateTemplate;

	int** dataNumeric_red; 
	int** dataNumericIdx_red;

	int* AllLevels_red;
	int* cnt_red;
	int* posArray_red;

} MemorySpace;


// struct ContainerMemory{ 
// 	int** sortedSample;
// 	int** Nxuiz;
// 	int* Ny;
// 	int* Nxui;
// 	int* Nx;
// 	int* Nxyuiz;
// 	int* Nyuiz;
// 	int* Nuiz;
// 	int* Nz;
// 	int* bridge;
// 	double* Pxyuiz;

// 	//continuous data
// 	int* samplesToEvaluate; 
// 	int* samplesToEvaluateTemplate;

// 	int** dataNumeric_red; 
// 	int** dataNumericIdx_red;

// 	int* AllLevels_red;
// 	int* cnt_red;
// 	int* posArray_red;

// };
#endif