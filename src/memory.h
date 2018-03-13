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
} MemorySpace;



struct ContainerMemory{ 
	int** sortedSample;
	int** Nxuiz;
	int* Ny;
	int* Nxui;
	int* Nx;
	int* Nxyuiz;
	int* Nyuiz;
	int* Nuiz;
	int* Nz;
	int* bridge;
};
#endif
