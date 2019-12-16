#ifndef IO_H_
#define IO_H_

void printDoubleMatrix(double** matrix, int x, int y, int process_rank);
void printLinearArray(double* array, int x, int y);
Map* readMapFile(const char* filename, int debug);
Parameters* readParameterFile(const char* filename);
void printUsefulMapInformation(Map* map);
void writeResultMatrix(char* filename, int xsize, int ysize,double** matrix, int debug);
void writeResultArray(char* filename, int xsize, int ysize,double* array, int debug);
void writeTestMap(char* filename, int debug);
void printGrid(Map* map);
void getFileNames(char* etaName, char* uName, char* vName, char* dir_name, unsigned int iteration);
int saveToDisk(double* etaTotal, double* uTotal, double* vTotal, unsigned int xSize, unsigned int ySize, unsigned int iteration, Parameters* params);
#endif  // IO_H_