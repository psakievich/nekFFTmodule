#ifndef PSVTKOUTPUT_
#define PSVTKOUTPUT_

extern "C" {
void cpphello_();
void dwritevtsc_(int* nPartition, int* nSize, int* nFlds, double* dPnts, double* dRFlds,double* dCFlds,char* chFileName, int nLenFN);
void dwritevts_(int* nPartition, int* nSize, int* nFlds, double* dPnts, double* dFlds,char* chFileName, int nLenFN);
}

#endif
