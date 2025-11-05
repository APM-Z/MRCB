#ifndef MRCB_SolHead_H
#define MRCB_SolHead_H

#include "MRCB.h"

#define NUMSYS      8
int creatFile(FILE** fp, const filopt_t* fopt, const prcopt_t* popt, const solopt_t* sopt, char* filename, char* profilename, char* logfilename);
#endif
