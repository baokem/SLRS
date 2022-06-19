#ifndef aligningFromBam_H_INCLUDED 
#define aligningFromBam_H_INCLUDED 
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "contig.h"

using namespace std;
using namespace BamTools;

typedef struct AligningResult{
	string barcodeid;
	int leftContigIndex;
	long int leftPosition;
	long int rightPosition;
}AligningResult;

typedef struct AligningResultHead {
	AligningResult * aligningResult;
	int aligningCount;
}AligningResultHead;

AligningResultHead* GetAligningResultHead(char* readFile, char* leftBamFile, char* rightBamFile, char* alignmentResultFile, ContigSetHead* contigSetHead, bool aligning, int contigLengthThreshold, int endlength);

long int MIN(long int a, long int b);

long int MAX(long int a, long int b);

int comp(const void* x, const void* y);


#endif