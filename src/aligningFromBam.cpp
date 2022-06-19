#ifndef aligningFromBam_CPP_INCLUDED 
#define aligningFromBam_CPP_INCLUDED 
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>
#include <bits/stdc++.h>
#include "aligningFromBam.h"

using namespace std;

long int MIN(long int a, long int b) {
	if (a > b) {
		return b;
	}
	else {
		return a;
	}
}

long int MAX(long int a, long int b) {
	if (a < b) {
		return b;
	}
	else {
		return a;
	}
}

AligningResultHead * GetAligningResultHead(char * readFile, char * leftBamFile, char * rightBamFile, char * alignmentResultFile,ContigSetHead * contigSetHead, bool aligning, int contigLengthThreshold, int endlength){
	
	FILE * fp;
	string readname,readname1,barcodeid;
	int aligningCount = 0;

	if (aligning == false) {
		if ((fp = fopen(alignmentResultFile, "w")) == NULL) {
			printf("%s, does not exist!", alignmentResultFile);
			exit(0);
		}
		BamReader bamReaderLeft;
		BamReader bamReaderRight;

		bamReaderLeft.Open(leftBamFile);
		bamReaderRight.Open(rightBamFile);

		BamAlignment alignmentLeft;
		BamAlignment alignmentRight;

		long int minScore = 30;
		long int leftRefID = -1;
		long int rightRefID = -1;
		long int leftPosition = -1;
		long int rightPosition = -1;
		long int leftOrientation = -1;
		long int rightOrientation = -1;
		while (bamReaderLeft.GetNextAlignment(alignmentLeft) && bamReaderRight.GetNextAlignment(alignmentRight)) {

			while ((alignmentLeft.AlignmentFlag & 0x900) != 0) {
				bamReaderLeft.GetNextAlignment(alignmentLeft);
				continue;
			}
			while ((alignmentRight.AlignmentFlag & 0x900) != 0) {
				bamReaderRight.GetNextAlignment(alignmentRight);
				continue;
			}
			long int ContigLength = contigSetHead->contigSet[alignmentLeft.RefID].contigLength;
			
			long int TailLength = MIN(endlength, ContigLength / 2);
			if ((alignmentLeft.IsMapped() && alignmentRight.IsReverseStrand() || alignmentRight.IsMapped() && alignmentLeft.IsReverseStrand()) && alignmentLeft.RefID == alignmentRight.RefID && alignmentLeft.MapQuality > minScore && alignmentRight.MapQuality > minScore && contigSetHead->contigSet[alignmentLeft.RefID].contigLength >= contigLengthThreshold && ((alignmentLeft.Position < TailLength && alignmentRight.Position < TailLength) || (alignmentLeft.Position > ContigLength - TailLength && alignmentRight.Position > ContigLength - TailLength))) {
				readname = alignmentLeft.Name;
				char* name = (char*)malloc(sizeof(char)*100);
				strcpy(name, readname.c_str());
				
				readname1 = strtok(name, "#");
				barcodeid = strtok(NULL, "#");
				if (barcodeid.compare("0_0_0")!=0) {
					fprintf(fp, "%s,%d,%d,%d\n", barcodeid.c_str(), alignmentLeft.RefID, alignmentLeft.Position, alignmentRight.Position);
				}
				free(name);
			}
		}
		fclose(fp);
	}
	
	if ((fp = fopen(alignmentResultFile, "r")) == NULL) {
        printf("%s, does not exist!", alignmentResultFile);
        exit(0);
    }
	char* p;
	const char * split = ","; 
	long int maxSize = 10000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	long int i = 0;

	while ((fgets(line, maxSize, fp)) != NULL) {
		aligningCount++;
	}
	
	fclose(fp);

	cout << "aligningCount:" << aligningCount << endl;

	AligningResultHead * aligningResultHead = (AligningResultHead *)malloc(sizeof(AligningResultHead));
	aligningResultHead->aligningCount = aligningCount;
	aligningResultHead->aligningResult = (AligningResult *)malloc(sizeof(AligningResult)*aligningCount);

	if ((fp = fopen(alignmentResultFile, "r")) == NULL) {
		printf("%s, does not exist!", alignmentResultFile);
		exit(0);
	}
	while((fgets(line, maxSize, fp)) != NULL){ 
		
		p = strtok(line,split);
		aligningResultHead->aligningResult[i].barcodeid = p;
		
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].leftContigIndex = atoi(p);
				
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].leftPosition = atoi(p);
		
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].rightPosition = atoi(p);		

		i++;
	}
	fclose(fp);
	//qsort(aligningResultHead->aligningResult, aligningResultHead->aligningCount, sizeof(aligningResultHead->aligningResult[0]), comp);
	return aligningResultHead;
}

int comp(const void* x, const void* y) {

	return strcmp((*(AligningResult*)x).barcodeid.c_str(), (*(AligningResult*)y).barcodeid.c_str());
}




#endif