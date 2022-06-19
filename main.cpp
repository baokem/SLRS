#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <malloc.h>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream> 
#include <string.h>
#include <cstdlib>
#include <set>
#include <bits/stdc++.h>

#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "scaffolding.h"

using namespace std;


int main(int argc,char** argv) {
	int maxSize = 2000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	char * resultOutPutDirectory = (char *)malloc(sizeof(char)*50);
	int contigLengthThreshold = 2000;
	int endlength = 30000;
	int readsnumperbarcode = 7;
	int endreadsnum = 3;
	char * alignmentResultFile = NULL;
	char * primarygraphFile = NULL;
	char * scaffoldGraphFile = NULL;
	char aligning = false;
	char * readFile = NULL;
	char * contigSetFile = NULL;
	char * leftReadBamFile = NULL;
	char * rightReadBamFile = NULL;
	int ch = 0;
	while ((ch = getopt(argc, argv, "c:r:l:o:t:e:n:f:a:")) != -1) {
		switch (ch) {
			case 'c': contigSetFile = (char *)(optarg); break;
			case 'r': rightReadBamFile = (char *)(optarg); break;
			case 'l': leftReadBamFile = (char *)optarg; break;
			case 'o': resultOutPutDirectory = (char *)optarg; break;
			case 't': contigLengthThreshold = atoi(optarg); break;
			case 'e': endlength = atoi(optarg); break;
			case 'n': readsnumperbarcode = atoi(optarg); break;
			case 'f': endreadsnum = atoi(optarg); break;
			case 'a': alignmentResultFile = (char*)(optarg); break;
			case 'b': readFile = (char*)(optarg); break;

			default: break; 
		}
	}
	
	if(opendir(resultOutPutDirectory) == NULL){  
		mkdir(resultOutPutDirectory, 0777);       
    }
	
	int fileNameLen = strlen(resultOutPutDirectory);

	ContigSetHead * contigSetHead = GetContigSet(contigSetFile, contigLengthThreshold);

	cout<<"contig count:"<<contigSetHead->contigCount<<endl;

	Trimend(contigSetHead);

	if(alignmentResultFile == NULL){
		alignmentResultFile = (char *)malloc(sizeof(char)* fileNameLen + 50);
		strcpy(alignmentResultFile, resultOutPutDirectory);
		strcat(alignmentResultFile, "/alignmentResultFile.fa");

	}else{
		aligning = true;
	}

	if (primarygraphFile == NULL) {
		primarygraphFile = (char*)malloc(sizeof(char) * fileNameLen + 50);
		strcpy(primarygraphFile, resultOutPutDirectory);
		strcat(primarygraphFile, "/primarygraphFile.fa");
	}

	if (scaffoldGraphFile == NULL) {
		scaffoldGraphFile = (char*)malloc(sizeof(char) * fileNameLen + 50);
		strcpy(scaffoldGraphFile, resultOutPutDirectory);
		strcat(scaffoldGraphFile, "/scaffoldGraphFile.fa");
	}

	if (aligning) {
		cout << "aligning is true" << endl;
	}
	cout << "starting GetAligningResultHead" << endl;                                                    
	AligningResultHead * aligningResultHead = GetAligningResultHead(readFile,leftReadBamFile, rightReadBamFile, alignmentResultFile, contigSetHead, aligning, contigLengthThreshold,endlength);
	cout << "starting GetPrimaryGraphHead" << endl;
	PrimaryGraphHead* primaryGraphHead = GetPrimaryGraphHead(contigSetHead, aligningResultHead, readsnumperbarcode, endreadsnum, endlength);
	
	Outputprimarygraph(primaryGraphHead, primarygraphFile);
	cout << "starting GetScaffoldGraphHeadFromPrimaryGraph" << endl;
	ScaffoldGraphHead * scaffoldGraphHead = GetScaffoldGraphHeadFromPrimaryGraph(contigSetHead, primaryGraphHead, contigLengthThreshold);

	OutputscaffoldGraphHead(scaffoldGraphHead, contigSetHead, scaffoldGraphFile);
	cout << "starting OptimizeScaffoldGraph" << endl;
	OptimizeScaffoldGraph(scaffoldGraphHead, contigSetHead);
	cout << "starting GetScaffoldingResult" << endl;
	ScaffoldSetHead * scaffoldSetHead = GetScaffoldingResult(scaffoldGraphHead, contigSetHead);

	char * scaffoldTagFileName = (char *)malloc(sizeof(char)*fileNameLen + 50);

	strcpy(scaffoldTagFileName, resultOutPutDirectory);

	strcat(scaffoldTagFileName, "/scaffold");

	OutPutScaffoldSet(scaffoldSetHead->scaffoldSet, contigSetHead, scaffoldTagFileName);
		
}
