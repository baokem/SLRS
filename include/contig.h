#ifndef CONTIG_H_INCLUDED 
#define CONTIG_H_INCLUDED 
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

using namespace std;


typedef struct Contig{
	char * contig;
	int contigLength;
}Contig;

typedef struct ContigSetHead{
	Contig * contigSet;
	long int contigCount;
	
}ContigSetHead;

ContigSetHead * GetContigSet(char * contigSetFile, int contigLengthThreshold);

char* ReverseComplement(char* temp);

void Trimend(ContigSetHead* contigSetHead);
#endif