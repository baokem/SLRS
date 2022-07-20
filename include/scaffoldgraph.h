#ifndef SCAFFOLDGRAPH_H_INCLUDED 
#define SCAFFOLDGRAPH_H_INCLUDED 
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>

#include "contig.h"
#include "aligningFromBam.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using namespace BamTools;

using namespace std;

typedef struct PrimaryGraph{
	int leftnode;
	int rightnode;
	int leftreadcount;
	int leftcount;
	int rightreadcount;
	int rightcount;
	int orientation;
	double weight;
}PrimaryGraph;

typedef struct PrimaryGraphHead {
	PrimaryGraph* pgraph;
	int barcodecount;
}PrimaryGraphHead;

typedef struct ScaffoldGraphEdge {
	int nodeIndex;
	long int sharedcount;
	long int allcount;
	long int* orientation;
	double gapDistance;
	double weight;
	int orientationType;
}ScaffoldGraphEdge;


typedef struct ScaffoldGraph {
	ScaffoldGraphEdge* edge;
	long int edgeCount;
	long int barcodeCount;
}ScaffoldGraph;


typedef struct ScaffoldGraphHead {
	ScaffoldGraph* scaffoldGraph;
	long int scaffoldGraphNodeCount;
}ScaffoldGraphHead;

PrimaryGraphHead* GetPrimaryGraphHead(ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead,int readsnumperbarcode, int endreadsnum, int endlength,int num);

ScaffoldGraphHead* GetScaffoldGraphHeadFromPrimaryGraph(ContigSetHead* contigSetHead, PrimaryGraphHead* primaryGraphHead, int contigLengthThreshold);

bool DetermineAddEdgeInScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, PrimaryGraphHead* primaryGraphHead, int startIndex, int endIndex, int contigLengthThreshold);

void SortEdgeWeight(ScaffoldGraphHead* scaffoldGraphHead);

void KeepEdgeWithMaxWeight(ScaffoldGraphHead* scaffoldGraphHead, long int count);

void OutputScaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, long int* contigPosition);

bool RecoverTwoEdgeInScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead);

bool RemoveEdgeInScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead, long int leftIndex, long int rightIndex);

bool InsertEdgeInScaffold(ScaffoldGraphHead* scaffoldGraphHead, long int leftIndex, long int edgeIndex, long int rightIndex);

void Outputprimarygraph(PrimaryGraphHead* primaryGraphHead, char* primarygraphFile);

void OutputscaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, char* scaffoldGraphFile);

int com(const void* a, const void* b);

#endif
