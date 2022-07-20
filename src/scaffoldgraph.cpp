#ifndef BUILDSCAFFOLDGRAPH_CPP_INCLUDED 
#define BUILDSCAFFOLDGRAPH_CPP_INCLUDED 
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <exception>
#include <math.h>
#include <set>
#include <bits/stdc++.h>
#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"

using namespace std;

PrimaryGraphHead* GetPrimaryGraphHead(ContigSetHead* contigSetHead, AligningResultHead* aligningResultHead, int readsnumperbarcode,int endreadsnum, int endlength,int num){
	string previousbarcode = "a";
	int max = 0;
	int second = 0;
	int third = 0;
	int maxcontigindex = -1;
	int secondcontigindex = -1;
	int thirdcontigindex = -1;
	long int barcodecount = 0;
	long int index = 0;
	long int start = 0;
	long int endindex = 0;
	long int temp = 0;
	int lefthead = 0;
	int lefttail = 0;
	int righthead = 0;
	int righttail = 0;
	long int ContigLength = 0;
	long int TailLength = 0;
	int secondtail = 0;
	int secondtailmaxposition = -1;
	int secondtailminposition = 100000000;
	int secondhead = 0;
	int secondheadmaxposition = -1;
	int secondheadminposition = 100000000;
	int thirdtail = 0;
	int thirdtailmaxposition = -1;
	int thirdtailminposition = 100000000;
	int thirdhead = 0;
	int thirdheadmaxposition = -1;
	int thirdheadminposition = 100000000;
	long long secondheadtotalposition = 0;
	long long secondtailtotalposition = 0;
	long long thirdheadtotalposition = 0;
	long long thirdtailtotalposition = 0;
	double secondvariation = 0;
	double thirdvariation = 0;
	double secondaverage = 0;
	double thirdaverage = 0;
	long long secondhe = 0;
	long long thirdhe = 0;
	PrimaryGraphHead* primaryGraphHead = (PrimaryGraphHead*)malloc(sizeof(PrimaryGraphHead));

	for (long int i = 0; i < aligningResultHead->aligningCount; i++) {
		if (aligningResultHead->aligningResult[i].barcodeid.compare(previousbarcode) != 0) {
			barcodecount++;
			previousbarcode = aligningResultHead->aligningResult[i].barcodeid;
		}
	}
	cout <<"barcodecount:"<< barcodecount << endl;
	primaryGraphHead->barcodecount = barcodecount;
	primaryGraphHead->pgraph = (PrimaryGraph*)malloc(sizeof(PrimaryGraph) * (primaryGraphHead->barcodecount));
	previousbarcode = "a";
	long int* Count = (long int*)malloc(sizeof(long int) * contigSetHead->contigCount);

	for (long int i = 0; i < contigSetHead->contigCount; i++) {
		Count[i] = 0;
	}
	for (long int i = 0; i < aligningResultHead->aligningCount; i++) {
		if (aligningResultHead->aligningResult[i].barcodeid.compare(previousbarcode) == 0 ) {
			endindex++;
		}
		else {
			if (endindex - start > readsnumperbarcode)
			{
				for (long int j = start; j < endindex; j++) {
					Count[aligningResultHead->aligningResult[j].leftContigIndex]++;
				}
				for (long int i = 0; i < contigSetHead->contigCount; i++) {
					if (Count[i] > max) {
						max = Count[i];
						maxcontigindex = i;
					}
				}
				primaryGraphHead->pgraph[index].leftnode = maxcontigindex;
				primaryGraphHead->pgraph[index].leftreadcount = max;
				for (long int i = 0; i < contigSetHead->contigCount; i++) {
					if (Count[i] > second && i != primaryGraphHead->pgraph[index].leftnode) {
						second = Count[i];
						secondcontigindex = i;
					}
				}

				if (second < endreadsnum || third > 0) {
					max = 0;
					second = 0;
					third = 0;
					maxcontigindex = -1;
					secondcontigindex = -1;
					thirdcontigindex = -1;
					lefthead = 0;
					lefttail = 0;
					righthead = 0;
					righttail = 0;
					secondtail = 0;
					secondtailmaxposition = -1;
					secondtailminposition = 100000000;
					secondhead = 0;
					secondheadmaxposition = -1;
					secondheadminposition = 100000000;
					thirdtail = 0;
					thirdtailmaxposition = -1;
					thirdtailminposition = 100000000;
					thirdhead = 0;
					thirdheadmaxposition = -1;
					thirdheadminposition = 100000000;
					secondheadtotalposition = 0;
					secondtailtotalposition = 0;
					thirdheadtotalposition = 0;
					thirdtailtotalposition = 0;
					secondvariation = 0;
					thirdvariation = 0;
					secondaverage = 0;
					thirdaverage = 0;
					secondhe = 0;
					thirdhe = 0;
					start = i;
					endindex = i + 1;
					previousbarcode = aligningResultHead->aligningResult[i].barcodeid;
					for (long int i = 0; i < contigSetHead->contigCount; i++) {
						Count[i] = 0;
					}
					continue;
				}
				
				for (long int i = 0; i < contigSetHead->contigCount; i++) {
					if (Count[i] > third && i != primaryGraphHead->pgraph[index].leftnode && i != secondcontigindex) {
						third = Count[i];
						thirdcontigindex = i;
					}
				}

				if (second - third < 2 && third > num) {
					for (long int j = start; j < endindex; j++) {
						ContigLength = contigSetHead->contigSet[aligningResultHead->aligningResult[j].leftContigIndex].contigLength;
						TailLength = MIN(endlength, ContigLength / 2);
						if (aligningResultHead->aligningResult[j].leftContigIndex == secondcontigindex) {
							if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
								secondheadtotalposition = aligningResultHead->aligningResult[j].leftPosition + secondheadtotalposition;
								secondhead++;
							}
							if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
								secondtailtotalposition = aligningResultHead->aligningResult[j].leftPosition + secondtailtotalposition;
								secondtail++;
							}
						}
						if (aligningResultHead->aligningResult[j].leftContigIndex == thirdcontigindex) {
							if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
								thirdheadtotalposition = aligningResultHead->aligningResult[j].leftPosition + thirdheadtotalposition;
								thirdhead++;
							}
							if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
								thirdtailtotalposition = aligningResultHead->aligningResult[j].leftPosition + thirdtailtotalposition;
								thirdtail++;
							}
						}
					}
					if (secondhead >= secondtail && secondhead != 0) {
						secondaverage = (double)((secondheadtotalposition - (secondheadminposition * secondhead)) / secondhead);
						for (long int j = start; j < endindex; j++) {
							ContigLength = contigSetHead->contigSet[aligningResultHead->aligningResult[j].leftContigIndex].contigLength;
							TailLength = MIN(endlength, ContigLength / 2);
							if (aligningResultHead->aligningResult[j].leftContigIndex == secondcontigindex) {
								if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
									secondhe = (long long)pow((double)aligningResultHead->aligningResult[j].leftPosition - secondheadminposition - secondaverage, (double)2) + secondhe;
								}
							}
						}
						secondvariation = sqrt(secondhe / secondhead) / secondaverage;

					}
					else {
						secondaverage = (double)((secondtailtotalposition - secondtailminposition * secondtail) / secondtail);
						for (long int j = start; j < endindex; j++) {
							ContigLength = contigSetHead->contigSet[aligningResultHead->aligningResult[j].leftContigIndex].contigLength;
							TailLength = MIN(endlength, ContigLength / 2);
							if (aligningResultHead->aligningResult[j].leftContigIndex == secondcontigindex) {
								if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
									secondhe = (long long)pow((double)aligningResultHead->aligningResult[j].leftPosition - secondtailminposition - secondaverage, (double)2) + secondhe;
								}
							}
						}
						secondvariation = sqrt(secondhe / secondtail) / secondaverage;
					}

					if (thirdhead >= thirdtail && thirdhead != 0) {
						thirdaverage = (double)((thirdheadtotalposition - thirdheadminposition * thirdhead) / thirdhead);
						for (long int j = start; j < endindex; j++) {
							ContigLength = contigSetHead->contigSet[aligningResultHead->aligningResult[j].leftContigIndex].contigLength;
							TailLength = MIN(endlength, ContigLength / 2);
							if (aligningResultHead->aligningResult[j].leftContigIndex == thirdcontigindex) {
								if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
									thirdhe = (long long)pow((double)aligningResultHead->aligningResult[j].leftPosition - thirdheadminposition - thirdaverage, (double)2) + thirdhe;
								}
							}
						}
						thirdvariation = sqrt(thirdhe / thirdhead) / thirdaverage;
					}
					else {
						thirdaverage = (double)((thirdtailtotalposition - thirdtailminposition * thirdtail) / thirdtail);
						for (long int j = start; j < endindex; j++) {
							ContigLength = contigSetHead->contigSet[aligningResultHead->aligningResult[j].leftContigIndex].contigLength;
							TailLength = MIN(endlength, ContigLength / 2);
							if (aligningResultHead->aligningResult[j].leftContigIndex == thirdcontigindex) {
								if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
									thirdhe = (long long)pow((double)aligningResultHead->aligningResult[j].leftPosition - thirdtailminposition - thirdaverage, (double)2) + thirdhe;
								}
							}
						}
						thirdvariation = sqrt(thirdhe / thirdtail) / thirdaverage;
					}
					if (((thirdhead >= secondhead && thirdhead >= secondtail) || (thirdtail >= secondhead && thirdtail >= secondtail)) && thirdvariation >= secondvariation) {
						primaryGraphHead->pgraph[index].rightnode = thirdcontigindex;
						primaryGraphHead->pgraph[index].rightreadcount = third;
					}
					else {
						primaryGraphHead->pgraph[index].rightnode = secondcontigindex;
						primaryGraphHead->pgraph[index].rightreadcount = second;
					}
					if (primaryGraphHead->pgraph[index].leftnode > primaryGraphHead->pgraph[index].rightnode) {
						temp = primaryGraphHead->pgraph[index].leftnode;
						primaryGraphHead->pgraph[index].leftnode = primaryGraphHead->pgraph[index].rightnode;
						primaryGraphHead->pgraph[index].rightnode = temp;

						temp = primaryGraphHead->pgraph[index].leftreadcount;
						primaryGraphHead->pgraph[index].leftreadcount = primaryGraphHead->pgraph[index].rightreadcount;
						primaryGraphHead->pgraph[index].rightreadcount = temp;
					}
					for (long int j = start; j < endindex; j++) {
						ContigLength = contigSetHead->contigSet[aligningResultHead->aligningResult[j].leftContigIndex].contigLength;
						TailLength = MIN(endlength, ContigLength / 2);
						if (aligningResultHead->aligningResult[j].leftContigIndex == primaryGraphHead->pgraph[index].leftnode) {
							if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
								lefthead++;
							}
							if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
								lefttail++;
							}
						}
						if (aligningResultHead->aligningResult[j].leftContigIndex == primaryGraphHead->pgraph[index].rightnode) {
							if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
								righthead++;
							}
							if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
								righttail++;
							}
						}
					}
					if (lefthead >= lefttail && righthead >= righttail) {//首首
						primaryGraphHead->pgraph[index].leftcount = lefthead;
						primaryGraphHead->pgraph[index].rightcount = righthead;
						primaryGraphHead->pgraph[index].orientation = 0;
					}
					if (lefthead <= lefttail && righthead <= righttail) {//尾尾
						primaryGraphHead->pgraph[index].leftcount = lefttail;
						primaryGraphHead->pgraph[index].rightcount = righttail;
						primaryGraphHead->pgraph[index].orientation = 1;
					}
					if (lefthead >= lefttail && righthead <= righttail) {//首尾
						primaryGraphHead->pgraph[index].leftcount = lefthead;
						primaryGraphHead->pgraph[index].rightcount = righttail;
						primaryGraphHead->pgraph[index].orientation = 2;
					}
					if (lefthead <= lefttail && righthead >= righttail) {//尾首
						primaryGraphHead->pgraph[index].leftcount = lefttail;
						primaryGraphHead->pgraph[index].rightcount = righthead;
						primaryGraphHead->pgraph[index].orientation = 3;
					}

					primaryGraphHead->pgraph[index].weight = (double)1000 * (primaryGraphHead->pgraph[index].leftcount + primaryGraphHead->pgraph[index].rightcount) / (MIN(endlength, contigSetHead->contigSet[primaryGraphHead->pgraph[index].leftnode].contigLength / 2) + MIN(endlength, contigSetHead->contigSet[primaryGraphHead->pgraph[index].rightnode].contigLength / 2));
					max = 0;
					second = 0;
					third = 0;
					maxcontigindex = -1;
					secondcontigindex = -1;
					thirdcontigindex = -1;
					lefthead = 0;
					lefttail = 0;
					righthead = 0;
					righttail = 0;
					secondtail = 0;
					secondtailmaxposition = -1;
					secondtailminposition = 100000000;
					secondhead = 0;
					secondheadmaxposition = -1;
					secondheadminposition = 100000000;
					thirdtail = 0;
					thirdtailmaxposition = -1;
					thirdtailminposition = 100000000;
					thirdhead = 0;
					thirdheadmaxposition = -1;
					thirdheadminposition = 100000000;
					secondheadtotalposition = 0;
					secondtailtotalposition = 0;
					thirdheadtotalposition = 0;
					thirdtailtotalposition = 0;
					secondvariation = 0;
					thirdvariation = 0;
					secondaverage = 0;
					thirdaverage = 0;
					secondhe = 0;
					thirdhe = 0;
					start = i;
					endindex = i + 1;
					previousbarcode = aligningResultHead->aligningResult[i].barcodeid;
					for (long int i = 0; i < contigSetHead->contigCount; i++) {
						Count[i] = 0;
					}
					index++;
					continue;
				}

				primaryGraphHead->pgraph[index].rightnode = secondcontigindex;
				primaryGraphHead->pgraph[index].rightreadcount = second;
				if (primaryGraphHead->pgraph[index].leftnode > primaryGraphHead->pgraph[index].rightnode) {
					temp = primaryGraphHead->pgraph[index].leftnode;
					primaryGraphHead->pgraph[index].leftnode = primaryGraphHead->pgraph[index].rightnode;
					primaryGraphHead->pgraph[index].rightnode = temp;

					temp = primaryGraphHead->pgraph[index].leftreadcount;
					primaryGraphHead->pgraph[index].leftreadcount = primaryGraphHead->pgraph[index].rightreadcount;
					primaryGraphHead->pgraph[index].rightreadcount = temp;
				}
				for (long int j = start; j < endindex; j++) {
					ContigLength = contigSetHead->contigSet[aligningResultHead->aligningResult[j].leftContigIndex].contigLength;
					TailLength = MIN(endlength, ContigLength / 2);
					if (aligningResultHead->aligningResult[j].leftContigIndex == primaryGraphHead->pgraph[index].leftnode) {
						if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
							lefthead++;
						}
						if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
							lefttail++;
						}
					}
					if (aligningResultHead->aligningResult[j].leftContigIndex == primaryGraphHead->pgraph[index].rightnode) {
						if (aligningResultHead->aligningResult[j].leftPosition < TailLength) {
							righthead++;
						}
						if (aligningResultHead->aligningResult[j].leftPosition > ContigLength - TailLength) {
							righttail++;
						}
					}
				}
				if (lefthead >= lefttail && righthead >= righttail) {//首首
					primaryGraphHead->pgraph[index].leftcount = lefthead;
					primaryGraphHead->pgraph[index].rightcount = righthead;
					primaryGraphHead->pgraph[index].orientation = 0;
				}
				if (lefthead <= lefttail && righthead <= righttail) {//尾尾
					primaryGraphHead->pgraph[index].leftcount = lefttail;
					primaryGraphHead->pgraph[index].rightcount = righttail;
					primaryGraphHead->pgraph[index].orientation = 1;
				}
				if (lefthead >= lefttail && righthead <= righttail) {//首尾
					primaryGraphHead->pgraph[index].leftcount = lefthead;
					primaryGraphHead->pgraph[index].rightcount = righttail;
					primaryGraphHead->pgraph[index].orientation = 2;
				}
				if (lefthead <= lefttail && righthead >= righttail) {//尾首
					primaryGraphHead->pgraph[index].leftcount = lefttail;
					primaryGraphHead->pgraph[index].rightcount = righthead;
					primaryGraphHead->pgraph[index].orientation = 3;
				}
				
				
				primaryGraphHead->pgraph[index].weight = (double)1000 * (primaryGraphHead->pgraph[index].leftcount + primaryGraphHead->pgraph[index].rightcount) / (MIN(endlength, contigSetHead->contigSet[primaryGraphHead->pgraph[index].leftnode].contigLength / 2) + MIN(endlength, contigSetHead->contigSet[primaryGraphHead->pgraph[index].rightnode].contigLength / 2));
				//primaryGraphHead->pgraph[index].weight = primaryGraphHead->pgraph[index].leftcount + primaryGraphHead->pgraph[index].rightcount;
				//primaryGraphHead->pgraph[index].weight = 1 + (primaryGraphHead->pgraph[index].leftcount + primaryGraphHead->pgraph[index].rightcount - 10) * 0.1;
				max = 0;
				second = 0;
				third = 0;
				maxcontigindex = -1;
				secondcontigindex = -1;
				thirdcontigindex = -1;
				lefthead = 0;
				lefttail = 0;
				righthead = 0;
				righttail = 0;
				secondtail = 0;
				secondtailmaxposition = -1;
				secondtailminposition = 100000000;
				secondhead = 0;
				secondheadmaxposition = -1;
				secondheadminposition = 100000000;
				thirdtail = 0;
				thirdtailmaxposition = -1;
				thirdtailminposition = 100000000;
				thirdhead = 0;
				thirdheadmaxposition = -1;
				thirdheadminposition = 100000000;
				secondheadtotalposition = 0;
				secondtailtotalposition = 0;
				thirdheadtotalposition = 0;
				thirdtailtotalposition = 0;
				secondvariation = 0;
				thirdvariation = 0;
				secondaverage = 0;
				thirdaverage = 0;
				secondhe = 0;
				thirdhe = 0;
				for (long int i = 0; i < contigSetHead->contigCount; i++) {
					Count[i] = 0;
				}
				index++;	
			}
			start = i;
			endindex = i + 1;
			previousbarcode = aligningResultHead->aligningResult[i].barcodeid;
		}	
	}
	primaryGraphHead->barcodecount = index;
	qsort(primaryGraphHead->pgraph, primaryGraphHead->barcodecount, sizeof(primaryGraphHead->pgraph[0]), com);
	return primaryGraphHead;
}

void Outputprimarygraph(PrimaryGraphHead* primaryGraphHead, char* primarygraphFile) {
	FILE* fp;
	if ((fp = fopen(primarygraphFile, "w")) == NULL) {
		printf("%s, does not exist!", primarygraphFile);
		exit(0);
	}
	for (long int i = 0; i < primaryGraphHead->barcodecount; i++) {
		fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%lf\n", primaryGraphHead->pgraph[i].leftnode, primaryGraphHead->pgraph[i].rightnode, primaryGraphHead->pgraph[i].leftreadcount, primaryGraphHead->pgraph[i].leftcount, primaryGraphHead->pgraph[i].rightreadcount, primaryGraphHead->pgraph[i].rightcount, primaryGraphHead->pgraph[i].orientation, primaryGraphHead->pgraph[i].weight);
	}
	fclose(fp);
}

int com(const void* a, const void* b) {
	if ((*(PrimaryGraph*)a).leftnode != (*(PrimaryGraph*)b).leftnode)
		return (*(PrimaryGraph*)a).leftnode - (*(PrimaryGraph*)b).leftnode;
	else return (*(PrimaryGraph*)a).rightnode - (*(PrimaryGraph*)b).rightnode;
}

ScaffoldGraphHead* GetScaffoldGraphHeadFromPrimaryGraph(ContigSetHead* contigSetHead, PrimaryGraphHead* primaryGraphHead, int contigLengthThreshold) {


	long int leftContigIndex = -1;
	long int rightContigIndex = -1;
	long int barcodecount = primaryGraphHead->barcodecount;
	cout << "barcodenumber:"<< barcodecount << endl;
	long int* edgeCount = (long int*)malloc(sizeof(long int) * contigSetHead->contigCount);
	long int* Count = (long int*)malloc(sizeof(long int) * contigSetHead->contigCount);
	for (long int i = 0; i < contigSetHead->contigCount; i++) {
		edgeCount[i] = 0;
		Count[i] = 0;
	}
	
	for (long int i = 0; i < barcodecount; i++) {

		Count[primaryGraphHead->pgraph[i].leftnode]++;
		Count[primaryGraphHead->pgraph[i].rightnode]++;
		if (!(primaryGraphHead->pgraph[i].leftnode == leftContigIndex && primaryGraphHead->pgraph[i].rightnode == rightContigIndex)) {
			leftContigIndex = primaryGraphHead->pgraph[i].leftnode;
			rightContigIndex = primaryGraphHead->pgraph[i].rightnode;
			
			edgeCount[leftContigIndex]++;
			edgeCount[rightContigIndex]++;
		}
	}
	
	ScaffoldGraphHead* scaffoldGraphHead = (ScaffoldGraphHead*)malloc(sizeof(ScaffoldGraphHead));
	scaffoldGraphHead->scaffoldGraphNodeCount = contigSetHead->contigCount;
	scaffoldGraphHead->scaffoldGraph = (ScaffoldGraph*)malloc(sizeof(ScaffoldGraph) * scaffoldGraphHead->scaffoldGraphNodeCount);

	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		scaffoldGraphHead->scaffoldGraph[i].edgeCount = edgeCount[i];
		scaffoldGraphHead->scaffoldGraph[i].barcodeCount = Count[i];
		scaffoldGraphHead->scaffoldGraph[i].edge = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge) * scaffoldGraphHead->scaffoldGraph[i].edgeCount);
		for (long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++) {
			scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex = -1;
			scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance = 0;
			scaffoldGraphHead->scaffoldGraph[i].edge[j].orientation = (long int*)malloc(sizeof(long int) * 4);
			for (long int p = 0; p < 4; p++) {
				scaffoldGraphHead->scaffoldGraph[i].edge[j].orientation[p] = 0;
			}
			scaffoldGraphHead->scaffoldGraph[i].edge[j].weight = 0;
			scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType = -1;
		}
		scaffoldGraphHead->scaffoldGraph[i].edgeCount = 0;
	}
	
	leftContigIndex = primaryGraphHead->pgraph[0].leftnode;
	rightContigIndex = primaryGraphHead->pgraph[0].rightnode;
	int startIndex = 0;
	int endIndex = 0;
	for (long int i = 0; i < barcodecount; i++) {
		if (primaryGraphHead->pgraph[i].leftnode == leftContigIndex && primaryGraphHead->pgraph[i].rightnode == rightContigIndex) {
			endIndex++;
			if (i == barcodecount - 1) {
				DetermineAddEdgeInScaffoldGraph(scaffoldGraphHead, contigSetHead, primaryGraphHead, startIndex, endIndex, contigLengthThreshold);
				
			}
		}
		else {
			DetermineAddEdgeInScaffoldGraph(scaffoldGraphHead, contigSetHead, primaryGraphHead, startIndex, endIndex, contigLengthThreshold);
			startIndex = i;
			endIndex = i+1;
			leftContigIndex = primaryGraphHead->pgraph[i].leftnode;
			rightContigIndex = primaryGraphHead->pgraph[i].rightnode;
			
		}
	}
	
	KeepEdgeWithMaxWeight(scaffoldGraphHead, 5);

	return scaffoldGraphHead;

}

bool DetermineAddEdgeInScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, PrimaryGraphHead* primaryGraphHead, int startIndex, int endIndex, int contigLengthThreshold) {

	long int leftContigIndex = primaryGraphHead->pgraph[startIndex].leftnode;
	long int rightContigIndex = primaryGraphHead->pgraph[startIndex].rightnode;
	long int sharedcount = endIndex - startIndex;
	long int leftContigLength = contigSetHead->contigSet[leftContigIndex].contigLength;
	long int rightContigLength = contigSetHead->contigSet[rightContigIndex].contigLength;

	if (leftContigLength < contigLengthThreshold || rightContigLength < contigLengthThreshold) {
		return false;
	}
	for (long int i = startIndex; i < endIndex; i++) {

		if (primaryGraphHead->pgraph[i].orientation == 0) {
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientation[0]++;
		}
		if (primaryGraphHead->pgraph[i].orientation == 1) {
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientation[1]++;
		}
		if (primaryGraphHead->pgraph[i].orientation == 2) {
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientation[2]++;
		}
		if (primaryGraphHead->pgraph[i].orientation == 3) {
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientation[3]++;
		}
	}

	long int maxOrientationCount = -1;
	long int maxOrientationIndex = -1;

	for (long int i = 0; i <= 3; i++) {
		if (scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientation[i] > maxOrientationCount) {
			maxOrientationCount = scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientation[i];
			maxOrientationIndex = i;
		}
	}

	for (long int i = startIndex; i < endIndex; i++) {
		if (primaryGraphHead->pgraph[i].orientation == maxOrientationIndex) {
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].weight = primaryGraphHead->pgraph[i].weight + scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].weight;
		}
	}
	int allcount = 0;
	for (long int i = 0; i < primaryGraphHead->barcodecount; i++) {
		if (primaryGraphHead->pgraph[i].leftnode == leftContigIndex) {
			if (maxOrientationIndex == 0 || maxOrientationIndex == 2) {
				if (primaryGraphHead->pgraph[i].orientation == 0 || primaryGraphHead->pgraph[i].orientation == 2) {
					allcount++;
				}
			}
			if (maxOrientationIndex == 1 || maxOrientationIndex == 3) {
				if (primaryGraphHead->pgraph[i].orientation == 1 || primaryGraphHead->pgraph[i].orientation == 3) {
					allcount++;
				}
			}
		}
		if (primaryGraphHead->pgraph[i].rightnode == leftContigIndex) {
			if (maxOrientationIndex == 0 || maxOrientationIndex == 2) {
				if (primaryGraphHead->pgraph[i].orientation == 0 || primaryGraphHead->pgraph[i].orientation == 3) {
					allcount++;
				}
			}
			if (maxOrientationIndex == 1 || maxOrientationIndex == 3) {
				if (primaryGraphHead->pgraph[i].orientation == 1 || primaryGraphHead->pgraph[i].orientation == 2) {
					allcount++;
				}
			}
		}
		if (primaryGraphHead->pgraph[i].leftnode == rightContigIndex) {
			if (maxOrientationIndex == 0 || maxOrientationIndex == 3) {
				if (primaryGraphHead->pgraph[i].orientation == 0 || primaryGraphHead->pgraph[i].orientation == 2) {
					allcount++;
				}
			}
			if (maxOrientationIndex == 1 || maxOrientationIndex == 2) {
				if (primaryGraphHead->pgraph[i].orientation == 1 || primaryGraphHead->pgraph[i].orientation == 3) {
					allcount++;
				}
			}
		}
		if (primaryGraphHead->pgraph[i].rightnode == rightContigIndex) {
			if (maxOrientationIndex == 0 || maxOrientationIndex == 3) {
				if (primaryGraphHead->pgraph[i].orientation == 0 || primaryGraphHead->pgraph[i].orientation == 3) {
					allcount++;
				}
			}
			if (maxOrientationIndex == 1 || maxOrientationIndex == 2) {
				if (primaryGraphHead->pgraph[i].orientation == 1 || primaryGraphHead->pgraph[i].orientation == 2) {
					allcount++;
				}
			}
		}
	}

	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientationType = maxOrientationIndex;

	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].nodeIndex = rightContigIndex;
	
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].allcount = allcount;

	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].sharedcount = sharedcount;
	
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].gapDistance = (double)scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].sharedcount / scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].allcount;
	
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount++;

	return true;

}

void SortEdgeWeight(ScaffoldGraphHead* scaffoldGraphHead) {
	long int temp = 0;
	double gap = 0;
	double weight = 0;
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		for (long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount - 1; j++) {
			for (long int p = j + 1; p < scaffoldGraphHead->scaffoldGraph[i].edgeCount; p++) {
				if (scaffoldGraphHead->scaffoldGraph[i].edge[j].weight < scaffoldGraphHead->scaffoldGraph[i].edge[p].weight) {
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].allcount;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].allcount = scaffoldGraphHead->scaffoldGraph[i].edge[p].allcount;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].allcount = temp;

					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].sharedcount;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].sharedcount = scaffoldGraphHead->scaffoldGraph[i].edge[p].sharedcount;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].sharedcount = temp;

					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex = scaffoldGraphHead->scaffoldGraph[i].edge[p].nodeIndex;;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].nodeIndex = temp;

					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType = scaffoldGraphHead->scaffoldGraph[i].edge[p].orientationType;;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].orientationType = temp;

					weight = scaffoldGraphHead->scaffoldGraph[i].edge[j].weight;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].weight = scaffoldGraphHead->scaffoldGraph[i].edge[p].weight;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].weight = weight;

					gap = scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance = scaffoldGraphHead->scaffoldGraph[i].edge[p].gapDistance;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].gapDistance = gap;
				}
			}
		}
	}
}

void KeepEdgeWithMaxWeight(ScaffoldGraphHead* scaffoldGraphHead, long int count) {

	SortEdgeWeight(scaffoldGraphHead);

	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		if (scaffoldGraphHead->scaffoldGraph[i].edgeCount > count) {
			scaffoldGraphHead->scaffoldGraph[i].edgeCount = count;
		}
	}
}

void OutputScaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, long int* contigPosition) {
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		for (long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++) {
			cout << i << "--" << scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex << "--" << scaffoldGraphHead->scaffoldGraph[i].edge[j].weight << "--" << scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance << "--type:" << scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType << endl;
			cout << "position:" << contigPosition[i] << "--" << contigPosition[scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex] << endl;
			cout << "length:" << contigSetHead->contigSet[i].contigLength << "--" << contigSetHead->contigSet[scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex].contigLength << endl;
			for (long int p = 0; p < 4; p++) {
				cout << scaffoldGraphHead->scaffoldGraph[i].edge[j].orientation[p] << ",";
			}
			cout << endl;
		}
		
	}

}

void OutputscaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, char* scaffoldGraphFile) {
	FILE* fp;
	if ((fp = fopen(scaffoldGraphFile, "w")) == NULL) {
		printf("%s, does not exist!", scaffoldGraphFile);
		exit(0);
	}
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		for (long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++) {
			fprintf(fp, "%ld--%d--%d--%lf--%ld--%lf\n", i, scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex, scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType, scaffoldGraphHead->scaffoldGraph[i].edge[j].weight, scaffoldGraphHead->scaffoldGraph[i].edge[j].sharedcount, scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance);
		}
	}
	fclose(fp);
}

bool RecoverTwoEdgeInScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead) {
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		for (long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++) {
			if (i < scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex) {
				InsertEdgeInScaffold(scaffoldGraphHead, i, j, scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex);
			}
		}
	}
	return true;
}

bool RemoveEdgeInScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead, long int leftIndex, long int rightIndex) {
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount; i++) {
		if (scaffoldGraphHead->scaffoldGraph[leftIndex].edge[i].nodeIndex == rightIndex) {
			for (long int j = i; j < scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount - 1; j++) {
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].nodeIndex = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].nodeIndex;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].weight = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].weight;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].orientationType = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].orientationType;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].orientation[0] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].orientation[0];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].orientation[1] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].orientation[1];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].orientation[2] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].orientation[2];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].orientation[3] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].orientation[3];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].allcount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].allcount;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].sharedcount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].sharedcount;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].gapDistance = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].gapDistance;
			}
			scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount--;
			return true;
		}
	}
	return false;
}

bool InsertEdgeInScaffold(ScaffoldGraphHead* scaffoldGraphHead, long int leftIndex, long int edgeIndex, long int rightIndex) {
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount; i++) {
		if (scaffoldGraphHead->scaffoldGraph[rightIndex].edge[i].nodeIndex == leftIndex) {
			return false;
		}
	}
	long int edgeCount = scaffoldGraphHead->scaffoldGraph[rightIndex].edgeCount;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].nodeIndex = leftIndex;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].weight = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].weight;
	if (scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType == 0 || scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType == 1) {
		scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientationType = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType;
	}
	else if (scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType == 2) {
		scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientationType = 3;
	}
	else {
		scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientationType = 2;
	}
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientation[0] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientation[0];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientation[1] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientation[1];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientation[2] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientation[3];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientation[3] = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientation[2];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].allcount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].allcount;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].sharedcount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].sharedcount;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].gapDistance = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].gapDistance;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edgeCount++;
	return true;

}
#endif
