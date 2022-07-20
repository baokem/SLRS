#ifndef SCAFFOLDING_CPP_INCLUDED 
#define SCAFFOLDING_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <typeinfo>
#include <vector>
#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "scaffolding.h"
#include "lp/lp_lib.h"


long int RemoveOrientationContradictionInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, long int contigCount, bool * contigOrientation){
	
	long int edgeNumber = 0;
    long int constraintNumber = 0;
	
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		edgeNumber = edgeNumber + scaffoldGraphHead->scaffoldGraph[i].edgeCount;
	}
	
	
	lprec * lp;
    int Ncol, * colno = NULL, ret = 0;
    REAL * row = NULL;
    
    Ncol = edgeNumber + contigCount;
    lp = make_lp(0, Ncol);
    if(lp == NULL){
        printf("couldn't construct a new linear programming model");
        exit(0);
    }
    int * weight = new int[edgeNumber]; 
    long int * edgeLeftNode = new long int[edgeNumber];
    long int * edgeRightNode = new long int[edgeNumber];
    
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL)){
        printf("couldn't new colno and row");
        exit(0);
    }
    
    set_add_rowmode(lp, TRUE);
    
    int ttt = 0;
	long int index = 0;
	long int p = 1;
	long int c = 10000;
    for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			if(scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType == 0 || scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType == 1){//首首或尾尾连接即左右端点一个正向一个反向
				index = 0;
				colno[index] = i + 1;
				row[index++] = 1;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = 1;
                        
				colno[index] = contigCount + p; 
				row[index++] = 1;
				if(!add_constraintex(lp, index, row, colno, LE, 2)){
					printf("couldn't add_constraintex");
					exit(0);
				}
				

				index = 0;
				colno[index] = i + 1;
				row[index++] = 1;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = 1;
                        
				colno[index] = contigCount + p; 
				row[index++] = -1;
				if(!add_constraintex(lp, index, row, colno, GE, 0)){
					printf("couldn't add_constraintex");
					exit(0);
				}
				
				ttt++;
				constraintNumber = constraintNumber + 2;
                        
            }else{
				index = 0;
				colno[index] = i + 1;
				row[index++] = 1;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = -1;
				
				colno[index] = contigCount + p;
				row[index++] = 1;
				if(!add_constraintex(lp, index, row, colno, LE, 1)){
					printf("couldn't add_constraintex");
					exit(0);
				}
				
				index = 0;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = 1;
				colno[index] = i + 1; 
				row[index++] = -1;
				
				colno[index] = contigCount + p;
				row[index++] = -1;
				if(!add_constraintex(lp, index, row, colno, GE, -1)){//>=-1
					printf("couldn't add_constraintex");
					exit(0);
				}
				
				constraintNumber = constraintNumber + 2;
				ttt++;
				
			}

			edgeLeftNode[p - 1] = i;
			edgeRightNode[p - 1] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
			p++;
			
			weight[p - 2] = scaffoldGraphHead->scaffoldGraph[i].edge[j].weight;
			
		}    
    }

    for(long int i = 0; i < contigCount + edgeNumber; i++){
        set_binary(lp, i + 1, TRUE);
    }
	
    p--;
    
    index = 0;
    for(long int i = 0; i < p; i++){
        colno[index] = contigCount + i + 1; 
        row[index] = weight[index];
        index++;
    }
    if(!set_obj_fnex(lp, index, row, colno)){
        printf("couldn't set_obj_fnex");
        exit(0);
    }
	set_add_rowmode(lp, FALSE);
	set_timeout(lp, 600);

    set_maxim(lp);
    set_scaling(lp, 128); 
    ret = solve(lp);

    if(!(ret == 0 || ret == 1)){
		cout<<"ee:"<<ret<<endl;
		return ret;
    }
    REAL * pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
    get_primal_solution(lp, pv);

    double temp = 1;
    long int result = 0;
    for(long int i = contigCount + constraintNumber + 1; i < constraintNumber + contigCount + p + 1; i++){
        
		if(pv[i] == 0){
           RemoveEdgeInScaffoldGraph(scaffoldGraphHead, edgeLeftNode[i-contigCount-constraintNumber-1], edgeRightNode[i-contigCount-constraintNumber-1]);
           result++;
        }
		
    }
	
	for(long int i = constraintNumber + 1; i < contigCount + constraintNumber + 1; i++){
        contigOrientation[i - constraintNumber - 1] = pv[i];
    }
    

    delete [] weight;
    delete [] edgeLeftNode;
    delete [] edgeRightNode;
    delete [] pv;
    
    delete_lp(lp);
    
    return result;
}

long int DetermineContigOrderInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, bool * contigOrientation, long int * contigPosition){
	
	long int contigCount = contigSetHead->contigCount;
	
	long int edgeNumber = 0;
    long int constraintNumber = 0;
	
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		edgeNumber = edgeNumber + scaffoldGraphHead->scaffoldGraph[i].edgeCount;
	}
	
	
	lprec * lp;
    int Ncol, * colno = NULL, ret = 0;
    REAL * row = NULL;
    
    Ncol = edgeNumber + contigCount;
    lp = make_lp(0, Ncol);
    if(lp == NULL){
        printf("couldn't construct a new linear programming model");
        exit(0);
    }
    int * weight = new int[edgeNumber]; 
    long int * edgeLeftNode = new long int[edgeNumber];
    long int * edgeRightNode = new long int[edgeNumber];
	long int * gapDistance = new long int[edgeNumber];
    
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL)){
        printf("couldn't new colno and row");
        exit(0);
    }
    
    set_add_rowmode(lp, TRUE);
    
    int ttt = 0;
	long int index = 0;
	long int p = 1;
	long int c = 10000;
    for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			if((scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType == 0 && contigOrientation[i] == 0) || (scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType == 1 && contigOrientation[i] == 1) || (scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType == 2 && contigOrientation[i] == 0) || (scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType == 3 && contigOrientation[i] == 1)){
				index = 0;
				colno[index] = i + 1;
				row[index++] = -1;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = 1;
				
				colno[index] = contigCount + p;
				row[index++] = c;
				
				if(!add_constraintex(lp, index, row, colno, GE, 1)){
					printf("couldn't add_constraintex");
					exit(0);
				}
				
				index = 0;
				colno[index] = i + 1;
				row[index++] = -1;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = 1;
				
				colno[index] = contigCount + p;
				row[index++] = -c;
				
				if(!add_constraintex(lp, index, row, colno, LE, 3)){
					printf("couldn't add_constraintex");
					exit(0);
				}
				
			}else{
				index = 0;
				colno[index] = i + 1;
				row[index++] = 1;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = -1;
				
				colno[index] = contigCount + p;
				row[index++] = c;
				
				if(!add_constraintex(lp, index, row, colno, GE, 1)){
					printf("couldn't add_constraintex");
					exit(0);
				}
				
				
				index = 0;
				colno[index] = i + 1;
				row[index++] = 1;
				colno[index] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex + 1; 
				row[index++] = -1;
				
				colno[index] = contigCount + p;
				row[index++] = -c;
				
				if(!add_constraintex(lp, index, row, colno, LE, 3)){
					printf("couldn't add_constraintex");
					exit(0);
				}
				
			}
			constraintNumber = constraintNumber + 2;
			
			edgeLeftNode[p - 1] = i;
			edgeRightNode[p - 1] = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
			p++;
			
			weight[p-2] = scaffoldGraphHead->scaffoldGraph[i].edge[j].weight;
			
		}
		index = 0;
		colno[index] = i + 1;
		row[index++] = 1;
		if(!add_constraintex(lp, index, row, colno, LE, contigCount + 1)){
			printf("couldn't add_constraintex");
			exit(0);
		}
		constraintNumber = constraintNumber + 1;
    }

    for(long int i = 0; i < contigCount; i++){
        set_int(lp, i + 1, TRUE);
		
		
    }
	for(long int i = contigCount; i < contigCount + edgeNumber; i++){
        set_binary(lp, i + 1, TRUE);
    }
	
    p--;
    
    index=0;
    for(long int i = 0; i < p; i++){
        colno[index] = contigCount + i + 1; 
        row[index] = weight[index];
        index++;
    }
    if(!set_obj_fnex(lp, index, row, colno)){
        printf("couldn't set_obj_fnex");
        exit(0);
    }
	set_add_rowmode(lp, FALSE);
	set_timeout(lp, 600);

    set_minim(lp); 

    ret = solve(lp);

    if(!(ret==0 || ret ==1)){
		cout<<"ee:"<<ret<<endl;
		return ret;
    }
    REAL * pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
    get_primal_solution(lp, pv);

    double temp = 1;
    long int result = 0;
    for(long int i = contigCount + constraintNumber + 1; i < constraintNumber + contigCount + p + 1; i++){
		if(pv[i] != 0){
           RemoveEdgeInScaffoldGraph(scaffoldGraphHead, edgeLeftNode[i-contigCount-constraintNumber-1], edgeRightNode[i-contigCount-constraintNumber-1]);
           result++;
        }
		
    }
	
	for(long int i = constraintNumber + 1; i < contigCount + constraintNumber + 1; i++){
		contigPosition[i - constraintNumber - 1] = pv[i];
		
    }
    
    delete [] weight;
    delete [] edgeLeftNode;
    delete [] edgeRightNode;
    delete [] pv;
    
    delete_lp(lp);
    return result;
}

void OptimizeScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead){
	
	bool * contigOrientation = (bool *)malloc(sizeof(bool) * scaffoldGraphHead->scaffoldGraphNodeCount);
	long int * contigPisition = (long int *)malloc(sizeof(long int) * scaffoldGraphHead->scaffoldGraphNodeCount);
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		contigOrientation[i] = false;
		contigPisition[i] = 1;
	}
	
	RemoveOrientationContradictionInScaffoldGraph(scaffoldGraphHead, contigSetHead, scaffoldGraphHead->scaffoldGraphNodeCount, contigOrientation);
	
	DetermineContigOrderInScaffoldGraph(scaffoldGraphHead, contigSetHead, contigOrientation, contigPisition);
	
	RecoverTwoEdgeInScaffoldGraph(scaffoldGraphHead);
	
	//OutputScaffoldGraphHead(scaffoldGraphHead, contigSetHead, contigPisition);
	
}

ScaffoldSetHead * GetScaffoldingResult(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead){
	
	ScaffoldSetHead * scaffoldSetHead = (ScaffoldSetHead *)malloc(sizeof(ScaffoldSetHead));
	scaffoldSetHead->scaffoldSet = NULL;
	
	bool * visited = (bool *)malloc(sizeof(bool)*contigSetHead->contigCount);
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		visited[i] = false;
	}
	
	long int startNodeIndex = -1;
	while(true){
		startNodeIndex = SelectNodeFromContigSet(contigSetHead, visited);
		if(startNodeIndex == -1){
			break;
		}
		GetSingleScaffold(scaffoldSetHead, scaffoldGraphHead, contigSetHead, startNodeIndex, visited);
	}
	return scaffoldSetHead;
	
}

void GetSingleScaffold(ScaffoldSetHead * scaffoldSetHead, ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, long int startNodeIndex, bool * visited){
	bool result = false;
	long int nodeIndex = startNodeIndex;
	
	ScaffoldSet * tempScaffoldSet = (ScaffoldSet *)malloc(sizeof(ScaffoldSet));
	tempScaffoldSet->next = NULL;
	if(scaffoldSetHead->scaffoldSet != NULL){
		tempScaffoldSet->next = scaffoldSetHead->scaffoldSet;
	}
	scaffoldSetHead->scaffoldSet = tempScaffoldSet;
	
	ContigSequence * tempContigSequence = (ContigSequence *)malloc(sizeof(ContigSequence));
	tempContigSequence->index = nodeIndex;
	tempContigSequence->orientation = 1;
	tempContigSequence->gapDistance = 0;
	tempContigSequence->next = NULL;
	scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequence;
	
	ContigSequence * tempContigSequenceNext = tempContigSequence;
	
	
	bool tail = true;
	
	while(true){
		tempContigSequenceNext = DeterminNextNode(scaffoldSetHead, scaffoldGraphHead, tempContigSequenceNext, visited, tail);
		if(tempContigSequenceNext == NULL){
			tail = false;
			tempContigSequenceNext = tempContigSequence;
			while(true){
				tempContigSequenceNext = DeterminNextNode(scaffoldSetHead, scaffoldGraphHead, tempContigSequenceNext, visited, tail);
				if(tempContigSequenceNext == NULL){
					visited[nodeIndex] = true;
					break;
				}
			}
			break;
		}
	}
	
}

ContigSequence *  DeterminNextNode(ScaffoldSetHead * scaffoldSetHead, ScaffoldGraphHead * scaffoldGraphHead, ContigSequence * tempContigSequence, bool * visited, bool & tail){
	long int nodeIndex = tempContigSequence->index;
	double maxWeight = 1000000;
	long int index = -1;
	
	bool type = false;
	
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraph[nodeIndex].edgeCount; i++){
		if(scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[i].orientationType == 0 || scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[i].orientationType == 2){
			type = false;
		}else{
			type = true;
		}
		if(maxWeight > scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[i].weight && type == tail && visited[scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[i].nodeIndex] != true ){
			maxWeight = scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[i].weight;
			index = i;
		}
	}
	if(index < 0){
		return NULL;
	}
	visited[scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].nodeIndex] = true;
	
	ContigSequence * tempContigSequenceNext = (ContigSequence *)malloc(sizeof(ContigSequence));
	tempContigSequenceNext->index = scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].nodeIndex;
	tempContigSequenceNext->gapDistance = 0;
	
	if(tempContigSequence->orientation == 1){
		if(tail == true && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 1){
			tempContigSequenceNext->orientation = 0;
			tempContigSequenceNext->next = NULL;
			tempContigSequence->next = tempContigSequenceNext;
			tail = false;
			tempContigSequence->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
		if(tail == true && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 3){
			tempContigSequenceNext->orientation = 1;
			tempContigSequenceNext->next = NULL;
			tempContigSequence->next = tempContigSequenceNext;
			tail = true;
			tempContigSequence->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
		if(tail == false && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 0){
			tempContigSequenceNext->orientation = 0;
			tempContigSequenceNext->next = tempContigSequence;
			scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequenceNext;
			tail = true;
			tempContigSequenceNext->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
		if(tail == false && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 2){
			tempContigSequenceNext->orientation = 1;
			tempContigSequenceNext->next = tempContigSequence;
			scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequenceNext;
			tail = false;
			tempContigSequenceNext->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
	}else{
		if(tail == true && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 1){
			tempContigSequenceNext->orientation = 1;
			tempContigSequenceNext->next = tempContigSequence;
			scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequenceNext;
			tail = false;
			tempContigSequenceNext->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
		if(tail == true && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 3){
			tempContigSequenceNext->orientation = 0;
			tempContigSequenceNext->next = tempContigSequence;
			scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequenceNext;
			tail = true;
			tempContigSequenceNext->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
		if(tail == false && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 0){
			tempContigSequenceNext->orientation = 1;
			tempContigSequenceNext->next = NULL;
			tempContigSequence->next = tempContigSequenceNext;
			tail = true;
			tempContigSequence->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
		if(tail == false && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].orientationType == 2){
			tempContigSequenceNext->orientation = 0;
			tempContigSequenceNext->next = NULL;
			tempContigSequence->next = tempContigSequenceNext;
			tail = true;
			tempContigSequence->gapDistance = GetGapDistance(index, nodeIndex, scaffoldGraphHead);
			//tempContigSequence->gapDistance = 500;
		}
	}
	
	return tempContigSequenceNext;
}

int GetGapDistance(long int index, long int nodeIndex, ScaffoldGraphHead* scaffoldGraphHead) {
	
	if (scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].gapDistance >= 0.4) {
		return 100;
	}
	else if (scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].gapDistance >= 0.2 && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].gapDistance < 0.4) {
		return 500;
	}
	else if (scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].gapDistance >= 0.1 && scaffoldGraphHead->scaffoldGraph[nodeIndex].edge[index].gapDistance < 0.2) {
		return 800;
	}
	else{
		return 1000;
	}
}

long int SelectNodeFromContigSet(ContigSetHead * contigSetHead, bool * visited){
	
	long int max = -1;
	long int index = -1;
	
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		if(visited[i] == true){
			continue;
		}
		if(contigSetHead->contigSet[i].contigLength > max){
			max = contigSetHead->contigSet[i].contigLength;
			index = i;
		}
	}
	
	return index;
	
	
}

void OutPutScaffoldSet(ScaffoldSet * scaffoldSet, ContigSetHead * contigSetHead, char * result){    
    long int i = 0;
    long int j = 0;

	Contig * contigSet = contigSetHead->contigSet;
	int contigCount = contigSetHead->contigCount;
	
	fflush(stdout);  
    setvbuf(stdout,NULL,_IONBF,0);
	
    bool * printContigIndex = new bool[contigCount];
    for(i = 0; i < contigCount; i++){
        printContigIndex[i] = false;
    }
	
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    scaffoldSet = tempScaffoldSet;
    
    char * scaffoldTagFileName = new char[1000];
    strcpy(scaffoldTagFileName, result);
    strcat(scaffoldTagFileName, "_tag.fa");
    ofstream ocoutTag;
    ocoutTag.open(scaffoldTagFileName);
    
    char * scaffoldSetFileName = new char[1000];   
    strcpy(scaffoldSetFileName, result);
    strcat(scaffoldSetFileName, "_set.fa");
	
	FILE * fp; 
    if((fp = fopen(scaffoldSetFileName, "w")) == NULL){
        printf("%s, does not exist!", scaffoldSetFileName);
        exit(0);
    }
	
	long int line = 200;
	char * tempLine = (char *)malloc(sizeof(char)*(line + 1));
    j = 0;
	int tempLength = 0;
	int tempGapDis = 0;
	
	
    while(tempScaffoldSet != NULL){
        
		ContigSequence * tempContigSequence = tempScaffoldSet->contigSequence;
        if(tempContigSequence == NULL){
            ocoutTag<<endl;
			tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        
        fprintf(fp, ">%ld\n", j);
        long int allLength = 0;
        tempLength = 0;
		tempGapDis = 0;
        while(tempContigSequence != NULL){
            
            if(printContigIndex[tempContigSequence->index] == true){
                ocoutTag<<"--";
            }
           	tempLength = tempGapDis + tempLength + strlen(contigSet[tempContigSequence->index].contig);
            printContigIndex[tempContigSequence->index] = true;
            ocoutTag<<tempContigSequence->index<<"("<<tempContigSequence->gapDistance<<"--"<<tempContigSequence->orientation<<"--"<<tempLength<<"),";
            tempGapDis = tempContigSequence->gapDistance;
            if(tempContigSequence->orientation==0){
                char * tempChar1 = ReverseComplement(contigSet[tempContigSequence->index].contig);
				if(tempContigSequence->gapDistance < 0 && tempContigSequence->next!=NULL){
					
					long int contigLength = strlen(tempChar1);
					if(contigLength + tempContigSequence->gapDistance < 0){
						tempContigSequence = tempContigSequence->next;
						continue;
					}
					
                	char * tempChar2 = (char *)malloc(sizeof(char)*(contigLength + tempContigSequence->gapDistance + 1));
					
					strncpy(tempChar2, tempChar1, contigLength + tempContigSequence->gapDistance);
					
					tempChar2[contigLength + tempContigSequence->gapDistance] = '\0';
					
                	free(tempChar1);
					
					tempChar1 = tempChar2;
            	}
                
				fprintf(fp, "%s",tempChar1);
				
				free(tempChar1);
            }else{
				if(tempContigSequence->gapDistance<0 && tempContigSequence->next!=NULL){
                	
					long int contigLength = strlen(contigSet[tempContigSequence->index].contig);
					
					if(contigLength + tempContigSequence->gapDistance < 0){
						tempContigSequence = tempContigSequence->next;
						continue;	
					}
					
					char * tempChar2 = (char *)malloc(sizeof(char)*(contigLength + tempContigSequence->gapDistance + 1));
					strncpy(tempChar2, contigSet[tempContigSequence->index].contig, contigLength + tempContigSequence->gapDistance);
					
					tempChar2[contigLength + tempContigSequence->gapDistance] = '\0';										
					fprintf(fp, "%s",tempChar2);					
					free(tempChar2);
            	}else{					
					fprintf(fp, "%s",contigSet[tempContigSequence->index].contig);
				}               
            }            
            if(tempContigSequence->gapDistance>=0 && tempContigSequence->next!=NULL){
                long int cc = tempContigSequence->gapDistance;
				for(long int tt = 0; tt<cc; tt++){
                     fprintf(fp, "N");
                }
            }
            tempContigSequence = tempContigSequence->next;
        }    
        ocoutTag<<endl;
		fprintf(fp, "\n");
        j++;
        tempScaffoldSet = tempScaffoldSet->next;
    }
    ocoutTag<<"----------------------------------------------------------------"<<endl;
    for(i = 0; i<contigCount; i++){
        if(printContigIndex[i] == false && contigSet[i].contig != NULL){
			fprintf(fp, ">%ld\n", j);
            fprintf(fp, "%s\n",contigSet[i].contig);
			ocoutTag<<i<<","<<endl;
            j++;
        }
    }
}

#endif
