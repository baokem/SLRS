#ifndef CONTIG_CPP_INCLUDED 
#define CONTIG_CPP_INCLUDED 
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include "contig.h"
#include <set>
using namespace std;

ContigSetHead * GetContigSet(char * contigSetFile, int contigLengthThreshold){
    
    ContigSetHead * contigSetHead = (ContigSetHead *)malloc(sizeof(ContigSetHead));
    contigSetHead->contigSet = NULL;
    contigSetHead->contigCount = 0;
    
    long int maxSize = 90000;
    char * contig = NULL;
    if(NULL == (contig = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }

    FILE * fp; 
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    while((fgets(contig, maxSize, fp)) != NULL){ 
       if(contig[0] == '>'){  
           contigSetHead->contigCount++; 
       }  
    }  
    fclose(fp);
    
    contigSetHead->contigSet = (Contig *)malloc(sizeof(Contig)*contigSetHead->contigCount);
	
    for(long int i = 0; i < contigSetHead->contigCount; i++){
        contigSetHead->contigSet[i].contig = NULL;
        contigSetHead->contigSet[i].contigLength = 0;
    }
    
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    
    long int allocateLength = 0;
    long int contigIndex = -1;
    while((fgets(contig, maxSize, fp)) != NULL){ 
       
       if(contig[0] == '>'){  
           
		   if(strlen(contig) == maxSize-1){              
               while((fgets(contig, maxSize, fp)) != NULL){
                   if(strlen(contig) != maxSize-1){
                       break;
                   }
               }        
           }

		   contigIndex++;
		   long int len = strlen(contig);
           continue;
       }
       
       
       long int extendLength = strlen(contig);
	   
	  
       if(contig[extendLength-1] == '\n' || contig[extendLength-1] == '\t' || contig[extendLength-1] == '\r'){
           extendLength--;
       }
	   if(contig[extendLength-1] == '\n' || contig[extendLength-1] == '\t' || contig[extendLength-1] == '\r'){
           extendLength--;
       }
	   
       long int contigLength = 0;
       char * tempContig = NULL;
       if(contigSetHead->contigSet[contigIndex].contig != NULL){
           if(contigSetHead->contigSet[contigIndex].contigLength + extendLength >= allocateLength){
               contigLength = contigSetHead->contigSet[contigIndex].contigLength;    
               contigSetHead->contigSet[contigIndex].contig = (char *)realloc(contigSetHead->contigSet[contigIndex].contig, allocateLength + maxSize + 1);

               allocateLength = allocateLength + maxSize + 1;
               
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigLength + extendLength] = '\0';    
               contigSetHead->contigSet[contigIndex].contigLength = contigLength + extendLength;
                       
           }else{
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigSetHead->contigSet[contigIndex].contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigSetHead->contigSet[contigIndex].contigLength + extendLength] = '\0';
               contigSetHead->contigSet[contigIndex].contigLength = contigSetHead->contigSet[contigIndex].contigLength + extendLength;
           }   
           
       }else{
           contigSetHead->contigSet[contigIndex].contig = (char *)malloc(sizeof(char)*(maxSize+1));
           strncpy(contigSetHead->contigSet[contigIndex].contig, contig, extendLength);
           contigSetHead->contigSet[contigIndex].contig[extendLength] = '\0';
           contigSetHead->contigSet[contigIndex].contigLength = extendLength;
           allocateLength = maxSize + 1;
       }  
    }  
	
    fflush(fp);
    fclose(fp);
	if(contigSetHead->contigCount <= 0){
		cout<<"Contig count is zero, please check the contig file!"<<endl;
		exit(0);
	}
	
    return contigSetHead;
}

void Trimend(ContigSetHead* contigSetHead) {
    for (int i = 0; i < contigSetHead->contigCount; i++) {
        if (contigSetHead->contigSet[i].contigLength > 25000) {
            contigSetHead->contigSet[i].contig[contigSetHead->contigSet[i].contigLength - 3000] = '\0';
            contigSetHead->contigSet[i].contig = contigSetHead->contigSet[i].contig + 3000;
            contigSetHead->contigSet[i].contigLength = contigSetHead->contigSet[i].contigLength - 6000;
        }
    }
}

char* ReverseComplement(char* temp) {
    long int len = strlen(temp);
    char* rcTemp = (char*)malloc(sizeof(char) * (len + 1));
    for (long int i = 0; i < len; i++) {
        if (temp[i] == 'A' || temp[i] == 'a') {
            rcTemp[len - 1 - i] = 'T';
        }
        else if (temp[i] == 'T' || temp[i] == 't') {
            rcTemp[len - 1 - i] = 'A';
        }
        else if (temp[i] == 'G' || temp[i] == 'g') {
            rcTemp[len - 1 - i] = 'C';
        }
        else if (temp[i] == 'C' || temp[i] == 'c') {
            rcTemp[len - 1 - i] = 'G';
        }
        else if (temp[i] == 'N' || temp[i] == 'n') {
            rcTemp[len - 1 - i] = 'N';
        }
    }
    rcTemp[len] = '\0';
    return rcTemp;
}

#endif