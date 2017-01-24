#ifndef _CMOIL_AMINOACIDS_H
#define _CMOIL_AMINOACIDS_H

#include <stdio.h>
#include <string.h>

namespace CMOIL {

const int NUM_LETT=25;     // number of alphabet used in AA

class AminoAcidNode
{
public:
  int index;            // index of AA[] array
  int type;             // type of AA_TYPE 
  int isAA;             // is an AminoAcid
  AminoAcidNode *kids[NUM_LETT];

  AminoAcidNode() {
    for (int i=0; i<NUM_LETT; i++)
      kids[i]=NULL;
  };
};

class AminoAcids 
{
public:
  enum AA_DISTRIBUTION {
  	MAX_isAA_HYD_index=10, MAX_isAA_POL_index=15, MAX_isAA_CHG_index=17, MAX_isAA_CHN_index=19,
  	MAX_isAA_UNK_index=23, MAX_noAA_POL_index=24, MAX_noAA_UNK_index=25  };
  enum { NUM_isAA=MAX_isAA_UNK_index+1,NUM_AA=MAX_noAA_UNK_index+1};   // NUM_AA:number of standard amino acids
  enum AA_TYPE { AA_HYD=0, AA_POL, AA_CHG, AA_CHN, AA_UNK, AA_UNI};    // amino acid types  

private:
  AminoAcidNode *all;
  static char AA[NUM_AA][4];

public:
  AminoAcids();
  ~AminoAcids();
  
  void PrintNames();
  int  GetType(char *aaName);
  int  IsAA(char *aaName);      // is amino acid?
};

const int NUM_AA_TYPE= AminoAcids::AA_UNI+1;

extern AminoAcids AminoAcidsObj;

}

#endif
