#include "cmoil_aminoacids.h"

namespace CMOIL {

char AminoAcids::AA[NUM_AA][4]= { "ALA", "CYS", "HIS", "ILE", "LEU", "MET", "PHE", "PRO", "TRP", "TYR", "VAL", // hydrophobic
                            "ASN", "GLN", "GLY", "SER", "THR",    // uncharged polar         AA_POL                     
                            "ARG", "LYS",                        // positively charged polar AA_CHG
                            "ASP", "GLU",                        // negatively charged polar AA_CHN
                            "GLX", "CYX", "HEM",  "HIP",         // unknown AA_UNK
                            // ------- for non-AAs  --------
                            "HOH",                                // AA_POL
                            // unknowns
                            "ACE"
};     // 21 amino acid name array

AminoAcids::AminoAcids() 
{
  all = new AminoAcidNode;

  AminoAcidNode *a, *b;
  for (int i=0; i<NUM_AA; i++)
  {
    a=all;
    for ( int j=0; j<3 && a!=NULL ; j++)    /* three letters for a AA name */
    {
      int k=AA[i][j]-'A';
      b = a->kids[k];
      if (b == NULL )
      {
        b= new AminoAcidNode;
        if ( b != NULL ) {
          a->kids[k]=b;
        } 
      } 
      a=b;
    }
    if ( a != NULL)
    {          
      if ( i<=MAX_isAA_HYD_index)
        a->type = AA_HYD;
      else if ( i<=MAX_isAA_POL_index)
        a->type = AA_POL; 
      else if ( i<=MAX_isAA_CHG_index)
        a->type = AA_CHG;
      else if ( i<=MAX_isAA_CHN_index )
        a->type = AA_CHN;
      else if (i<=MAX_isAA_UNK_index )
        a->type = AA_UNK;
      else if ( i<=MAX_noAA_POL_index)
        a->type = AA_POL;
      else
        a->type = AA_UNK;
     
      if (i<=MAX_isAA_UNK_index)
        a->isAA=1;      // is an amino acid
      else
        a->isAA=0;
    }
  }
}

AminoAcids::~AminoAcids() 
{ 
  if (all != NULL )
  {
    AminoAcidNode *a, *b, *c;
    for (int i=0; i<NUM_LETT; i++)
    {
      a = all->kids[i];
      if (a != NULL) 
      {
        for (int j=0; j<NUM_LETT; j++)
        {
          b=a->kids[j];
          if ( b!=NULL)
          {
            for (int k=0; k<NUM_LETT; k++)
            {
              c=b->kids[k];
              if ( c != NULL)
                delete c;
            }
            delete b;
          }
        }
        delete a;
      }
    }
  }
}

void AminoAcids::PrintNames()
{
  if ( all==NULL)
    return;

  AminoAcidNode *a, *b, *c;
  for (int i=0; i<NUM_LETT; i++)
  {
    a = all->kids[i];
    if (a != NULL) 
    {
      for (int j=0; j<NUM_LETT; j++)
      {
        b=a->kids[j];
        if ( b!=NULL)
        {
          for (int k=0; k<NUM_LETT; k++)
          {
            c=b->kids[k];
            if ( c != NULL)
              printf("%c%c%c:  %d\n", i+'A', j+'A', k+'A', c->type);
          }
        }
      }
    }
  }
}

/* get AA pol type */
int AminoAcids::GetType(char *aaName)
{
  AminoAcidNode *a=all;
  for (int i=0; i<3 && a!=NULL ; i++)
    a = a->kids[aaName[i]-'A'];
  return (a!=NULL)?a->type:AA_UNK; 
}

int AminoAcids::IsAA(char *aaName)
{
  AminoAcidNode *a=all;
  for (int i=0; i<3 && a!=NULL ; i++)
    a = a->kids[aaName[i]-'A'];
  return (a!=NULL)?a->isAA:0; 
}

AminoAcids AminoAcidsObj;
}
 

