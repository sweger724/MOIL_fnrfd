//*************************************************************************************************
//*  Filename:   cmoil_ce.cpp
//*
//*  Description: 
//*    overlap with ce program.
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  09-25-2004	Baohua Wang	Initial Development
//*
//*************************************************************************************************
#include "cmoil.h"
#include "cmoil_rms.h"
#include "cmoil_globals.h"

namespace CMOIL {

CeAlign::CeAlign() {
  alignType='B';
  numAlignedMonos=0;
  numGaps=0;
  for (int i=0; i<ALIGN_MAXLEN; i++)
    align[i]=0;
}

// use alignment array to get match directory, bypass align string presentation
// align in int[4] group, i.e. int[4] is an alignment segment.
// [4i+0][4i+1] is the start,stop of the first aligned sequence segment
// [4i+2][4i+3] is the start,stop of the sencond aligned sequence segment
int CeAlign::alignAA2Match(AtomArray &aa, COORDINATES &match, int seq_1or2)
{
  Atom *aptr=aa.atomList;
  int  *alignStart = (seq_1or2==1)?align:align+2;  // seq#1 crd start at 0, while seq#2 start at 2
  int  n=0;

  for (int i=0, rn=-1; (i<aa.numAtom); i++, aptr++)
  {
    rn = aptr->rn+1;    // change rn from index to order for comparing with align
    if ( rn >= *alignStart )
    {
      while ( (*alignStart)>0 && rn> *(alignStart+1) )     // move to the aligned segment [*align, *(align+1)]
        alignStart += 4;

      if ( *alignStart <=0 )
        break;
    
      if ( rn >= *alignStart  &&                           // in the align interval 
           ( !strncmp(aptr->atomName, "CA", 4 ) || 
             alignType=='B' && (                           // backbone
             !strncmp(aptr->atomName, "C", 4 )  ||
             !strncmp(aptr->atomName, "O", 4 )  ||
             !strncmp(aptr->atomName, "N", 4 )  )   )  )
      {
        match.x[n]=aptr->orgX;
        match.y[n]=aptr->orgY;
        match.z[n]=aptr->orgZ;
        n++;
        //printf("%d:%d:%s:%s\t", n,aptr->rn, aa.resName[aptr->rn], aptr->atomName);
      }
    }
  }
  //printf("\n\n");
  return n;
}

int CeAlign::overlap(AtomArray &aa1, AtomArray &aa2)
{
   // align two chains by running CE.exe
   //
   char aline[LINE_MAXLEN+LINE_MAXLEN]="", outfile[LINE_MAXLEN];
   double mm[3][4];

   outfilename(aa1, aa2, outfile); 

   char chain1, chain2;
   get_tmppdb(aa1, NULL, &chain1);
   get_tmppdb(aa2, NULL, &chain2);
   sprintf(aline, "%s%s %s %c %s %c %s > %s", ExeDir, WRAP_CE, aa1.InParm.crdName, chain1, aa2.InParm.crdName, chain2, ExeDir, outfile); 
   system(aline);
   FILE *outf = fopen(outfile, "r");
   if ( outf == NULL) 
     return 0;

   printf("CE output:\n--------------------\n");
   int i=0;
   char *cptr;
   while ( fgets(aline, LINE_MAXLEN, outf ) != NULL )
   {
      printf("%s", aline);
      if ( (cptr=strstr(aline, "Alignment length ="))!=NULL ) 
         sscanf(cptr, "Alignment length = %d Rmsd = %gA", &numAlignedMonos, &rmsd);
      else if (i<3 && sscanf(aline,"%*s = (%lf)*X1 + (%lf)*Y1 + (%lf)*Z1 + (%lf)", &(mm[i][0]),&(mm[i][1]),&(mm[i][2]),&(mm[i][3]) )==4 )
         i++;
   }
   fclose(outf);

   if (i < 3 )
     return 0;

   Atom *atom=aa2.atomList;
   double x,y,z;
   for (i=0; i<aa2.numAtom; i++, atom++)
   {
     x=atom->orgX*mm[0][0] + atom->orgY*mm[0][1]+atom->orgZ*mm[0][2]+mm[0][3];
     y=atom->orgX*mm[1][0] + atom->orgY*mm[1][1]+atom->orgZ*mm[1][2]+mm[1][3];
     z=atom->orgX*mm[2][0] + atom->orgY*mm[2][1]+atom->orgZ*mm[2][2]+mm[2][3];
     atom->orgX=x;
     atom->orgY=y;
     atom->orgZ=z;
   }
   aa1.ResetCrd();
   aa2.ResetCrd();
   aa1.SubtractCoMass();
   aa2.SubtractCoMass(aa1.CoMX, aa1.CoMY, aa1.CoMZ);
   return aa2.numAtom;
}

// return the number of continuous aligned segments
// return a two dim array of seq1_seg1[start,stop], seq2_seg[1][start,stop], seq2_seg2.... 
// of the segments, the array ends with {0,0} elements
// preq:  (char*)ExeDir has been set 
//
int CeAlign::align_2pdbs(char *pdb1, char chain1, char *pdb2, char chain2, char *outfile )
{
   const char thisSeqGap='-', otherSeqGap='.';  // indicates mono align type

   // align two chains by running CE.exe
   //
   char aline[LINE_MAXLEN]="", seqs[2][LINE_MAXLEN]={"",""}, *cptr=NULL, seqC;
   int  jChain=0,k, m, p ;
   unsigned int i=0;
   int  monolen[2]={0,0}, inGaps[2]={0,0},   // two chain must start with aligned mono
        nsegs[2]={0,0},   segStarts[2]={0,0}; 

   sprintf(aline, "%s%s %s %c %s %c > %s", ExeDir, WRAP_CE, pdb1, chain1, pdb2, chain2, outfile); 
   system(aline);
   FILE *outf = fopen(outfile, "r");
   if ( outf == NULL) 
     return 0;
   
   // parse CE.exe output
   align[0]=align[2]=0;
   printf("CE output:\n--------------------\n");
   while ( fgets(aline, LINE_MAXLEN, outf ) != NULL )
   {
      printf("%s", aline);
      if ( strncmp("Alignment length =", aline, 18) == 0 )
      {
         numAlignedMonos=atoi(aline+18); 
         numGaps=atoi(strstr(aline+18, "Gaps =")+6);
         jChain=0;                                // starting from first aligned sequence
      }
      else if ( sscanf(aline, "Chain%*d:%d%s", segStarts+jChain, seqs[jChain]) == 2 )    // is an alignment line?
      {                                           // segStarts & seqs are assigned  
        monolen[jChain]=0;                        // number of actual mono length from the segStarts 
        if ( align[jChain*2]==0 )     
        {            
          align[jChain*2]=segStarts[jChain];      // init align start: seq#0_align[0] and seq#1_align[2];
          inGaps[jChain]=0;                       // read in the first alined mono
        }
                         
        if ( jChain==1 )                          // the second alignment been read, process two aligned sequences
        {
          for (k=0, m=0; k<2; k++ )               // marked gaps in both chains
          {
            m=(k+1)%2;                            // k: this sequence;  m:the index of other sequence
            for (i=0; i< strlen(seqs[k]); i++ )
            {
              if ( seqs[k][i]==thisSeqGap)   
                seqs[m][i]=otherSeqGap;           // mark other alignment with '.'
              if ( seqs[k][i]=='X')               // X: unknown residue
              {
                seqs[k][i]=otherSeqGap;
                if ( seqs[m][i] != thisSeqGap ) 
                  seqs[m][i]=otherSeqGap;          
              }
            }
          }
          //printf("\n--CMOIL: Chain1: %s\n--CMOIL: Chain2: %s\n", seqs[0], seqs[1]);

          for (k=0; k<2; k++ )                    // for two structures
          {
            for (i=0; i<strlen(seqs[k]); i++ )    // processing one aligned line at a time
            {
              seqC=seqs[k][i];                    // the ith mono of this alignment
              p = 4*nsegs[k]+2*k;                 // position of current (start,stop) in the array 
              if ( seqC==thisSeqGap || seqC==otherSeqGap)
              {
                 if (inGaps[k]==0) {              // change from non-gap to gaps, output alignment
                   segStarts[k] += monolen[k];    // start of following un-aligned segment
                   monolen[k]=0;                  // m starts counting number of unalinged monos
                   align[p+1]=segStarts[k]-1;     // end of this alinged segment
                   nsegs[k]++;                    // increase the number of aligned segments
                   inGaps[k]=1;
                   //printf("%d: %d,%d\n", k, align[p], align[p+1]);
                 }   // else do nothing for staying in gaps
              } else if ( inGaps[k]==1) {         // change from gaps to non-gap
                 segStarts[k] += monolen[k];      // set no-gap start
                 monolen[k]=0;                    // m starts counting aligned monos
                 align[p]  = segStarts[k];    
                 inGaps[k]=0;
              }
              if (  seqC != thisSeqGap )
                 monolen[k]++;                               // m: length of the aligned/gaped sequence
            }
          }
        }
        jChain = (jChain+1) % 2;                  // siwtch between two chainID 1/0
     }
   }
   fclose(outf);

   for (k=0; k<2; k++)
   {
      if (inGaps[k]==0)
      {
         p = 4*nsegs[k]+2*k;            // position of current start,stop in the array 
         segStarts[k] += monolen[k]-1;                //  set non-gap stop
         align[p+1]=segStarts[k];
         nsegs[k]++;                   
      }
   }
   p++;
   align[++p]=0;                      // terminate the array
   align[++p]=0;
   align[++p]=0;
   align[++p]=0;

   printf("--------------------\nCMOIL: Aligned by CE: %d segment(s), saved in %s\n", nsegs[0],outfile);
   printf("Chain1\t  <->  Chain2\n");
   int m1=0, m2=0;
   for (p=0; p<nsegs[0]*4; p+=4)
   {
      printf("%3d,%3d\t  <->  %3d,%3d\n", align[p], align[p+1], align[p+2], align[p+3]);
      m1 += align[p+1]-align[p]+1;
      m2 += align[p+3]-align[p+2]+1;
   }
   printf("Aligned #residures:  %d, %d\n", m1, m2);
   return numAlignedMonos;
}

// auto alignment by running CE
void CeAlign::get_tmppdb(AtomArray &aa, char *tmppdb, char *chain)
{
  int chainid=-1;
  FILE *fptr=NULL;

  if ( tmppdb != NULL)
  {
    strcpy(tmppdb, aa.InParm.crdName);
    strcat(tmppdb, ".tmp");
    fptr= fopen(tmppdb, "w");
    if ( fptr==NULL)
    { 
      *tmppdb='\0';
      return;
    }
  }

  if ( aa.chains==NULL || aa.displayChainid==-1 )   
  {
    *chain='-';
    chainid=-1;
  }
  else 
  {
    *chain= (aa.chains+aa.displayChainid)->id;   // align selected chain only
     chainid = aa.displayChainid ; 
  }
  
  if (*chain==' ') 
    *chain='-';       // ce uses - for ' ' chain
  if ( tmppdb != NULL) 
  {
    aa.OutputCrdInPDB(fptr, chainid);
    fclose(fptr);
  }

}

void CeAlign::outfilename(AtomArray &aa1, AtomArray &aa2, char *outfile)
{
  strncpy(outfile, aa1.InParm.crdName, LINE_MAXLEN); 
  outfile[LINE_MAXLEN-1]='\0';
  char *c=strrchr(outfile, DIRDELIM);     // file basename
  if ( c == NULL ) 
    c=outfile;
  c=strchr(c, '.');
  if ( c != NULL)
    *c = '\0';
  strcat(outfile, "_");
  strcat(outfile, aa2.moleculeName);
  strcat(outfile, ".align");
}

void CeAlign::matches(AtomArray &aa1, AtomArray &aa2, COORDINATES &match1, COORDINATES &match2, char matchType)
{
   alignType=matchType;

   //  write to two tmp pdb files as CE requires
   //  for pdb files using current display chainID to align
   char pdb1[LINE_MAXLEN]="", pdb2[LINE_MAXLEN]="", chain1='-', chain2='-';
   char outfile[LINE_MAXLEN]="";
   get_tmppdb(aa1, pdb1, &chain1);
   get_tmppdb(aa2, pdb2, &chain2);
   outfilename(aa1, aa2, outfile); 

   // align two chains by running CE.exe
   align_2pdbs(pdb1, chain1, pdb2, chain2, outfile);
   UNLINK(pdb1);
   UNLINK(pdb2);

   // compose the alignment string
   if (align[0] != 0)
   {
      int n1=alignAA2Match(aa1, match1, 1);
      int n2=alignAA2Match(aa2, match2, 2);        // align+2:  first aligned segment of chain2
      if ( n1 < n2 )
        numAlignedAtoms=-1;
      else if ( n2<n1)
        numAlignedAtoms=-2;
      else
        numAlignedAtoms=n1;
      printf("#AutoAlignedAtoms: %d,%d\n", n1, n2);
   } else
      numAlignedAtoms=0; 
}

}

