/*
 * gifmerge.c
 *
 * Copyright (C) 1990,1991,1992,1993 by Mark Podlipec. 
 * All rights reserved.
 *
 * This software may be freely copied, modified and redistributed
 * without fee provided that this copyright notice is preserved 
 * intact on all copies and modified copies.
 * 
 * There is no warranty or other guarantee of fitness of this software.
 * It is provided solely "as is". The author(s) disclaim(s) all
 * responsibility and liability with respect to this software's usage
 * or its effect upon hardware or computer systems.
 *
 */
 /*
  * Description:
  *
  * This program reads a gif91(see XAnim docs) file and merges the listed
  * gif files into one gif file. 
  *
  * Eventually, I'd like to have this program compare the current image
  * with the previous image and check to see if only a small section of
  * the screen changed from the previous image. Worth a shot.
  */

 /*
  * Rev 1.00	23Jul91	Mark Podlipec
  *	creation
  * Rev 1.01	08Jan92	Mark Podlipec
  *     use all colormaps, not just 1st.
  *
  * Rev 1.2    20Dec95  Rene Mueller
  *     command-line input (no longer txtfile needed)
  *
  * Rev 1.3    05Feb96  Rene Mueller (kiwi@iis.ee.ethz.ch)
  *     GIF89a transparency, and "Netscape2.0" application extension (looping)
  * Rev 1.31   14May96  Rene Mueller (kiwi@iis.ee.ethz.ch)
  *     disposal selectable
  * Rev 1.32   16Jul96  Rene Mueller (kiwi@iis.ee.ethz.ch)
  *     logical position per image manipulating
  * Rev 1.33   22Jul96  Rene Mueller (kiwi@iis.ee.ethz.ch)
  *     -notransp and -nopos added
  * Rev 1.34   29Jul96  Rene Mueller (kiwi@iis.ee.ethz.ch)
  *     -extract, extract a merged gif into single GIFs,
  *      ideal to "steal" animation and customize them ;-)
  */
/***
  * 02Oct03   Baohua Wang (bwang@tc.cornell.edu)
  *     Modified to merge gif files for MOIL-cmoil with VC++. 
  *     Added namespace and following functions:
  *        GIF_Merge_Start(), GIF_Merge(), and GIF_Merge_End().  
  *     Delete main(). 
***/

#define DA_REV 1.33
#include "gifmerge.h"
#include <iostream>
#include <stdlib.h>

namespace GifMerge {

using namespace std;

#define MAXVAL  4100            /* maxval of lzw coding size */
#define MAXVALP 4200

int debug_flag = 0;  /* make these options */
int verbose = 0;  /* make these options */
int imagex = 0;
int imagey = 0;
int imagec = 0;

GIF_Color gif_cmap[256];

ULONG GIF_Get_Code(FILE *, FILE *);
void GIF_Decompress(FILE *fp, FILE *fout);
int  GIF_Get_Short(FILE *fp,FILE *fout,int first_time);
int  GIF_Put_Short(FILE *fout,unsigned int data);
int  GIF_Get_Next_Entry(FILE *fp);
void GIF_Add_To_Table(register ULONG body,register ULONG next,register ULONG index);
void GIF_Send_Data(register int index);
void GIF_Clear_Table();
void GIF_Screen_Header(FILE *fp,FILE *fout,int first_time);
int GIF_Get_Short(FILE *fp,FILE *fout,int first_time);
void GIF_Image_Header(FILE *fp,FILE *fout,int first_time);

GIF_Screen_Hdr gifscrn;
GIF_Image_Hdr gifimage;
GIF_Table table[MAXVALP];

ULONG root_code_size,code_size,CLEAR,EOI,INCSIZE;
ULONG nextab;
ULONG gif_mask[16] = {1,1,3,7,15,31,63,127,255,511,1023,2047,4095,8191,0,0};
ULONG gif_ptwo[16] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,0,0};
UBYTE gif_buff[MAXVALP];
ULONG gif_block_size;
int num_bits,bits;

int pic_i;
char gif_file_name[256];
int screen_was_last;

int disposal = 2;
int repeats = -1;
int delay = 50;
int transp = -1;
int pos_set = 0;
int xpos = 0, ypos = 0;

FILE *MergedFptr=NULL;

int TheEnd()
{
 return -1;
}

int TheEnd(char *p)
{
 fprintf(stderr,"%s",p);
 return TheEnd();
}

void GIF_Decompress(FILE *fp, FILE *fout)
{
 register ULONG code,old;

 pic_i = 0;
 bits=0;
 num_bits=0;
 gif_block_size=0;
    /* starting code size of LZW */
 root_code_size=(fgetc(fp) & 0xff); fputc(root_code_size,fout);
 GIF_Clear_Table();                /* clear decoding symbol table */

 code=GIF_Get_Code(fp,fout);

 if (code==CLEAR) 
 {
  GIF_Clear_Table(); 
  code=GIF_Get_Code(fp,fout);
 }
 /* write code(or what it currently stands for) to file */
 GIF_Send_Data(code);   
 old=code;
 code=GIF_Get_Code(fp,fout);
 do
 {
  if (table[code].valid==1)    /* if known code */
  {
       /* send it's associated string to file */
    GIF_Send_Data(code);
    GIF_Get_Next_Entry(fp);       /* get next table entry (nextab) */
    GIF_Add_To_Table(old,code,nextab);  /* add old+code to table */
    old=code;
  }
  else      /* code doesn't exist */
  {
    GIF_Add_To_Table(old,old,code);   /* add old+old to table */
    GIF_Send_Data(code);
    old=code;
  }
  code=GIF_Get_Code(fp,fout);
  if (code==CLEAR)
  { 
   GIF_Clear_Table();
   code=GIF_Get_Code(fp,fout);
   GIF_Send_Data(code);
   old=code;
   code=GIF_Get_Code(fp,fout);
  }
 } while(code!=EOI);
}

int GIF_Get_Next_Entry(FILE *fp)
{
   /* table walk to empty spot */
 while(  (table[nextab].valid==1)
       &&(nextab<MAXVAL)
      ) nextab++;
 /* 
  * Ran out of space??!?  Something's roached 
  */
 if (nextab>=MAXVAL)    
 { 
  fprintf(stderr,"Error: GetNext nextab=%ld\n",nextab);
  fclose(fp);
  return TheEnd();
 }
 if (nextab==INCSIZE)   /* go to next table size (and LZW code size ) */
 {
   /* fprintf(stderr,"GetNext INCSIZE was %ld ",nextab); */
   code_size++; INCSIZE=(INCSIZE*2)+1;
   if (code_size>=12) code_size=12;
/*   fprintf(stderr,"<%ld>",INCSIZE); */
 }
 return 0;
}
/*  body is associated string
    next is code to add to that string to form associated string for
    index
*/     

void GIF_Add_To_Table(register ULONG body,register ULONG next,register ULONG index)
{
 if (index>MAXVAL)
 { 
  fprintf(stderr,"Error index=%ld\n",index);
 }
 else
 {
  table[index].valid=1;
  table[index].data=table[next].first;
  table[index].first=table[body].first;
  table[index].last=body;
 }
}

void GIF_Send_Data(register int index)
{
 register int i,j;
 i=0;
 do         /* table walk to retrieve string associated with index */
 { 
  gif_buff[i]=table[index].data; 
  i++;
  index=table[index].last;
  if (i>MAXVAL)
  { 
   fprintf(stderr,"Error: Sending i=%ld index=%ld\n",i,index);
   TheEnd();
  }
 } while(index>=0);

 /* now invert that string since we retreived it backwards */
 i--;
 for(j=i;j>=0;j--)
 {
  /*pic[pic_i] = gif_buff[j] | gif_pix_offset;*/
  pic_i++;
 }
}


/* 
 * initialize string table 
 */
void GIF_Init_Table()       
{
 register int maxi,i;

if (debug_flag) fprintf(stderr,"Initing Table...");
 maxi=gif_ptwo[root_code_size];
 for(i=0; i<maxi; i++)
 {
  table[i].data=i;   
  table[i].first=i;
  table[i].valid=1;  
  table[i].last = -1;
 }
 CLEAR=maxi; 
 EOI=maxi+1; 
 nextab=maxi+2;
 INCSIZE = (2*maxi)-1;
 code_size=root_code_size+1;
}


/* 
 * clear table 
 */
void GIF_Clear_Table()   
{
 register int i;
if (debug_flag) fprintf(stderr,"Clearing Table...\n");
 for(i=0;i<MAXVAL;i++) table[i].valid=0;
 GIF_Init_Table();
}

/*CODE*/
ULONG GIF_Get_Code(FILE *fp, FILE *fout) /* get code depending of current LZW code size */
{
 ULONG code;
 register int tmp;

 while(num_bits < (int)code_size)
 {
  /**** if at end of a block, start new block */
  if (gif_block_size==0) 
  {
   tmp = fgetc(fp);
   if (tmp >= 0 )
   {
    fputc(tmp,fout);
    gif_block_size=(ULONG)(tmp);
   }
   else TheEnd("EOF in data stream\n");
  }

  tmp = fgetc(fp);   gif_block_size--;
  if (tmp >= 0)
  {
   fputc(tmp,fout);
   bits |= ( ((ULONG)(tmp) & 0xff) << num_bits );
   num_bits+=8;
  }
  else TheEnd("EOF in data stream\n");
 }
  
 code = bits & gif_mask[code_size];
 bits >>= code_size;
 num_bits -= code_size; 


 if (code>MAXVAL)
 { 
  fprintf(stderr,"\nError! in stream=%lx \n",code); 
  fprintf(stderr,"CLEAR=%lx INCSIZE=%lx EOI=%lx code_size=%lx \n",
                                           CLEAR,INCSIZE,EOI,code_size); 
  code=EOI;
 }

 if (code==INCSIZE)
 {
  if (code_size<12)
  {
   code_size++; INCSIZE=(INCSIZE*2)+1;
  }
  else if (debug_flag) fprintf(stderr,"<13?>"); 
 }

 return(code);
}


/* 
 * read GIF header 
 */
void GIF_Screen_Header(FILE *fp,FILE *fout,int first_time)
{
 int temp,i;

 if (fp == NULL)
   return;

 for(i=0;i<6;i++) {
  temp = fgetc(fp);
/*   if (first_time==TRUE) fputc(temp,fout); */
 }
 if(first_time) 
   fputs("GIF89a",fout);
 gifscrn.width  = GIF_Get_Short(fp,fout,first_time);
 gifscrn.height = GIF_Get_Short(fp,fout,first_time);
 temp=fgetc(fp);		 if (first_time==TRUE) fputc(temp,fout);
 gifscrn.m       =  temp & 0x80;
 gifscrn.cres    = (temp & 0x70) >> 4;
 gifscrn.pixbits =  temp & 0x07;
 gifscrn.bc  = fgetc(fp);	 if (first_time==TRUE) fputc(gifscrn.bc,fout);
 temp=fgetc(fp);		 if (first_time==TRUE) fputc(temp,fout);
 imagec=gif_ptwo[(1+gifscrn.pixbits)];

 if (verbose)
  fprintf(stderr,"Screen: %ldx%ldx%ld m=%ld cres=%ld bkgnd=%ld pix=%ld\n",
    gifscrn.width,gifscrn.height,imagec,gifscrn.m,gifscrn.cres,
    gifscrn.bc,gifscrn.pixbits);

 if (gifscrn.m)
 {
  for(i=0;i<imagec;i++)
  {
   gif_cmap[i].cmap.red   = temp = fgetc(fp); 
           if (first_time==TRUE) fputc(temp,fout);
   gif_cmap[i].cmap.green = temp = fgetc(fp); 
           if (first_time==TRUE) fputc(temp,fout);
   gif_cmap[i].cmap.blue  = temp = fgetc(fp); 
           if (first_time==TRUE) fputc(temp,fout);
  }
 }
 screen_was_last = TRUE;
 if(gifscrn.m&&(transp>=0||delay>=0)) {
   int ix = 0, max_dist = 3*256;
   if(transp>=0) {
      for(i=0; i<imagec; i++) {
         int dist = 
            abs(gif_cmap[i].cmap.red-(transp&255))+
            abs(gif_cmap[i].cmap.green-((transp>>8)&255))+
            abs(gif_cmap[i].cmap.blue-(transp>>16));
         if(dist<max_dist) 
            ix = i, max_dist = dist;
      } 
      if(max_dist==0)   /* info at http://www.iis.ee.ethz.ch/~kiwi/GIFMerge/gifspecs.txt */   
         ;
/*          fprintf(stderr,"Transparent color matched fully\n"); */
      else
         fprintf(stderr,"Transparent not matched fully, col #%d (%d,%d,%d) used now\n",ix,
            gif_cmap[ix].cmap.red,gif_cmap[ix].cmap.green,gif_cmap[ix].cmap.blue);
   }
   fputc(0x21,fout);
   fputc(0xF9,fout);
   fputc(0x04,fout);
   fputc((transp>=0?0x01:0x00)|(disposal<<2),fout);
   fputc(delay&255,fout);
   fputc((unsigned)delay>>8,fout);
   fputc(ix,fout);
   fputc(0x00,fout);
 }
 if(first_time&&repeats>=0) { /* documentation look at */
   fputc(0x21,fout);          /* http://www.reiworld.com/royalef/gifabout.htm */
   fputc(0xFF,fout);
   fputc(0x0B,fout);
   fputs("NETSCAPE2.0",fout);
   fputc(0x03,fout);
   fputc(0x01,fout);
   fputc(repeats&255,fout);
   fputc((unsigned)repeats>>8,fout);
   fputc(0x00,fout);
 }
}

void GIF_Image_Header(FILE *fp,FILE *fout,int first_time)
{
 int temp,tnum,i,tmp;

 tmp = GIF_Get_Short(fp,fout,0); if(!pos_set) xpos = tmp;
 gifimage.left   = xpos; GIF_Put_Short(fout,xpos);
 tmp = GIF_Get_Short(fp,fout,0); if(!pos_set) ypos = tmp;
 gifimage.top    = ypos; GIF_Put_Short(fout,ypos);
 gifimage.width  = GIF_Get_Short(fp,fout,1);
 gifimage.height = GIF_Get_Short(fp,fout,1);
 temp=fgetc(fp); 

 gifimage.m        = temp & 0x80;
 gifimage.i        = temp & 0x40;
 gifimage.pixbits  = temp & 0x07;

 if (screen_was_last && (first_time==FALSE)) temp |= 0x80;
 temp &= 0xf8;
 temp |= gifscrn.pixbits;
 fputc(temp,fout);

 imagex=gifimage.width;
 imagey=gifimage.height;
 tnum=gif_ptwo[(1+gifimage.pixbits)];
 if (verbose)
  fprintf(stderr,"Image: %ldx%ldx%ld m=%ld i=%ld pix=%ld \n",
    imagex,imagey,tnum,gifimage.m,gifimage.i,gifimage.pixbits);

 /* if there is an image cmap, then use it */
 if (gifimage.m)
 {
  for(i=0;i<tnum;i++)
  {
   gif_cmap[i].cmap.red   = temp = fgetc(fp); fputc(temp,fout);
   gif_cmap[i].cmap.green = temp = fgetc(fp); fputc(temp,fout);
   gif_cmap[i].cmap.blue  = temp = fgetc(fp); fputc(temp,fout);
  }
 }  /* else if screen was last not 1st time */
 else if (screen_was_last && (first_time==FALSE))
 {
  for(i=0;i<imagec;i++)
  {
   fputc(gif_cmap[i].cmap.red  ,fout);
   fputc(gif_cmap[i].cmap.green,fout);
   fputc(gif_cmap[i].cmap.blue ,fout);
  }
 }
 screen_was_last = FALSE; 
}


/*
 *
 */
int GIF_Get_Short(FILE *fp,FILE *fout,int first_time)
{
 if ( fp != NULL)
 {
   register int temp,tmp1;
   temp=fgetc(fp);	 if (first_time==TRUE) fputc(temp,fout);
   tmp1=fgetc(fp);	 if (first_time==TRUE) fputc(tmp1,fout);
   return(temp|( (tmp1) << 8 ));
 } else
   return 0; 
}


/*
 *
 */
int GIF_Put_Short(FILE *fout,unsigned int data)
{
 if ( fout != NULL)
 {
   fputc((unsigned char)(data&255),fout);      /* lo */
   fputc((unsigned char)((data>>8)&255),fout); /* hi */
 }
 return data;
}

/*
 * Read a GIF file, outputting to fname as we go.
 * It would be faster to read and write the individual blocks,
 * but eventually we'd like to optimize based on changes from
 * previous images(ie only a small section of the image changed.
 */
void GIF_Read_File(char *fname, int first_image)
{
 FILE *fp;
 int  i;

 if ( (fp=fopen(fname,"rb"))==NULL)
 { 
  fprintf(stderr,"Can't open %s for reading.\n",fname); 
  TheEnd();
  return;
 }

 GIF_Screen_Header(fp,MergedFptr,first_image);

 /*** read until  ,  separator */
 do
 {
  i=fgetc(fp);
  if ( (i<0) && feof(fp))
  {
   fclose(fp);
   TheEnd("GIF_Read_Header: Unexpected End of File\n");
  }
 } while(i != ',');

 fputc(',',MergedFptr); /* image separator */

 GIF_Image_Header(fp,MergedFptr,first_image);

 /*** Setup ACTION for IMAGE */

 GIF_Decompress(fp,MergedFptr);
 fputc(0,MergedFptr);  /* block count of zero */

 fclose(fp);
} // GIF_Read_File()

// added by B. Wang
// mergedGif:  the file for merged GIF stream

static int NextIsFirstGIF=1;
void GIF_Merge_Start(char *mergedGif, int delayI, int repeatsI, int disposalI)
{
  if (NextIsFirstGIF && mergedGif!=NULL) 
  {
    MergedFptr=fopen(mergedGif, "w+b");
    delay=delayI;
    repeats=repeatsI;
    disposal=disposalI;
  } 
}

// singleGif:   the gif file going to be merged into mergedGif
int GIF_Merge(char *singleGif)
{
  if (MergedFptr==NULL)
    return -1;

  GIF_Read_File(singleGif, NextIsFirstGIF);
  if (NextIsFirstGIF!= 0 )
    NextIsFirstGIF=0;
  return 0;
}

// calling sequence:  [GIF_Merge_Start()]->[GIF_Merge() loop]->[GIF_Merge_End()]
void GIF_Merge_End()
{
  if (MergedFptr != NULL )
  {
    fputc(0x3B,MergedFptr);        /* GIF Trailer */
    fclose(MergedFptr);
    MergedFptr=NULL;
    NextIsFirstGIF=1;
  }
}

} //end namespace GifMerge
