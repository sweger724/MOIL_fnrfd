//*************************************************************************************************
//*  Filename:   cmoil_keysim.cpp
//*
//*  Description: 
//*    Simulate keystrokes
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Feb 5. 2002	Baohua Wang	Initial Development
//*
//*************************************************************************************************
#include <windows.h>
#include "stdio.h"

int main (int argc, char **argv)
{
  if (argc < 2 )
    return -1;
  else
  {
    HWND cmoil=FindWindow(NULL,  "Molecular Viewer 0.1a");
    if ( cmoil != NULL ) 
    {
      SetForegroundWindow(cmoil);
      for (unsigned int i=0; i < strlen(argv[1]); i++ ) 
        PostMessage(cmoil, WM_KEYDOWN, argv[1][i], WM_KEYUP);
    }
  }
  return 0;
}

//PressAButton(Button2->Handle); 

