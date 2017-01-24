#ifndef _CMOIL_MSGBOARD_H
#define _CMOIL_MSGBOARD_H
//****************************************************************************************************
//*  Filename:   msgboard.hpp
//*
//*  Description: class for display text in the bottom of the window
//*
//*  History:
//*  Date		Developer	Description
//*  ------------ ----------- ------------------------------------------------------------------------
//*  June 2001	Baohua Wang	Initial Development
//*
//****************************************************************************************************
#include "cmoil_const.h"

namespace CMOIL {
class MessageBoard {
private:
  char text[LINE_MAXLEN];
  int  displayonce;                // display at lease once, will not erase

public:
  MessageBoard();
  MessageBoard(char *);
  void set(char *);   			// copy message to string "text" for display
  void append(char *);			// append message for the display
  void erase();				// erase all messages
  void print(int inBlackBg=1, int lowerRight=0 );	// if inBlackBg=true, text is printed on black background
  void once();                      // display once only
  char *GetText() { return text; };
	
};

}

#endif

