
#ifndef MODELLIST_H
#define MODELLIST_H

#include "zuilist2.h"

enum ZmoilItemType { ZIT_Invalid=-1, ZIT_AtomArray, ZIT_Model, ZIT_Chain };

void dataListClear();
void dataListUpdate( int listScroll=0 );

namespace CMOIL {
	class AtomArray;
	class Model;
	class Chain;
}

CMOIL::AtomArray* modelListGetSelected( char* listName, CMOIL::Model** model=0, CMOIL::Chain** chain=0 );
	// output vars model, chain receive selected items if non NULL, or are NULL on return
	// returns the (parent) AtomArray selected if anything was found selected.

void setAllListItems( char *listname, int state );


#endif

