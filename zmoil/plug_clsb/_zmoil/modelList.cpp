// @ZBS {
//		*MODULE_DEPENDS zuilist2.h zuilist2.cpp
// }

// module
#include "_zmoil.h"
#include "modelList.h"
#include "mainutil.h"
#include "cmoil.h"

// zbslib
#include "tbutil.h"
#include "ztmpstr.h"
#include "zuilist2.h"
#include "ztime.h"
#include "zwildcard.h"
#include "zmathtools.h"

extern char *escapeQuotes( char* string );

// stdlib/os
#ifdef WIN32
#include "malloc.h"
#else
#include "memory.h"
#endif

#include "string.h"
#include "ctype.h"
#include "float.h"
#include "math.h"
#include "wingl.h"
#include "GL/gl.h"
#include "GL/glu.h"

#ifdef __USE_GNU
#define stricmp strcasecmp
#endif

#define TRACE_HICOLOR 0x80C099FF

using namespace CMOIL;
extern CMOIL::AtomArray * manualAlign;

// List GUI Create
//----------------------------------------------------------------------------------

ZUI * dataZuiCreate( ZUIList *list, ListItem *li, ZUI *panel ) {
	
	switch( li->itemType ) {
		case ZIT_Chain: {
			Chain *chain = (Chain*)li->item;
			AtomArray *aa = (AtomArray*)li->parent->parent->item;

			
			ZUI *zname = ZUI::factory( 0, "ZUIText" );
			//zname->putS( "layout_forceH", "20" );
			zname->putS( "text", ZTmpStr( "Chain %c", chain->id ) );
			zname->attachTo( panel );

			list->setCheckboxMessages( li, ZTmpStr( "type=Zmoil_ChainShow show=1 model=%ld chainId=%d", long(aa), chain->id ),
				ZTmpStr( "type=Zmoil_ChainShow show=0 model=%ld chainId=%d", long(aa), chain->id ), chain->display );
		}
	    break;

		case ZIT_Model: {
			Model *model = (Model*)li->item;
			AtomArray *aa = (AtomArray*)li->parent->item;
			
			panel->putS( "pack_side", "t" );
		
			ZUI *name = ZUI::factory( 0, "ZUIText" );
			name->putS( "layout_cellFill", "wh" );
			name->attachTo( panel );
			name->putS( "text", ZTmpStr( "Model %d", model->modelName )  );
		
			ZUI *zdesc = ZUI::factory( 0, "ZUIText" );
			zdesc->putS( "layout_cellFill", "wh" );
			zdesc->attachTo( panel );
			zdesc->putS( "text", "some description here" );

			list->setCheckboxMessages( li, ZTmpStr( "type=Zmoil_ModelShow show=1 model=%ld modelName=%d", long(aa), model->modelName ),
									       ZTmpStr( "type=Zmoil_ModelShow show=0 model=%ld modelName=%d", long(aa), model->modelName ) );
		}
		break;

		case ZIT_AtomArray: {
			AtomArray *aa = (AtomArray*)li->item;

			list->setCheckboxMessages( li, ZTmpStr( "type=Zmoil_ModelDraw draw=1 model=%ld", long(aa) ),
										   ZTmpStr( "type=Zmoil_ModelDraw draw=0 model=%ld", long(aa) ) );

			int colorIndex = aa->properties.getI( "colorIndex", -1 );
			if( colorIndex >= 0 ) {
				list->setItemColorTag( li, tRGBA_to_Int( zmoilColors[colorIndex] ) );
			}

			panel->putS( "pack_side", "t" );
			
			ZUI *zsrc = ZUI::factory( 0, "ZUIText" );
			zsrc->putS( "layout_cellFill", "wh" );
			zsrc->attachTo( panel );
			zsrc->putS( "text", getTail( aa->InParm.fileName, 30 ) );
			
			ZUI *desc = ZUI::factory( 0, "ZUIText" );
			desc->attachTo( panel );
			ZTmpStr def( "%d Residues, %d Atoms", aa->numRes, aa->numAtom );
			desc->putS( "text", def.s );

			ZUI *subpanel = ZUI::factory( 0, "ZUIPanel" );
			subpanel->putS( "layout_cellFill", "wh" );
			subpanel->putS( "pack_side", "l" );
			subpanel->putS( "panelColor", "0" );
			subpanel->attachTo( panel );
			
			ZUI *z = ZUI::factory( 0, "ZUICheck" );
			z->putS( "text", "Blink   " );
			z->putI( "selected", aa->bBlink );
			z->putS( "sendMsgOnSelect", ZTmpStr( "type=Zmoil_ModelBlink blink=1 model=%ld", long(aa) )  );
			z->putS( "sendMsgOnUnselect", ZTmpStr( "type=Zmoil_ModelBlink blink=0 model=%ld", long(aa) ) );
			z->attachTo( subpanel );
			
			// cycle toggle
			z = ZUI::factory( 0, "ZUICheck" );
			z->putS( "text", "Cycle   " );
			z->putI( "selected", 0 );
			z->putI( "disabled", aa->numModels < 2 );
			z->putS( "sendMsgOnSelect", ZTmpStr( "type=Zmoil_ModelCycle cycle=1 model=%ld", long(aa) )  );
			z->putS( "sendMsgOnUnselect", ZTmpStr( "type=Zmoil_ModelCycle cycle=0 model=%ld", long(aa) ) );
			z->attachTo( subpanel );
			
			// opacity control for this model
			z = ZUI::factory( 0, "ZUIText" );
			z->putS( "text", "Opacity:" );
			z->attachTo( subpanel );
			
			z = ZUI::factory( 0, "ZUIVar" );
			z->putS( "style", "styleCmoilTwiddler" );
			z->putS( "type", "int" );
			z->putF( "layout_forceW", 40 );
			z->putI( "val", 255 );
			z->putP( "valPtr", &(aa->alpha) );
			z->putI( "rangeLow", 0 );
			z->putI( "rangeHigh", 255 );
			z->attachTo( subpanel );
		}
			break;
	}
	return panel;
}

void dataItemRender( ListItem *li, float _x, float _y, float _w, float _h, int selected ) {
	// this is called *before* the zui for the data item renders, allowing us to paint the 
	// background to indicate selection, etc...
	assert( li->itemZui );

	if( selected ) {
		// custom rendering, before itemZui drawn.
		int color = ZUI::colorPaletteHash.getU( "dataFrameSelected" );
		glBegin( GL_QUADS );
			glColor4ub( (color&0xFF000000)>>24, (color&0x00FF0000)>>16, (color&0x0000FF00)>>8, color&0x000000FF );
			glVertex2f( _x, _y );
			glVertex2f( _x+_w, _y );
			glVertex2f( _x+_w, _y+_h );
			glVertex2f( _x, _y+_h );
		glEnd();
	}
}

int addDataToList( AtomArray **aaptrs, int count, ZUIList *list, int setSelected=0 ) {
	// Add the collections and traces in pd to the list.
	// If seriesCount is given, then rules governing which data
	// should be added and which should be selectable are applied.

	int countAdded = 0;
	for( int i=0; i<count; i++ ) {
		AtomArray *aa = Main::AAptrs[ i ];
		ListItem *liAA = list->addItem( aa, ZIT_AtomArray );
		countAdded++;

		/*
		if( setSelected && tc->properties.getI( "highlight", 0 ) ) {
			list->setSelectedListItem( liCollection );
		}
	    */

		int mcount = aa->numModels;
		for( int m=0; m<mcount; m++ ) {
			if( m==49 ) {
				int d=1;
			}
			ListItem *liModel = list->addItem( aa->models + m, ZIT_Model, liAA ); 
			countAdded++;
			
			for( int c=0; c<aa->numChains; c++ ) {
				if( aa->chains[ c ].modelIdx == m ) {
					ListItem *liChain = list->addItem( aa->chains + c, ZIT_Chain, liModel );
					countAdded++;
				}
			}
		}
	}
	return countAdded;
}

void dataListClear() {
	ZUIList *list = (ZUIList*)ZUI::zuiFindByName( "dataList" );
	if( list ) {
		list->reset();
		list->dirty();
	}
}

void dataListUpdate( int listScroll ) {
	// Populate the raw data list with any data that the current
	// project has.  

	ZUIList *list = (ZUIList*)ZUI::zuiFindByName( "dataList" );
	if( list ) {
		list->reset();
		list->dirty();
		list->putI( "scrollYStep", 40 );
		list->itemZuiCreate = dataZuiCreate;
		list->itemRender = dataItemRender;
		list->putS( "clickMsg", "type=dataItemClick" );
		list->putF( "scrollY", (float)listScroll );
		list->putI( "treeView", 1 );
		list->putI( "itemPanelCheckbox", 1 );

		addDataToList( Main::AAptrs, Main::NumImages, list, 1 );
		list->expandItem( -1, 0 );
			// collapse all items by default
		int itemsVisible = list->showChildren( -1, 0 );
		if( list->selectedId >= 0 ) {
			itemsVisible = list->showChildren( list->selectedId, 0 );
		}
	}
}

ZMSG_HANDLER( dataItemClick ) {
	ListItem *li = (ListItem*)STRING_TO_PTR( msg, listItem );
	ListItem *lastSelected = (ListItem*)STRING_TO_PTR( msg, lastSelected );
	CMOIL::AtomArray *aa = 0;
	CMOIL::Model *model = 0;
	CMOIL::Chain *chain = 0;
	if( li ) {
		if( li->itemType == ZIT_AtomArray ) {
			aa = (AtomArray*)li->item;
		}
		else if( li->itemType == ZIT_Model ) {
			model = (Model*)li->item;
			aa = (AtomArray*)li->parent->item;
		}
		else if( li->itemType == ZIT_Chain ) {
			chain = (Chain*)li->item;
			model = (Model*)li->parent->item;
			aa    = (AtomArray*)li->parent->parent->item;
		}
		else {
			assert( !"Bad itemType in modelList" );
		}
	}
	if( aa && manualAlign ) {
		manualAlign = aa;
	}
	if( li != lastSelected && aa ) {
		updateIndexSlider( 0 );
	}
}

//-------------------------------------------------------------

CMOIL::AtomArray * modelListGetSelected( char* listName, CMOIL::Model** _model, CMOIL::Chain** _chain ) {
	// helper fn to get app-specific items from ZUIList
	
	CMOIL::AtomArray *aa = 0;
	CMOIL::Model *model = 0;
	CMOIL::Chain *chain = 0;
	
	ZUIList *list = (ZUIList*)ZUI::zuiFindByName( listName );
	if( list ) {
		ListItem *li = list->getSelectedListItem();
		if( li ) {
			if( li->itemType == ZIT_AtomArray ) {
				aa = (AtomArray*)li->item;
			}
			else if( li->itemType == ZIT_Model ) {
				model = (Model*)li->item;
				aa = (AtomArray*)li->parent->item;
			}
			else if( li->itemType == ZIT_Chain ) {
				chain = (Chain*)li->item;
				model = (Model*)li->parent->item;
				aa    = (AtomArray*)li->parent->parent->item;
			}
			else {
				assert( !"Bad itemType in modelList" );
			}
		}
	}
	if( _model ) { *_model = model; }
	if( _chain ) { *_chain = chain; }
	
	return aa;
}

void setAllListItems( char *listName, int state ) {
	ZUIList *list = (ZUIList*)ZUI::zuiFindByName( listName );
	if( list ) {
		for( int i=0; i<list->items.count; i++ ) {
			list->setCheckboxState( &list->items[i], state );
		}
	}
}

//-------------------------------------------------------------
	

