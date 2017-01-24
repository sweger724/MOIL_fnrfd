// @ZBS {
//		*MODULE_OWNER_NAME zstacktrace
// }

#ifndef ZSTACKTRACE_H
#define ZSTACKTRACE_H

int zStackTrace( int size, unsigned int *ips, unsigned int *args );
	// This function traces up the stack and creates a useful
	// stack dump good for printing after an assert of fatal.
	// Extracts instrution pointers and argumetn bases into the the two
	// passed-in arrays, not to exceed 'size'
	// Returns the number extracted, useful for passing to stackTraceGetString

void zStackTraceBuildSymbolTable();
	// Call this function before any calls to stackTraceGetString to build
	// the symbol table so it can be searched quickly.
	// This funtion mallocs a big chunk of memory for the tables
	// so it is best if you do it early instead of
	// waiting until fatal time

char *zStackTraceGetString( int size, unsigned int *ips, unsigned int *args );
	// Given a list of instruction pointers and arguments generted from
	// a call to stack trace, this function
	// parses the symbol tables inside of the exe file
	// in order to print a useful report.  However, it
	// only understands COFF format debugging.  You must enable
	// COFF in the Project Settings Linker Options.  You may,
	// if you choose, enable both COFF and MS formats.
	//   If COFF data is not available, it will print out using
	// a hex dump, and may print more information than is valid.
	// The traversal of the stack is predicated on assumptions
	// about non-optimized MSVC stack frames.  Thus, if optimizations
	// are on, results may not be predicatble.
	//   If the skipAssert flag is true, it will not dump
	// and functions with the word assert in them.  Useful
	// for ignoring the call to assert which is not helpful.

char *zStackTraceAndGetString();
	// Calls stackTrace first, then stackTraceGetString with the results

int zStackTraceOpenFiles();
void zStackTraceCloseFiles();
	// Use these two functions around stackTrace calls as an optimzation
	// if you are going to be doing a lot of them in a row to cut down on repeative file opens
	// If you don't call this, the stack trace will open and close automatically


#ifdef WIN32
char *zStackTraceDemangleMSDEVFuncName(
	char *mangledName, 
	char funcName[], int funcNameMax,
	char args[], int maxArgs, 
	int &numArgs
);
#endif
	// This function is used by the stackTrace function in
	// order to deceipher the mangled C++ function names
	// generated by MSDEV in order to determine the argument list.  
	// The parser for this is very complicated and unfortunately
	// may not be able to parse all names.  In these cases,
	// stackTrace should resort to printing a hex dump.

#endif