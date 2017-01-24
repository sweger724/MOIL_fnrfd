REM *************************************************************
REM This batch file builds MOIL programs on win32 using the minGW
REM environment which includes mingw32-make.  g77 can also be 
REM installed as part of this package.
REM
REM 		http://www.mingw.org/download.shtml
REM
REM Templates for the required makefiles exist in this folder
REM with the suffix '.template'.  Copy these files to the 
REM makefile of the same name, omitting the .template suffix.
REM e.g. copy Makefile.template Makefile
REM
REM You may also use the batch file copyTemplates.bat to copy
REM all relevant templates to respective makefiles.    
REM 
REM Then edit the copied files (not the templates) to specify
REM the compiler you are using.  Entries for g77 and ifort   
REM exist for windows-based builds.  Note also the entries
REM near the top that define macros for MV and RM which should
REM be commented IN for windows, and OUT for linux etc.
REM 
REM See the readme file in this folder for links to further 
REM build documentation.
REM **************************************************************

REM 1
REM ************************** make everything with LENGTH-CHMIN.BLOCK using Makefile.chmin'
mingw32-make clean
copy /Y COMMON\LENGTH-CHMIN.BLOCK COMMON\LENGTH.BLOCK
mingw32-make everything -f Makefile.chmin

REM 2
REM ************************** make exe/sdel exe/sdp with LENGTH-SDEL.BLOCK using Makefile.paral'
REM *** ( parallel programs not built on win32 )
REM mingw32-make clean
REM copy /Y COMMON\LENGTH-SDEL.BLOCK COMMON\LENGTH.BLOCK
REM mingw32-make -f Makefile.sdel exe/sdel exe/sdp 

REM 3
REM ************************** make exe/dynapt exe/coupledDyna exe/dyna_prl with LENGTH_SAVE.BLOCK using Makefile.paral'
REM *** ( parallel programs not built on win32 )
REM mingw32-make clean
REM copy /Y COMMON\LENGTH_SAVE.BLOCK COMMON\LENGTH.BLOCK
REM mingw32-make -f Makefile.sdel exe/dynapt exe/coupledDyna exe/dyna_prl

REM 4
REM ************************** make with LENGTH-SECOND.BLOCK using Makefile.second'
mingw32-make clean
copy /Y COMMON\LENGTH-SECOND.BLOCK COMMON\LENGTH.BLOCK
mingw32-make -f Makefile.second

REM 5
REM ************************** make with LENGTH_SAVE.BLOCK using Makefile'
mingw32-make clean
copy /Y COMMON\LENGTH_SAVE.BLOCK COMMON\LENGTH.BLOCK
mingw32-make

REM 6
REM ************************** make sdelS with LENGTH-SDEL.BLOCK using Makefile'
mingw32-make clean
copy /Y COMMON\LENGTH-SDEL.BLOCK COMMON\LENGTH.BLOCK
mingw32-make exe/sdelS

REM 7
REM ************************** make exe/dynapress with LENGTH_SAVE_PRESS.BLOCK using Makefile'
mingw32-make clean
copy /Y COMMON\LENGTH_SAVE_PRESS.BLOCK COMMON\LENGTH.BLOCK
mingw32-make exe/dynapress

REM 8
echo '************************** make Directional Milestoning executables with LENGTH-DIM.BLOCK using Makefile'
mingw32-make clean
copy /Y  COMMON\LENGTH-DIM.BLOCK COMMON\LENGTH.BLOCK
mingw32-make exe/dim_prepare exe/dim_run

REM 9
echo '************************** make Directional Milestoning parallel executables with LENGTH-DIM.BLOCK using Makefile.paral'
echo ' ( parallel programs not built on win32 )'
REM mingw32-make clean
REM copy /Y  COMMON\LENGTH-DIM.BLOCK COMMON\LENGTH.BLOCK
REM mingw32-make -f Makefile.paral exe/dim_sample

REM 10
REM ************************* cleanup
copy /Y COMMON\LENGTH_SAVE.BLOCK COMMON\LENGTH.BLOCK
mingw32-make clean

REM 11 -- windows only
REM ************************************************
REM Copying files to moil/moil.exe for win32
move /Y exe\*.exe ..\moil.exe
