#!/bin/sh 
rm -f $3.ntf
$2 < $3 > $4
echo $1\#$3\#$4 > $3.ntf
