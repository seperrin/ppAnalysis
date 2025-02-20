#!/bin/bash
##################################################
validateout=`dirname $0`
validatetime=`date`
validated="0";
error=0
if [ -z $validateout ]
then
    validateout="."
fi

cd $validateout;
validateworkdir=`pwd`;

echo "*******************************************************" >> stdout
echo "* Automatically generated validation script           *" >> stdout

echo "* Time:    $validatetime " >> stdout
echo "* Dir:     $validateout" >> stdout
echo "* Workdir: $validateworkdir" >> stdout
echo "* ----------------------------------------------------*" >> stdout
ls -la ./ >> stdout
echo "* ----------------------------------------------------*" >> stdout

##################################################

if [ ! -f stderr ] ; then
   error=1
   echo "* ########## Job not validated - no stderr  ###"  >> stdout
   echo "Error = $error"  >> stdout
fi
parArch=`grep -Ei "Cannot Build the PAR Archive" stderr`
segViol=`grep -Ei "Segmentation violation" stderr`
segFault=`grep -Ei "Segmentation fault" stderr`
glibcErr=`grep -Ei '\*\*\* glibc detected \*\*\*' stderr`

if [ "$parArch" != "" ] ; then
   error=1
   echo "* ########## Job not validated - PAR archive not built  ###"  >> stdout
   echo "$parArch"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$segViol" != "" ] ; then
   error=1
   echo "* ########## Job not validated - Segment. violation  ###"  >> stdout
   echo "$segViol"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$segFault" != "" ] ; then
   error=1
   echo "* ########## Job not validated - Segment. fault  ###"  >> stdout
   echo "$segFault"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$glibcErr" != "" ] ; then
   error=1
   echo "* ########## Job not validated - *** glibc detected ***  ###"  >> stdout
   echo "$glibcErr"  >> stdout
   echo "Error = $error"  >> stdout
fi
if ! [ -f EventStat_temp.root ] ; then
   error=1
   echo "Output file EventStat_temp.root not found. Job FAILED !" >> stdout
   echo "Output file EventStat_temp.root not found. Job FAILED !" >> stderr
fi
if ! [ -f AnalysisResults.root ] ; then
   error=1
   echo "Output file AnalysisResults.root not found. Job FAILED !" >> stdout
   echo "Output file AnalysisResults.root not found. Job FAILED !" >> stderr
fi
if ! [ -f muonGrid.root ] ; then
   error=1
   echo "Output file muonGrid.root not found. Job FAILED !" >> stdout
   echo "Output file muonGrid.root not found. Job FAILED !" >> stderr
fi
if ! [ -f outputs_valid ] ; then
   error=1
   echo "Output files were not validated by the analysis manager" >> stdout
   echo "Output files were not validated by the analysis manager" >> stderr
fi
if [ $error = 0 ] ; then
   echo "* ----------------   Job Validated  ------------------*" >> stdout
   echo "* === Logs std* will be deleted === "
   rm -f std*
fi
echo "* ----------------------------------------------------*"
echo "*******************************************************"
cd -
exit $error
