#!/bin/tcsh

#
# Simple script to cycle through all the sub-directories and compare 
# Makefile, *.f, and *.INC files in two different MPI distributions of MPI_CMFGEN.
#

setenv RED "\033[31m"
setenv GREEN "\033[32m"
setenv YELLOW "\033[33m"
setenv BLUE "\033[34m"
setenv RESET "\033[0m"

echo "${RED}" > /dev/tty
echo "  This program compares all Makefiles, Fortran files, and *.INC files" > /dev/tty
echo "  in 2 MPI_CMFGEN directory structures. If one directory arguments is supplied," > /dev/tty
echo "  the comparison directory is taken as the pwd. Two directory arguments may" > /dev/tty
echo "  also be supplied. Output is to Diff_sum. Diff_output is corrupted." > /dev/tty
echo " " > /dev/tty
echo "  IT IS NOT TO BE USED FOR COMPARING CMFGEN DISTRIBUTIONS ${RESET}" > /dev/tty
echo " " > /dev/tty

if ($1 == "")then
  echo " Need to supply at least one directory  argument" > /dev/tty
  exit
endif

rm -f Diff_sum

if ($2 == "")then
  set X1="."
  set X2=$1
else
  set X1=$1
  set X2=$2
endif

echo "${BLUE}  Current    directory is" $PWD > /dev/tty
echo "  Comparison directory is" $X2 > /dev/tty
echo "${RESET} " > /dev/tty

echo -n " Please enter any character to continue: " > /dev/tty
set jnk_char = $<
echo " "  > /dev/tty

echo " " > Diff_sum
echo "Current directory is:" >> Diff_sum
pwd >> Diff_sum
echo " " >> Diff_sum

$cmfdist/com/com_diff.sh $X1/com $X2/com
cat Diff_output >> Diff_sum
 
$cmfdist/com/main_diff.sh $X1/blas $X2/blas
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/lpack $X2/lpack
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/main $X2/main
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main $X2/new_main
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/mod_subs $X2/new_main/mod_subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/mpi_output $X2/mpi_output
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/gam_transport $X2/new_main/gam_transport
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs $X2/new_main/subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/auto $X2/new_main/subs/auto
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/chg $X2/new_main/subs/chg
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/two $X2/new_main/subs/two
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/non_therm $X2/new_main/subs/non_therm
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/newsubs $X2/newsubs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/plane $X2/plane
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/subs $X2/subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/subs/chg $X2/subs/chg
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/subs/two $X2/subs/two
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/stark $X2/stark
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/txt_files $X2/txt_files
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/tools $X2/tools
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/unix  $X2/unix
cat Diff_output >> Diff_sum
