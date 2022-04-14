#!/bin/sh

dest="${TMP}\\ChemChartsGlysadePython"
df_dir="c:\\gareth\\packages\\data_functions"
src_dir="c:\\gareth\\src\\pylib"

echo dest is $dest

if [ -d $dest ]; then
	echo "removing current distribution"
	rm -r $dest
fi

mkdir $dest
mkdir $dest\\jobs
touch $dest\\jobs\\empty_file
echo "Copying Python"
xcopy.exe $df_dir\\Python39 $dest\\Python /s /e /i /q

echo "Copying binaries"
xcopy.exe $df_dir\\bin $dest\\bin /s /e /i /q

echo "Copying Pylib"
mkdir $dest\\pylib
xcopy.exe $src_dir\\ruse $dest\\pylib\\ruse /s /e /i /q
xcopy.exe $src_dir\\df $dest\\pylib\\df /s /e /i /q

