#!/bin/sh

echo 'using "convert" from ImageMagick to do convert to PNG'
for INP in *.ps
do
    newname=`basename $INP .ps`
    convert -geometry 100% $INP $newname.png
    convert -rotate 90 $newname.png $newname.png
    echo "convert $INP to $newname.png completely"
done
echo "process ended, please check your graphical files"
