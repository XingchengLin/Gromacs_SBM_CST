TDIR=..

for f in `find src`
do

  if [ -d $f ]
  then
    echo "Ignoring $f"
  else
    echo "Copying $f"
    cp $f $TDIR/$f
  fi

done
