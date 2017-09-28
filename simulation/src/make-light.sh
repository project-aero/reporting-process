#!/usr/bin/env bash

cd $outputdir
../src/get-res.R
cd ..
lightname="${outputdir}-light"
zip -r $lightname ${outputdir}/ --exclude \*.rda
