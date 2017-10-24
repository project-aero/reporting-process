#!/usr/bin/env bash

lightname="${outputdir}-light"
zip -r $lightname ${outputdir}/ --exclude \*.rda
