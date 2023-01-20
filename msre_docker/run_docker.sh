#!/usr/bin/env bash
PWD=`pwd`

mountdir=$1
if [ "a$1" == "a" ]; then
	mountdir=notebooks
fi
if [ ! -d $mountdir ]; then
	mkdir $mountdir
fi

docker run -it -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -e JUPYTER_TOKEN=docker -v ${PWD}/${mountdir}:/home/usr/notebooks copenhagenatomics/msre:0.1.1
