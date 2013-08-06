IRTK_BUILD := ../../../build

TMP_DIR := ${IRTK_BUILD}/tmp
SHELL := /bin/bash

all:
	mkdir -p ${TMP_DIR}
	python setup.py build_ext --build-temp ${TMP_DIR} --build-lib ${IRTK_BUILD}/lib/irtk --pyrex-c-in-temp #--force
	mv ${IRTK_BUILD}/lib/irtk/*.so ${IRTK_BUILD}/lib/irtk/ext/
	mv ${IRTK_BUILD}/lib/irtk/ext/_irtk.so ${IRTK_BUILD}/lib/irtk/

