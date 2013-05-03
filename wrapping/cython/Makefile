IRTK_BUILD := ../../build

TMP_DIR := ${IRTK_BUILD}/tmp
SHELL := /bin/bash

all:
	if [ -e ${TMP_DIR} ]; then rm -r ${TMP_DIR}; fi
	mkdir -p ${TMP_DIR}
	python generate_templates.py ${TMP_DIR}
	cp -r irtk ${IRTK_BUILD}/lib/
	python setup.py build_ext --build-temp ${TMP_DIR} --build-lib ${IRTK_BUILD}/lib/irtk --force --pyrex-c-in-temp 
	mv ${IRTK_BUILD}/lib/irtk/_template.so ${IRTK_BUILD}/lib/irtk/ext/
	mv ${IRTK_BUILD}/lib/irtk/_slic.so ${IRTK_BUILD}/lib/irtk/ext/
	#chmod +x examples/*

clean:
	rm -r ${IRTK_BUILD}/lib/irtk
