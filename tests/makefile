VPATH=../repos/DAZZ_DB
PATH:=../repos/DAZZ_DB:${PATH}
export PATH

default: test-open
test-open: DBrm
	nim c -r -o:open.exe open
DBrm:
	cd ../repos/DAZZ_DB; make -j
	which DBrm
	which fasta2DB
