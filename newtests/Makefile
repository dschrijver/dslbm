TEST_CASES = $(wildcard */)

all: .PHONY

.PHONY: $(TEST_CASES:/=)

.SILENT:
$(TEST_CASES:/=): 
	echo $@
	cd $@;\
	cp params.h ../../params.h;\
	cp definitions.h ../../definitions.h
	cd ../;\
	make cleandata && make && make run_local n=$(n)
	cd $@;\
	python3 get_result.py
	cd ../;\
	sleep 2

all: $(TEST_CASES:/=)

cleandata:
	cd ../;\
	make cleandata
