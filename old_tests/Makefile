TEST_CASES = $(wildcard */)

all: .PHONY

.PHONY: $(TEST_CASES:/=)

.SILENT:
$(TEST_CASES:/=): 
	echo $@
	cd ../;\
	mv params.h params_old.h;\
	mv definitions.h definitions_old.h
	cd $@;\
	cp params.h ../../params.h;\
	cp definitions.h ../../definitions.h
	cd ../;\
	make cleandata && make && make run_local n=$(n)
	cd $@;\
	python3 get_result.py
	cd ../;\
	mv params_old.h params.h;\
	mv definitions_old.h definitions.h;\
	sleep 2

all: $(TEST_CASES:/=)
