TEST_CASES = $(wildcard */)

all: .PHONY

.PHONY: $(TEST_CASES:/=) $(TEST_CASES:/=_cluster)

.SILENT:
$(TEST_CASES:/=): 
	echo $@
	cd $@;\
	cp params.h ../../params.h;\
	cp definitions.h ../../definitions.h
	cd ../;\
	make cleandata && make && make run_local n=$(n);\
	python3 tests/$@/get_result.py

$(TEST_CASES:/=_cluster): %_cluster: 
	echo $@
	cd $*;\
	cp params.h ../../params.h;\
	cp definitions.h ../../definitions.h
	cd ../;\
	make run name=$(name)
	cd $*;\
	cp get_result.py ../../runs/$(name)

cleandata:
	cd ../;\
	make cleandata
