TEST_CASES = $(wildcard */)

all: .PHONY

.PHONY: $(TEST_CASES:/=)

.SILENT:
$(TEST_CASES:/=): 
	echo $@
	@cd $@;\
	rm -f NOT_PASSED PASSED;\
	cp params.h ../../params.h;\
	cp definitions.h ../../definitions.h
	@cd ../;\
	make cleandata > std.out;\
	echo "Compiling...";\
	make > std.out;\
	echo "Running...";\
	make run_local n=$(n) > std.out
	@cd $@;\
	./compare.sh
