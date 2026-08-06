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
	make dslbm > std.out;\
	echo "Running...";\
	mpirun -np $(n) ./dslbm > std.out
	@cd $@;\
	./compare.sh
