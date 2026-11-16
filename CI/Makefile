TEST_CASES = $(wildcard */)

all: .PHONY

.PHONY: $(TEST_CASES:/=)

.SILENT:
$(TEST_CASES:/=): 
	echo $@
	@cd $@;\
	rm -f NOT_PASSED PASSED;\
	mv ../../params.h ../../params_old.h;\
	mv ../../definitions.h ../../definitions_old.h;\
	cp params.h ../../params.h;\
	cp definitions.h ../../definitions.h
	@cd ../;\
	make cleandata > std.out;\
	echo "Compiling...";\
	make dslbm > std.out;\
	echo "Running...";\
	mpirun -np $(n) ./dslbm > std.out
	@cd $@;\
	./compare.sh;\
	mv ../../params_old.h ../../params.h;\
	mv ../../definitions_old.h ../../definitions.h