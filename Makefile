all: sketo

clean:
	rm -f matrix_sketo

nproc:
	@awk -F ":" -f nproc.awk Machinefile

sketo:
	sketocxx matrix_sketo.cpp -o matrix_sketo

.PHONY: all clean nproc sketo

