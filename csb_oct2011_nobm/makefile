CILK = /opt/intel/composerxe-2011.3.174/compiler
INCADD = -I$(CILK)/include -I$(CILK)/examples/include
LIBADD = -L$(CILK)/lib/intel64

GCCOPT = -O2 -fno-rtti -fno-exceptions # -ftree-vectorize
INTELOPT = -O2 -fno-rtti -fno-exceptions -parallel -restrict -ipo -xHOST -no-prec-div -fno-inline-functions #-vec_report3
DEB = -g -DNOBM

seqsym: sym_spmv_test.cpp csbsym.cpp csbsym.h utility.h friends.h SSEspmv.o
	icpc -cilk-serialize $(INCADD) $(INTELOPT) -o seqsym sym_spmv_test.cpp SSEspmv.o

parsym: sym_spmv_test.cpp csbsym.cpp csbsym.h utility.h friends.h SSEspmv.o
	icpc $(INCADD) $(INTELOPT) -o parsym sym_spmv_test.cpp SSEspmv.o 

symanal: sym_spmv_test.cpp csbsym.cpp csbsym.h utility.h friends.h SSEspmv.o
	icpc -DSTATS $(INCADD) $(INTELOPT) -o symanal sym_spmv_test.cpp SSEspmv.o -lcilkutil

seqspmv: csb_spmv_test.cpp bicsb.cpp bicsb.h bmcsb.cpp bmcsb.h friends.h utility.h SSEspmv.o
	icpc -cilk-serialize $(INCADD) $(INTELOPT) -o seqspmv csb_spmv_test.cpp SSEspmv.o

parspmv: csb_spmv_test.cpp bicsb.cpp bicsb.h bmcsb.cpp bmcsb.h friends.h utility.h SSEspmv.o 
	icpc $(INCADD) $(INTELOPT) -o parspmv csb_spmv_test.cpp SSEspmv.o

parspmv_nobm: csb_spmv_test.cpp bicsb.cpp bicsb.h friends.h utility.h
	icpc $(INCADD) $(INTELOPT) -DNOBM -o parspmv_nobm csb_spmv_test.cpp

parspmvt: csb_spmvt_test.cpp bicsb.cpp bicsb.h utility.h
	icpc $(INCADD) $(INTELOPT) -o parspmvt csb_spmvt_test.cpp

SSEspmv.o: SSEspmv.cpp
	icpc -DICC -fast -c SSEspmv.cpp	

clean:	
	rm -f seqspmv
	rm -f seqsym
	rm -f parspmv
	rm -f parsym 
	rm -f parspmvt
	rm -f parspmv_nobm
	rm -f *.o
