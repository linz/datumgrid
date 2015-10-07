cl /EHsc /Fe..\ms\datumgrid.exe grid.cpp datumgrid.cpp lineqn.cpp smplarry.cpp obseq.cpp grdintrp.cpp progress.cpp bltmatrx.cpp grdobseq.cpp readcpf.cpp symmatrx.cpp contrlpt.cpp
del *.obj
copy datumgrid.help ..\ms
