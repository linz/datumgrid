LIBFILES= \
    bltmatrx.cpp \
    contrlpt.cpp \
    grdintrp.cpp \
    grdobseq.cpp \
    grid.cpp \
    lineqn.cpp \
    obseq.cpp \
    progress.cpp \
    readcpf.cpp \
    smplarry.cpp \
    symmatrx.cpp \
#    token.cpp \
#    tokenstr.cpp

LIBO=$(LIBFILES:.cpp=.o)

CPPFLAGS=-g

all: ../datumgrd ../griddef

libdatumgrid.a: $(LIBO)
	ar -rc libdatumgrid.a $(LIBO)
	ranlib libdatumgrid.a


../datumgrd: datumgrd.o libdatumgrid.a
	g++ $(CPP_FLAGS) datumgrd.o libdatumgrid.a -lm -o ../datumgrd
../griddef: griddef.o libdatumgrid.a
	g++ $(CPP_FLAGS) griddef.o libdatumgrid.a -lm -o ../griddef

clean:
	rm -f $(LIBO) datumgrd.o  griddef.o libdatumgrd.a
# 
# classtst: classtst.o libdatumgrid.a
# 	g++ $(CPP_FLAGS) classtst.o libdatumgrid.a -lm -o classtst
# tokentst: tokentst.o libdatumgrid.a
# 	g++ $(CPP_FLAGS) tokentst.o libdatumgrid.a -lm -o tokentst
# satest: satest.o libdatumgrid.a
# 	g++ $(CPP_FLAGS) satest.o libdatumgrid.a -lm -o satest
# testle: testle.o libdatumgrid.a
# 	g++ $(CPP_FLAGS) testle.o libdatumgrid.a -lm -o testle
