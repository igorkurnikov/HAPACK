CXX = g++
CXXFLAGS = -fopenmp -O3 -std=c++11 -c -fpic -I/usr/include/python3.6 -I../tinyxml -DEXIT_SUCCESS=0

SWIG = swig
SWIGFLAGS = -python -py3 -shadow -c++

LDFLAGS = -shared -fopenmp -std=c++11

_pnpsll.so: buildworldni.o contworld.o mapio.o nernstplanksolver.o pbwithljsolver.o \
           pnpdebug.o pnpinterfaces.o pnps.o pnpsapp.o pnpstructs.o \
           pnputil.o poissonboltzmannsolver.o poissonnernstplancksolver.o poissonsolver.o \
           sissgausps.o pnpsll_wrap.o tinyxml.o  tinyxmlerror.o  tinyxmlparser.o -lz
	$(CXX) $(LDFLAGS) $^  -o $@
	mkdir -p ../pnpsll
	cp _pnpsll.so pnpsll.py ../pnpsll

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $? -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $? -o $@

%.o: ../tinyxml/%.cpp
	$(CXX) $(CXXFLAGS) $? -o $@

%_wrap.cxx: %.i
	$(SWIG) $(SWIGFLAGS) $<

clean:
	rm -f *.o pnpsll_wrap.cxx *.pyc pnpsll.py _pnpsll.so ../pnpsll/_pnpsll.so ../pnpsll/pnpsll.py \
	    ../pnpsll/pnpsll.pyc

