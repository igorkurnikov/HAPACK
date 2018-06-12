CXX = g++
CXXFLAGS = -O3 -c -fpic -I/usr/include/python2.7 -I../tinyxml -DEXIT_SUCCESS=0

SWIG = swig
SWIGFLAGS = -c++ -python -shadow

LDFLAGS = -shared

_llpnps.so: buildworldni.o contworld.o mapio.o nernstplanksolver.o pbwithljsolver.o \
           pnpdebug.o pnpinterfaces.o pnps.o pnpsapp.o pnpstructs.o \
           pnputil.o poissonboltzmannsolver.o poissonnernstplancksolver.o poissonsolver.o \
           sissgausps.o llpnps_wrap.o tinyxml.o  tinyxmlerror.o  tinyxmlparser.o -lz
	$(CXX) $(LDFLAGS) $^  -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $? -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $? -o $@

%.o: ../tinyxml/%.cpp
	$(CXX) $(CXXFLAGS) $? -o $@

%_wrap.cxx: %.i
	$(SWIG) $(SWIGFLAGS) $<

clean:
	rm -f *.o llpnps_wrap.cxx *.pyc llpnps.py _llpnps.so

