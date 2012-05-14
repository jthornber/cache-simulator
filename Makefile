SOURCE=\
	main.cc

OBJECTS=$(subst .cc,.o,$(SOURCE))

CXXFLAGS=\
	--std=c++0x \
	-Wall \
	-g \
	-I.

.SUFFIXES: .o .cc

%.o: %.cc
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) -o $@ $<

cache_sim: $(OBJECTS)
	g++ $(CXXFLAGS) -o $@ $+