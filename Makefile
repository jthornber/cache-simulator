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

pdf_approx: $(OBJECTS)
	g++ $(CXXFLAGS) -o $@ $+