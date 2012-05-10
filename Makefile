SOURCE=\
	main.cc

OBJECTS=$(subst .cc,.o,$(SOURCE))

.SUFFIXES: .o .cc

%.o: %.cc
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) -o $@ $<

pdf_approx: $(OBJECTS)
	g++ $(CXXFLAGS) -o $@ $+