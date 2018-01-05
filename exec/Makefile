override HOME := ./..
include $(HOME)/Makefile

%.exe : %.cpp Makefile $(OBJS)
	rm -rf solution*.vtk
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) $< -o $@ 
	
