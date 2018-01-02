override HOME := ./..
include $(HOME)/Makefile 

%.exe : %.cpp Makefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) $< -o $@ 
	