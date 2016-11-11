$(OBJDIR)/IntPtsFFT.o	:IntPtsFFT.f;		$(F77) -c $(FL4) $< -o $@
$(OBJDIR)/psVtkOutput.o	:$(NEKFFTSRC)/psVtkOutput.cpp;		mpic++ -I$(VTK_INC) -c $(cFL2) $< -o $@
