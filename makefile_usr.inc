$(OBJDIR)/IntPtsFFT.o	:IntPtsFFT.f;		$(F77) -c $(FL4) $< -o $@
$(OBJDIR)/psVtkOutput.o	:psVtkOutput.cpp;		mpic++ -I$(VTK_INC) -c $(cFL2) $< -o $@
