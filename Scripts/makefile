# Makefile for uboonecode v06_70_01

CPPFLAGS=-I $(BOOST_INC) \
		     -I $(CANVAS_ROOT_IO_INC) \
		     -I $(CANVAS_INC) \
		     -I $(CETLIB_INC) \
		     -I $(CETLIB_EXCEPT_INC) \
         -I $(FHICLCPP_INC) \
		     -I $(GALLERY_INC) \
		     -I $(LARCOREOBJ_INC) \
		     -I $(LARCOREALG_INC) \
		     -I $(LARDATAOBJ_INC) \
		     -I $(NUSIMDATA_INC) \
		     -I $(ROOT_INC)

CXXFLAGS=-std=c++17 #-Wall #-pedantic #-Werror
CXX=g++
LDFLAGS=$$(root-config --libs) \
		      -L $(CANVAS_ROOT_IO_LIB) -l canvas_root_io\
					-L $(CANVAS_LIB) -l canvas \
					-L $(CETLIB_LIB) -l cetlib \
          -L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
					-L $(GALLERY_LIB) -l gallery \
			    -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
			    -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
			    -L $(LARDATAOBJ_LIB) -l lardataobj_Simulation -l lardataobj_RecoBase -l lardataobj_MCBase -l lardataobj_RawData -l lardataobj_OpticalDetectorData -l lardataobj_AnalysisBase

PlotDiffusion: PlotDiffusion.cxx
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

PlotDiffusionMultiPlane: PlotDiffusionMultiPlane.cxx
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

all: PlotDiffusion PlotDiffusionMultiPlane

clean:
	rm *.o PlotDiffusion
	rm *.o PlotDiffusionMultiPlane
