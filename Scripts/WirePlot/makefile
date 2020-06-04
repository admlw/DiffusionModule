CPPFLAGS=-I $(BOOST_INC) \
				 -I $(CANVAS_ROOT_IO_INC) \
				 -I $(CANVAS_INC) \
				 -I $(CETLIB_INC) \
				 -I $(CETLIB_EXCEPT_INC) \
				 -I $(FHICLCPP_INC) \
				 -I $(GALLERY_INC) \
				 -I $(LARCOREOBJ_INC) \
				 -I $(LARDATAOBJ_INC) \
				 -I $(NUSIMDATA_INC) \
				 -I $(HEP_CONCURRENCY_INC)\
				 -I $(LARCOREALG_INC)\
				 -I $(ROOT_INC)

CXXFLAGS=-std=c++17 #-Wall #-pedantic #-Werror
CXX=g++
LDFLAGS=$$(root-config --libs) \
				-L $(CANVAS_ROOT_IO_LIB) -l canvas_root_io\
				-L $(CANVAS_LIB) -l canvas\
				-L $(CETLIB_LIB) -l cetlib \
				-L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
				-L $(GALLERY_LIB) -l gallery \
				-L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
				-L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
				-L $(LARDATAOBJ_LIB) -l lardataobj_Simulation -l lardataobj_RecoBase -l lardataobj_MCBase -l lardataobj_RawData -l lardataobj_OpticalDetectorData -l lardataobj_AnalysisBase \
				-L $(HEP_CONCURRENCY_LIB) -l hep_concurrency\
				-L $(LARCOREALG_LIB) -l larcorealg_Geometry -l larcorealg_GeoAlgo\

GetWaveform: GetWaveform.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

GetSimChan: GetSimChan.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

GetWire: GetWire.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

all: GetWaveform GetSimChan GetWire

clean:
	rm *.o GetWaveform
	rm *.o GetSimChan
	rm *.o GetWire

