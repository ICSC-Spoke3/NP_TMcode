include Makefile.in

ifeq ($(BUILDFORTRAN), yes)
	FORTRAN_OBJS=../src/cluster/clu.o ../src/cluster/edfb_clu.o ../src/inclusion/edfb_inclu.o ../src/inclusion/inclu.o ../src/sphere/edfb_sph.o ../src/sphere/sph.o ../src/trapping/frfme.o ../src/trapping/lffft.o
	FORTRAN_BINS=cluster/clu cluster/edfb_clu inclusion/edfb_inclu inclusion/inclu sphere/edfb_sph sphere/sph trapping/frfme trapping/lffft
else
	FORTRAN_OBJS=
	FORTRAN_BINS=
endif # BUILDFORTRAN
ifeq ($(LIBMODE), shared)
	NPTM_LIB=libnptm/libnptm.so
	CXXFLAGSLIB=-fPIC -DPIC
else
	NPTM_LIB=libnptm/libnptm.a
	CXXFLAGSLIB=-static
endif # LIBMODE
NP_DOC_SRCS=../doc/src/config.dox
NP_LIBNPTM_SRCS=../src/libnptm/algebraic.cpp ../src/libnptm/clu_subs.cpp ../src/libnptm/Commons.cpp ../src/libnptm/Configuration.cpp ../src/libnptm/cublas_calls.cpp ../src/libnptm/file_io.cpp ../src/libnptm/inclu_subs.cpp ../src/libnptm/logging.cpp ../src/libnptm/lapack_calls.cpp ../src/libnptm/magma_calls.cpp ../src/libnptm/outputs.cpp ../src/libnptm/Parsers.cpp ../src/libnptm/sph_subs.cpp ../src/libnptm/tfrfme.cpp ../src/libnptm/tra_subs.cpp ../src/libnptm/TransitionMatrix.cpp ../src/libnptm/utils.cpp
NP_LIBNPTM_OBJS=../src/libnptm/algebraic.o ../src/libnptm/clu_subs.o ../src/libnptm/Commons.o ../src/libnptm/Configuration.o ../src/libnptm/cublas_calls.o ../src/libnptm/file_io.o ../src/libnptm/inclu_subs.o ../src/libnptm/logging.o ../src/libnptm/lapack_calls.o ../src/libnptm/magma_calls.o ../src/libnptm/outputs.o ../src/libnptm/Parsers.o ../src/libnptm/sph_subs.o ../src/libnptm/tfrfme.o ../src/libnptm/tra_subs.o ../src/libnptm/TransitionMatrix.o ../src/libnptm/utils.o
NP_CLUSTER_SRCS=../src/cluster/np_cluster.cpp ../src/cluster/cluster.cpp
NP_CLUSTER_OBJS=../src/cluster/np_cluster.o ../src/cluster/cluster.o
NP_CLUSTER_BINS=cluster/np_cluster
NP_INCLUSION_SRCS=../src/inclusion/np_inclusion.cpp ../src/inclusion/inclusion.cpp
NP_INCLUSION_OBJS=../src/inclusion/np_inclusion.o ../src/inclusion/inclusion.o
NP_INCLUSION_BINS=inclusion/np_inclusion
NP_SPHERE_SRCS=../src/sphere/np_sphere.cpp ../src/sphere/sphere.cpp
NP_SPHERE_OBJS=../src/sphere/np_sphere.o ../src/sphere/sphere.o
NP_SPHERE_BINS=sphere/np_sphere
NP_TRAPPING_SRCS=../src/trapping/np_trapping.cpp ../src/trapping/cfrfme.cpp ../src/trapping/clffft.cpp
NP_TRAPPING_OBJS=../src/trapping/np_trapping.o ../src/trapping/cfrfme.o ../src/trapping/clffft.o
NP_TRAPPING_BINS=trapping/np_trapping
NP_TESTING_OBJS=../src/testing/test_cluster_outputs.o ../src/testing/test_inclusion_outputs.o ../src/testing/test_sphere_outputs.o ../src/testing/test_ParticleDescriptor.o ../src/testing/test_TEDF.o ../src/testing/test_TTMS.o
NP_TESTING_BINS=testing/test_cluster_outputs testing/test_inclusion_outputs testing/test_sphere_outputs testing/test_ParticleDescriptor testing/test_TEDF testing/test_TTMS

all: $(NPTM_LIB) $(FORTRAN_BINS) $(NP_CLUSTER_BINS) $(NP_INCLUSION_BINS) $(NP_SPHERE_BINS) $(NP_TRAPPING_BINS) $(NP_TESTING_BINS)

docs: $(NP_LIBNPTM_SRCS) $(NP_CLUSTER_SRCS) $(NP_DOC_SRCS) $(NP_INCLUSION_SRCS) $(NP_SPHERE_SRCS) $(NP_TRAPPING_SRCS)
	cd ../doc/src; doxygen config.dox

libnptm/libnptm.a: $(NP_LIBNPTM_OBJS)
	$(AR) $(ARFLAGS) $@ $(NP_LIBNPTM_OBJS)

libnptm/libnptm.so: $(NP_LIBNPTM_OBJS)
	$(CXX) $(CXXFLAGS) $(NP_LIBNPTM_OBJS) -shared -o $@ $(CXXLDFLAGS)
	chmod ugo-x $@

cluster/clu: ../src/cluster/clu.o
	$(FC) $(FCFLAGS) $? -o $@

cluster/edfb_clu: ../src/cluster/edfb_clu.o
	$(FC) $(FCFLAGS) $? -o $@

cluster/np_cluster: $(NPTM_LIB) $(NP_CLUSTER_OBJS)
	$(CXX) $(CXXFLAGS) $(NP_CLUSTER_OBJS) -o $@ $(CXXLDFLAGS)

inclusion/edfb_inclu: ../src/inclusion/edfb_inclu.o
	$(FC) $(FCFLAGS) $? -o $@

inclusion/inclu: ../src/inclusion/inclu.o
	$(FC) $(FCFLAGS) $? -o $@

inclusion/np_inclusion: $(NPTM_LIB) $(NP_INCLUSION_OBJS)
	$(CXX) $(CXXFLAGS) $(NP_INCLUSION_OBJS) -o $@ $(CXXLDFLAGS)

sphere/edfb_sph: ../src/sphere/edfb_sph.o
	$(FC) $(FCFLAGS) $? -o $@

sphere/sph: ../src/sphere/sph.o
	$(FC) $(FCFLAGS) $? -o $@

sphere/np_sphere: $(NPTM_LIB) $(NP_SPHERE_OBJS)
	$(CXX) $(CXXFLAGS) $(NP_SPHERE_OBJS) -o $@ $(CXXLDFLAGS)

trapping/frfme: ../src/trapping/frfme.o
	$(FC) $(FCFLAGS) $? -o $@

trapping/lffft: ../src/trapping/lffft.o
	$(FC) $(FCFLAGS) $? -o $@

trapping/np_trapping: $(NPTM_LIB) $(NP_TRAPPING_OBJS)
	$(CXX) $(CXXFLAGS) $(NP_TRAPPING_OBJS) -o $@ $(CXXLDFLAGS)

testing/test_ParticleDescriptor: $(NPTM_LIB) ../src/testing/test_ParticleDescriptor.o
	$(CXX) $(CXXFLAGS) ../src/testing/test_ParticleDescriptor.o -o $@ $(CXXLDFLAGS)

testing/test_cluster_outputs: $(NPTM_LIB) ../src/testing/test_cluster_outputs.o
	$(CXX) $(CXXFLAGS) ../src/testing/test_cluster_outputs.o -o $@ $(CXXLDFLAGS)

testing/test_inclusion_outputs: $(NPTM_LIB) ../src/testing/test_inclusion_outputs.o
	$(CXX) $(CXXFLAGS) ../src/testing/test_inclusion_outputs.o -o $@ $(CXXLDFLAGS)

testing/test_sphere_outputs: $(NPTM_LIB) ../src/testing/test_sphere_outputs.o
	$(CXX) $(CXXFLAGS) ../src/testing/test_sphere_outputs.o -o $@ $(CXXLDFLAGS)

testing/test_TEDF: $(NPTM_LIB) ../src/testing/test_TEDF.o
	$(CXX) $(CXXFLAGS) ../src/testing/test_TEDF.o -o $@ $(CXXLDFLAGS)

testing/test_TTMS: $(NPTM_LIB) ../src/testing/test_TTMS.o
	$(CXX) $(CXXFLAGS) ../src/testing/test_TTMS.o -o $@ $(CXXLDFLAGS)

../src/cluster/clu.o: ../src/cluster/clu.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/cluster/edfb_clu.o: ../src/cluster/edfb_clu.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/cluster/np_cluster.o: ../src/cluster/np_cluster.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/cluster/cluster.o: ../src/cluster/cluster.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/inclusion/edfb_inclu.o: ../src/inclusion/edfb_inclu.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/inclusion/inclu.o: ../src/inclusion/inclu.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/inclusion/np_inclusion.o: ../src/inclusion/np_inclusion.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/inclusion/inclusion.o: ../src/inclusion/inclusion.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/sphere/edfb_sph.o: ../src/sphere/edfb_sph.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/sphere/sph.o: ../src/sphere/sph.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/sphere/np_sphere.o: ../src/sphere/np_sphere.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/sphere/sphere.o: ../src/sphere/sphere.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/trapping/frfme.o: ../src/trapping/frfme.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/trapping/lffft.o: ../src/trapping/lffft.f
	$(FC) $(FCFLAGS) -c $? -o $@

../src/trapping/np_trapping.o: ../src/trapping/np_trapping.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/trapping/cfrfme.o: ../src/trapping/cfrfme.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/trapping/clffft.o: ../src/trapping/clffft.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

../src/libnptm/%.o: ../src/libnptm/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGSLIB) -c $< -o $@

../src/testing/%.o: ../src/testing/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean wipe

clean:
	rm -rf $(FORTRAN_OBJS) $(NP_LIBNPTM_OBJS) $(NP_CLUSTER_OBJS) $(NP_INCLUSION_OBJS) $(NP_SPHERE_OBJS) $(NP_TRAPPING_OBJS) $(NP_TESTING_OBJS)

wipe:
	rm -rf $(FORTRAN_OBJS) $(NP_LIBNPTM_OBJS) $(NP_CLUSTER_OBJS) $(NP_INCLUSION_OBJS) $(NP_SPHERE_OBJS) $(NP_TRAPPING_OBJS) $(NP_TESTING_OBJS)
	rm -rf $(FORTRAN_BINS) libnptm/libnptm.a libnptm/libnptm.so $(NP_CLUSTER_BINS) $(NP_INCLUSION_BINS) $(NP_SPHERE_BINS) $(NP_TRAPPING_BINS) $(NP_TESTING_BINS)
	rm -rf ../doc/build
