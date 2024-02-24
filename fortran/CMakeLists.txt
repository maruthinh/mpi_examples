cmake_minimum_required(VERSION 3.15)
project(mpi_examples C CXX Fortran)

set(openmpi_dir $ENV{HOME}/softwares/openmpi-4.1.2/openmpi-4.1.2)

set(MPI_C_COMPILER ${openmpi_dir}/bin/mpicc)
set(MPI_CXX_COMPILER ${openmpi_dir}/bin/mpicxx)
set(MPI_Fortran_COMPILER ${openmpi_dir}/bin/mpifort)

set(MPI_LIBRARIES $ENV{HOME}/softwares/openmpi-4.1.2/openmpi-4.1.2/lib/libmpi.so 
	$ENV{HOME}/softwares/openmpi-4.1.2/openmpi-4.1.2/lib/libmpi_mpifh.so)

add_library(MPI::MPI_Fortran IMPORTED INTERFACE)
set_property(TARGET MPI::MPI_Fortran PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_LIBRARIES})

add_executable(test_scatter.out scatterv4d_arr_with_mpi_data_type.f90)
target_include_directories(test_scatter.out PRIVATE ${openmpi_dir}/lib)
target_link_libraries(test_scatter.out PRIVATE MPI::MPI_Fortran)
