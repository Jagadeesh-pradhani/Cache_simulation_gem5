###############################################################################
# SPECtools2000
#   Version 0.3
#   This makefile should NOT be changed other than to customize for a 
#   particular benchmark.  Users may overide variables by placing them
#   in the Makefile.spec file.  Or by using the "runspec" program to 
#   build the Makefile.spec file for them.
###############################################################################
#	BENCHMARK	-- standard definitions for this benchmark
all:
	f95 -fopenmp -c -o fma3d.o                fma3d.f90
	f95 -fopenmp -c -o beam_.o                beam_.f90
	f95 -fopenmp -c -o include_file_.o                include_file_.f90
	f95 -fopenmp -c -o penta_.o                penta_.f90
	f95 -fopenmp -c -o segment_set_.o                segment_set_.f90
	f95 -fopenmp -c -o body_force_.o                body_force_.f90
	f95 -fopenmp -c -o indx_.o                indx_.f90
	f95 -fopenmp -c -o periodic_bc_.o                periodic_bc_.f90
	f95 -fopenmp -c -o sliding_interface_.o                sliding_interface_.f90
	f95 -fopenmp -c -o constrained_node_.o                constrained_node_.f90
	f95 -fopenmp -c -o layering_.o                layering_.f90
	f95 -fopenmp -c -o plate_pair_.o                plate_pair_.f90
	f95 -fopenmp -c -o sliding_node_.o                sliding_node_.f90
	f95 -fopenmp -c -o contact_node_.o                contact_node_.f90
	f95 -fopenmp -c -o location_.o                location_.f90
	f95 -fopenmp -c -o platq_.o                platq_.f90
	f95 -fopenmp -c -o spot_weld_.o                spot_weld_.f90
	f95 -fopenmp -c -o contact_surface_.o                contact_surface_.f90
	f95 -fopenmp -c -o lsold_.o                lsold_.f90
	f95 -fopenmp -c -o platt_.o                platt_.f90
	f95 -fopenmp -c -o spring_.o                spring_.f90
	f95 -fopenmp -c -o coord_.o                coord_.f90
	f95 -fopenmp -c -o massprop_.o                massprop_.f90
	f95 -fopenmp -c -o pressure_bc_.o                pressure_bc_.f90
	f95 -fopenmp -c -o spring_bc_.o                spring_bc_.f90
	f95 -fopenmp -c -o damper_.o                damper_.f90
	f95 -fopenmp -c -o material_.o                material_.f90
	f95 -fopenmp -c -o property_.o                property_.f90
	f95 -fopenmp -c -o state_variables_.o                state_variables_.f90
	f95 -fopenmp -c -o damper_bc_.o                damper_bc_.f90
	f95 -fopenmp -c -o mean_stress_.o                mean_stress_.f90
	f95 -fopenmp -c -o shared_common_data.o                shared_common_data.f90
	f95 -fopenmp -c -o stress_.o                stress_.f90
	f95 -fopenmp -c -o displacement_bc_.o                displacement_bc_.f90
	f95 -fopenmp -c -o membq_.o                membq_.f90
	f95 -fopenmp -c -o qa_record_.o                qa_record_.f90
	f95 -fopenmp -c -o tabulated_function_.o                tabulated_function_.f90
	f95 -fopenmp -c -o element_set_.o                element_set_.f90
	f95 -fopenmp -c -o membt_.o                membt_.f90
	f95 -fopenmp -c -o relink_scratch_.o                relink_scratch_.f90
	f95 -fopenmp -c -o tetra_.o                tetra_.f90
	f95 -fopenmp -c -o enumerated_sets_.o                enumerated_sets_.f90
	f95 -fopenmp -c -o motion_.o                motion_.f90
	f95 -fopenmp -c -o results_.o                results_.f90
	f95 -fopenmp -c -o tied_bc_.o                tied_bc_.f90
	f95 -fopenmp -c -o force_.o                force_.f90
	f95 -fopenmp -c -o nodal_point_mass_.o                nodal_point_mass_.f90
	f95 -fopenmp -c -o rigid_body_.o                rigid_body_.f90
	f95 -fopenmp -c -o truss_.o                truss_.f90
	f95 -fopenmp -c -o force_bc_.o                force_bc_.f90
	f95 -fopenmp -c -o node_.o                node_.f90
	f95 -fopenmp -c -o rigid_body_mass_.o                rigid_body_mass_.f90
	f95 -fopenmp -c -o value_.o                value_.f90
	f95 -fopenmp -c -o gauge1d_.o                gauge1d_.f90
	f95 -fopenmp -c -o node_set_.o                node_set_.f90
	f95 -fopenmp -c -o rigid_wall_bc_.o                rigid_wall_bc_.f90
	f95 -fopenmp -c -o velocity_ic_.o                velocity_ic_.f90
	f95 -fopenmp -c -o gauge2d_.o                gauge2d_.f90
	f95 -fopenmp -c -o nonreflecting_bc_.o                nonreflecting_bc_.f90
	f95 -fopenmp -c -o section_1d_.o                section_1d_.f90
	f95 -fopenmp -c -o gauge3d_.o                gauge3d_.f90
	f95 -fopenmp -c -o nrbc_data_.o                nrbc_data_.f90
	f95 -fopenmp -c -o section_2d_.o                section_2d_.f90
	f95 -fopenmp -c -o hexah_.o                hexah_.f90
	f95 -fopenmp -c -o output_.o                output_.f90
	f95 -fopenmp -c -o segment_.o                segment_.f90
	f95 -fopenmp -c -o lsold.o                lsold.f90
	f95 -fopenmp -c -o damper.o                damper.f90
	f95 -fopenmp -c -o spring.o                spring.f90
	f95 -fopenmp -c -o material_00.o                material_00.f90
	f95 -fopenmp -c -o material_10.o                material_10.f90
	f95 -fopenmp -c -o material_11.o                material_11.f90
	f95 -fopenmp -c -o material_17.o                material_17.f90
	f95 -fopenmp -c -o material_22.o                material_22.f90
	f95 -fopenmp -c -o material_25.o                material_25.f90
	f95 -fopenmp -c -o material_32.o                material_32.f90
	f95 -fopenmp -c -o material_33.o                material_33.f90
	f95 -fopenmp -c -o material_34a.o                material_34a.f90
	f95 -fopenmp -c -o material_36.o                material_36.f90
	f95 -fopenmp -c -o material_38.o                material_38.f90
	f95 -fopenmp -c -o material_dm.o                material_dm.f90
	f95 -fopenmp -c -o material_sp.o                material_sp.f90
	f95 -fopenmp -c -o sort.o                sort.f90
	f95 -fopenmp -c -o pdb.o                pdb.f90
	f95 -fopenmp -c -o beam.o                beam.f90
	f95 -fopenmp -c -o membq.o                membq.f90
	f95 -fopenmp -c -o membt.o                membt.f90
	f95 -fopenmp -c -o penta.o                penta.f90
	f95 -fopenmp -c -o tetra.o                tetra.f90
	f95 -fopenmp -c -o hexah.o                hexah.f90
	f95 -fopenmp -c -o platq.o                platq.f90
	f95 -fopenmp -c -o truss.o                truss.f90
	f95 -fopenmp -c -o platt.o                platt.f90
	f95 -fopenmp -c -o fma1.o                fma1.f90
	f95 -fopenmp -c -o getirv.o                getirv.f90
	f95 -fopenmp -c -o relink.o                relink.f90
	f95 -fopenmp -c -o output.o                output.f90
	f95 -fopenmp -c -o fma2.o                fma2.f90
	f95 -fopenmp -c -o partition.o                partition.f90
	f95 -fopenmp -c -o strain.o                strain.f90
	f95 -fopenmp -c -o slide.o                slide.f90
	f95 -fopenmp         fma3d.o beam_.o include_file_.o penta_.o segment_set_.o body_force_.o indx_.o periodic_bc_.o sliding_interface_.o constrained_node_.o layering_.o plate_pair_.o sliding_node_.o contact_node_.o location_.o platq_.o spot_weld_.o contact_surface_.o lsold_.o platt_.o spring_.o coord_.o massprop_.o pressure_bc_.o spring_bc_.o damper_.o material_.o property_.o state_variables_.o damper_bc_.o mean_stress_.o shared_common_data.o stress_.o displacement_bc_.o membq_.o qa_record_.o tabulated_function_.o element_set_.o membt_.o relink_scratch_.o tetra_.o enumerated_sets_.o motion_.o results_.o tied_bc_.o force_.o nodal_point_mass_.o rigid_body_.o truss_.o force_bc_.o node_.o rigid_body_mass_.o value_.o gauge1d_.o node_set_.o rigid_wall_bc_.o velocity_ic_.o gauge2d_.o nonreflecting_bc_.o section_1d_.o gauge3d_.o nrbc_data_.o section_2d_.o hexah_.o output_.o segment_.o lsold.o damper.o spring.o material_00.o material_10.o material_11.o material_17.o material_22.o material_25.o material_32.o material_33.o material_34a.o material_36.o material_38.o material_dm.o material_sp.o sort.o pdb.o beam.o membq.o membt.o penta.o tetra.o hexah.o platq.o truss.o platt.o fma1.o getirv.o relink.o output.o fma2.o partition.o strain.o slide.o     -o fma3d

run :
	./fma3d < fma3d.in	

