/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2017, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"
// put custom code modules here!
#include "./custom_modules/custom.h"
std::ofstream* Genotype::genome_file;

using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	bool XML_status = false;
	if( argc > 1 )
	{ XML_status = load_PhysiCell_config_file( argv[1] ); }
	else
	{ XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" ); }
	if( !XML_status )
	{ exit(-1); }

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);

	// PNRG setup
	SeedRandom( );

	// time setup
	std::string time_units = "min";

	/* Microenvironment setup */

	setup_microenvironment(); // modify this in the custom code

	/* PhysiCell setup */

	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 20;
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

	/* Users typically start modifying here. START USERMODS */

	create_cell_types();

 // setup_tissue(); //comment this out because using the code below
  setup_tissue("./input_cells_fibr_example2.txt");
  //  input file : input_cells...txt needs to be space detelimited (save with excel, mv command in terminal from .prn to .txt)


	/* Users typically stop modifying here. END USERMODS */

	double therapy_activation_time = 60 * 24 * 2; // start therapy at day 1 or 5
	double therapy_deactivation_time = 60 * 18 * 1; // start therapy at hour 18


	// set MultiCellDS save options

	set_save_biofvm_mesh_as_matlab( true ); //change this to save microenvironment
	set_save_biofvm_data_as_matlab( true );
	set_save_biofvm_cell_data( true );
	set_save_biofvm_cell_data_as_custom_matlab( true );

	// set the performance timers
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();

	// save a simulation snapshot

	char dfilename[1024];
	sprintf( dfilename , "%s/initial" , PhysiCell_settings.folder.c_str() );
	save_PhysiCell_to_MultiCellDS_xml_pugi( dfilename , microenvironment , PhysiCell_globals.current_time );

	// save a quick SVG cross section through z = 0, after setting its
	// length bar to 200 microns

	PhysiCell_SVG_options.length_bar = 200;

	// for simplicity, set a pathology coloring function
// Yannis: here define the coloring function in the custom.ccp
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
	//std::vector<std::string> (*cell_coloring_function)(Cell*) = simple_cell_coloring;

	sprintf( dfilename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() );
	SVG_plot( dfilename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );


	// create of simulation report file
	std::ofstream report_file;
	char filename[1024];
	sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() );
	report_file.open(filename); 	// create the data log file
	report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;

	// create clone report file
	std::ofstream clone_file;
	char clone_filename[1024];
	sprintf( clone_filename, "%s/clone_report.txt" , PhysiCell_settings.folder.c_str() );
	clone_file.open(clone_filename);
	clone_file<<"birth time\tgenotype\tbirth rate\tdeath rate"<<std::endl;
	// ancestor clone
	Genotype::genome_file = &clone_file;

	// main loop
	try
	{
		// Begin by outputting ancestor clone information to genome_file
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			*(Genotype::genome_file) << PhysiCell_globals.current_time << "\t" <<
			(*all_cells)[i]->genotype.genotype_id << "\t" <<
			(*all_cells)[i]->phenotype.cycle.data.transition_rates[0][0] << "\t" <<
			(*all_cells)[i]->phenotype.cycle.data.transition_rates[0][1] << std::endl;
		}

		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
//				// include treatment at time therapy_activation_time
			if( fabs(PhysiCell_globals.current_time - therapy_activation_time)<= 0.01*diffusion_dt )
	 {
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			// switch all currently living cells to under treatment types
			if( (*all_cells)[i]->type_name == "tumorcell" )
			{
				(*all_cells)[i]->convert_to_cell_definition(sensitive_cell_trmt);
			}
	//		else if( (*all_cells)[i]->type_name == "resistant" )
	//		{
	//			(*all_cells)[i]->convert_to_cell_definition(resistant_trmt);
			}
		}

		//if( fabs(PhysiCell_globals.current_time - therapy_deactivation_time)<= 0.01*diffusion_dt )
		//{
		//for( int i=0; i < (*all_cells).size(); i++ )
		//{
		// switch all currently living cells to under treatment types
		//if( (*all_cells)[i]->type_name == "sensitive cell under treatment" )
		//{
		//	(*all_cells)[i]->convert_to_cell_definition(motile_cell);
		//}
		//		else if( (*all_cells)[i]->type_name == "resistant" )
		//		{
		//			(*all_cells)[i]->convert_to_cell_definition(resistant_trmt);
		//}
		//}



			// save data if it's time.
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout );
				log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);

//				char filename[1024];
			//	sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );
			if( PhysiCell_settings.enable_full_saves == true )
			{
				sprintf( dfilename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );

				save_PhysiCell_to_MultiCellDS_xml_pugi( dfilename , microenvironment , PhysiCell_globals.current_time );
			}

				PhysiCell_globals.full_output_index++;
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}

			// save SVG files
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{
					sprintf( dfilename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index );

					SVG_plot( dfilename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

					PhysiCell_globals.SVG_output_index++;
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// update the microenvironment (enable this for microenviroment simulation)
			microenvironment.simulate_diffusion_decay( diffusion_dt );


			// run PhysiCell
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			PhysiCell_globals.current_time += diffusion_dt;

		}

		log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
		report_file.close();
		clone_file.close();
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}

	// save a final simulation snapshot

	// timer

	std::cout << std::endl << "Total simulation runtime: " << std::endl;
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() );
	std::cout << std::endl;

	return 0;
}
