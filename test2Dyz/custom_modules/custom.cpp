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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "./custom.h"

// macros for customizable cell parameters and scenarios
// parameters for growth: BT474
#define BIRTH_RATE 0.00025    //  BR=0.00025/ DR=0.000002: results in 5-fold increase in cell numbers
#define DEATH_RATE 0.000002  //  BR=0.00025/ DR=0.000002: results in 5-fold increase in cell numbers

// parameters for therapy: BT474

#define BIRTH_RATE_TRMT 0.000002
#define DEATH_RATE_TRMT 0.0003

#define TUMORFIBR_RATE 2 //ratio of tumor / fibroblasts growth rates

#define BIRTH_RATE_FIBR 0.00005
#define DEATH_RATE_FIBR 0.000002

#define BIRTH_RATE_FIBR_TRMT 0.01
#define DEATH_RATE_FIBR_TRMT 0.002


// copy_number_alteration_model or cna_fitness_model
#define CNA_MODEL cna_fitness_model
// other approach is to set this function to empty_function
// #define CNA_MODEL empty_cna_function;

#define CNA_PROB 0


// declare cell definitions here

Cell_Definition motile_cell;
Cell_Definition fibroblast_cell;
Cell_Definition sensitive_cell_trmt;
//Cell_Definition fibroblast_cell_trmt;

// Yannis: here are my definitions for fibroblasts
/*

fibroblast_cell = cell_defaults;
fibroblast_cell.name = "fibroblasts"
fibroblast_cell.type = 2;
//*/


//

Cycle_Model birth_death;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs


	// housekeeping
	create_birthdeath_model();
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// Name the default cell type

	cell_defaults.type = 0;
	cell_defaults.name = "tumor cell";

	// End adding new live phase info for cell death/removal
	cell_defaults.functions.cycle_model = birth_death; // NOTE: Uses stochastic splitting according to birth death process
	cell_defaults.genotype.genotype_id="0";
	cell_defaults.genotype.alter_genotype = CNA_MODEL;
	cell_defaults.genotype.cna_prob = CNA_PROB;

	// set default_cell_functions;

	//yannis: will probably need
	// cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.update_phenotype = NULL;
	// cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_therapy;
	cell_defaults.functions.volume_update_function = empty_function;

	// only needed for a 2-D simulation:
	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;
	//*/

	// make sure the defaults are self-consistent.
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions );

	// Set cell size to have radius of 10 - note radius is updated based on volume
	double radius = 10.0;
	cell_defaults.phenotype.geometry.radius = radius;
	static double four_thirds_pi = 4.188790204786391;
	cell_defaults.phenotype.volume.total = pow(radius, 3) * four_thirds_pi;

	// first find index for a few key variables.
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );

	// initially no necrosis
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0;

	// set oxygen uptake / secretion parameters for the default cell type
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10;
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38;

// soluble factor secretion decay_rates

cell_defaults.phenotype.secretion.saturation_densities[1] = 1;  //yannis added this for

	// add custom data here, if any



	// Now, let's define another cell type.
	// It's best to just copy the default and modify it.

	// make this cell type randomly motile, less adhesive, greater survival,
	// and less proliferative

	motile_cell = cell_defaults;
	motile_cell.type = 1;
	motile_cell.name = "tumorcell";

	// make sure the new cell type has its own reference phenotyhpe

	motile_cell.parameters.pReference_live_phenotype = &( motile_cell.phenotype );

	// enable random motility
	motile_cell.phenotype.motility.is_motile = false;//true;
	motile_cell.phenotype.motility.persistence_time = 15.0; // 15 minutes
	motile_cell.phenotype.motility.migration_speed = 0.25; // 0.25 micron/minute
	motile_cell.phenotype.motility.migration_bias = 0.0;// completely random

	// Set cell-cell adhesion to 5% of other cells
	motile_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 0.05;

	// Set apoptosis to zero
	motile_cell.phenotype.death.rates[apoptosis_model_index] = 0.0;

// define birth and death for cells
motile_cell.phenotype.cycle.data.transition_rate(0,0) = BIRTH_RATE;
motile_cell.phenotype.cycle.data.transition_rate(0, 1) = DEATH_RATE;

// should this rule apply for motile cells?
// motile_cell.functions.update_phenotype = tumor_fibroblast_cell_rule;

sensitive_cell_trmt = motile_cell;
sensitive_cell_trmt.type = 21;
sensitive_cell_trmt.name = "sensitive cell under treatment";
sensitive_cell_trmt.phenotype.cycle.data.transition_rate(0,0) = BIRTH_RATE_TRMT;
sensitive_cell_trmt.phenotype.cycle.data.transition_rate(0, 1) = DEATH_RATE_TRMT;

// change the function below to allow interactions
sensitive_cell_trmt.functions.update_phenotype = NULL;
sensitive_cell_trmt.functions.update_phenotype = tumor_fibroblast_cell_rule_global;


	fibroblast_cell = cell_defaults;
	fibroblast_cell.type = 2;
	fibroblast_cell.name = "fibroblast";

// define birth and death for cells
	fibroblast_cell.phenotype.cycle.data.transition_rate(0,0) = BIRTH_RATE_FIBR;
	fibroblast_cell.phenotype.cycle.data.transition_rate(0, 1) = DEATH_RATE_FIBR;

// define fibroblast secretion of solubble factor
fibroblast_cell.phenotype.secretion.secretion_rates[1] = 10;

//  fibroblast_cell_trmt = motile_cell;
//	fibroblast_cell_trmt.type = 22;
//	fibroblast_cell_trmt.name = "fibroblast cell under treatment";
//	fibroblast_cell_trmt.phenotype.cycle.data.transition_rate(0,0) = BIRTH_RATE_fibr_TRMT;
//	fibroblast_cell_trmt.phenotype.cycle.data.transition_rate(0, 1) = DEATH_RATE_fibr_TRMT;
	// make sure the new cell type has its own reference phenotyhpe

	fibroblast_cell.parameters.pReference_live_phenotype = &( fibroblast_cell.phenotype );

	// enable random motility
	//fibroblast_cell.phenotype.motility.is_motile = false;//true;
	//fibroblast_cell.phenotype.motility.persistence_time = 15.0; // 15 minutes
	//fibroblast_cell.phenotype.motility.migration_speed = 0.25; // 0.25 micron/minute
	//fibroblast_cell.phenotype.motility.migration_bias = 0.0;// completely random

	// Set cell-cell adhesion to 5% of other cells
	fibroblast_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 0.05;

	// Set apoptosis to zero
	fibroblast_cell.phenotype.death.rates[apoptosis_model_index] = 0.0;

	// NOTE: To alter the transition rate of a specific types
	//motile_cell.phenotype.cycle.data.transition_rate(0,0) = 1.0;

	return;
}

void setup_microenvironment( void )
{
	// set domain parameters

	default_microenvironment_options.X_range = {-300, 300};
	default_microenvironment_options.Y_range = {-300, 300};
	default_microenvironment_options.Z_range = {-30, 30};
	default_microenvironment_options.simulate_2D = true; // 3D!

	// no gradients need for this example

	default_microenvironment_options.calculate_gradients = true;

	// add fibroblast secreted factor

	microenvironment.add_density( "fibroblast_secreted_factor", "dimensionless" );
	microenvironment.diffusion_coefficients[1] = 1e3;
	microenvironment.decay_rates[1] = .1;

	// let BioFVM use oxygen as the default

	default_microenvironment_options.use_oxygen_as_first_field = true;

	// set Dirichlet conditions

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	// define specific conditions
	default_microenvironment_options.Dirichlet_condition_vector[0] = 38; // physioxic conditions
	default_microenvironment_options.Dirichlet_condition_vector[1] = 0;

	default_microenvironment_options.Dirichlet_activation_vector[1] = false;  // no Dirichlet for the fibroblast_secreted_factor

	// if there are more substrates, resize accordingly
	// Yannis commented out this
// std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
//	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

//void introduce_therapy( void )
//{

//	motile_cell.phenotype.death.rates[apoptosis_model_index] = 1.0;

//	double BIRTH_RATE = 0.001;
//	double DEATH_RATE = 0.1;

	// idea: we'll introduce therapy by changing birth and death rate of cells

	//double worker_fraction = 0.10;
//	double worker_fraction = 0; //yannis changed to 0 from 0.10 (see above)
//	int number_of_injected_cells = 500;
//
//	double left_coordinate = 600.0;
//	double right_cooridnate = 700.0;
//
//	double bottom_coordinate = -700;
//	double top_coordinate = 700;
//
//	for( int i=0 ;i < number_of_injected_cells ; i++ )
//	{
//		std::vector<double> position = {0,0,0};
//		position[0] = left_coordinate + (right_cooridnate-left_coordinate)*UniformRandom();
//		position[1] = bottom_coordinate + (top_coordinate-bottom_coordinate)*UniformRandom();
//
//		Cell* pCell;
//		if( UniformRandom() <= worker_fraction )
//		{ pCell = create_cell( worker_cell ); }
//		else
//		{ pCell = create_cell( cargo_cell ); }
//		pCell->assign_position( position );
//	}
//
//	return;
//}


void setup_tissue( void )
{
	// create some cells near the origin

	Cell* pC;

	pC = create_cell( motile_cell );
	pC->assign_position( 0.0, 0.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( -40.0, 0.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( -80.0, 0.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( -80.0, 40.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( -80.0, 80.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( 40.0, 00.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( 40.0, 40.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( 0.0, -40.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( -40.0, -40.0, 0.0 );
	pC->genotype.genotype_id="0";

	pC = create_cell( motile_cell );
	pC->assign_position( -80.0, -40.0, 0.0 );
	pC->genotype.genotype_id="0";


	pC = create_cell( fibroblast_cell );
	pC->assign_position( 150.0, -8.0, 0.0 );

	pC = create_cell( fibroblast_cell );
	pC->assign_position( -100.0, -20.0, 0.0 );

	pC = create_cell( fibroblast_cell );
	pC->assign_position( -150.0, -120.0, 0.0 );

  /* Commented out to initialize with single cell
	pC = create_cell();
	pC->assign_position( -2.0, -6.0, 1.0 );

	pC = create_cell();
	pC->assign_position( -5.0, 8.0, -7.0 );

	// now create a motile cell

	pC = create_cell( motile_cell );
	pC->assign_position( 5.0, -8.0, 3.0 );
  */

	return;
}

// function to import from a file
void setup_tissue( std::string input_file )
{
	// create some cells near the origin
	std::ifstream infile(input_file);
	Cell* pC;
	int iter = 1;

	while (infile)
  {
    std::string s;
    if (!getline( infile, s )) break;

    std::istringstream ss( s );
		int x_coord;
		int y_coord;
		int z_coord;
		std::string cell_type;

		ss >> x_coord;
		ss >> y_coord;
		ss >> z_coord;
		ss >> cell_type;
		if(cell_type == "tumorcell")
		{
			pC = create_cell( motile_cell );
		}
		else if(cell_type == "fibroblast")
		{
			pC = create_cell( fibroblast_cell );
		}
//		else
//		{
//			pC = create_cell( motile_cell);
//		}

		pC->assign_position( x_coord, y_coord, z_coord );
		pC->genotype.genotype_id="0"+ std::to_string(iter);
  }
	return;
}



std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with the Ki67 coloring

	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

// if the cell is a tumor cell and not dead, paint it red

if( pCell->type == 1 )
{
	 output[0] = "green"; // this is the cytoplasm
	 output[2] = "black"; // this is the nucles;
	 output[3] = "black"; // black nuclear outline color

}

	// if the cell is a fibroblast and not dead, paint it red

	if( pCell->type == 2 )
	{
		 output[0] = "blue";
		 output[2] = "black";
		 output[3] = "black"; // black nuclear outline color

	}

	return output;
}


void create_birthdeath_model( void )
{
	// Create a 2-phase model for live and "dead" cells
	// When live cells die, they transition to dead cells but are
	// removed in the process - trick to get around issues with
	// default exit upon division options for live cell types
	birth_death.code = PhysiCell_constants::live_cells_cycle_model;
	birth_death.name = "BirthDeath";

	birth_death.data.time_units = "min";

	birth_death.add_phase( PhysiCell_constants::live , "Live" );
	birth_death.add_phase(PhysiCell_constants::custom_phase , "Dead");

	birth_death.add_phase_link( 0 , 0 , NULL );
	birth_death.transition_rate(0, 0) = BIRTH_RATE;
	birth_death.add_phase_link( 0 , 1 , NULL );
	birth_death.transition_rate(0, 1) = DEATH_RATE;

	birth_death.phases[0].entry_function = standard_live_phase_entry_function;

	birth_death.phases[0].division_at_phase_exit = false;
	birth_death.phases[1].division_at_phase_exit = false;


	birth_death.phase_links[0][0].exit_function = phase_link_division;
	birth_death.phase_links[0][1].exit_function = phase_link_death;

	return;
}

void phase_link_death( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.flagged_for_removal = true;
	return;
}

void phase_link_division( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.flagged_for_division = true;
	return;
}


void copy_number_alteration_model( Cell* pCell, Genotype& genotype)
{
	if ( genotype.cna_prob > 0.0)
	{
		double cna_runif = UniformRandom();
		if (cna_runif < genotype.cna_prob)
		{
			genotype.copynumber_change();

			(*genotype.genome_file) <<
				PhysiCell_globals.current_time << "\t" <<
				genotype.genotype_id << "\t" <<
				pCell->phenotype.cycle.data.transition_rates[0][0] << "\t" <<
				pCell->phenotype.cycle.data.transition_rates[0][1] << std::endl;
		}
	}
	return;
}

// Create model with increased fitness for new mutations from exponential RNG
void cna_fitness_model(Cell* pCell, Genotype& genotype)
{
	if ( genotype.cna_prob > 0.0)
	{
		double cna_runif = UniformRandom();
		if (cna_runif < genotype.cna_prob)
		{
			// Add an additional fitness that is an exponential increase with rate 30
			double addl_fitness = generate_double_exponential_rv(50.0);
			pCell->phenotype.cycle.data.transition_rates[0][0] = fmax(pCell->phenotype.cycle.data.transition_rates[0][0]+addl_fitness, 0);

			genotype.copynumber_change();

			(*genotype.genome_file) <<
				PhysiCell_globals.current_time << "\t" <<
				genotype.genotype_id << "\t" <<
				pCell->phenotype.cycle.data.transition_rates[0][0] << "\t" <<
				pCell->phenotype.cycle.data.transition_rates[0][1] << std::endl;
		}
	}
	return;
}

double generate_double_exponential_rv(double rate)
{
	double fitness = -log(1 - UniformRandom()) / rate;
	if(UniformRandom() <= 0.5)
	{
		fitness *= -1;
	}
	return fitness;
}

// Create tumor_fibroblast rule //based on worker_cell_motility function in cancer_biorobbots
void tumor_fibroblast_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int signal_index = microenvironment.find_density_index( "fibroblast_secreted_factor" );
	static double detection_threshold = 0.3;  // i need to figure out this threshold

// if tumor cell is dead don't bother
if( phenotype.death.dead == true )
{
	// the cell death functions don't automatically turn off custom functions,
	// since those are part of mechanics.

	// Let's just fully disable now.
	pCell->functions.custom_cell_rule = NULL;
	return;
}

// find if concentration for tumor cell is at specific range
if( pCell->nearest_density_vector()[signal_index] > detection_threshold )
{
	//phenotype.motility.is_motile = false;
	//pCell->functions.update_migration_bias = NULL;
	phenotype.cycle.data.transition_rate(0,0) = BIRTH_RATE_FIBR_TRMT;
	phenotype.cycle.data.transition_rate(0, 1) = DEATH_RATE_FIBR_TRMT;
  return;
}
return;

}

// Create tumor_fibroblast rule //based on worker_cell_motility function in cancer_biorobbots
void tumor_fibroblast_cell_rule_global( Cell* pCell, Phenotype& phenotype, double dt )
{
// if tumor cell is dead don't bother
if( phenotype.death.dead == true )
{
	pCell->functions.custom_cell_rule = NULL;
	return;
}

int nfibroblasts=0;
// find the number of fibroblasts
for( int i=0; i < (*all_cells).size(); i++ )
{
	// switch all currently living cells to under treatment types
	if( (*all_cells)[i]->type_name == "fibroblast" )
	{
		nfibroblasts=nfibroblasts+1;
	}
}

//alpha_death-beta_death = rate of death without fibroblasts
// death_rate = alpha_death + (beta_death / (beta_death + N_fibroblasts))
//
// assume that for infinite fibroblast number the tumor cell death rate is 10times less
double Nfold_protect= 0.01 ; // this is the fold-reduction in tumor cell death rate is decreased when infinite fibroblasts are present
double Nfold_growth=130; //was 120 for many simulations to match BT474-AR22 experiments at 100nM lap
double alpha_death=Nfold_protect*DEATH_RATE_TRMT; // this is the death rate for infinite number of fibr
double beta_death=(Nfold_protect-1)*DEATH_RATE_TRMT;  //
double alpha_birth=Nfold_growth*BIRTH_RATE_TRMT; // this is the maximum  birth rate increase through fibroblasts
//double beta_death=(Nfold_protect-1)*DEATH_RATE_TRMT;  //

// change birth_rate only to get tumor cells to grow in the presence of therapy
phenotype.cycle.data.transition_rate(0,0) = BIRTH_RATE_TRMT + (alpha_birth*nfibroblasts/(nfibroblasts+1)) ; // BIRTH_RATE

//DEATH RATE was a function of fibroblasts, but to get growth, change birth-rate only - fix death rate
phenotype.cycle.data.transition_rate(0, 1) = DEATH_RATE_TRMT;
//phenotype.cycle.data.transition_rate(0, 1) = alpha_death + (beta_death/(beta_death+nfibroblasts));

return;

}
