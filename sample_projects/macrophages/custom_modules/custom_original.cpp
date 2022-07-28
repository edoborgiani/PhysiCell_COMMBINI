/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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
#include <string> 
#include <unistd.h>

#include "./custom.h"

// declare cell definitions here 

Cell_Definition macro0; 
Cell_Definition macro1;
Cell_Definition macro2;
Cell_Definition deb;
Cell_Definition neutro;
Cell_Definition MSC;

double ecf=1e12; //exponential correction factor due to small concentration values

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 

	// housekeeping 
	
	initialize_default_cell_definition();
	//create_standard_cell_death_models();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "cell"; 

	int debris_index = microenvironment.find_density_index( "debris" );
	int TNF_index = microenvironment.find_density_index( "TNFa" );
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" );
	int IFN_index = microenvironment.find_density_index( "IFNg" );

	// set default cell cycle model 
	
	cell_defaults.functions.cycle_model = live; 
	cell_defaults.phenotype.cycle.data.transition_rate(0,0) = 0.0; //parameters.doubles("macro_proliferation_ratio"); //proliferation ratio
	cell_defaults.phenotype.volume.fluid_change_rate=0.0;
	//cell_defaults.phenotype.molecular.internalized_total_substrates[debris_index]=1.0;

	// set default_cell_functions; 

	//cell_defaults.functions.update_phenotype = default_function; 
	
	//cell_defaults.functions.update_phenotype = macrophage0_function;
	//cell_defaults.functions.update_migration_bias = macrophage_chemotaxis; 

	// needed for a 2-D simulation: 
	
	/* grab code from heterogeneity */ 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 

	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	//cell_defaults.phenotype.sync_to_functions(macro0.functions ); 
	//cell_deafults.phenotype.secretion.secretion_rates[...] =0;

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int napoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Neutrophil Apoptosis" );
	//std::cout << "NAP" << " " << napoptosis_model_index << std::endl; 

	//std::cout<<"cytoindex: "<<cyto_index<<std::endl;

	//int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	//int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	// initially no death 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[napoptosis_model_index] = 0.0; 

	// initially no birth 
	//cell_defaults.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.0 ; 

	// not motile 
	cell_defaults.phenotype.motility.is_motile = true; 

	// set virus uptake / secretion parameters for the default (epithelial) cell type 
	//cell_defaults.phenotype.secretion.uptake_rates[virus_index] = 
	//	parameters.doubles("viral_internalization_rate"); 
	//cell_defaults.phenotype.secretion.secretion_rates[cyto_index] = 5; 
	//cell_defaults.phenotype.secretion.saturation_densities[virus_index] = 0; 

	// add custom data here, if any 
	
	// release virus at death or when eaten 
	cell_defaults.phenotype.molecular.fraction_released_at_death[debris_index]= 1.0; //100% of debris released into microenvironment
	//cell_defaults.phenotype.molecular.fraction_transferred_when_ingested[ virus_index ]= 1.0; 

	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it.

	macro0 = cell_defaults; 
	macro0.type = 1; 
	macro0.name = "macrophage M0";

	macro1 = cell_defaults; 
	macro1.type = 2; 
	macro1.name = "macrophage M1";

	macro2 = cell_defaults; 
	macro2.type = 3; 
	macro2.name = "macrophage M2";

	neutro = cell_defaults; 
	neutro.type = 4; 
	neutro.name = "neutrophils PMN";

	MSC = cell_defaults; 
	MSC.type = 5; 
	MSC.name = "MSC";

	//deb = cell_defaults;
	//deb.type = 4;
	//deb.name = "debris0";
	
	//M0 PROPERTIES

	macro0.functions.update_phenotype = macrophage0_function; 
	macro0.phenotype.sync_to_functions( cell_defaults.functions );
	double mvol=parameters.doubles("macro_volume");
	macro0.phenotype.volume.multiply_by_ratio(mvol/macro0.phenotype.volume.total); // cell volume = 1000 um^3
	//macro0.parameters.pReference_live_phenotype = &( macro0.phenotype ); 

	macro0.functions.cycle_model = live; 
	macro0.phenotype.cycle.data.transition_rate(0,0) = parameters.doubles("k_P0"); //proliferation ratio
	macro0.phenotype.volume.fluid_change_rate=0.0;
	macro0.phenotype.molecular.fraction_released_at_death[debris_index]= 1.0; //100% of debris released into microenvironment
	//cell_defaults.phenotype.molecular.internalized_total_substrates[debris_index]=1.0;

	//macro0.phenotype.cycle.data.transition_rate(0,0)=0.0;

	macro0.functions.update_migration_bias = macrophage_chemotaxis; 

	// enable random motility 
	macro0.phenotype.motility.is_motile = true; 
	macro0.phenotype.motility.persistence_time = 1.0; //parameters.doubles( "macrophage_persistence_time" );
	macro0.phenotype.motility.migration_speed = parameters.doubles( "m_ms" ); 
	//macro0.phenotype.motility.migration_bias = 0.5;
	
	macro0.phenotype.death.rates[apoptosis_model_index] = parameters.doubles("k_A0");
	//macrophage.phenotype.motility.migration_bias_direction={-1,0,0};

	// Set cell-cell adhesion relative to other cells 
	//macrophage.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	// no birth 
	//macrophage.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.05;
	
	// SECRETION + UPTAKE
	//macro0.phenotype.secretion.uptake_rates[debris_index] = 0.0;
	macro0.phenotype.secretion.secretion_rates[TNF_index] = parameters.doubles("k_T0"); // same as M1?
	macro0.phenotype.secretion.saturation_densities[TNF_index] = 1; // ng/um^3
	macro0.phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG0"); // same as M1?
	macro0.phenotype.secretion.saturation_densities[TGF_index] = 1; // ng/um^3
	macro0.phenotype.secretion.secretion_rates[IL10_index] = parameters.doubles("k_I0"); // same as M1?
	macro0.phenotype.secretion.saturation_densities[IL10_index] = 1; // ng/um^3

	macro0.phenotype.secretion.secretion_rates[IFN_index] = parameters.doubles("k_IF0"); // not sure
	macro0.phenotype.secretion.saturation_densities[IFN_index] = 1; // ng/um^3

	//	parameters.doubles("viral_internalization_rate"); 

	// MACROPHAGES M1
	macro1=macro0;

	macro1.functions.update_phenotype = macrophage1_function; 
	macro1.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	macro1.parameters.pReference_live_phenotype = &( macro1.phenotype ); 
	macro1.phenotype.secretion.secretion_rates[TNF_index] = parameters.doubles("k_T1");
	macro1.phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG1");
	macro1.phenotype.secretion.secretion_rates[IL10_index] = parameters.doubles("k_I1");
	macro1.phenotype.secretion.secretion_rates[IFN_index] = parameters.doubles("k_IF1");

	// MACROPHAGES M2
	macro2=macro0;

	macro2.functions.update_phenotype = macrophage2_function; 
	macro2.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	macro2.parameters.pReference_live_phenotype = &( macro2.phenotype ); 
	macro2.phenotype.secretion.secretion_rates[TNF_index] = parameters.doubles("k_T2");
	macro2.phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG2");
	macro2.phenotype.secretion.secretion_rates[IL10_index] = parameters.doubles("k_I2");

	//cell mechanics

/* 	macro0.phenotype.mechanics.cell_cell_adhesion_strength*=0.0;
	macro1.phenotype.mechanics.cell_cell_adhesion_strength*=0.0;
	macro2.phenotype.mechanics.cell_cell_adhesion_strength*=0.0; */

/* 	macro0.phenotype.mechanics.cell_cell_repulsion_strength*=3.0;
	macro1.phenotype.mechanics.cell_cell_repulsion_strength*=3.0;
	macro2.phenotype.mechanics.cell_cell_repulsion_strength*=3.0; */

	// NEUTROPHILS PMN

	neutro.functions.update_phenotype = neutrophils_function;
	neutro.phenotype.sync_to_functions( cell_defaults.functions );
	double nvol=parameters.doubles("neutro_volume");
	neutro.phenotype.volume.multiply_by_ratio(nvol/neutro.phenotype.volume.total); // cell volume = 600 um^3 (diam= 10 um)
	//macro0.parameters.pReference_live_phenotype = &( macro0.phenotype ); 

	neutro.functions.cycle_model = live; 
	neutro.phenotype.cycle.data.transition_rate(0,0) = 0.0; //proliferation ratio
	neutro.phenotype.volume.fluid_change_rate=0.0;
	neutro.phenotype.molecular.internalized_total_substrates[debris_index]=1.6e-3;
	neutro.phenotype.molecular.fraction_released_at_death[debris_index]= 1.0; //100% of debris released into microenvironment
	//cell_defaults.phenotype.molecular.internalized_total_substrates[debris_index]=1.0;

	//macro0.phenotype.cycle.data.transition_rate(0,0)=0.0;

	neutro.functions.update_migration_bias = neutrophils_chemotaxis;

	// enable random motility 
	neutro.phenotype.motility.is_motile = true; 
	neutro.phenotype.motility.persistence_time = 1.0; //parameters.doubles( "macrophage_persistence_time" );
	neutro.phenotype.motility.migration_speed = parameters.doubles( "n_ms" ); 
	//macro0.phenotype.motility.migration_bias = 0.5;
	
	neutro.phenotype.death.rates[apoptosis_model_index] = 0.0;
	neutro.phenotype.death.rates[napoptosis_model_index] = parameters.doubles("k_AN"); 
	//std::cout<< "test apoptosis " << std::endl;
	//neutro.phenotype.cycle.data.transition_rate(100,100) = 1.0;
	//macrophage.phenotype.motility.migration_bias_direction={-1,0,0};

	// Set cell-cell adhesion relative to other cells 
	//macrophage.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	// no birth 
	//macrophage.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.05;
	
	// SECRETION + UPTAKE
	//macro0.phenotype.secretion.uptake_rates[debris_index] = 0.0;
	neutro.phenotype.secretion.secretion_rates[TNF_index] = parameters.doubles("k_TN"); 
	neutro.phenotype.secretion.saturation_densities[TNF_index] = 1; // ng/um^3
	neutro.phenotype.secretion.secretion_rates[IFN_index] = parameters.doubles("k_IFN"); 
	neutro.phenotype.secretion.saturation_densities[IFN_index] = 1; // ng/um^3

	// MESENCHIMAL STEM CELLS MSC

	//MSC.functions.update_phenotype = mesenchimal_function;
	MSC.phenotype.sync_to_functions( cell_defaults.functions );
	double MSCvol=parameters.doubles("mesenchimal_volume");
	MSC.phenotype.volume.multiply_by_ratio(MSCvol/MSC.phenotype.volume.total); // cell volume = 600 um^3 (diam= 10 um)
	//macro0.parameters.pReference_live_phenotype = &( macro0.phenotype ); 

	MSC.functions.cycle_model = live; 
	MSC.phenotype.cycle.data.transition_rate(0,0) = parameters.doubles("mesenchimal_proliferation_ratio"); //proliferation ratio
	MSC.phenotype.volume.fluid_change_rate=0.0;
	MSC.phenotype.molecular.internalized_total_substrates[debris_index]=1.6e-3;
	MSC.phenotype.molecular.fraction_released_at_death[debris_index]= 0.0; //100% of debris released into microenvironment
	//cell_defaults.phenotype.molecular.internalized_total_substrates[debris_index]=1.0;

	//macro0.phenotype.cycle.data.transition_rate(0,0)=0.0;

	//MSC.functions.update_migration_bias = mesenchimal_chemotaxis;

	// enable random motility 
	MSC.phenotype.motility.is_motile = true; 
	MSC.phenotype.motility.persistence_time = 1.0; //parameters.doubles( "macrophage_persistence_time" );
	MSC.phenotype.motility.migration_speed = parameters.doubles( "mesenchimal_migration_speed" ); 
	//macro0.phenotype.motility.migration_bias = 0.5;
	
	MSC.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "mesenchimal_migration_speed" );
	//std::cout<< "test apoptosis " << std::endl;
	//neutro.phenotype.cycle.data.transition_rate(100,100) = 1.0;
	//macrophage.phenotype.motility.migration_bias_direction={-1,0,0};

	// Set cell-cell adhesion relative to other cells 
	//macrophage.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	// no birth 
	//macrophage.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.05;

	build_cell_definitions_maps(); 
	//display_cell_definitions( std::cout ); 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
// now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = true; 
	
	//set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 5 , 0.0); // 5% o2
	std::vector<bool> bc_activ(5, true); // USER
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = bc_vector; 

	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 

	//NUOVA PARTE
/* 	for( unsigned int k=0 ; k < default_microenvironment_options.mesh.z_coordinates.size() ; k++ )
	{
		int I = 10; 
		// set Dirichlet conditions along the xmin outer edges 
		for( unsigned int j=0 ; j < default_microenvironment_options.mesh.y_coordinates.size() ; j++ )
		{
			//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
			//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
			// set the value 
			default_microenvironment_options.add_dirichlet_node( default_microenvironment_options.voxel_index(I,j,k) , bc_vector );
			
			// set the activation 
			//microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(I,j,k) , 
			//default_microenvironment_options.Dirichlet_xmin );
			
		}
	} */
	
	// initialize BioFVM 
	
	initialize_microenvironment();

/* 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 150; i<200; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=145 ; j < 255 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,k) , bc_activ); 
				
			}
		}
	} 

 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 0; i<90; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=145 ; j < 255 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,k) , bc_activ); 
				
			}
		}
	}	

 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 170; i<300; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=255 ; j < 310 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,k) , bc_activ); 
				
			}
		}
	} 

 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 170; i<300; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=90 ; j < 145 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,k) , bc_activ); 
				
			}
		}
	}	

 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 0; i<130; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=255 ; j < 310 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,k) , bc_activ); 
				
			}
		}
	} 

 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 0; i<130; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=90 ; j < 145 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,k) , bc_activ); 
				
			}
		}
	} */

/*	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int j = 175; j<200; j++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int i=30 ; i < 90 ; i++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
			}
		}
	} 

	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for (int i = 90; i<145; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=150 ; j < 200 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				
			}
		}
	} 

	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int j = 150; j<200; j++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int i=90 ; i < 115 ; i++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				
			}
		}
	} 

	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 115; i<145; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j=1 ; j < 150 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
			}
		}
	}  */

/* 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int j = 100; j<=135; j++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int i= 50 ; i < 115 ; i++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
			}
		}
	}  */

/* 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 50; i<145; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j= 80 ; j < 135 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
			}
		}
	}  */

/* 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		int J = 80; 
		// set Dirichlet conditions along the xmin outer edges 
		for( unsigned int i= 50 ; i < 90 ; i++ )
		{
			//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
			//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
			// set the value 
			microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,J,k) , bc_vector );
			
		}
	} */

/* 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int i = 90; i<145; i++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int j= 1 ; j < 80 ; j++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
				
			}
		}
	}  */

/* 	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
	{
		for(int j = 1; j<25; j++)
		{
			// set Dirichlet conditions along the xmin outer edges 
			for( unsigned int i= 30 ; i < 90 ; i++ )
			{
				//std::cout<< "test " << microenvironment.mesh.y_coordinates.size() << std::endl;
				//std::cout<< "test " << microenvironment.voxel_index(I,j,k) << std::endl;
				// set the value 
				microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,k) , bc_vector );
			}
		}
	} */

	//microenvironment.name = "tissue"; 
	
	int debris_ID = microenvironment.find_density_index( "debris" ); 
	
	//microenvironment.diffusion_coefficients[debris_ID] = 0.0;
	//parameters.doubles("cyto_diffusion_coefficient");  
	//microenvironment.decay_rates[debris_ID] = 0;  
	
/* 	microenvironment.add_dirichlet_node(200,default_microenvironment_options.Dirichlet_xmax_values );
	microenvironment.set_substrate_dirichlet_activation(200,default_microenvironment_options.Dirichlet_xmax); 

	std::cout << "test " << microenvironment.voxel_index(450,450,0) << std::endl; */

	// set a new initial condition for debris

	double xDebris;
	double yDebris;

	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		// x coordinate of the nth voxel's center
		xDebris = microenvironment.mesh.voxels[n].center[0];
		yDebris = microenvironment.mesh.voxels[n].center[1];
		// access kth substrate of the nth voxel
		//if((xDebris<400.0 && yDebris<500.0 && yDebris>=350.0 && xDebris>=150.0)||(xDebris<150.0 && yDebris<750 && yDebris>=350 && xDebris>=-250)||(xDebris<150.0 && yDebris<-200 && yDebris>=-750 && xDebris>=-250)||(xDebris<-250 && yDebris<750 && yDebris>=-750 && xDebris>-450))
/* 		if((xDebris>-450.0 && yDebris>350 && yDebris<=(750.0-400.0*(((xDebris+450.0)*(xDebris+450.0))/(850.0*850.0))))||(xDebris<150.0 && yDebris<-200 && yDebris>=-750 && xDebris>=-250)||(xDebris<-250 && yDebris<350 && yDebris>=-750 && xDebris>-450))
		{
			microenvironment(n)[debris_ID] = 1e-4; //*(1-0.001*abs(400.0-xDebris)); //	particles/micron^3
		} */ // USER - 1/4 callus

		// USER - full callus 
 		if((yDebris>650.0 && yDebris<=(750.0-100.0*(((xDebris)*(xDebris))/(750.0*750.0))))||(xDebris<750.0 && yDebris<500.0 && yDebris>=0.0 && xDebris>-750.0)||(xDebris<375.0 && yDebris<650.0 && yDebris>=500.0 && xDebris>-375.0))
		{
			microenvironment(n)[debris_ID] = 1e-4; //*(1-0.001*abs(400.0-xDebris)); //	particles/micron^3
		}
 		if((yDebris<-650.0 && yDebris>=(-750.0+100.0*(((xDebris)*(xDebris))/(750.0*750.0))))||(xDebris<750.0 && yDebris>=-500.0 && yDebris<0.0 && xDebris>-750.0)||(xDebris<375.0 && yDebris>=-650.0 && yDebris<-500.0 && xDebris>-375.0))
		{
			microenvironment(n)[debris_ID] = 1e-4; //*(1-0.001*abs(400.0-xDebris)); //	particles/micron^3
		}
	}

	
	// display the microenvironment again 
	
	microenvironment.display_information( std::cout ); 
	
	return; 
}

void release_new_macrophages( void)
{
	std::cout << "--> Macrophages releasing";
	
	Cell* pC;

	double volume_total=4.6875e6*4; //um^3 volume within the region of interest - calculated by hand (remove *4 if 1/4 domain used)
	
	int cell_total=0;
 	for(int i=1; i<(*all_cells).size();i++)
	{
		if(((*all_cells)[i]->type==1)||((*all_cells)[i]->type==2)||((*all_cells)[i]->type==3))
		{
			cell_total++;
		}
	}
	double macro_concentration=cell_total/volume_total; // cell/um^3
	
	int debris_index = microenvironment.find_density_index( "debris" );
	double debris_particles=0;
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		debris_particles+=microenvironment(n)[debris_index]*1000; //particles, 1000 um^3 volume of 1 voxel
	}
	double debris_density=debris_particles/volume_total; // particles/um^3 

	//std::cout<< "Macrophages concentration = " << macro_concentration << std::endl;
	//std::cout<< "Debris density = " << debris_density << std::endl;

	double Mmax=parameters.doubles("M0_max"); // cell/micron^3
	double kmax=parameters.doubles("k_Rmax"); //h^-1  
	double RM=kmax*(1-(macro_concentration/Mmax))*debris_density*volume_total;

	double fraction_marrow=parameters.doubles("frac_marrow"); // -

	//std::cout<< "RMacro = " << RM << std::endl;

	int perimRelease_tissues;
	int perimRelease_marrow;
	int perimPosition;

	int released_M0=ceil(RM);

	for(int z=0;z<released_M0;z++)
	{
		perimRelease_tissues=750.0;
		perimRelease_marrow=500.0;
		//std::cout<< "released = " << released_M0 << std::endl;
		int zone=rand()%4;
		int side=rand()%100;
		if(side < (1-fraction_marrow)*100)
		{
			pC=create_cell(macro0);
			perimPosition=rand()%perimRelease_tissues;
			switch(zone){
			case 0:
				pC->assign_position(perimPosition,749.0-100.0*((perimPosition)*(perimPosition)/(750.0*750.0)),0.0);
				break;
			case 1:
				pC->assign_position(-perimPosition,749.0-100.0*((perimPosition)*(perimPosition)/(750.0*750.0)),0.0);
				break;
			case 2:
				pC->assign_position(-perimPosition,-749.0+100.0*((perimPosition)*(perimPosition)/(750.0*750.0)),0.0);
				break;
			case 3:
				pC->assign_position(perimPosition,-749.0+100.0*((perimPosition)*(perimPosition)/(750.0*750.0)),0.0);
				break;
			}
		}
		else 
		{
			pC=create_cell(macro0);
			perimPosition=rand()%perimRelease_marrow;
			switch(zone){
			case 0:
				pC->assign_position(749.0,-501.0+perimPosition,0.0);
				break;
			case 1:
				pC->assign_position(-749.0,-501.0+perimPosition,0.0);
				break;
			case 2:
				pC->assign_position(-749.0,501.0-perimPosition,0.0);
				break;
			case 3:
				pC->assign_position(749.0,501.0-perimPosition,0.0);
				break;
			}
		}
	}

	std::cout << " DONE" << std::endl;
	std::cout << released_M0 << " macrophages released! <--" << std::endl;

/* 	for(int z=0;z<released_M0;z++)
	{
		perimRelease=1797.0;
		perimPosition=rand()%perimRelease;
		//std::cout<< "released = " << released_M0 << std::endl;
		if(perimPosition < 599.0)
		{
			pC=create_cell(macro0);
			pC->assign_position(perimPosition-449.0,749.0,0.0);
		}
		else if(perimPosition < 850.0)
		{
			pC=create_cell(macro0);
			perimPosition -= 599.0;
			pC->assign_position(149.0,749.0-perimPosition,0.0);
		}
		else if(perimPosition < 1100.0)
		{
			pC=create_cell(macro0);
			perimPosition -= 850.0;
			pC->assign_position(149.0+perimPosition,499.0,0.0);
		}
		else if(perimPosition < 1249.0)
		{
			pC=create_cell(macro0);
			perimPosition -= 1100.0;
			pC->assign_position(399.0,499.0-perimPosition,0.0);
		}
		else 
		{
			pC=create_cell(macro0);
			perimPosition -= 1249.0;
			pC->assign_position(149.0,-201.0-perimPosition,0.0);
		}
	} */
	
	return; 
}

void release_new_neutrophils( void )
{
	std::cout << "--> Neutrophils releasing";

	Cell* pC;

	double volume_total=4.6875e6*4; //um^3 volume within the region of interest - calculated by hand (remove *4 if 1/4 domain used)

	int cell_total=0;
 	for(int i=1; i<(*all_cells).size();i++)
	{
		if(((*all_cells)[i]->type==4))
		{
			cell_total++;
		}
	}
	double neutro_concentration=cell_total/volume_total; // cell/um^3
	
	int debris_index = microenvironment.find_density_index( "debris" );
	double debris_particles=0;
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		debris_particles+=microenvironment(n)[debris_index]*1000; //particles, 1000 um^3 volume of 1 voxel
	}
	double debris_density=debris_particles/volume_total; // particles/um^3 

	//std::cout<< "Macrophages concentration = " << macro_concentration << std::endl;
	//std::cout<< "Debris density = " << debris_density << std::endl;

	double PMNmax=parameters.doubles("PMN_max"); // cell/micron^3
	double kmax=parameters.doubles("k_RNmax"); //h^-1 
	double RM=kmax*(1-(neutro_concentration/PMNmax))*debris_density*volume_total;

	//std::cout<< "Rneutro = " << RM << std::endl;

	int perimRelease;
	int perimPosition;

	int released_PMN=ceil(RM); 

	for(int z=0;z<released_PMN;z++)
	{
		perimRelease=498.0;
		perimPosition=rand()%perimRelease;
		int zone=rand()%4;
		pC=create_cell(neutro);
		switch (zone){
		case 0:
			pC->assign_position(749.0,-501.0+perimPosition,0.0);
			break;
		case 1:
			pC->assign_position(-749.0,-501.0+perimPosition,0.0);
			break;
		case 2:
			pC->assign_position(-749.0,501.0-perimPosition,0.0);
			break;
		case 3:
			pC->assign_position(749.0,501.0-perimPosition,0.0);
			break;
		}
	}

	std::cout << " DONE" << std::endl;
	std::cout << released_PMN << " neutrophils released! <--" << std::endl;

	return; 
}

void setup_tissue( void )
{

	std::cout << "--> Tissue Setup: " << std::endl;
	int debris_index = microenvironment.find_density_index( "debris" ); 
	// create some cells near the origin
	
	double length_x = microenvironment.mesh.bounding_box[3] - 
		microenvironment.mesh.bounding_box[0]; 
		
	double length_y = microenvironment.mesh.bounding_box[4] - 
		microenvironment.mesh.bounding_box[1]; 

	// create some cells near the origin
	
	Cell* pC;
	int initial_M0=parameters.doubles("initial_M0"); // cell/mm^2
	int initial_PMN=parameters.doubles("initial_PMN"); // cell/mm^2

	int num_res_M0=ceil(initial_M0/0.5); // this model is 0.5 mm^2
	int num_PMN=ceil(initial_PMN/0.5); // this model is 0.5 mm^2 

	/*int num_debris=200;

	for(int z=0;z<num_debris;z++)
	{
		bool debrisPositioned=false;

		float xDebris=NormalRandom(400,400);
		float yDebris=((float)rand()/(float)(RAND_MAX/1500.0))-750.0;
		
		//while(debrisPositioned==false)
		//{
			if((xDebris<400.0 && yDebris<500.0 && yDebris>350.0 && xDebris>=150.0)||(xDebris<150.0 && yDebris<750 && yDebris>350 && xDebris>=-250)||(xDebris<150.0 && yDebris<-200 && yDebris>-750 && xDebris>=-250)||(xDebris<-250 && yDebris<750 && yDebris>-750 && xDebris>-450))
			{
				pC=create_cell(deb);
				pC->set_total_volume(200.0);
				//pC->phenotype.volume.rupture_volume=1e20;
				pC->phenotype.molecular.internalized_total_substrates[nPCyto]=1.0;
				//std::cout<< "test " << xDebris << " " << yDebris << std::endl;
				pC->assign_position(xDebris,yDebris,0.0);
			//	debrisPositioned=true;
			}
			else
			{
				z--;
			}
		//}
	}*/

/* 	int perimRelease;
	int perimPosition;

	for(int z=0;z<num_extrav_M0;z++)
	{
		perimRelease=1398.0;
		perimPosition=rand()%perimRelease;
		//std::cout<< "released = " << released_M0 << std::endl;
		if(perimPosition < 850.0)
		{
			pC=create_cell(macro0);
			perimPosition-=450.0;
			pC->assign_position(perimPosition,749.0-400.0*((perimPosition+450.0)*(perimPosition+450.0)/(850.0*850.0)),0.0);
		}
		else 
		{
			pC=create_cell(macro0);
			perimPosition -= 849.0;
			pC->assign_position(149.0,-201.0-perimPosition,0.0);
		}
	} */

	std::cout << "-> Macrophages initialization";
	int macrox;
	int macroy;

	for(int z=0;z<num_res_M0;z++)
	{
		macrox=rand()%750;
		macroy=rand()%750;
		//if((neux>600 && neuy>1250) || (neux>600 && neuy<1100) || (neux>200 && neuy<1100 && neuy>550))
		if((macroy>750.0-100.0*((macrox*macrox)/(750.0*750.0))) || (macrox>750.0 && macroy<750.0) || (macrox>375.0 && macroy<650.0 && macroy>500.0))
		{
			z--;
		} 
		else 
		{
			pC=create_cell(macro0);
			int zone=rand()%4;
			switch(zone){
			case 0:
				pC->assign_position(macrox,macroy,0.0);
				break;
			case 1:
				pC->assign_position(-macrox,macroy,0.0);
				break;
			case 2:
				pC->assign_position(-macrox,-macroy,0.0);
				break;
			case 3:
				pC->assign_position(macrox,-macroy,0.0);
				break;
			}
			//pC->phenotype.death.models[0]->transition_rate(0,1) = 1.0;
		}
	}

	std::cout << " DONE" << std::endl;
	std::cout << num_res_M0 << " M0 macrophages released! <-" << std::endl;

/* 	for(int z=0;z<num_M0;z++)
	{
		perimRelease=1797;
		perimPosition=rand()%perimRelease;
		//std::cout<< "Perimeter = " << perimPosition << std::endl;
		if(perimPosition < 599)
		{
			pC=create_cell(macro0);
			pC->assign_position(perimPosition-449.0,749.0,0.0);
		}
		else if(perimPosition < 850)
		{
			pC=create_cell(macro0);
			perimPosition -= 599;
			pC->assign_position(149.0,749.0-perimPosition,0.0);
		}
		else if(perimPosition < 1100)
		{
			pC=create_cell(macro0);
			perimPosition -= 850;
			pC->assign_position(149.0+perimPosition,499.0,0.0);
		}
		else if(perimPosition < 1249)
		{
			pC=create_cell(macro0);
			perimPosition -= 1100;
			pC->assign_position(399.0,499.0-perimPosition,0.0);
		}
		else 
		{
			pC=create_cell(macro0);
			perimPosition -= 1249;
			pC->assign_position(149.0,-201.0-perimPosition,0.0);
		}
	} */

	std::cout << "-> Neutrophils initialization";
	int neux;
	int neuy;

	for(int z=0;z<num_PMN;z++)
	{
		neux=rand()%750;
		neuy=rand()%750;
		//if((neux>600 && neuy>1250) || (neux>600 && neuy<1100) || (neux>200 && neuy<1100 && neuy>550))
		if((neuy>750.0-100.0*((neux*neux)/(750.0*750.0))) || (neux>750.0 && neuy<750.0) || (neux>375.0 && neuy<650.0 && neuy>500.0))
		{
			z--;
		} 
		else 
		{
			pC=create_cell(neutro);
			int zone=rand()%4;
			switch(zone){
			case 0:
				pC->assign_position(neux,neuy,0.0);
				break;
			case 1:
				pC->assign_position(-neux,neuy,0.0);
				break;
			case 2:
				pC->assign_position(-neux,-neuy,0.0);
				break;
			case 3:
				pC->assign_position(neux,-neuy,0.0);
				break;
			}
			//pC->phenotype.death.models[0]->transition_rate(0,1) = 1.0;
		}
	}

	std::cout << " DONE" << std::endl;
	std::cout << num_PMN << " neutrophils released! <-" << std::endl;

/* 	int perimRelease;
	int perimPosition;

	for(int z=0;z<num_M0;z++)
	{
		pC=create_cell(macro0);
		
		perimRelease=1800;
		perimPosition=rand() % perimRelease;
		if(perimPosition < 600)
		{
			pC->assign_position(perimPosition-450,750.0,0.0);
		}
		else if(perimPosition < 850)
		{
			perimPosition -= 600;
			pC->assign_position(150.0,750.0-perimPosition,0.0);
		}
		else if(perimPosition < 1100)
		{
			perimPosition -= 850;
			pC->assign_position(150.0+perimPosition,500.0,0.0);
		}
		else if(perimPosition < 1250)
		{
			perimPosition -= 1100;
			pC->assign_position(400.0,500.0-perimPosition,0.0);
		}
		/* else if(perimPosition < 1800)
		{
			perimPosition -= 1250;
			pC->assign_position(400.0-perimPosition,350.0,0.0);
		}
		else if(perimPosition < 2350)
		{
			perimPosition -= 1800;
			pC->assign_position(-250.0,350.0-perimPosition,0.0);
		}
		else if(perimPosition < 2750)
		{
			perimPosition -= 2350;
			pC->assign_position(perimPosition-250.0,-200.0,0.0);
		} 
		else 
		{
			perimPosition -= 1250;
			pC->assign_position(150.0,-200.0-perimPosition,0.0);
		}
		//pC->assign_position(NormalRandom(0, 100),NormalRandom(0, 100),0.0);
		//pC->phenotype.molecular.internalized_total_substrates[nCyto]=1000;
	} */
	
	std::cout << " Tissue Setup COMPLETE <--" << std::endl;
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output =  false_cell_coloring_live_dead(pCell); 

	if(pCell->type == 2)
	{
		 output[0] = "red";  
	}

	if(pCell->type == 3)
	{
		 output[0] = "blue";  
	}

	output[1]=output[0];
		
	if(pCell->type == 4) // type==1 - test01
	{
		 output[0] = "yellow"; // "green" - test01
		 output[1] = "black"; // "green" - test01
	}
	
	return output; 
}

std::vector<std::string> immunefluorescence_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output =  false_cell_coloring_live_dead(pCell); 

	output[2]="rgb(49,89,242)";
	output[3]=output[2];

	if(pCell->type == 1)
	{
		 output[0] = "rgb(102,255,0)";  
	}

	if(pCell->type == 2)
	{
		 output[0] = "rgb(235,236,240)";  
	}

	if(pCell->type == 3)
	{
		 output[0] = "rgb(255,22,12)";  
	}
		
	if(pCell->type == 4)
	{
		 output[0] = "rgb(49,89,242)"; 
	}
	
	output[1]=output[0];
	return output; 
}


void macrophage0_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	//Microenvironment* p_density_vectors;
	// bookkeeping 
	
	//static int nPCyto = microenvironment.find_density_index( "cytokine" ); 
	/* int x->pCell.position(1);
	int y->pCell.position(2); */

	int debris_index = microenvironment.find_density_index("debris");	
	int TNF_index = microenvironment.find_density_index("TNFa");	
	int IL10_index = microenvironment.find_density_index("IL10");

	pCell->phenotype.molecular.internalized_total_substrates[debris_index]=1e-3; //5e-5;
	//pCell->phenotype.death.models[0]->transition_rate(0,1) = 1.0 / (8.6 * 60.0);
	
	// digest virus particles inside me 
	
	//static double implicit_Euler_constant = 
	//	(1.0 + dt * parameters.doubles("virus_digestion_rate") );
	//phenotype.molecular.internalized_total_substrates[nVirus] /= implicit_Euler_constant; 
	
	// check for contact with a cell
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// SURVIVAL

	int nsize=0;

	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		std::vector<double> displacement = pTestCell->position;
		displacement -= pCell->position;
		double distance = norm( displacement ); 
		
		double max_distance = pCell->phenotype.geometry.radius + pTestCell->phenotype.geometry.radius; 
		max_distance *= 1.1; 

		if(distance < max_distance)
			{
				nsize++; 
			} 
	}
	//std::cout << "survival " << nsize << std::endl;
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
 	pCell->phenotype.death.rates[apoptosis_model_index] = (nsize+1) * parameters.doubles("k_A0");
	pCell->phenotype.cycle.data.transition_rate(0,0)=parameters.doubles("k_P0")/(nsize+1);  
/* 	if(pCell->phenotype.cycle.data.transition_rate(0,0)>2e-4){
		std::cout << "fact: " << pCell->phenotype.cycle.data.transition_rate(0,0) << std::endl; 
	} */
	
	// INGEST CELL
	
 	//for( int n=0; n < pCell->cells_in_my_container().size() ; n++ )
	
	/*for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		// if it is not me and not a macrophage 
  		 if( pCell->nearest_density_vector()[nCyto] > threshold )
		{
			double fact=1+(pCell->nearest_density_vector()[nCyto]/100);
			//std::cout << "fact: " << fact << std::endl; 
			pCell->phenotype.volume.total*=fact;
			//std::cout << "volume: " << pCell->phenotype.volume.total << std::endl;
		}   
			// calculate distance to the cell 
			std::vector<double> displacement = pTestCell->position;
			displacement -= pCell->position;
			double distance = norm( displacement ); 
			
			double max_distance = pCell->phenotype.geometry.radius + 
				pTestCell->phenotype.geometry.radius; 
			max_distance *= 1.1; 
			
			// if it is not a macrophage, test for viral load 
			// if high viral load, eat it. 
		
			if( pTestCell->type == 4 && distance < max_distance )
			{
				//std::cout << "\t\tnom nom nom" << std::endl; 
				pCell->ingest_cell( pTestCell ); 
			} 
	}*/

	//std::cout << "test:" << microenvironment.mesh.dx << std::endl;

	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dx/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE

	double ke1=parameters.doubles("k_e1");
	double aed=parameters.doubles("a_ed");
	double RD=microenvironment(PositionFinal)[debris_index];
	RD/=(aed+microenvironment(PositionFinal)[debris_index]);

	pCell->phenotype.secretion.uptake_rates[debris_index] = ke1*RD;

	// PROLIFERATION

	int IFN_index = microenvironment.find_density_index( "IFNg" );
	pCell->phenotype.cycle.data.transition_rate(0,0)*=(parameters.doubles("k_0")+((1-parameters.doubles("k_0"))/(1.0+exp(ecf*(microenvironment(PositionFinal)[IFN_index]-parameters.doubles("a_0")))))); //proliferation ratio
/* 	if(pCell->phenotype.cycle.data.transition_rate(0,0)>2e-4){
		std::cout << "fact2: " << pCell->phenotype.cycle.data.transition_rate(0,0) << std::endl; 
	} */

	// CYTOKINE SECRETION

	//int TNF_index = microenvironment.find_density_index("TNFa");
	//DYNAMICS TO CHECK
	pCell->phenotype.secretion.secretion_rates[TNF_index]=parameters.doubles("k_T0");
	pCell->phenotype.secretion.secretion_rates[TNF_index]*=(1+(parameters.doubles("k_TI")/(1+1.*exp(ecf*(parameters.doubles("a_TI")-microenvironment(PositionFinal)[IFN_index])))));
	pCell->phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("k_IF0")*exp(-parameters.doubles("a_IF0")*microenvironment(PositionFinal)[TNF_index]);
	
	// DIFFERENTIATE INTO M1 OR M2

 	bool polarized=false;

	double k01=parameters.doubles("k_01");
	double a01=parameters.doubles("a_01");
	double g1=k01*microenvironment(PositionFinal)[TNF_index];
	g1/=(a01+microenvironment(PositionFinal)[TNF_index]);

	double gtest1=(double)rand() / RAND_MAX;

	//std::cout << "g01" << g1 << std::endl; //NEW
	//std::cout << "gtest " << gtest1 << std::endl; //NEW

	if(gtest1<=g1) 
	{
		pCell->type=2;
		polarized=true;
		pCell->functions.update_phenotype = macrophage1_function; 

		//std::cout << "M0 -> M1!" << std::endl; 
	}

	double k02=parameters.doubles("k_02");
	double a02=parameters.doubles("a_02");
	double g2=k02*microenvironment(PositionFinal)[IL10_index];
	g2/=(a02+microenvironment(PositionFinal)[IL10_index]);

	double gtest2=(double)rand() / RAND_MAX;

	//std::cout << "g02 " << g2 << std::endl; //NEW
	//std::cout << "gtest " << gtest2 << std::endl; //NEW

	if(gtest2<=g2 && polarized==false) 
	{
		pCell->type=3;
		pCell->functions.update_phenotype = macrophage2_function; 
		//std::cout << "M0 -> M2!" << std::endl; 
	} 

	boundaries(pCell);
	return; 
}

void boundaries(Cell* pCell)
{	
	Cell* pC = NULL; 

	std::vector<double> CellPosition = pCell->position;

	int debris_index = microenvironment.find_density_index("debris");

/*  	if (CellPosition[0]>-450 && CellPosition[0]<150 && CellPosition[1]>750)
	{
		//pCell->assign_position(CellPosition[0],1500-CellPosition[1],CellPosition[2]); 
		pCell->die();
	} 
  	else if (CellPosition[0]>150 && CellPosition[1]>500 && CellPosition[0]-150 <= CellPosition[1]-500)
	{
		//pCell->assign_position(300-CellPosition[0],CellPosition[1],CellPosition[2]); 
		pCell->die();
	}
 	else if (CellPosition[0]>150 && CellPosition[1]>500)
	{
		//pCell->assign_position(300-CellPosition[0],CellPosition[1],CellPosition[2]); 
		pCell->die();
	}
 	else if (CellPosition[0]>400 && CellPosition[1]>350 && CellPosition[1]<500)
	{
		//pCell->assign_position(800-CellPosition[0],CellPosition[1],CellPosition[2]); 
		pCell->die();
	} */
	if (CellPosition[1]>=650.0 && CellPosition[1]>750.0-100.0*((CellPosition[0])*(CellPosition[0])/(750.0*750.0)))
	{
		//pCell->assign_position(CellPosition[0],1500-CellPosition[1],CellPosition[2]); 
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
	if (CellPosition[1]<=-650.0 && CellPosition[1]<-750.0+100.0*((CellPosition[0])*(CellPosition[0])/(750.0*750.0)))
	{
		//pCell->assign_position(CellPosition[0],1500-CellPosition[1],CellPosition[2]); 
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
/*  	else if (CellPosition[0]>200.0 && CellPosition[1]<1100.0 && CellPosition[1]>-201 && CellPosition[0]+251 >= -(CellPosition[1]-349))
	{
		//pCell->assign_position(CellPosition[0],700-CellPosition[1],CellPosition[2]);
		pCell->die();
	}
 	else if (CellPosition[0]>-251 && CellPosition[1]<349 && CellPosition[1]>-201 && CellPosition[0]+251 >= (CellPosition[1]+201))
	{
		//pCell->assign_position(CellPosition[0],-400-CellPosition[1],CellPosition[2]);
		pCell->die();
	} */
 	else if (CellPosition[0]>=375.0 && CellPosition[1]<=650.0 && CellPosition[1]>=500.0)
	{
		//pCell->assign_position(-500.0-CellPosition[0],CellPosition[1],CellPosition[2]);
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
 	else if (CellPosition[0]<=-375.0 && CellPosition[1]<=650.0 && CellPosition[1]>=500.0)
	{
		//pCell->assign_position(-500.0-CellPosition[0],CellPosition[1],CellPosition[2]);
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
 	else if (CellPosition[0]>=375.0 && CellPosition[1]>=-650.0 && CellPosition[1]<=-500.0)
	{
		//pCell->assign_position(-500.0-CellPosition[0],CellPosition[1],CellPosition[2]);
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
 	else if (CellPosition[0]<=-375.0 && CellPosition[1]>=-650.0 && CellPosition[1]<=-500.0)
	{
		//pCell->assign_position(-500.0-CellPosition[0],CellPosition[1],CellPosition[2]);
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
 	else if (CellPosition[0]>=750.0 && CellPosition[1]<=-500.0 && CellPosition[1]>=500.0)
	{
		//pCell->assign_position(300-CellPosition[0],CellPosition[1],CellPosition[2]);
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
 	else if (CellPosition[0]<=-750.0 && CellPosition[1]<=-500.0 && CellPosition[1]>=500.0)
	{
		//pCell->assign_position(300-CellPosition[0],CellPosition[1],CellPosition[2]);
		pCell->phenotype.molecular.internalized_total_substrates[debris_index]=0.0;
		pCell->die();
	}
/*  	else if (CellPosition[1]<-750)
	{
		//pCell->assign_position(CellPosition[0],-1500-CellPosition[1],CellPosition[2]);
		pC=create_cell();
		pC->copy_data(pCell);	
		pC->copy_function_pointers(pCell);
		pC->assign_position(CellPosition[0],-1500-CellPosition[1],CellPosition[2]);
		pC->set_total_volume(pCell->phenotype.volume.total); 
		//set_total_volume(pCell->phenotype.volume.total); 
		pCell->die();
	}
 	else if (CellPosition[0]<-450.0)
	{
		//pCell->assign_position(-900-CellPosition[0],CellPosition[1],CellPosition[2]);
		pC=create_cell();
		pC->assign_position(-900.0-CellPosition[0],CellPosition[1],CellPosition[2]);		
		std::cout << "cell " << pCell->ID << " type " << pCell->type << " position " << CellPosition[0] << " " << CellPosition[1] << " "<< CellPosition[2] << std::endl;
		pC->ID=pCell->ID;
		pC->set_total_volume(pCell->phenotype.volume.total);
		pC->copy_data(pCell);	
		pC->copy_function_pointers(pCell);
		std::cout << "cell " << pC->ID << " type " << pC->type << " position " << -900.0-CellPosition[0] << " " << CellPosition[1] << " "<< CellPosition[2] << std::endl;
		std::vector<double> NewcellPosition = pC->position;
		std::cout << "cell " << pC->ID << " type " << pC->type << " position " << NewcellPosition[0] << " " << NewcellPosition[1] << " "<< NewcellPosition[2] << std::endl;
		pCell->die();		
		//sleep(10);
		//set_total_volume(pCell->phenotype.volume.total);
	}  */

 	/* if (CellPosition[0]<=-50)
	{
		pCell->assign_position(-19,CellPosition[1],CellPosition[2]);
	}  */
	
	return; 
}

void macrophage1_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	int debris_index = microenvironment.find_density_index("debris");	
	pCell->phenotype.molecular.internalized_total_substrates[debris_index]=1e-3; //5e-5;
	//pCell->phenotype.death.models[0]->transition_rate(0,1) = 1.0 / (8.6 * 60.0);
	//pCell->phenotype.volume.total=10000.0;
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// SURVIVAL

	int nsize=0;

	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		std::vector<double> displacement = pTestCell->position;
		displacement -= pCell->position;
		double distance = norm( displacement ); 
		
		double max_distance = pCell->phenotype.geometry.radius + pTestCell->phenotype.geometry.radius; 
		max_distance *= 1.1; 

		if(distance < max_distance)
			{
				nsize++; 
			} 
	}
	//std::cout << "survival " << nsize << std::endl;
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	pCell->phenotype.death.rates[apoptosis_model_index] = (nsize+1) * parameters.doubles("k_A1");
	pCell->phenotype.cycle.data.transition_rate(0,0)=parameters.doubles("k_P1")/(nsize+1);

	// INGEST DEBRIS - readapt
	
 	//for( int n=0; n < pCell->cells_in_my_container().size() ; n++ )
	/*for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		// if it is not me and not a macrophage 
/*  		 if( pCell->nearest_density_vector()[nCyto] > threshold )
		{
			double fact=1+(pCell->nearest_density_vector()[nCyto]/100);
			//std::cout << "fact: " << fact << std::endl; 
			pCell->phenotype.volume.total*=fact;
			//std::cout << "volume: " << pCell->phenotype.volume.total << std::endl;
		}   */
			// calculate distance to the cell 
			/*std::vector<double> displacement = pTestCell->position;
			displacement -= pCell->position;
			double distance = norm( displacement ); 
			
			double max_distance = pCell->phenotype.geometry.radius + 
				pTestCell->phenotype.geometry.radius; 
			max_distance *= 1.1; 
			
			// if it is not a macrophage, test for viral load 
			// if high viral load, eat it. 
		
			if( pTestCell->type == 4 && distance < max_distance )
			{
				//std::cout << "\t\tnom nom nom" << std::endl; 
				pCell->ingest_cell( pTestCell ); 
			} 
	}*/
	//boundaries(pCell);
	
	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dx/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE

	double ke1=parameters.doubles("k_e1");
	double aed=parameters.doubles("a_ed");
	double RD=microenvironment(PositionFinal)[debris_index];
	RD/=(aed+microenvironment(PositionFinal)[debris_index]);

	pCell->phenotype.secretion.uptake_rates[debris_index] = ke1*RD;

	// Cytokine secretion influence 
	
	int TNF_index = microenvironment.find_density_index("TNFa");
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" ); 
	int IFN_index = microenvironment.find_density_index( "IFNg" );
	//DYNAMICS TO CHECK
	pCell->phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("k_IF1")*exp(-parameters.doubles("a_IF1")*microenvironment(PositionFinal)[TNF_index]); // not sure!!!
	pCell->phenotype.secretion.secretion_rates[TNF_index]=parameters.doubles("k_T1")*(parameters.doubles("k_TIL")*exp(-parameters.doubles("a_TIL")*(microenvironment(PositionFinal)[IL10_index]))+parameters.doubles("b_TIL"));
	pCell->phenotype.secretion.secretion_rates[TNF_index]*=(parameters.doubles("k_TTGF")*exp(-parameters.doubles("a_TTGF")*(microenvironment(PositionFinal)[TGF_index]))+parameters.doubles("b_TTGF"));
	pCell->phenotype.secretion.secretion_rates[TNF_index]*=(1+(parameters.doubles("k_TI")/(1+1.*exp(ecf*(parameters.doubles("a_TI")-microenvironment(PositionFinal)[IFN_index])))));
	pCell->phenotype.secretion.secretion_rates[IL10_index]=parameters.doubles("k_I1")*(1+((parameters.doubles("a_I1")*microenvironment(PositionFinal)[TGF_index])/(1+(parameters.doubles("a_I1")*microenvironment(PositionFinal)[TGF_index]))));
	pCell->phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG1");

	//PROLIFERATION
	pCell->phenotype.cycle.data.transition_rate(0,0)*=(parameters.doubles("k_1")+((1-parameters.doubles("k_1"))/(1+exp(ecf*(microenvironment(PositionFinal)[IFN_index]-parameters.doubles("a_1"))))));; //proliferation ratio

	// DIFFERENTIATE INTO M2

/* 	double k01=parameters.doubles("k01");
	double a01=parameters.doubles("a01");
	double g1=k01*microenvironment(PositionFinal)[TNF_index];
	g1/=(a01+microenvironment(PositionFinal)[TNF_index]);

	double gtest1=(double)rand() / RAND_MAX;

	//std::cout << "value: " << g1*6 << std::endl;

	if(gtest1>g1*6) 
	{ */

	double k12=parameters.doubles("k_12");
	double a12=parameters.doubles("a_12");
	double g2=k12*microenvironment(PositionFinal)[IL10_index];
	g2/=(a12+microenvironment(PositionFinal)[IL10_index]);

	double gtest2=(double)rand() / RAND_MAX;

	//std::cout << "g12 " << g2 << std::endl; //NEW
	//std::cout << "gtest " << gtest2 << std::endl; //NEW

	//std::cout << "pre M1 -> M2!" << std::endl; 

	if(gtest2<=g2)
	{
		pCell->type=3;
		//std::cout << "M1 -> M2!" << std::endl; 
		pCell->functions.update_phenotype = macrophage2_function; 
	}

	//}

	boundaries(pCell);
	return;
}

void macrophage2_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	int debris_index = microenvironment.find_density_index("debris");	
	pCell->phenotype.molecular.internalized_total_substrates[debris_index]=1e-3; //5e-5;
	//pCell->phenotype.death.models[0]->transition_rate(0,1) = 1.0 / (8.6 * 60.0);
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// SURVIVAL

	int nsize=0;

	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		std::vector<double> displacement = pTestCell->position;
		displacement -= pCell->position;
		double distance = norm( displacement ); 
		
		double max_distance = pCell->phenotype.geometry.radius + pTestCell->phenotype.geometry.radius; 
		max_distance *= 1.1; 

		if(distance < max_distance)
			{
				nsize++; 
			} 
	}
	//std::cout << "survival " << nsize << std::endl;
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	pCell->phenotype.death.rates[apoptosis_model_index] = (nsize+1) * parameters.doubles("k_A2");
	pCell->phenotype.cycle.data.transition_rate(0,0) = parameters.doubles("k_P2")/(nsize+1);

	// INGEST DEBRIS - readapt
	
 	//for( int n=0; n < pCell->cells_in_my_container().size() ; n++ )
	/*for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		// if it is not me and not a macrophage 
/*  		 if( pCell->nearest_density_vector()[nCyto] > threshold )
		{
			double fact=1+(pCell->nearest_density_vector()[nCyto]/100);
			//std::cout << "fact: " << fact << std::endl; 
			pCell->phenotype.volume.total*=fact;
			//std::cout << "volume: " << pCell->phenotype.volume.total << std::endl;
		}   */
			// calculate distance to the cell 
			/*std::vector<double> displacement = pTestCell->position;
			displacement -= pCell->position;
			double distance = norm( displacement ); 
			
			double max_distance = pCell->phenotype.geometry.radius + 
				pTestCell->phenotype.geometry.radius; 
			max_distance *= 1.1; 
			
			// if it is not a macrophage, test for viral load 
			// if high viral load, eat it. 
		
			if( pTestCell->type == 4 && distance < max_distance )
			{
				//std::cout << "\t\tnom nom nom" << std::endl; 
				pCell->ingest_cell( pTestCell ); 
			} 
	}*/
	//boundaries(pCell);

	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dx/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE

	double ke2=parameters.doubles("k_e2");
	double aed=parameters.doubles("a_ed");
	double RD=microenvironment(PositionFinal)[debris_index];
	RD/=(aed+microenvironment(PositionFinal)[debris_index]);

	pCell->phenotype.secretion.uptake_rates[debris_index] = ke2*RD;

	//Cytokine secretion
	//DYNAMICS TO CHECK
	int TNF_index = microenvironment.find_density_index("TNFa");
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" );
	int IFN_index = microenvironment.find_density_index( "IFNg" );
	pCell->phenotype.secretion.secretion_rates[TNF_index]=parameters.doubles("k_T2");
	pCell->phenotype.secretion.secretion_rates[TNF_index]*=(1+(parameters.doubles("k_TI")/(1+1.*exp(ecf*(parameters.doubles("a_TI")-microenvironment(PositionFinal)[IFN_index])))));
	pCell->phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG2");
	pCell->phenotype.secretion.secretion_rates[IL10_index] = parameters.doubles("k_I2");
	pCell->phenotype.secretion.secretion_rates[IFN_index] = 0;

	// PROLIFERATION
	pCell->phenotype.cycle.data.transition_rate(0,0)*=(parameters.doubles("k_2")+((1-parameters.doubles("k_2"))/(1+exp(ecf*(microenvironment(PositionFinal)[IFN_index]-parameters.doubles("a_2"))))));; //proliferation ratio

	//DIFFERENTIATE INTO M1

/* 	double k02=parameters.doubles("k02");
	double a02=parameters.doubles("a02");
	double g2=k02*microenvironment(PositionFinal)[IL10_index];
	g2/=(a02+microenvironment(PositionFinal)[IL10_index]);

	double gtest1=(double)rand() / RAND_MAX;

	if(gtest1>g2*6) 
	{ */

	double k21=parameters.doubles("k_21");
	double a21=parameters.doubles("a_21");
	double g1=k21*microenvironment(PositionFinal)[TNF_index];
	g1/=(a21+microenvironment(PositionFinal)[TNF_index]);

	double gtest2=(double)rand() / RAND_MAX;

	//std::cout << "g21 " << g1 << std::endl; //NEW
	//std::cout << "gtest " << gtest2 << std::endl; //NEW


	if(gtest2<=g1)  
	{
		pCell->type=2;
		//std::cout << "M2 -> M1!" << std::endl; 
		pCell->functions.update_phenotype = macrophage1_function; 
	}

	//}

	boundaries(pCell);
	return;
}

void neutrophils_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	int debris_index = microenvironment.find_density_index("debris");	
	int TNF_index = microenvironment.find_density_index("TNFa");	
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" );
	int IFN_index = microenvironment.find_density_index( "IFNg" );

	//pCell->phenotype.molecular.internalized_total_substrates[debris_index]=6e-4; //5e-5;
	//pCell->phenotype.death.models[0]->transition_rate(0,1) = 1.0;
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// More debris it eats, faster it dies	

	int napoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Neutrophil Apoptosis" );
	pCell->phenotype.death.rates[napoptosis_model_index] = parameters.doubles("k_AN"); 
	
	double neutrodeb=pCell->phenotype.molecular.internalized_total_substrates[debris_index];
	//std::cout << "NAP " << pCell->phenotype.molecular.internalized_total_substrates[debris_index] <<std::endl; 

	if (neutrodeb-1.6e-3>0)
	{
		pCell->phenotype.death.rates[napoptosis_model_index] *= (1+(neutrodeb-1.6e-3));
		//std::cout << "NAP " << pCell->phenotype.molecular.internalized_total_substrates[debris_index] <<std::endl; 
	}

/* 	if (pCell->phenotype.cycle.model().code>5)
	{
		std::cout << "NAP " << pCell->phenotype.cycle.model().code <<std::endl; 
	} */

	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dx/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE // work on this

	double ken=parameters.doubles("k_en");
	double aedn=parameters.doubles("a_edn");
	double RD=microenvironment(PositionFinal)[debris_index];
	RD/=(aedn+microenvironment(PositionFinal)[debris_index]);

	pCell->phenotype.secretion.uptake_rates[debris_index] = ken*RD;
	
	// CYTOKINE SECRETION 

	//int TNF_index = microenvironment.find_density_index("TNFa");
	//int IFN_index = microenvironment.find_density_index( "IFNg" );
	//int IL12_index = microenvironment.find_density_index( "IL12" ); 

	pCell->phenotype.secretion.secretion_rates[TNF_index]=parameters.doubles("k_TN");
	pCell->phenotype.secretion.secretion_rates[IL10_index]=0;
	pCell->phenotype.secretion.secretion_rates[TGF_index]=0;
	pCell->phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("k_IFN"); //  /(1+exp(2.5e-12-microenvironment(PositionFinal)[IFN_index]));

	boundaries(pCell);
	return; 
}

/* void debris_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	//boundaries(pCell);
	return;
} */

/* void epithelial_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping
	
	static int nVirus = microenvironment.find_density_index( "virus" ); 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	
	// compare against viral load. Should I commit apoptosis? 
	
	double virus = phenotype.molecular.internalized_total_substrates[nVirus]; 
	if( virus >= parameters.doubles("burst_virion_count") )
	{
		std::cout << "\t\tburst!" << std::endl; 
		pCell->lyse_cell(); // start_death( apoptosis_model_index );
		pCell->functions.update_phenotype = NULL; 
		return; 
	}

	// replicate virus particles inside me 
	
	if( virus >= parameters.doubles("min_virion_count") ) 
	{
		double new_virus = parameters.doubles( "viral_replication_rate" ); 
		new_virus *= dt;
		phenotype.molecular.internalized_total_substrates[nVirus] += new_virus; 
	}
//	static double implicit_Euler_constant = 
//		(1.0 + dt * parameters.doubles("virus_digestion_rate") );
//	phenotype.molecular.internalized_total_substrates[nVirus] /= implicit_Euler_constant; 
	
	
	// if I have too many 
	// if I have too many 

	return; 
}  */

void macrophage_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double bias = parameters.doubles("m_mb");
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	
	phenotype.motility.migration_bias = bias; 
	
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( debris_index ); 
	double denominator =  norm( phenotype.motility.migration_bias_direction ) + 1e-17; 
	
	phenotype.motility.migration_bias_direction /= denominator; 
	
	return; 
}

void neutrophils_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double bias = parameters.doubles("n_mb");
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	
	phenotype.motility.migration_bias = bias; 
	
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( debris_index ); 
	double denominator =  norm( phenotype.motility.migration_bias_direction ) + 1e-17; 
	
	phenotype.motility.migration_bias_direction /= denominator; 
	
	return; 
}

std::vector<double> integrate_total_substrates( void )
{
	// start with 0 vector 
	std::vector<double> out( microenvironment.number_of_densities() , 0.0 ); 

	// integrate extracellular substrates 
	for( unsigned int n = 0; n < microenvironment.number_of_voxels() ; n++ )
	{
		// out = out + microenvironment(n) * dV(n) 
		axpy( &out , microenvironment.mesh.voxels[n].volume , microenvironment(n) ); 
	}

	// inte
	for( unsigned int n=0; n < (*all_cells).size(); n++ )
	{
		Cell* pC = (*all_cells)[n];
		out += pC->phenotype.molecular.internalized_total_substrates;
	}
	
	return out; 
}

std::vector<Cell*> get_possible_neighbors( Cell* pCell )
{
	std::vector<Cell*> neighbors = {}; 

	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end =
		pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ neighbors.push_back( *neighbor ); }

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{ neighbors.push_back( *neighbor ); }
	}
	
	return neighbors; 
}