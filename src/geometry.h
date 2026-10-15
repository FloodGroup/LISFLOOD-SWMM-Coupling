#pragma once
#include "lisflood.h"

typedef struct Geometry
{
	int xsz;
	int ysz;
	NUMERIC_TYPE blx;
	NUMERIC_TYPE bly;
	NUMERIC_TYPE tly;
	NUMERIC_TYPE dx;
	NUMERIC_TYPE dy;
} Geometry;


typedef struct vconnect
{
	int  count_index; 
	int* node_Vindex; 
	int* dem_Vindex; 
	int* dem_count;
	int* node_dem_slot; // compact cell slot for each mapped inlet

	double* node_surfelev; 

	double* node_c0; 
	double* node_area; 
	double* node_qoutmax; 

	double* node_cw; 
	double* node_w; 
	double* node_qinmax; 

	double* node_c1;
	double* node_area1;

	double* flow_1D_node; 
	double* flow_2D_node; 
	double* flow_2D_dem; 

	int* flow_neg_count; 

	int count_dem_index; 
	int* dem_Vindex_uni; 
	int* dem_xi;         
	int* dem_yi;         
	int* flow_type;      

	// Outfalls whose
	// coordinates fall on a valid LISFLOOD DEM cell are stored here.
	int outfall_count;
	int* outfall_node_Vindex;
	int* outfall_dem_Vindex;
	int outfall_dem_count;
	int* outfall_dem_Vindex_uni;
	unsigned char* outfall_dynamic_stage; // 1 for mapped COUPLED outfalls
	int dynamic_outfall_count;
}Vconnect;
