#include "msu_hbt/hbt.h"

Chbt_cell::Chbt_cell(){
	int ix,iy,iz;
	neighbor.resize(3);
	for(ix=0;ix<3;ix++){
		neighbor[ix].resize(3);
		for(iy=0;iy<3;iy++){
			neighbor[ix][iy].resize(3);
			for(iz=0;iz<3;iz++){
				neighbor[ix][iy][iz]=NULL;
			}
		}
	}
	partlist_a.clear();
	partlist_b.clear();
}