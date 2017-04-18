#ifndef BUNDLES_H
#define BUNDLES_H
/*
 *  bundles.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/6/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <boost/bind.hpp>
#include <vector>
#include "hits.h"

class FragmentFactory
{
public:
	FragmentFactory(SAMHitFactory& fac, FILE* hfile)
	: sam_hit_fac(fac), hit_file(hfile), _next_line_num(0) {}

	bool next_fragment(MateHit& mate_out);
	
	SAMHitFactory& hit_factory() { return sam_hit_fac; } 
	
	void reset() 
	{ 
		rewind(hit_file); 
		_next_line_num = 0;
	}
	
private:
	SAMHitFactory sam_hit_fac;
	FILE* hit_file;

    int _next_line_num;
};

#endif
