/*
 *  bundles.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/6/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <list>
#include <map>

#include "common.h"
#include "bundles.h"

using namespace std;

bool FragmentFactory::next_fragment(MateHit& fragment_out)
{	
	char left_bwt_buf[2048];
    char right_bwt_buf[2048];
    
    shared_ptr<ReadHit> read_1(new ReadHit());
    shared_ptr<ReadHit> read_2(new ReadHit());
    int refid = 0;
    bool treated = false;
    
    while (true && !feof(hit_file))
    {
        fgets(left_bwt_buf, 2048, hit_file);
        char* nl = strrchr(left_bwt_buf, '\n');
        if (nl) *nl = 0;
        
        if (*left_bwt_buf == '@')
            continue;
        
        
        
        bool left_read_good = sam_hit_fac.get_hit_from_buf(_next_line_num, 
                                                           left_bwt_buf, 
                                                           *read_1, 
                                                           false);
        char* left_tab = strchr(left_bwt_buf, '\t');
        string left_name;
        if (left_tab)
        {
            *left_tab = 0;
            left_name = left_bwt_buf;
        }
        
        
        fgets(right_bwt_buf, 2048, hit_file);
        nl = strrchr(right_bwt_buf, '\n');
        if (nl) *nl = 0;
        
        bool right_read_good = sam_hit_fac.get_hit_from_buf(_next_line_num, 
                                                           right_bwt_buf, 
                                                           *read_2, 
                                                           false);
        
        if (!left_read_good && !feof(hit_file))
            continue;
        
        if (feof(hit_file))
            return false;
        
        
        char* right_tab = strchr(right_bwt_buf, '\t');
        string right_name;
        if (right_tab)
        {
            *right_tab = 0;
            right_name = right_bwt_buf;
        }
        
        if (!right_read_good && !feof(hit_file))
            continue;
        
        if (feof(hit_file))
            return false;
        
        // 84696373 is "*" under FNV hash -> SAM unaligned read record
        if (read_1->ref_id() == 84696373 || read_2->ref_id() == 84696373)   
            continue;
        
        if (read_1->ref_id() != read_2->ref_id())
		{
			int a = 4;
			continue;
        }
		
		assert (read_1->insert_id() == read_2->insert_id());
        assert (read_1->insert_id() != 0);
		
		if (read_1->antisense_align() == read_2->antisense_align())
            continue;
		
		//fprintf(stderr, "processing fragment %llu\n", read_1->insert_id());
        
        if (left_name != right_name)
            continue;
        
        assert (read_1->insert_id() == read_2->insert_id());
        
        //if (read_1->antisense_align() == read_2->antisense_align())
        //    continue;
        
        //fprintf(stderr, "Treating %s and %s as pairs from a fragment\n",
        //        left_name.c_str(), right_name.c_str());
        
        refid = read_1->ref_id(); 
		
		//fprintf(stderr, "%s: %s\n", left_name.c_str(), ++left_tab);
		//fprintf(stderr, "%s: %s\n", right_name.c_str(), ++right_tab);
		
		int a = 43;
		
        break;
    }
    
//    if (treated)
//    {
//        if (!read_2->first_in_pair()) 
//        {
//            treated = true;
//        }
//        else
//        {
//            treated = false;
//        }
//    }
//    else
//    {
//        if (!read_2->first_in_pair()) 
//        {
//            treated = false;
//        }
//        else
//        {
//            treated = true;
//        }   
//    }
    
    
    fragment_out = MateHit(refid,
                           treated,
                           read_1, 
                           read_2,
                           0, 
                           0);
    
    return true;
}
