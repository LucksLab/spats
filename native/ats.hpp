
#ifndef __SPATS_ATS_HPP_INCLUDED__
#define __SPATS_ATS_HPP_INCLUDED__


#include <assert.h>
#include <iostream>

#if 1

#define __FILENAME__ (strrchr(__FILE__,'/') ? strrchr(__FILE__,'/')+1 : __FILE__)
#define ___ATS_ASSERT(cond,failureBlock) do { if (!(cond)) { failureBlock; assert(0); } } while (0)
#define ATS_WHERE_BLOCK(msg) { std::cerr << "ASSERT " << msg << " at " << __FILENAME__ << ":" << __FUNCTION__ << ":" << __LINE__ << std::endl; }
#define ATS_ASSERT_NOT_REACHED() ___ATS_ASSERT(0, ATS_WHERE_BLOCK("NOTREACHED"))
#define ATS_ASSERT(cond) ___ATS_ASSERT((cond), ATS_WHERE_BLOCK("failure"))

#define ATS_PRINT(...) do { printf(__VA_ARGS__); printf("\n"); } while (0)

#else

#define ATS_ASSERT_NOT_REACHED()
#define ATS_ASSERT(cond)
#define ATS_PRINT(...)

#endif

#define ATS_IGNORE(...)

#define ATS_WARN ATS_PRINT
#define ATS_DEBUG ATS_PRINT
#define ATS_VERBOSE ATS_IGNORE

#define ATS_THROW_IF(theCond) do { if (theCond) throw std::exception(); } while(0)
#define ATS_RET_IF(theCond) do { if (theCond) { ATS_ASSERT_NOT_REACHED(); return; } } while(0)
#define ATS_RETVAL_IF(theCond,theVal) do { if (theCond) { ATS_ASSERT_NOT_REACHED(); return (theVal); } } while(0)


#endif // __SPATS_ATS_HPP_INCLUDED__
