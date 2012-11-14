/* 
 * File:   LogError.h
 * Author: Ben
 *
 * Created on April 2, 2010, 12:16 PM
 */

#undef trace
#define trace(format, ...)

#undef debug
#define debug(format, ...)

#undef info
#define info(format, ...)

#undef warn
#define warn(format, ...) 

#undef error
#define error(format, ...) fprintf(procFile, "ERROR: %s %d: ", __FILE__, __LINE__); fprintf(procFile, format, __VA_ARGS__); fflush(procFile);


