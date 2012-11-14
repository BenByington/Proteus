/* 
 * File:   LogTrace.h
 * Author: Ben
 *
 * Created on April 2, 2010, 12:18 PM
 */

//TODO Add timing functionality to the other log levels
#undef trace
#define trace(format, ...) fprintf(procFile, "TRACE: %s %d - %f: ", __FILE__, __LINE__, (double)clock()); fprintf(procFile, format, __VA_ARGS__); fflush(procFile);

#undef debug
#define debug(format, ...) fprintf(procFile, "DEBUG: %s %d - %f: ", __FILE__, __LINE__, (double)clock()); fprintf(procFile, format, __VA_ARGS__); fflush(procFile);

#undef info
#define info(format, ...) fprintf(procFile, "INFO:  %s %d - %f: ", __FILE__, __LINE__, (double)clock()); fprintf(procFile, format, __VA_ARGS__); fflush(procFile);

#undef warn
#define warn(format, ...) fprintf(procFile, "WARN:  %s %d - %f: ", __FILE__, __LINE__, (double)clock()); fprintf(procFile, format, __VA_ARGS__); fflush(procFile);

#undef error
#define error(format, ...) fprintf(procFile, "ERROR: %s %d - %f: ", __FILE__, __LINE__, (double)clock()); fprintf(procFile, format, __VA_ARGS__); fflush(procFile);


