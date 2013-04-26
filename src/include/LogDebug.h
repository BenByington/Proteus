/* 
 * File:   LogDebug.h
 * Author: Ben
 *
 * Created on April 2, 2010, 11:58 AM
 */


#undef trace
#define trace(format, ...)

#undef debug
#define debug(format, ...) fprintf(procFile, "DEBUG: %s %d: ", __FILE__, __LINE__); fprintf(procFile, format, ##__VA_ARGS__); fflush(procFile);

#undef info
#define info(format, ...) fprintf(procFile, "INFO:  %s %d: ", __FILE__, __LINE__); fprintf(procFile, format, ##__VA_ARGS__); fflush(procFile);

#undef warn
#define warn(format, ...) fprintf(procFile, "WARN:  %s %d: ", __FILE__, __LINE__); fprintf(procFile, format, ##__VA_ARGS__); fflush(procFile);

#undef error
#define error(format, ...) fprintf(procFile, "ERROR: %s %d: ", __FILE__, __LINE__); fprintf(procFile, format, ###__VA_ARGS__); fflush(procFile);

