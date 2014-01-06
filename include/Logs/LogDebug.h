/*
 * Copywrite 2013 Benjamin Byington
 *
 * This file is part of the IMHD software package
 * 
 * IMHD is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public Liscence as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * IMHD is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with IMHD.  If not, see <http://www.gnu.org/licenses/>
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
#define error(format, ...) fprintf(procFile, "ERROR: %s %d: ", __FILE__, __LINE__); fprintf(procFile, format, ##__VA_ARGS__); fflush(procFile);

