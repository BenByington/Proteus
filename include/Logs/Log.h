/*
 * Copywrite 2013 Benjamin Byington
 *
 * This file is part of the Proteus software package
 * 
 * Proteus is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Proteus is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with Proteus.  If not, see <http://www.gnu.org/licenses/>
 */

/*********
 * This is a light weight logging system that is configured at compile time.
 * There are five defined levels of logging macros, which in order of priority,
 * are error(), warn(), info(), debug() and trace().  These macros have the
 * same arguments as printf (indeed the arguments are passed TO printf), so they
 * take a formatted string with optional values for that string.  
 * 
 * The level of debugging is controlled through the use of including the 
 * specific header file for that level of logging, where the name of the logging
 * file indicates the most detailed level of logging to be written.  The 
 * specific version included in this file is the default for the entire program,
 * and can be changed locally through additional include statements.  Logging 
 * statements of a lower priority than is enabled (i.e. debug() statements when 
 * only Info and above are in use) are defined as empty macros and are removed 
 * by the compiler.
 * 
 * For example:  
 *   - If you wish to globally suppress trace and debug statements, then at the
 *     bottom of this file you place #include LogInfo.h
 *   - If you are debugging a section of code and locally wish more detailed 
 *     messages then you braked the region of interest with the relevant
 *     includes:
 *     #include LogTrace.h
 *       <code snippet>
 *     #include Log.h
 * 
 * Note: 
 *   - If locally changing the level of logging, always remember to again
 *     include Log.h afterwards, so that regions of code later in the file
 *     remain unaffected.
 *   - Trace creates output messages during some of the tight loops.  It is
 *     NOT recommended to enable this as the program default.
 *********/

#ifndef _LOG_H
#define	_LOG_H

#include "Precision.h"
#include <stdio.h>
#include <time.h>

extern FILE * procFile;

/*
 * Initialization routine for the io stream used for logging.  
 * 
 * Only call once during execution, and before any logging macros are used!!!
 */
void initLogging();

/*
 * Cleans up logging after program execution.
 */
void endLogging();

#endif	/* _LOG_H */

//Include one of the log levels here.  The program will default to this level
//of verbosity.
#include "Logs/LogInfo.h"

