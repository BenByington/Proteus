/* 
 * File:   IO.h
 * Author: Ben
 *
 * Created on March 31, 2010, 3:19 PM
 */

#ifndef _IO_H
#define	_IO_H

#include <string.h>
#include "Field.h"

void testIO();
void initIO();

void writeSpatial(p_field f, char * name);
void readSpatial(p_field f, char * name);

void writeCheckpoint();
void readCheckpoint();

void performOutput();

#endif	/* _IO_H */

