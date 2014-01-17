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

#ifndef BASETEST_H
#define	BASETEST_H

class Test
{
protected:
    Test(){};
    
public:
    virtual ~Test(){}
    
    virtual void setup() = 0;
    virtual bool run() = 0;
    virtual void cleanup() = 0;
};

#endif	/* BASETEST_H */

