/*
 * Copyright 2013 Benjamin Byington
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

#include "Test.h"
#include "TestRunner.h"

#include "TopDags/TopDagsTest.h"

void runTests()
{
    int numTests = 0;
    Test ** tests = registerTests(numTests);
    
    bool * results = new bool[numTests];
    
    for(int i = 0; i < numTests; i++)
    {
        tests[i]->setup();
        results[i] = tests[i]->run();
        tests[i]->cleanup();
    }
    
    for(int i = 0; i < numTests; i++)
    {
        delete tests[i];
    }
    delete [] tests;
    delete [] results;
    
}

Test ** registerTests(int & numTests)
{
    numTests = 1;
    
    Test ** ret = new Test*[numTests];
    
    ret[0] = new TopDagsTest();
    
    return ret;
}