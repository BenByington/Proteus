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

#include "Properties.h"
#include "Environment.h"
#include "Logs/Log.h"
#include "Physics.h"

#include <fstream>
#include <string.h>
#include <ctype.h>
#include <cstdlib>
#include <algorithm>

using namespace std;

#define PROBLEM_SIZE "ProblemSize"
#define IO "IO"
#define INITIAL_CONDITIONS "InitialConditions"
#define PHYSICS "Physics"
#define FORCINGS "Forcings"
#define INTEGRATION "Integration"

const string on("on");
const string off("off");

/*
 * The parameters file is divided into several different "topics", and each 
 * topic is handled by a different, but similar routine (to reduce clutter).
 * Only one of these will receive comments (below) since they all follow
 * the same pattern.
 */
void parseProblemSize(iostream & in);
void parseIO(iostream & in);
void parseIC(iostream & in);
void parsePhysics(iostream & in);
void parseForcings(iostream & in);
void parseIntegration(iostream & in);

//Currently does nothing...
void init();

/*
 * Open up the configuration file, and call the sub-parsers to grab environment 
 * variables out of it.
 */
void loadPrefs(char * loc)
{
    info("Loading Preferences from file %s\n",loc);
    const char starts = '[';
    const char ends = ']';
    
    string line;
    string section;
    fstream input(loc, ios_base::in);

    init();

    int ione;
    int itwo;
    while(!getline(input, line).eof())
    {
        //Burn through lines until we come to the start of a section
        ione = line.find_first_of(starts);
        itwo = line.find_first_of(ends);

        if(ione == -1 || itwo == -1 || itwo < ione)
            continue;

        //Figure out which section we stumbled upon, and hand it to the
        //appropriate parser.
        section = line.substr(ione+1, itwo - ione - 1);

        if(section == PROBLEM_SIZE)
        {
            parseProblemSize(input);
        }
        else if(section == IO)
        {
            parseIO(input);
        }
        else if(section == INITIAL_CONDITIONS)
        {
            parseIC(input);
        }
        else if(section == PHYSICS)
        {
            parsePhysics(input);
        }
        else if(section == FORCINGS)
        {
            parseForcings(input);
        }
        else if(section == INTEGRATION)
        {
            parseIntegration(input);
        }
    }
}

void init()
{
}

void parseProblemSize(iostream & in)
{
    /*
     * These are all of the valid strings we recognize as belonging to this 
     * section.  Any other strings present will result in a warning message
     * being placed in the logs.
     */
    const string snx("nx");
    const string sny("ny");
    const string snz("nz");
    const string sxmx("xmx");
    const string symx("ymx");
    const string szmx("zmx");
    const string shdiv("hdiv");
    const string svdiv("vdiv");

    string line;
    string one;
    string two;
    int index;

    //Read the file line by line
    debug("Loading PROBLEM_SIZE\n",0);
    while(!getline(in, line).eof())
    {
        //If the line does not contain an assigment, it is garbage.  If it
        //contains a '[' then we have hit the end of our section and we must
        //surrender control back to the calling routine.
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if((int)line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;
        
        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        //Stupid long if/else if to see which of the expected inputs we have now
        //Any input found will be parsed, stored into the appropriate variable,
        //and then logged.
        if((int)one.find(snx) != -1)
        {
            nx = atoi(two.c_str());
            debug("nx = %d\n", nx);
        }
        else if((int)one.find(sny) != -1)
        {
            ny = atoi(two.c_str());
            debug("ny = %d\n", ny);
        }
        else if((int)one.find(snz) != -1)
        {
            nz = atoi(two.c_str());
            debug("nz = %d\n", nz);
        }
        else if((int)one.find(sxmx) != -1)
        {
            xmx = atof(two.c_str());
            debug("xmx = %g\n", xmx);
        }
        else if((int)one.find(symx) != -1)
        {
            ymx = atof(two.c_str());
            debug("ymx = %g\n", ymx);
        }
        else if((int)one.find(szmx) != -1)
        {
            zmx = atof(two.c_str());
            debug("zmx = %g\n", zmx);
        }
        else if((int)one.find(shdiv) != -1)
        {
            hdiv = atoi(two.c_str());
            debug("hdiv = %d\n", hdiv);
        }
        else if((int)one.find(svdiv) != -1)
        {
            vdiv = atoi(two.c_str());
            debug("vdiv = %d\n", vdiv);
        }
        else
        {
           warn("Found unknown value in properties file!!  %s\n", line.c_str());
        }
    }

}

void parseIO(iostream & in)
{
    const string nip("N_IO_PROCS");
    const string statusr("statusRate");
    const string spatialr("spatialRate");
    const string scalarr("scalarRate");
    const string scalarpf("scalarPerF");
    const string sCheckRate("checkRate");

    string line;
    string one;
    string two;
    int index;

    debug("Loading IO\n",0)
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if((int)line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(nip) != -1)
        {
            n_io_nodes = atoi(two.c_str());
            debug("n_io_nodes = %d\n", n_io_nodes);
        }
        else if((int)one.find(statusr) != -1)
        {
            statusRate = atoi(two.c_str());
            debug("statusRate = %d\n", statusRate);
        }
        else if((int)one.find(spatialr) != -1)
        {
            spatialRate = atoi(two.c_str());
            debug("spatialRate = %d\n", spatialRate);
        }
        else if((int)one.find(scalarr) != -1)
        {
            scalarRate = atoi(two.c_str());
            debug("scalarRate = %d\n", scalarRate);
        }
        else if((int)one.find(scalarpf) != -1)
        {
            scalarPerF = atoi(two.c_str());
            debug("scalarPerF = %d\n", scalarPerF);
        }
        else if((int)one.find(sCheckRate) != -1)
        {
            checkRate = atoi(two.c_str());
            debug("checkpoint frequency = %d\n", checkRate);
        }
        else
        {
            warn("Found unknown value!!:  %s\n", line.c_str());
        }
    }
}

void parseIC(iostream & in)
{
    const string st("startType");
    const string sStartDir("startDir");

    string line;
    string one;
    string two;
    int index;

    debug("Loading Initial Conditions\n",0)
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if((int)line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(st) != -1)
        {
            transform(two.begin(), two.end(), two.begin(), ::tolower);
            trace("Lower case version of string: %s\n", two.c_str());
            if(two.find(string("scratch")) != two.npos)
            {
                startFlag = SCRATCH;
            }
            else if(two.find(string("spatial")) != two.npos)
            {
                startFlag = SPATIAL;
            }
            else if(two.find(string("checkpoint")) != two.npos)
            {
                startFlag = CHECKPOINT;
            }
            else
            {
                error("Start Type %s is invalid\n", two.c_str());
            }
            debug("startFlag = %d\n", startFlag);
        }
        else if((int)one.find(sStartDir) != -1)
        {
            int len = two.length()+1;
            startDir = (char*)malloc(len);
            strcpy(startDir, two.c_str());

            debug("Start Directory is %s\n", startDir);
        }
        else
        {
            warn("Found unknown value!!:  %s\n", line.c_str());
        }
    }
}

void parsePhysics(iostream & in)
{
    const string sMomentum("momentumEQ");
    const string sTemperature("temperatureEQ");
    const string sMagnetic("magneticEQ");
    const string sMomAdvect("momAdvection");
    const string sViscosity("viscosity");
    const string sBuoyancy("buoyancy");
    const string sMagBuoy("magBuoy");
    const string sLorentz("lorentz");
    const string sTDiff("tdiff");
    const string sTempAdvection("tempAdvection");
    const string sTempBackground("tempBackground");
    const string sMagDiff("magdiff");
    const string sMagAdvect("magAdvect");
    const string sPR("Pr");
    const string sRA("Ra");
    const string sPM("Pm");
    const string sAlpha("alpha");
    const string sMBuoyScale("mBuoyScale");
    
    string line;
    string one;
    string two;
    int index;

    debug("Loading Physics Parameters\n",0)
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if((int)line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        
        if((int)one.find(sMomentum) != -1)
        {
            if((int)two.find(on) != -1)
                momEquation = 1;
            else if((int)two.find(off) != -1)
                momEquation = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            debug("Momentum equation flag: %d\n", momEquation);
        }
        else if((int)one.find(sTemperature) != -1)
        {
            if((int)two.find(on) != -1)
                tEquation = 1;
            else if((int)two.find(off) != -1)
                tEquation = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            debug("Temperature equation flag: %d\n", tEquation);
        }
        else if((int)one.find(sMagnetic) != -1)
        {
            if((int)two.find(on) != -1)
                magEquation = 1;
            else if((int)two.find(off) != -1)
                magEquation = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            debug("Magnetic equation flag: %d\n", magEquation);
        }
        else if((int)one.find(sViscosity) != -1)
        {
            if((int)two.find(on) != -1)
                viscosity = 1;
            else if((int)two.find(off) != -1)
                viscosity = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }
            debug("viscosity term flag: %d\n", viscosity);
        }
        else if((int)one.find(sMomAdvect) != -1)
        {
            if((int)two.find(on) != -1)
                momAdvection = 1;
            else if((int)two.find(off) != -1)
                momAdvection = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }
            debug("momentum advection flag: %d\n", momAdvection);
        }
        else if((int)one.find(sBuoyancy) != -1)
        {
            if((int)two.find(on) != -1)
                buoyancy = 1;
            else if((int)two.find(off) != -1)
                buoyancy = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("buoyancy flag: %d\n", buoyancy);
        }
        else if((int)one.find(sMagBuoy) != -1)
        {
            if((int)two.find(on) != -1)
                magBuoy = 1;
            else if((int)two.find(off) != -1)
                magBuoy = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("magBuoy flag: %d\n", buoyancy);
        }
        else if((int)one.find(sLorentz) != -1)
        {
            if((int)two.find(on) != -1)
                lorentz = 1;
            else if((int)two.find(off) != -1)
                lorentz = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("lorentz flag: %d\n", lorentz);
        }
        else if((int)one.find(sTDiff) != -1)
        {
            if((int)two.find(on) != -1)
                tDiff = 1;
            else if((int)two.find(off) != -1)
                tDiff = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Temperature Diffusion flag: %d\n", tDiff);
        }
        else if((int)one.find(sTempAdvection) != -1)
        {
            if((int)two.find(on) != -1)
                tempAdvection = 1;
            else if((int)two.find(off) != -1)
                tempAdvection = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Temperature Advection flag: %d\n", tempAdvection);
        }
        else if((int)one.find(sTempBackground) != -1)
        {
            if((int)two.find(on) != -1)
                tempBackground = 1;
            else if((int)two.find(off) != -1)
                tempBackground = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Temperature Background flag: %d\n", tempAdvection);
        }
        else if((int)one.find(sMagDiff) != -1)
        {
            if((int)two.find(on) != -1)
                magDiff = 1;
            else if((int)two.find(off) != -1)
                magDiff = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Magnetic Diffusion flag: %d\n", magDiff);
        }
        else if((int)one.find(sMagAdvect) != -1)
        {
            if((int)two.find(on) != -1)
                magAdvect = 1;
            else if((int)two.find(off) != -1)
                magAdvect = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Magnetic \"Advection\" flag: %d\n", magAdvect);
        }
        else if((int)one.find(sPR) != -1)
        {
            Pr = atof(two.c_str());
            debug("Pr set to %g\n", Pr);
        }
        else if((int)one.find(sRA) != -1)
        {
            Ra = atof(two.c_str());
            debug("Ra set to %g\n", Ra);
        }
        else if((int)one.find(sPM) != -1)
        {
            Pm = atof(two.c_str());
            debug("Magnetic Pr set to %g\n", Pm);
        }
        else if((int)one.find(sAlpha) != -1)
        {
            alpha = atof(two.c_str());
            debug("alpha set to %g\n", alpha);
        }
        else if((int)one.find(sMBuoyScale) != -1)
        {
            magBuoyScale = atof(two.c_str());
            debug("mBuoyScale set to %g\n", alpha);
        }
        else
        {
            warn("Found unknown value!!:  %s %s\n", one.c_str(), two.c_str());
        }
    }
}

void parseForcings(iostream & in)
{
    const string sMomStaticForcing("momStaticForcing");
    const string sMagStaticForcing("magStaticForcing");
    const string sMomTimeForcing("momTimeForcing");
    const string sForcingFile("forcingFile");
    const string sMagForcingFile("magForcingFile");
    const string sMagTimeForcing("magTimeForcing");
    const string sKinematic("kinematic");
    const string sMomOmega("momOmega");
    const string sMomEps("momEps");
    const string sMagK("magK");
    const string sMagW("magW");
    const string sMagB0("magB0");
    const string sRecentering("recentering");
    const string sDefineHalt("defineHalt");
    const string sSanitize("sanitizeBoundaries");

    string line;
    string one;
    string two;
    int index;

    debug("Loading Forcing Parameters\n",0)
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());

        index = line.find_first_of('=');
        if((int)line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(sForcingFile) != -1)
        {
            int len = two.length()+1;
            forceFile = (char*)malloc(len);
            strcpy(forceFile, two.c_str());

            debug("Forcing file is = %s\n", forceFile);
        }
        else if((int)one.find(sMagForcingFile) != -1)
        {
            int len = two.length()+1;
            magForceFile = (char*)malloc(len);
            strcpy(forceFile, two.c_str());

            debug("Magnetic Forcing file is = %s\n", magForceFile);
        }
        else if((int)one.find(sMomStaticForcing) != -1)
        {
            if((int)two.find(on) != -1)
                momStaticForcing = 1;
            else if((int)two.find(off) != -1)
                momStaticForcing = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            debug("Momentum static forcing flag: %d\n", momStaticForcing);
        }
        else if((int)one.find(sMagStaticForcing) != -1)
        {
            if((int)two.find(on) != -1)
                magStaticForcing = 1;
            else if((int)two.find(off) != -1)
                magStaticForcing = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            debug("Magnetic static forcing flag: %d\n", momStaticForcing);
        }
        else if((int)one.find(sMomTimeForcing) != -1)
        {
            if((int)two.find(on) != -1)
                momTimeForcing = 1;
            else if((int)two.find(off) != -1)
                momTimeForcing = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            debug("Momentum time forcing flag: %d\n", momTimeForcing);
        }
        else if((int)one.find(sKinematic) != -1)
        {
            if((int)two.find(on) != -1)
                kinematic = 1;
            else if((int)two.find(off) != -1)
                kinematic = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Kinematic problem flag: %d\n", kinematic);
        }
        else if((int)one.find(sMagTimeForcing) != -1)
        {
            if((int)two.find(on) != -1)
                magTimeForcing = 1;
            else if((int)two.find(off) != -1)
                magTimeForcing = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Magnetic time forcing flag: %d\n", magTimeForcing);
        }
        else if((int)one.find(sMomOmega) != -1)
        {
            momOmega = atof(two.c_str());
            debug("momOmega = %g\n", momOmega);
        }
        else if((int)one.find(sMomEps) != -1)
        {
            momEps = atof(two.c_str());
            debug("momEps = %g\n", momEps);
        }
        else if((int)one.find(sMagK) != -1)
        {
            magK = atof(two.c_str());
            debug("magK = %g\n", magK);
        }
        else if((int)one.find(sMagW) != -1)
        {
            magW = atof(two.c_str());
            debug("magW = %g\n", magW);
        }
        else if((int)one.find(sMagB0) != -1)
        {
            magB0 = atof(two.c_str());
            debug("magB0 = %g\n", magB0);
        }
        else if((int)one.find(sRecentering) != -1)
        {
            const string sMaxByCenter("maxByCenter");
            if((int)two.find(off) != -1)
                recentering = NOCENTERING;
            else if((int)two.find(sMaxByCenter) != -1)
                recentering = BYMAXCENTER;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Recentering flag: %d\n", recentering);
        }
        else if((int)one.find(sDefineHalt) != -1)
        {
            if((int)two.find(off) != -1)
                recenterTerminate = 0;
            //else if((int)two.find(off) != -1)
            //    magTimeForcing = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("recenterTerminate Function: %s\n", two);
        }
        else if((int)one.find(sSanitize) != -1)
        {
            if((int)two.find(on) != -1)
                sanitize = 1;
            else if((int)two.find(off) != -1)
                sanitize = 0;
            else
            {
                warn("unrecognized option %s for %s\n", two.c_str(), one.c_str());
            }
            debug("Sanitize boundaries flag: %d\n", sanitize);
        }
        else
        {
            warn("Found unknown value!!:  %s %s\n", one.c_str(), two.c_str());
        }
    }
}

void parseIntegration(iostream & in)
{
    const string sSafety("safetyFactor");
    const string sMaxSteps("maxSteps");
    const string sMaxTime("maxTime");

    string line;
    string one;
    string two;
    int index;

    debug("Loading Integration Parameters\n",0)
    while(!getline(in, line).eof())
    {
        debug("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if((int)line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(sSafety) != -1)
        {
            safetyFactor = atof(two.c_str());
            debug("Safety factor = %g\n", safetyFactor);
        }
        else if((int)one.find(sMaxSteps) != -1)
        {
            maxSteps = atoi(two.c_str());
            debug("maxSteps = %d\n", maxSteps);
        }
        else if((int)one.find(sMaxTime) != -1)
        {
            maxTime = atof(two.c_str());
            debug("maxTime = %f\n", maxTime);
        }
        else
        {
            warn("Found unknown value!!:  %s\n", line.c_str());
        }
    }
}


