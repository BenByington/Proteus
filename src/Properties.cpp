/* 
 * File:   Properties.cpp
 * Author: Ben
 * 
 * Created on March 11, 2010, 10:02 PM
 */

#include "Properties.h"
#include "Environment.h"
#include "Log.h"

#include <fstream>
#include <string>

using namespace std;

#define PROBLEM_SIZE "ProblemSize"
#define IO "IO"
#define INITIAL_CONDITIONS "InitialConditions"
#define PHYSICS "Physics"
#define INTEGRATION "Integration"

const string on("on");
const string off("off");

void parseProblemSize(iostream & in);
void parseIO(iostream & in);
void parseIC(iostream & in);
void parsePhysics(iostream & in);
void parseIntegration(iostream & in);
void init();

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
        ione = line.find_first_of(starts);
        itwo = line.find_first_of(ends);

        if(ione == -1 || itwo == -1 || itwo < ione)
            continue;

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

    debug("Loading PROBLEM_SIZE\n",0);
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if(line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;
        
        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(snx) != -1)
        {
            nx = atoi(two.c_str());
            trace("nx = %d\n", nx);
        }
        else if((int)one.find(sny) != -1)
        {
            ny = atoi(two.c_str());
            trace("ny = %d\n", ny);
        }
        else if((int)one.find(snz) != -1)
        {
            nz = atoi(two.c_str());
            trace("nz = %d\n", nz);
        }
        else if((int)one.find(sxmx) != -1)
        {
            xmx = atof(two.c_str());
            trace("xmx = %g\n", xmx);
        }
        else if((int)one.find(symx) != -1)
        {
            ymx = atof(two.c_str());
            trace("ymx = %g\n", ymx);
        }
        else if((int)one.find(szmx) != -1)
        {
            zmx = atof(two.c_str());
            trace("zmx = %g\n", zmx);
        }
        else if((int)one.find(shdiv) != -1)
        {
            hdiv = atoi(two.c_str());
            trace("hdiv = %d\n", hdiv);
        }
        else if((int)one.find(svdiv) != -1)
        {
            vdiv = atoi(two.c_str());
            trace("vdiv = %d\n", vdiv);
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

    string line;
    string one;
    string two;
    int index;

    debug("Loading IO\n",0)
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if(line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(nip) != -1)
        {
            n_io_nodes = atoi(two.c_str());
            trace("n_io_nodes = %d\n", n_io_nodes);
        }
        else if((int)one.find(statusr) != -1)
        {
            statusRate = atoi(two.c_str());
            trace("statusRate = %d\n", statusRate);
        }
        else if((int)one.find(spatialr) != -1)
        {
            spatialRate = atoi(two.c_str());
            trace("spatialRate = %d\n", spatialRate);
        }
        else if((int)one.find(scalarr) != -1)
        {
            scalarRate = atoi(two.c_str());
            trace("scalarRate = %d\n", scalarRate);
        }
        else if((int)one.find(scalarpf) != -1)
        {
            scalarPerF = atoi(two.c_str());
            trace("scalarPerF = %d\n", scalarPerF);
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
        if(line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(st) != -1)
        {
            startType = atoi(two.c_str());
            trace("startType = %d\n", startType);
        }
        else if((int)one.find(sStartDir) != -1)
        {
            int len = two.length()+1;
            startDir = (char*)malloc(len);
            strcpy(startDir, two.c_str());

            trace("Start Directory is %s\n", startDir);
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
    const string sMomStaticForcing("momStaticForcing");
    const string sMomTimeForcing("momTimeForcing");
    const string sForcingFile("forcingFile");
    const string sBuoyancy("buoyancy");
    const string sLorentz("lorentz");
    const string sTDiff("tdiff");
    const string sTempAdvection("tempAdvection");
    const string sMagDiff("magdiff");
    const string sMagAdvect("magAdvect");
    const string sKinematic("kinematic");
    const string sMagTimeForcing("magTimeForcing");
    const string RE("Re");
    const string RA("Ra");
    const string PRM("Prm");
    const string sAlpha("alpha");

    string line;
    string one;
    string two;
    int index;

    debug("Loading Physics Parameters\n",0)
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if(line.find_first_of("[") != -1)
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

            trace("Forcing directory = %s\n", forceFile);
        }
        if((int)one.find(sMomStaticForcing) != -1)
        {
            if((int)two.find(on) != -1)
                momStaticForcing = 1;
            else if((int)two.find(off) != -1)
                momStaticForcing = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            trace("Momentum static forcing flag: %d\n", momStaticForcing);
        }
        if((int)one.find(sMomTimeForcing) != -1)
        {
            if((int)two.find(on) != -1)
                momTimeForcing = 1;
            else if((int)two.find(off) != -1)
                momTimeForcing = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            trace("Momentum time forcing flag: %d\n", momTimeForcing);
        }
        else if((int)one.find(sMomentum) != -1)
        {
            if((int)two.find(on) != -1)
                momEquation = 1;
            else if((int)two.find(off) != -1)
                momEquation = 0;
            else
            {
                warn("unrecognized option %s for %s", two.c_str(), one.c_str());
            }

            trace("Momentum equation flag: %d\n", momEquation);
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

            trace("Temperature equation flag: %d\n", tEquation);
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

            trace("Magnetic equation flag: %d\n", magEquation);
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
            trace("viscosity term flag: %d\n", viscosity);
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
            trace("momentum advection flag: %d\n", momAdvection);
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
            trace("buoyancy flag: %d\n", buoyancy);
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
            trace("lorentz flag: %d\n", lorentz);
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
            trace("Temperature Diffusion flag: %d\n", tDiff);
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
            trace("Temperature Advection flag: %d\n", tempAdvection);
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
            trace("Magnetic Diffusion flag: %d\n", magDiff);
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
            trace("Magnetic \"Advection\" flag: %d\n", magAdvect);
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
            trace("Kinematic problem flag: %d\n", kinematic);
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
            trace("Magnetic time forcing flag: %d\n", magAdvect);
        }
        else if((int)one.find(RE) != -1)
        {
            Re = atof(two.c_str());
            trace("Re set to %g\n", Re);
        }
        else if((int)one.find(RA) != -1)
        {
            Ra = atof(two.c_str());
            trace("Ra set to %g\n", Ra);
        }
        else if((int)one.find(PRM) != -1)
        {
            Prm = atof(two.c_str());
            trace("Magnetic Pr set to %g\n", Prm);
        }
        else if((int)one.find(sAlpha) != -1)
        {
            alpha = atof(two.c_str());
            trace("alpha set to %g\n", alpha);
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
    const string sNSteps("nSteps");

    string line;
    string one;
    string two;
    int index;

    debug("Loading Integration Parameters\n",0)
    while(!getline(in, line).eof())
    {
        trace("Reading line %s\n", line.c_str());
        index = line.find_first_of('=');
        if(line.find_first_of("[") != -1)
            return;
        if(index == -1)
            continue;

        one = line.substr(0, index);
        two = line.substr(index+1, line.size()-1);

        if((int)one.find(sSafety) != -1)
        {
            safetyFactor = atof(two.c_str());
            trace("Safety factor = %g\n", safetyFactor);
        }
        else if((int)one.find(sNSteps) != 1)
        {
            nSteps = atoi(two.c_str());
            trace("nSteps = %d\n", nSteps);
        }
        else
        {
            warn("Found unknown value!!:  %s\n", line.c_str());
        }
    }
}


