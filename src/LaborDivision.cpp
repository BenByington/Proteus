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

#include <stdlib.h>
#include <stdio.h>

#include "Environment.h"
#include "LaborDivision.h"
#include "Logs/Log.h"

/*
 * This routine is in charge of things related to the size of the problem.
 * Most notably it calculates the number of wavemodes, which are derived from
 * the resolution passed in from the configuration file
 */
void lab_initGeometry()
{
    info("Initializing Problem Geometry\n");

    //physical extent of each grid cell in length units
    dx = xmx / nx;
    dy = ymx / ny;
    dz = zmx / nz;
    trace("dx dy dz: %g %g %g\n", dx, dy, dz);

    //An FFT produces the same number of wave modes as grid points.  However,
    //because this is a real to complex fft, the first fft is only have the
    //size because of the symmetries involved.
    nkx = nx;
    nky = ny/2 + 1;    //Integer division is intentional.
    nkz = nz;
    trace("nkx nky nkz: %d %d %d\n", nkx, nky, nkz);

    //Due to the nature of the pseudo-spectral algorithm, it is necessary to
    //eliminate the middle third of wave numbers (which are the smallest length
    //scales) in order to remove aliasing effects from the calculation of
    //nonlinear terms.
    dealias_kx.min = (nx/2) * 2 / 3 + 1;
    dealias_ky.min = (ny/2) * 2 / 3 + 1;
    dealias_kz.min = (nz/2) * 2 / 3 + 1;
    trace("Smallest dialised kx ky kz: %d %d %d\n", dealias_kx.min, dealias_ky.min, dealias_kz.min);

    dealias_kx.max = nkx  - dealias_kx.min;
    dealias_ky.max = nky - 1;      //This is the first dimension FFT'ed
    dealias_kz.max = nkz  - dealias_kz.min;
    trace("Largest dialised kx ky kz: %d %d %d\n", dealias_kx.max, dealias_ky.max, dealias_kz.max);

    dealias_kx.width = dealias_kx.max - dealias_kx.min + 1;
    dealias_ky.width = dealias_ky.max - dealias_ky.min + 1;
    dealias_kz.width = dealias_kz.max - dealias_kz.min + 1;
    trace("Width of dealiasing kx ky kz: %d %d %d\n", dealias_kx.width, dealias_ky.width, dealias_kz.width)

    //stores the number of wavemodes that exist after de-aliasing
    ndkx = nkx - dealias_kx.width;
    ndky = nky - dealias_ky.width;
    ndkz = nkz - dealias_kz.width;
    trace("Number of retained waves %d %d %d\n", ndkx, ndky, ndkz);

    info("Problem Geometry Done\n");
}


/*
 * Here we must set up the MPI communication groups.  There are 5:
 * 
 * 1. Compute nodes - a logical 2D array of processes excluding IO nodes.
 * 2. Horizontal groups - Communication along rows of compute nodes.
 * 3. Vertical groups - Communication along columns of compute nodes.
 * 4. IO groups - 1 IO node and an integer number of compute rows.
 * 5. File group - Between all IO nodes for parallel file writes.
 * 
 * Nodes are laid out in the following order:
 * 1. The first nxm nodes are compute nodes, were there are n rows and m columns
 * 2. The compute nodes are laid out in row major order
 * 3. The last X nodes are io nodes, where there are 1 < X <= n io nodes.
 * 4. Simple layout example:
 *      c01 c02 c03 c04    I17
 *      c05 c06 c07 c08    I18
 *      c09 c10 c11 c12    I19
 *      c13 c14 c15 c16    I20
 */
void lab_initGroups()
{
    info("Initializing Communication Groups\n");

    int i;

    //Divide rows (contiguously) among io nodes as evenly as possible.  
    //If they cannot be divided evenly, the extra go to the first set.
    int div = vdiv / n_io_nodes;
    int rem = vdiv % n_io_nodes;
    io_layers = (indexes*)malloc(n_io_nodes*sizeof(indexes));

    for(i = 0; i < rem; i++)
    {
        io_layers[i].width = div+1;
        trace("io_layers[%d].width = %d\n", i, io_layers[i].width);
    }
    for(i = rem; i < n_io_nodes; i++)
    {
        io_layers[i].width = div;
        trace("io_layers[%d].width = %d\n", i, io_layers[i].width);
    }

    //we have widths, now explicitly record the first and last index for rows
    //owned by each io node
    io_layers[0].min = 0;
    io_layers[0].max = io_layers[0].width-1;
    trace("io_layers[%d].min = %d\n", 0, io_layers[0].min);
    trace("io_layers[%d].max = %d\n", 0, io_layers[0].max);
    for(i = 1; i < n_io_nodes; i++)
    {
        io_layers[i].min = io_layers[i-1].max + 1;
        io_layers[i].max = io_layers[i].min + io_layers[i].width - 1;
        trace("io_layers[%d].min = %d\n", i, io_layers[i].min);
        trace("io_layers[%d].max = %d\n", i, io_layers[i].max);
    }

    //Sanity check...
    if(io_layers[n_io_nodes-1].max != vdiv-1)
    {
        error("Logic for setting IO layers is broken!!!  %d layers calculated but %d layers in group\n", io_layers[n_io_nodes-1].max, vsize-1);
        abort();
    }

    my_io_layer = -1;

    //layout for the pencils.  A row is defined by consecutive processors
    //in the global group

    if(grank < hdiv * vdiv)
    {
        compute_node = 1;
        io_node = 0;
    }
    else if(grank < hdiv*vdiv + n_io_nodes)
    {
        compute_node = 0;
        io_node = 1;
    }
    //this is just in case someone set mpi to run with more nodes than were
    //needed.  Its wasteful, but at least we can still run gracefully.
    else
    {
        compute_node = 0;
        io_node = 0;
        warn("I am a useless node!  More nodes provided than were requested in parameters file!\n");
    }

    info("compute_node flag set to %d\n", compute_node);
    info("io_node flag set to %d\n", io_node);

    //grab the global group
    MPI_Group global;
    MPI_Comm_group(MPI_COMM_WORLD, &global);

    
    //logical compute grid is set up so that rows are contiguous in the global
    //group, and compute nodes come before the IO nodes.
    if(compute_node)
    {
        debug("Setting up hcomm and vcomm\n");
        //set up the groups for inernal transposes and such
        int row = grank / hdiv;
        int col = grank % hdiv;
        debug("row col = %d %d\n", row, col);

        int ytrip[1][3];
        int ztrip[1][3];

        ytrip[0][0] = row * hdiv;             //start index of our row
        ytrip[0][1] = (row+1)*hdiv - 1;       //final index of our row
        ytrip[0][2] = 1;                      //stride between elements of our row
        trace("tripplet for horizontal group: %d %d %d\n", ytrip[0][0], ytrip[0][1], ytrip[0][2]);

        ztrip[0][0] = col;                    //start index of our column
        ztrip[0][1] = (vdiv-1)*hdiv + col;    //final index of our column
        ztrip[0][2] = hdiv;                   //stride between elements of our column
        trace("tripplet for vertical group: %d %d %d\n", ztrip[0][0], ztrip[0][1], ztrip[0][2]);

        trace("Creating slab groups\n");
        //create groups for our horizontal and vertical associations
       MPI_Group hgroup;
        MPI_Group vgroup;
        MPI_Group_range_incl(global, 1, ytrip, &hgroup);
        MPI_Group_range_incl(global, 1, ztrip, &vgroup);

        trace("Creating communicators for groups\n");
        //get the communicators for our groups
        MPI_Comm_create(MPI_COMM_WORLD, hgroup, &hcomm);
        MPI_Comm_create(MPI_COMM_WORLD, vgroup, &vcomm);

        trace("Getting rank and size for groups\n");
        MPI_Comm_rank(hcomm, &hrank);
        MPI_Comm_rank(vcomm, &vrank);
        MPI_Comm_size(hcomm, &hsize);
        MPI_Comm_size(vcomm, &vsize);

        //figure out which io layer we are in
        my_io_layer = 0;
        while(vrank > io_layers[my_io_layer].max)
            my_io_layer++;
        debug("We belong to IO layer %d\n", my_io_layer);
    }
    else
    {
        trace("Doing dummy call for slab group creation\n");
        //collective routine.  Our IO nodes need to at least check in!
        MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, &hcomm);
        MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, &vcomm);
    }

    //prep computational comm
    if(compute_node)
    {
        debug("Setting up compute nodes group\n");
        int ctrip[1][3];

        ctrip[0][0] = 0;                 //start index of our row
        ctrip[0][1] = hdiv*vdiv-1;       //final index of our row
        ctrip[0][2] = 1;                 //stride between elements of our row
        trace("Tripplet for compute node group: %d %d %d\n", ctrip[0][0], ctrip[0][1], ctrip[0][2]);

        trace("Creating group\n");
        MPI_Group cgroup;
        MPI_Group_range_incl(global, 1, ctrip, &cgroup);

        trace("Getting communicator\n");
        //get the communicators for our groups
        MPI_Comm_create(MPI_COMM_WORLD, cgroup, &ccomm);

        trace("Getting rank of size\n");
        MPI_Comm_rank(ccomm, &crank);
        MPI_Comm_size(ccomm, &csize);
    }
    else
    {
        trace("Dummy call for compute node group creation\n");
        //collective routine.  Our IO nodes need to at least check in!
        MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, &ccomm);
    }


    if(io_node)
    {
        debug("Create group for parallel IO\n");
        int iotrip[1][3];

        iotrip[0][0] = vdiv * hdiv;                     //start index of our row
        iotrip[0][1] = iotrip[0][0] + n_io_nodes - 1;   //final index of our row
        iotrip[0][2] = 1;                               //stride between elements of our row
        trace("Tripplet for group creation: %d %d %d\n", iotrip[0][0], iotrip[0][1], iotrip[0][2]);

        //create groups for our horizontal and vertical associations
        trace("Creating Group\n");
        MPI_Group fgroup;
        MPI_Group_range_incl(global, 1, iotrip, &fgroup);

        trace("Getting communicator\n");
        //get the communicators for our groups
        MPI_Comm_create(MPI_COMM_WORLD, fgroup, &fcomm);

        trace("Getting rank and size\n");
        MPI_Comm_rank(fcomm, &frank);
        MPI_Comm_size(fcomm, &fsize);

        my_io_layer = frank;
        info("We belong to IO layer %d\n", my_io_layer);
       
    }
    else
    {
        trace("Dummy call for parallel IO group creation\n");
        //again the compute nodes need to check in for the collective operation
        MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, &fcomm);
    }

    if(compute_node || io_node)
    {
        debug("Creating group for moving data to IO nodes\n");
        //create the comm group for consolidating data to IO nodes
        int iotrip[2][3];

        //We want the IO node to be root in this group.
        iotrip[0][0] = vdiv * hdiv + my_io_layer;
        iotrip[0][1] = vdiv * hdiv + my_io_layer;
        iotrip[0][2] = 1;

        //Compute nodes next.
        iotrip[1][0] = hdiv * io_layers[my_io_layer].min;   //start index of our row
        iotrip[1][1] = hdiv * io_layers[my_io_layer].max + hdiv - 1;   //final index of our row
        iotrip[1][2] = 1;                                    //stride between elements of our row
        trace("Tripplets for group: %d %d %d %d %d %d\n", iotrip[0][0], iotrip[0][1], iotrip[0][2], iotrip[1][0], iotrip[1][1], iotrip[1][2])

        trace("Creating group\n");
        MPI_Group iogroup;
        MPI_Group_range_incl(global, 2, iotrip, &iogroup);

        trace("Creating Comm\n");
        //get the communicators for our groups
        MPI_Comm_create(MPI_COMM_WORLD, iogroup, &iocomm);

        trace("Getting rank and size\n");
        MPI_Comm_rank(iocomm, &iorank);
        MPI_Comm_size(iocomm, &iosize);
    }
    else
    {
        trace("Dummy call for parallel IO group creation\n");
        //again the compute nodes need to check in for the collective operation
        MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, &iocomm);
    }

    info("Communication Groups Done\n");

}

/*
 * Here we divide the computational domain among processors.  We define the
 * mapping of our 2D logical grid of compute nodes to our 3D grid of data, so
 * that each processor knows how it's local data arrays fit into the full 3D
 * array
 * 
 * For all of these, we will give each processor a roughly equal division of
 * contiguous data layers, and any dimensions that do not divide equally will
 * distribute the remainder among the first layers.
 * 
 * Note:  Spatial arrays are stored as [z][x][y] with y in processor. 
 *        Spectral arrays are stored as [kx][ky][kz] with kz in processor
 */
void lab_initDistributions()
{
    info("Initializing Work Distributions\n");
    int i;
    int d,r;

    //These all follow the same pattern so I'll only comment the first section
    
    //Find how many x go in each layer
    d = nx / hdiv;
    r = nx % hdiv;
    debug("Dividing x: d r = %d %d\n", d,r);
    all_x = (indexes*)malloc((hdiv+1) * sizeof(indexes)); 
    max_x = all_x + hdiv;               //the max is stored in the final element
    max_x->width = (r == 0 ? d : d+1);  //we only care about width for the max
    //set widths for each layer.  First layers get any extra that exist.
    for(i = 0; i < r; i++)
    {
        all_x[i].width = d+1;
        trace("all_x[%d].width = %d\n", i, all_x[i].width);
    }
    for(i = r; i < hdiv; i++)
    {
        all_x[i].width = d;
        trace("all_x[%d].width = %d\n", i, all_x[i].width);
    }
    //calculate min and max index for each layer in the logical global array.
    all_x[0].min = 0;
    all_x[0].max = all_x[0].width-1;
    trace("all_x[%d].min = %d\n", 0, all_x[0].min);
    trace("all_x[%d].max = %d\n", 0, all_x[0].max);
    for(i = 1; i < hdiv; i++)
    {
        all_x[i].min = all_x[i-1].max + 1;
        all_x[i].max = all_x[i-1].max + all_x[i].width;
        trace("all_x[%d].min = %d\n", i, all_x[i].min);
        trace("all_x[%d].max = %d\n", i, all_x[i].max);
    }

    d = nz / vdiv;
    r = nz % vdiv;
    debug("Dividing z: d r = %d %d\n", d,r);
    all_z = (indexes*)malloc((vdiv+1) * sizeof(indexes));
    max_z = all_z + vdiv;
    max_z->width = (r == 0 ? d : d+1);
    for(i = 0; i < r; i++)
    {
        all_z[i].width = d+1;
        trace("all_z[%d].width = %d\n", i, all_z[i].width);
    }
    for(i = r; i < vdiv; i++)
    {
        all_z[i].width = d;
        trace("all_z[%d].width = %d\n", i, all_z[i].width);
    }
    all_z[0].min = 0;
    all_z[0].max = all_z[0].width-1;
    trace("all_z[%d].min = %d\n", 0, all_z[0].min);
    trace("all_z[%d].max = %d\n", 0, all_z[0].max);
    for(i = 1; i < vdiv; i++)
    {
        all_z[i].min = all_z[i-1].max + 1;
        all_z[i].max = all_z[i-1].max + all_z[i].width;
        trace("all_z[%d].min = %d\n", i, all_z[i].min);
        trace("all_z[%d].max = %d\n", i, all_z[i].max);
    }

    d = ndky / hdiv;
    r = ndky % hdiv;
    debug("Dividing ky: d r = %d %d\n", d,r);
    all_ky = (indexes*)malloc((hdiv+1) * sizeof(indexes));
    max_ky = all_ky + hdiv;
    max_ky->width = (r == 0 ? d : d+1);
    for(i = 0; i < r; i++)
    {
        all_ky[i].width = d+1;
        trace("all_ky[%d].width = %d\n", i, all_ky[i].width);
    }
    for(i = r; i < hdiv; i++)
    {
        all_ky[i].width = d;
        trace("all_ky[%d].width = %d\n", i, all_ky[i].width);
    }
    all_ky[0].min = 0;
    all_ky[0].max = all_ky[0].width-1;
    trace("all_ky[%d].min = %d\n", 0, all_ky[0].min);
    trace("all_ky[%d].max = %d\n", 0, all_ky[0].max);
    for(i = 1; i < hdiv; i++)
    {
        all_ky[i].min = all_ky[i-1].max + 1;
        all_ky[i].max = all_ky[i-1].max + all_ky[i].width;
        trace("all_ky[%d].min = %d\n", i, all_ky[i].min);
        trace("all_ky[%d].max = %d\n", i, all_ky[i].max);
    }

    d = ndkx / vdiv;
    r = ndkx % vdiv;
    debug("Dividing kx: d r = %d %d\n", d,r);
    all_kx = (indexes*)malloc((vdiv+1) * sizeof(indexes));
    max_kx = all_kx + vdiv;
    max_kx->width = (r == 0 ? d : d+1);
    for(i = 0; i < r; i++)
    {
        all_kx[i].width = d+1;
        trace("all_kx[%d].width = %d\n", i, all_kx[i].width);
    }
    for(i = r; i < vdiv; i++)
    {
        all_kx[i].width = d;
        trace("all_kx[%d].width = %d\n", i, all_kx[i].width);
    }
    all_kx[0].min = 0;
    all_kx[0].max = all_kx[0].width-1;
    trace("all_kx[%d].min = %d\n", 0, all_kx[0].min);
    trace("all_kx[%d].max = %d\n", 0, all_kx[0].max);
    for(i = 1; i < vdiv; i++)
    {
        all_kx[i].min = all_kx[i-1].max + 1;
        all_kx[i].max = all_kx[i-1].max + all_kx[i].width;
        trace("all_kx[%d].min = %d\n", i, all_kx[i].min);
        trace("all_kx[%d].max = %d\n", i, all_kx[i].max);
    }

    //Set up an alias for our processors personal information inside the full
    //array
    if(compute_node)
    {
        my_x = all_x + hrank;
        my_z = all_z + vrank;
        my_ky = all_ky + hrank;
        my_kx = all_kx + vrank;
    }

    //TODO:  Find a better place for this.  It REALLY doesn't belong here...
    //Find number of z layers in a current IO communcate group.  Note that a 
    //given io group may contain multiple rows of compute nodes.
    nz_layers = 0;
    for(i = io_layers[my_io_layer].min; i <= io_layers[my_io_layer].max; i++)
        nz_layers += all_z[i].width;

    //set up convenience variables telling us how large our main local arrays
    //will be.
    if(compute_node)
    {
        spatialCount = my_x->width * my_z->width * ny;
        spectralCount = my_kx->width * my_ky->width * ndkz;
        debug("Size of spatial array: %d\n", spatialCount);
        debug("Size of spectral array: %d\n", spectralCount);
    }

    info("Work Distrubution Done\n");
}

/*
 * Take out the trash!
 * 
 * Calling this before computations are finished will of course make everything
 * die a fiery death.
 */
void lab_finalize()
{
    free(all_x);
    free(all_z);
    free(all_kx);
    free(all_ky);
    free(io_layers);
}

