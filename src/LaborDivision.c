/* 
 * File:   LaborDivision.cpp
 * Author: Ben
 * 
 * Created on March 8, 2010, 3:57 PM
 */

#include <stdlib.h>
#include <stdio.h>

#include "Environment.h"
#include "LaborDivision.h"
#include "Log.h"

void lab_initGeometry()
{
    info("Initializing Problem Geometry\n");

    dx = xmx / nx;
    dy = ymx / ny;
    dz = zmx / nz;
    trace("dx dy dz: %g %g %g\n", dx, dy, dz);

    nkx = nx;
    nky = ny/2 + 1;
    nkz = nz;
    trace("nkx nky nkz: %d %d %d\n", nkx, nky, nkz);

    dealias_kx.min = (nx/2) * 2 / 3 + 1;
    dealias_ky.min = (ny/2) * 2 / 3 + 1;
    dealias_kz.min = (nz/2) * 2 / 3 + 1;
    trace("Smallest dialised kx ky kz: %d %d %d\n", dealias_kx.min, dealias_ky.min, dealias_kz.min);

    dealias_kx.max = nkx  - dealias_kx.min;
    dealias_ky.max = nky - 1;      //This is the first dimension transposed
    dealias_kz.max = nkz  - dealias_kz.min;
    trace("Largest dialised kx ky kz: %d %d %d\n", dealias_kx.max, dealias_ky.max, dealias_kz.max);

    dealias_kx.width = dealias_kx.max - dealias_kx.min + 1;
    dealias_ky.width = dealias_ky.max - dealias_ky.min + 1;
    dealias_kz.width = dealias_kz.max - dealias_kz.min + 1;
    trace("Width of dealiasing kx ky kz: %d %d %d\n", dealias_kx.width, dealias_ky.width, dealias_kz.width)

    ndkx = nkx - dealias_kx.width;
    ndky = nky - dealias_ky.width;
    ndkz = nkz - dealias_kz.width;
    trace("Number of retained waves %d %d %d\n", ndkx, ndky, ndkz);

    info("Problem Geometry Done\n");
}

void lab_initGroups()
{
    info("Initializing Communication Groups\n");

    int i;

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
    else
    {
        compute_node = 0;
        io_node = 0;
    }

    info("compute_node flag set to %d\n", compute_node);
    info("io_node flag set to %d\n", io_node);

    //grab the global group
    MPI_Group global;
    MPI_Comm_group(MPI_COMM_WORLD, &global);

    
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

        iotrip[0][0] = vdiv * hdiv + my_io_layer;
        iotrip[0][1] = vdiv * hdiv + my_io_layer;
        iotrip[0][2] = 1;

        iotrip[1][0] = hdiv * io_layers[my_io_layer].min;   //start index of our row
        iotrip[1][1] = hdiv * io_layers[my_io_layer].max + hdiv - 1;   //final index of our row
        iotrip[1][2] = 1;                                    //stride between elements of our row
        trace("Tripplets for group: %d %d %d %d %d %d\n", iotrip[0][0], iotrip[0][1], iotrip[0][2], iotrip[1][0], iotrip[1][1], iotrip[1][2])

        trace("Creating group\n");
        //create groups for our horizontal and vertical associations
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

void lab_initDistributions()
{
    info("Initializing Work Distributions\n");
    int i;
    int d,r;

    d = nx / hdiv;
    r = nx % hdiv;
    debug("Dividing x: d r = %d %d\n", d,r);
    all_x = (indexes*)malloc((hdiv+1) * sizeof(indexes));
    max_x = all_x + hdiv;
    max_x->width = (r == 0 ? d : d+1);
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

    if(compute_node)
    {
        my_x = all_x + hrank;
        my_z = all_z + vrank;
        my_ky = all_ky + hrank;
        my_kx = all_kx + vrank;
    }

    //TODO:  Find a better place for this.  It REALLY doesn't belong here...
    nz_layers = 0;
    for(i = io_layers[my_io_layer].min; i <= io_layers[my_io_layer].max; i++)
        nz_layers += all_z[i].width;

    if(compute_node)
    {
        spatialCount = my_x->width * my_z->width * ny;
        spectralCount = my_kx->width * my_ky->width * ndkz;
        debug("Size of spatial array: %d\n", spatialCount);
        debug("Size of spectral array: %d\n", spectralCount);
    }

    info("Work Distrubution Done\n");
}

void lab_finalize()
{
    free(all_x);
    free(all_z);
    free(all_kx);
    free(all_ky);
    free(io_layers);
}

