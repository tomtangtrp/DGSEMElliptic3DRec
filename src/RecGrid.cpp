#include <iostream>
#include "RecGrid.hpp"
#include "Connection.hpp"
#include "ConnectionTable.hpp"

// Custom Constructor on [0,1]x[0,1]
RecGrid::RecGrid(int Nel_x, int Nel_y, int Nel_z)
{
    axL = 0.0;
    bxL = 1.0;
    ayL = 0.0;
    byL = 1.0;
    azL = 0.0;
    bzL = 1.0;

    nx=Nel_x;
    ny=Nel_y;
    nz=Nel_z;

}



//Custom Constructor on [ax, bx]x[ay, by]
RecGrid::RecGrid(double ax, double bx, double ay, double by, double az, double bz, int Nel_x, int Nel_y, int Nel_z)
{
    axL = ax;
    bxL = bx;
    ayL = ay;
    byL = by;
    azL = az;
    bzL = bz;

    nx=Nel_x;
    ny=Nel_y;
    nz=Nel_z;
}

ConnectionTable RecGrid::getRecConnectionTable() 
{
    ConnectionTable RecConnectionTable;

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iz = 0; iz < nz; ++iz) {
                int e1 = eid_from_ixiyiz(ix, iy, iz);

                // +z neighbor (FRONT - BACK)
                if (iz < nz - 1) {
                    int e2 = eid_from_ixiyiz(ix, iy, iz + 1);
                    Connection connect;
                    connect.ElmConnect = std::make_tuple(e1, e2);
                    connect.FaceConnect = std::make_tuple(1, 0); // FRONT (1), BACK (0)
                    RecConnectionTable.add_connection(connect);
                }

                // +x neighbor (RIGHT - LEFT)
                if (ix < nx - 1) {
                    int e2 = eid_from_ixiyiz(ix + 1, iy, iz);
                    Connection connect;
                    connect.ElmConnect = std::make_tuple(e1, e2);
                    connect.FaceConnect = std::make_tuple(3, 2); // RIGHT (3), LEFT (2)
                    RecConnectionTable.add_connection(connect);
                }

                // +y neighbor (TOP - BOTTOM)
                if (iy < ny - 1) {
                    int e2 = eid_from_ixiyiz(ix, iy + 1, iz);
                    Connection connect;
                    connect.ElmConnect = std::make_tuple(e1, e2);
                    connect.FaceConnect = std::make_tuple(5, 4); // TOP (5), BOTTOM (4)
                    RecConnectionTable.add_connection(connect);
                }
            }
        }
    }

    return RecConnectionTable;
}

