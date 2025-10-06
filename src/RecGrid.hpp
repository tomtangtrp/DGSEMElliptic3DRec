#ifndef RECGRIDHEADERDEF
#define RECGRIDHEADERDEF

#include "Connection.hpp"
#include "ConnectionTable.hpp"

class RecGrid
{
    private:
        // SquareGrid of [axL, bxL] x [ayL, byL] x [azL, bzL]
        double axL;
        double bxL;
        double ayL;
        double byL;
        double azL;
        double bzL;

        // 
        int nx; // Nel_x = nCol;
        int ny; // Nel_y = nRow;
        int nz; // Nel_y = nRow;
        // Connection table for the structured square grid
        ConnectionTable RecConnectionTable; 
    public:
        // Custom constructor for [0,1]x[0,1]
        RecGrid(int Nel_x, int Nel_y, int Nel_z);
        // Custom constructor
        RecGrid(double ax, double bx, double ay, double by, double az, double bz, int Nel_x, int Nel_y, int Nel_z);
        int eid_from_ixiyiz(int ix, int iy, int iz) const {return iz + nz * (ix + nx * iy);}
    
        ConnectionTable getRecConnectionTable();
    

};
#endif
