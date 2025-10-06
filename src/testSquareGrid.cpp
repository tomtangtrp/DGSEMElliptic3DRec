#include <iostream>
#include "RecGrid.hpp"
#include "ConnectionTable.hpp"

int main(int argc, char* argv[])
{
    int Nel_x = 2;
    int Nel_y = 2;
    int Nel_z = 2;
    RecGrid mRecGrid(Nel_x, Nel_y, Nel_z);
    ConnectionTable mSqaureConnectionTable = mRecGrid.getRecConnectionTable();
    mSqaureConnectionTable.printConnectionTable();

    return 0;
}