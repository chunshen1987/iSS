// Ver 1.1

#include <vector>
#include <ostream>
#include <string>
#include <iostream>
#include "Table.h"

#ifndef TableFunction_h
#define TableFunction_h


class TableFunction
{
  private:
    Table* mappingTable;
  public:
    int interpolation_model; // controls which model to use in interpolation
    TableFunction();
    TableFunction(string filename);
    TableFunction(Table&);
    ~TableFunction();
    bool isMappingTableLoaded();
    void deleteMappingTable();
    void loadMappingTableFromFile(string);
    void setMappingTable(long, double, double);
    void resetMappingTable(long, double defaultValue=0.0);
    double map(double);
    double invert(double);
    std::vector<double>* getX();
    std::vector<double>* getY();
    double getXMin();
    double getXMax();
    double getYMin();
    double getYMax();
    long getNumberOfRows();
    void printFunction(ostream& os=std::cout);
};

#endif

/*-----------------------------------------------------------------------
 Change logs:

 02-07-2012:
 -- Ver 1.0:
    First version.
 03-15-2012:
 -- Ver 1.1:
    Several changes for convenience: char* -> string; printFunction
    now accepts a stream argument; resetMappingTable "reset" the table.
-----------------------------------------------------------------------*/
