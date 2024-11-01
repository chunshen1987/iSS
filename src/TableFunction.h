// Ver 1.1.1

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include "Table.h"

#ifndef TableFunction_h
#define TableFunction_h

class TableFunction {
  private:
    Table* mappingTable;

  public:
    int interpolation_model;  // controls which model to use in interpolation
    TableFunction();
    TableFunction(std::string filename);
    TableFunction(Table&);
    ~TableFunction();
    bool isMappingTableLoaded();
    void deleteMappingTable();
    void loadMappingTableFromFile(std::string);
    void setMappingTable(long, double, double);
    void resetMappingTable(long, double defaultValue = 0.0);
    double map(double);
    double invert(double);
    std::vector<double>* getX();
    std::vector<double>* getY();
    double getXMin();
    double getXMax();
    double getYMin();
    double getYMax();
    long getNumberOfRows();
    void printFunction(std::ostream& os = std::cout);
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
 04-24-2012:
 -- Ver 1.1.1:
    The default parameter for interpolation is changed to 5 from 6 for
    better compatibility.
-----------------------------------------------------------------------*/
