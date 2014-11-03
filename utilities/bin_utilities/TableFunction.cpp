// Ver. 1.0.1
// Once a mapping table is loaded, use (or overload) "map" and "invert" to get
// function values and their inverses.
// The mapping table can be loaded from a file using "loadMappingTableFromFile"
// function or manually constructed using "setMappingTable" function.
// The X and Y vector can be returned using getX and getY funtions.

#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include "arsenal.h"
#include "TableFunction.h"

using namespace std;

//----------------------------------------------------------------------
// Constructors:
//----------------------------------------------------------------------

//----------------------------------------------------------------------
TableFunction::TableFunction() {
  mappingTable = NULL;
  deleteMappingTable();
  interpolation_model = 6;
};

//----------------------------------------------------------------------
TableFunction::TableFunction(string filename) {
  mappingTable = NULL;
  loadMappingTableFromFile(filename);
  interpolation_model = 6;
};

//----------------------------------------------------------------------
TableFunction::TableFunction(Table& tab) {
  mappingTable = new Table(tab);
  interpolation_model = 6;
};

//----------------------------------------------------------------------
TableFunction::~TableFunction() { deleteMappingTable(); };

//----------------------------------------------------------------------
bool TableFunction::isMappingTableLoaded() {return mappingTable;};

//----------------------------------------------------------------------
void TableFunction::deleteMappingTable()
// Release memory and set variables to initial values.
{
  if (isMappingTableLoaded()) delete mappingTable;
  mappingTable = NULL;
}


//----------------------------------------------------------------------
void TableFunction::loadMappingTableFromFile(string mappingTable_filename)
// The TableFunction mappingTable file (mappingTable_filename) is assumed to be a n-column file:
{
  // delete old mappingTable first
  deleteMappingTable();
  // get new ones
  mappingTable = new Table;
  mappingTable->loadTableFromFile(mappingTable_filename);
};



//----------------------------------------------------------------------
void TableFunction::setMappingTable(long i, double x, double y)
// The pair (x,y) is set to the i-th row of the table. i starts from 1.
{
  if(!isMappingTableLoaded()) mappingTable = new Table;
  mappingTable->set(1,i,x);
  mappingTable->set(2,i,y);
}


//----------------------------------------------------------------------
void TableFunction::resetMappingTable(long i, double defaultValue)
// Reset the table to have i rows, and all filled with defaultValue.
{
  deleteMappingTable();
  mappingTable = new Table(2,i,defaultValue);
}

//----------------------------------------------------------------------
double TableFunction::map(double x)
// Return the value of function at x
// -- x: input value
// -- mode: used in interpolation. 1: linear direct, 2: linear mono
//   5: cubic direct, 6: cubic mono
{
  if (isMappingTableLoaded())
  {
    return mappingTable->interp(1,2,x,interpolation_model);
  }
  else
  {
    cout << "TableFunction::map error: no mappingTable loaded." << endl;
  }
  return 0;
}


//----------------------------------------------------------------------
double TableFunction::invert(double y)
// Return the inverse value of y
//
{
  if (isMappingTableLoaded())
  {
    return mappingTable->invert(1,2,y,interpolation_model);
  }
  else
  {
    cout << "TableFunction::invert error: no mappingTable loaded." << endl;
  }
  return 0;
}


//----------------------------------------------------------------------
// Return X and Y vectors
//----------------------------------------------------------------------
vector<double>* TableFunction::getX() {return mappingTable->getColumn(1);}
vector<double>* TableFunction::getY() {return mappingTable->getColumn(2);}


//----------------------------------------------------------------------
void TableFunction::printFunction(ostream& os) { mappingTable->printTable(os); }


//----------------------------------------------------------------------
// Return info
double TableFunction::getXMin() { return mappingTable->getFirst(1); }
double TableFunction::getXMax() { return mappingTable->getLast(1); }
double TableFunction::getYMin() { return mappingTable->getFirst(2); }
double TableFunction::getYMax() { return mappingTable->getLast(2); }
long TableFunction::getNumberOfRows() { return mappingTable->getNumberOfRows(); }
