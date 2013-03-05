// Ver. 1.3.1
// The Table data file is a n-column double data file (see loadTableFromFile).
// Note the (column,row) convention used in this class.
// Note that all indices start with 1.

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <vector>
#include <string>
#include "stdlib.h"
#include <cmath>
#include "arsenal.h"
#include "Table.h"

#define INIT_GUESS 1
#define ACCURACY 1e-4

// Used to invert the table
Table *zq_global_table;
long zq_global_colX, zq_global_colY;
int zq_global_mode;

using namespace std;

//----------------------------------------------------------------------
// Constructors:
//----------------------------------------------------------------------

//----------------------------------------------------------------------
Table::Table() { data = NULL; numberOfCols = 0; numberOfRows = 0; };

//----------------------------------------------------------------------
Table::Table(string filename)
// This version load table from file.
{ 
  data = NULL;
  loadTableFromFile(filename); 
}

//----------------------------------------------------------------------
Table::Table(Table& tab_in)
// This version COPY content from tab_in to table.
{
  // get info
  long cols = tab_in.getNumberOfCols();
  long rows = tab_in.getNumberOfRows();
  // allocate memory
  data = NULL; numberOfCols = 0; numberOfRows = 0;
  extendTable(cols, rows);
  // copy
  for (long i=1; i<=numberOfCols; i++)
  for (long j=1; j<=numberOfRows; j++)
    set(i, j, tab_in.get(i,j));
}

//----------------------------------------------------------------------
Table::Table(long i, long j, double defaultValue) 
// This version allocates a table of i-col and j-rows.
{ 
  data = NULL; numberOfCols = 0; numberOfRows = 0;
  extendTable(i, j, defaultValue);
}

//----------------------------------------------------------------------
Table::Table(double** data_in, long size1, long size2) 
// This version COPY content from data_in to the internal table.
// The 1st index of the data_in should go from 0 to size1 and the
// 2nd index of data_in should go from 0 to size2.
{
  data = NULL; numberOfCols = 0; numberOfRows = 0;
  loadTableFromDoubleArray(data_in, size1, size2);
}

//----------------------------------------------------------------------
Table::~Table() { releaseBlockData(data); };


//----------------------------------------------------------------------
void Table::deleteTable()
// Release memory and set variables to initial values.
{
  releaseBlockData(data);
  data = NULL; numberOfCols = 0; numberOfRows = 0;
}

//----------------------------------------------------------------------
void Table::extendTable(long i, long j, double defaultValue)
// Extend the table to i columns and j rows (block), the empty places
// are filled with defaultValue.
// Note that this function uses the internal variables
// numberOfXXXs to check dimension, instead of inquiring it again for
// efficiency.
{
  if (!data) data = new vector< vector<double>* >;
  // extend columns first
  if (i>numberOfCols)
  {
    for (long ii=numberOfCols+1; ii<=i; ii++)
    {
      data->push_back(new vector<double>(numberOfRows,defaultValue));
    }
    numberOfCols = i;
  }
  // extend rows first
  if (j>numberOfRows)
  {
    for (long ii=0; ii<numberOfCols; ii++) (*data)[ii]->resize(j, defaultValue);
    // for (long jj=numberOfRows+1; jj<=j; jj++)
    // {
    //   (*data)[ii]->push_back(defaultValue);
    // }
    numberOfRows = j;
  }
}

//----------------------------------------------------------------------
double Table::get(long i, long j)
// Return the (i-col,j-row) element of the table, starting from (1,1)
{
  if (data) return (*(*data)[i-1])[j-1];
  else
  {
    cout << "Table::get error: no data table loaded." << endl;
    return 0.0;
  }
}

//----------------------------------------------------------------------
void Table::set(long i, long j, double value, double defaultValue)
// Set the (i-col,j-row) element of the table, starting from (1,1) to
// value. defaultValue is used when extending the table.
{
  if (i>numberOfCols || j>numberOfRows) extendTable(i,j,defaultValue);
  (*(*data)[i-1])[j-1] = value;
}


//----------------------------------------------------------------------
vector<double>* Table::getColumn(long colX)
// Return the vector with column index colX (starting from 1)
{
  return (*data)[colX-1];
}


//----------------------------------------------------------------------
void Table::setAll(double defaultValue)
// Set the all element to defaultValue.
{
  long i=numberOfCols, j=numberOfRows;
  deleteTable();
  extendTable(i,j,defaultValue);
}


//----------------------------------------------------------------------
// Return size info
long Table::getNumberOfCols() {return numberOfCols;}
long Table::getNumberOfRows() {return numberOfRows;}


//----------------------------------------------------------------------
void Table::printTable(ostream& os)
{
  for (long j=1; j<=numberOfRows; j++)
  {
    for (long i=1; i<=numberOfCols; i++)
    {
      os << scientific <<  setw(15) << setprecision(8) << get(i,j) << "  ";
    }
    os << endl;
  }
}


//----------------------------------------------------------------------
void Table::loadTableFromFile(string data_filename)
// The Table data file (data_filename) is assumed to be a n-column file:
{
  // delete old data first
  deleteTable();
  // get new ones
  fstream fs(data_filename.c_str());
  if (fs.is_open()==false)
  {
    cout << "Table::loadTableFromFile error: the data file " << data_filename << " cannot be opened." << endl;
    exit(-1);
  }
  data = readBlockData(fs);
  numberOfCols = (*data).size();
  numberOfRows = (*(*data)[0]).size();
  fs.close();
};


//----------------------------------------------------------------------
void Table::loadTableFromDoubleArray(double** data_in, long size1, long size2)
// This version COPY content from data_in to the internal table.
// The 1st index of the data_in should go from 0 to size1 and the
// 2nd index of data_in should go from 0 to size2.
{
  deleteTable();
  extendTable(size1, size2);
  for (long i=0; i<size1; i++)
  for (long j=0; j<size2; j++)
    set(i+1,j+1, data_in[i][j]);
}


//----------------------------------------------------------------------
double Table::getFirst(long col)
// Return the first element in column "col"
{
  return get(col, 1);
}

//----------------------------------------------------------------------
double Table::getLast(long col)
// Return the last element in column "col"
{
  return get(col, numberOfRows);
}


//----------------------------------------------------------------------
// For all the following functions, column index starts with 1
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double Table::interp(long colX, long colY, double xx, int mode)
// Return yy interpolated from xx from the table given by colX and colY.
// The mode parameter is the one controls the way interpolation is done.
// Mode: 1: linear direct, 2: linear mono,
//       5: cubic direct, 6: cubic mono,
//       10: nearest direct, 11: nearest mono
{
  switch (mode)
  {
    case 1:
      return interpLinearDirect((*data)[colX-1], (*data)[colY-1], xx);
      break;
    case 2:
      return interpLinearMono((*data)[colX-1], (*data)[colY-1], xx);
      break;
    case 5:
      return interpCubicDirect((*data)[colX-1], (*data)[colY-1], xx);
      break;
    case 6:
      return interpCubicMono((*data)[colX-1], (*data)[colY-1], xx);
      break;
    case 10:
      return interpNearestDirect((*data)[colX-1], (*data)[colY-1], xx);
      break;
    case 11:
      return interpNearestMono((*data)[colX-1], (*data)[colY-1], xx);
      break;
    default:
      cout << "Table::interp error: mode parameter " << mode << " not recogonized." << endl;
      return 0.0;
  }
}

// used in root search to invert tables
//----------------------------------------------------------------------
double zq_global_invert_hook(double xx) { return zq_global_table->interp(zq_global_colX, zq_global_colY, xx, zq_global_mode); }
//----------------------------------------------------------------------
double Table::invert(long colX, long colY, double yy, int mode)
// Return the x value corresponding to yy from table given by colX and colY.
// This function assumes that data in colX are monotonically increasing.
{
  zq_global_table = this;
  zq_global_colX = colX;
  zq_global_colY = colY;
  zq_global_mode = mode;
  return invertFunc(&zq_global_invert_hook, yy, (*(*data)[colX-1])[0], (*(*data)[colX-1])[numberOfRows-1], (*(*data)[colX-1])[1]-(*(*data)[colX-1])[0], (*(*data)[colX-1])[INIT_GUESS], ACCURACY);
}
