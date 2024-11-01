// Ver 1.6.0

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#ifndef Table_h
#define Table_h

double zq_global_invert_hook(double);
class Table {
  private:
    std::vector<std::vector<double>*>* data;
    long numberOfCols, numberOfRows;

  public:
    Table();
    Table(std::string);
    Table(Table&);
    Table(long, long, double defaultValue = 0.0);
    Table(double**, long, long);
    ~Table();
    void deleteTable();
    void loadTableFromFile(std::string);
    void loadTableFromDoubleArray(double**, long, long);
    void extendTable(long, long, double defaultValue = 0);
    double get(long, long);
    void set(long, long, double, double defaultValue = 0);
    std::vector<double>* getColumn(long);
    void setAll(double defaultValue = 0.0);
    long getNumberOfRows();
    long getNumberOfCols();
    long getSizeDim1() { return getNumberOfCols(); };
    long getSizeDim2() { return getNumberOfRows(); };
    void printTable(std::ostream& os = std::cout);
    double getFirst(long);
    double getLast(long);
    double interp(
        long, long, double, int mode = 6, bool allowExtrapolation = false);
    double invert(
        long, long, double, int mode = 6, bool allowExtrapolation = false);
    double interp2(double, double, int mode = 1);
};

#endif

/*-----------------------------------------------------------------------
 Change logs:

 02-03-2012:
    First version. Note that the macro INIT_GUESS needs to be customized
    for different tables for efficiency.
 02-04-2012:
    The line "zq_global_table = this;" is moved to Table::invertMono
    to support multi-instances of the Table class.
 02-06-2012:
 -- Ver 1.1:
    Several functions added: set, get, extendTable, getNumberOfXXXs, printTable.
    Note that the "set" function will extend the table if necessary.
    Several interpolation functions have been combined into one "interp"
    function, which accepts a mode parameter to determine the actual
    method for performing interpolation. So is the "invert" function.
03-14-2012:
 -- Ver 1.1.1:
    Functions getSizeDim1/2, printTable(ostream), Table(long,long), setAll
    added for convenience.
    Declarations using char* are replaced by string.
 -- Ver 1.3:
    Self copy "Table(Table &)" enabled.
    Added a check for indices in "set" function for efficiency.
 03-21-2012:
 -- Ver 1.3.1:
    The extendTable function now calls resize() function from vector
    class for efficiency.
 04-11-2012:
 -- Ver 1.4.0:
    Functions added: interp2.
 07-30-2012:
 -- Ver 1.4.1:
    Bug fix: in interp2 the error criteria should be (col_idx_int>=numberOfCols
|| col_idx_int<1 || row_idx_int>=numberOfRows || row_idx_int<1) instead of
(col_idx_int>numberOfCols || col_idx_int<1 || row_idx_int>numberOfRows ||
row_idx_int<1). Also boundary safty control added. 08-01-2012:
 -- Ver 1.5.0:
    A bilinear extroplation method is added to the interp2 function.
 08-02-2012:
 -- Ver 1.6.0:
    Extrapolations for interp and invert are now supported.
    A cubic interpolation / extroplation method is added to the interp2
function.
-----------------------------------------------------------------------*/
