/***********************************************************************
ParameterReader class is used to simplify the process of reading parameters from
an input file and/or from the command line. Version 1.02 (03-14-2012) Zhi Qiu
***********************************************************************/

#ifndef _ParameterReaderHeader
#define _ParameterReaderHeader

#include <string>
#include <vector>

using std::string;

class ParameterReader {
  private:
    std::vector<string>* names;
    std::vector<double>* values;  // store all parameter names and values

    // all substring after "symbol" in "str" will be removed
    string _removeComments(string str, string commentSymbol);

    // phrase an equation like "x=1", assume string has no comments
    void _phraseEquationWithoutComments(string equation);

    // give the index of parameter with "name", or -1 if it does not exist
    long _find(std::string name);

  public:
    ParameterReader();
    ~ParameterReader();

    // read and phrase one setting string like "x=1"
    void phraseOneLine(
        string str, string commentSymbol = static_cast<string>("#"));

    // read in parameters from a file
    void readFromFile(
        string filename, string commentSymbol = static_cast<string>("#"));

    // read in parameter from argument list.
    // The process starts with index="start_from".
    void readFromArguments(
        long argc, char* argv[],
        string commentSymbol = static_cast<string>("#"), long start_from = 1);

    // check if parameter with "name" exists
    bool exist(string name);

    // set the parameter with "name" to value "value"
    void setVal(string name, double value);

    // return the value for parameter with "name"
    double getVal(string name);
    double getVal(string name, double defaultValue);

    // print out all parameters to the screen
    void echo();
};

#endif

/***********************************************************************
Changelog:
09-20-2011: Ver1.01
 -- Bug fix: If the parameter file that is passed to the readFromFile
    function does not exist, the program stops instead of going into
    infinite loops.
03-13-2012: Ver1.02
 -- Bug fix: The commentSymbol parameter was not passed to the
    phraseOneLine function in readFromFile.
***********************************************************************/
