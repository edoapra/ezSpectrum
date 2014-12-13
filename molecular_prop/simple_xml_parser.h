#ifndef _simple_xml_parser_h
#define _simple_xml_parser_h

/*! \file simple_xml_parser.h
\brief A simple class for SAX XML parsing using eXpat library   (vm 2006-09)
*/

#include <expat.h>

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>


class simpleXMLparser
{
 public:

  simpleXMLparser();
  simpleXMLparser(const simpleXMLparser& other);
  simpleXMLparser& operator= (const simpleXMLparser& other);
  ~simpleXMLparser();
  //! assign file
  void assignFile(const char* fileName); // dirty game to assign the file not when class is declared
  //! copies <input.xml> to std::cout
  void printInputFile();

  //-------------------------------------------------------------------------
  //---- Handlers to set up a path through the tree-structure of the xml file:

  //! cleans the "path":
  simpleXMLparser& reset();

  //! adds one node to "the path" 
  simpleXMLparser& node(const char* elementName, int elementNumber);
  //! adds one node to "the path"
  simpleXMLparser& node(const char* elementName) { return node(elementName,1); }

  //! removes 'nSteps' nodes from the end of the current "path" 
  simpleXMLparser& stepBack(int nSteps);
  //! removes one node from the end of the current "path" 
  simpleXMLparser& stepBack() { return stepBack(1); };

  //! Checks if current path exists:
  bool Check();
  bool CheckSubNode(const char* subNodeName)
    {bool returnBool; node(subNodeName); returnBool=Check(); stepBack(); return returnBool;};

  //! Get string value of the current node
  const char* getTextValue();
  //! Get string value of the attrName attribute of the current node
  const char* getAttrValue(const char* attrName);
  //! Get string value of the current node
  const char* value() {return getTextValue(); } 
  //! Get string value of the attrName attribute of the current node
  const char* value(const char* attrName) { return getAttrValue(attrName); } 

  //---- Get numeric values:
  //! Get double numeric value of the current node
  double getDoubleValue();
  //! Get double numeric value of the attrName attribute of the current node
  double getDoubleValue(const char* attrName);
  //! Get double numeric value of the subnode ("subNodeName"). But keep the position in the node-tree unchanged;
  double getSubNodeDbl(const char* subNodeName)
           {double returnDbl=node(subNodeName).getDoubleValue(); stepBack(); return returnDbl;};
  //! Get int numeric value of the current node
  int getIntValue();
  //! Get int numeric value of the attrName attribute of the current node
  int getIntValue(const char* attrName);
  //! Get bool value of the current node
  bool getBoolValue();
  //! Get bool value of the attrName attribute of the current node
  bool getBoolValue(const char* attrName);
  //! Get first word (skip non-letters/numbers) of the current node
  void getWordValue(std::string& returnWord);
  //! Get first word (skip non-letters/numbers) of the attrName attribute of the current node
  void getWordValue(std::string& returnWord, const char* attrName);

  //! Get int numeric value of the subnode ("subNodeName"). But keep the position in the node-tree unchanged;
  int getSubNodeInt(const char* subNodeName)
           {int returnInt=node(subNodeName).getIntValue(); stepBack(); return returnInt;};
  

  //-------------------------------------------------------------------------
  // Error "handlers":

  //! terminate program if true; prints out an error message:
  void exitOnFormatError(bool ifExit); 

  //! returns current path:
  const char* getPath();


  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

 protected:
  //! .xml file name for parsing
  const char* xmlFileName;             

 private:
  bool ifFileAssigned;
  
  //! Parser itself (expat)
  XML_Parser XMLp;                   
  
  //! Struct: tag name and number(among the same tags at the particular depth)
  struct XMLlevel{                   
    std::string XMLtagName;                              
    int XMLtagNumber;
  };
  
  //! current tag for search;
  XMLlevel searchTag;                 
  //! previous found tag (at previous XML depth);
  XMLlevel prevTag;           
  
  //! 'path' list to particular tag  formed from ("name1",number1)("name2",number2)... by "index" operator();
  std::list<XMLlevel> originalXMLpathList;    
  std::list<XMLlevel> tXMLpathList;    
 
  //! If 'path' list was processed
  bool ifTriedToPosition;             

  //! current text value -- updated every time by "text-handler"
  std::string currentTextStr;         

  //! return value -- updated only when the correct position is fund; return "NULL" if nothing found
  std::string returnTextStr;          

  int currentTagNumber;

  //! current Depth in xml (+1 for each opening tag, -1 for each closing)
  int Depth;                          
  //! depth of previous oppening tag from the List
  int prevDepth;                      

  //! value of the depth where a tag should located
  int searchDepth;                    

  //! Flag -- ifStop -- stop searching through xml
  bool ifStop;                       

  struct tagAttribute{                
    std::string attrName;                              
    std::string attrValueStr;
  };

  tagAttribute tagAttrTmp;

  //!  All attributes of the current tag
  std::list<tagAttribute> currentAttr;
  //!  All attributes of the returned tag
  std::list<tagAttribute> returnAttr; 


  
  //! remove this:
  #define BUFFSIZE   8192             

  //! read from .xml file by BUFFSIZE portions
  char Buff[BUFFSIZE];                

  std::ostringstream errS; //

  std::string currAtrName;
  
  std::string tmpS;
  std::ostringstream tmpStr;


  ///////////////////////////////////////////////////////
  // Headers of the Handlers for 'eXpat' XML parser:
  ///////////////////////////////////////////////////////
  
  //! Handler for an opening tag:
  static void XMLCALL openingHndl (void *data, const char *el, const char **attr);
  //! Handler for a closing tag:
  static void XMLCALL closingHndl(void *data, const char *el);
  //! Handler for text in between an open and closing tag:
  static void XMLCALL textHndl(void *data, const char *str, const int lng);


  /*! Resets the parser and searches for the tag described by 'path' 
   list (XMLpathList); when found, handlers (described above) put 
   value of the Text to "returnTextStr" and Atrributes to "returnAttr":
  */

  void CherchezLaTag();

  const char* formString(int iNumber);

};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

class My_istringstream
{

  std::istringstream iStr;
  bool ifLetterOrNumber(char Ch)
    {
      bool return_bool;

      if ( ((int(Ch)<=int('Z')) and (int(Ch)>=int('A')))  
	   or ((int(Ch)<=int('z')) and (int(Ch)>=int('a'))) 
	   or ((int(Ch)<=int('9')) and (int(Ch)>=int('0')))  )
	return_bool=true;
      else 
	return_bool=false;
      
  return return_bool;
    };

 public:
  My_istringstream(const char* str){ iStr.str(str); iStr.clear(); };
  My_istringstream(const My_istringstream& other){ iStr.str(other.iStr.str()); iStr.clear(); };

  std::string str(){return iStr.str(); };
  bool fail(){return iStr.fail(); };

  int getNextInt(){int next; iStr>>next; return next; };
  double getNextDouble(){double next; iStr>>next; return next; };

  void getNextWord(std::string& next)
    {
      char tmp_char;

      //skip spaces (non-letters/numbers):
      tmp_char=' ';
      while ( not(ifLetterOrNumber(tmp_char) ) and not(iStr.fail()) )
	iStr.get(tmp_char);  

      //read the word:
      next="";
      if (not(iStr.fail()))
	do {  
	  next+=tmp_char; 
	  iStr.get( tmp_char );  
	} while( ifLetterOrNumber(tmp_char) and not(iStr.fail()) );
    };
};


#endif


