/*! \file simple_xml_parser.C
\ingroup (DATA_CLASSES)
*/

#include "simple_xml_parser.h"

//----------------------
simpleXMLparser::simpleXMLparser()
{
  ifFileAssigned=false;
}

//----------------------
simpleXMLparser::simpleXMLparser(const simpleXMLparser& other)
{
  ifFileAssigned=false;
  assignFile(other.xmlFileName);
}

//----------------------
simpleXMLparser& simpleXMLparser::operator=(const simpleXMLparser& other)
{
  ifFileAssigned=false;
  assignFile(other.xmlFileName);
  return *this;
}


//------------------------
void simpleXMLparser::assignFile(const char* fileName)
{
  xmlFileName = fileName;

  //try to open <input.xml> file
  std::ifstream xmlFile(xmlFileName);
  if  ( xmlFile.fail() ){
    std::cout << "Error: File \"" << xmlFileName << "\" does not exist\n" ;
    exit(-1);
  }
  xmlFile.close();


  ifTriedToPosition=false;
  if (! ifFileAssigned )
    {
      XMLp = XML_ParserCreate(NULL);

      if (! XMLp) {
	std::cout << "\nError: couldn't allocate memory for parser\n\n";
	exit(-1);
      }
    }
  simpleXMLparser::reset();
  errS.clear();
  
}


//------------------------
void simpleXMLparser::printInputFile()
{
  std::ifstream xmlFile(xmlFileName, std::ios::binary);
  std::cout  << xmlFile.rdbuf();
  xmlFile.close();
}


//----------------------
simpleXMLparser::~simpleXMLparser()
{
  XML_ParserFree(XMLp);
}

//----------------------
simpleXMLparser& simpleXMLparser::node(const char* elementName, const int elementNumber)
{
  XMLlevel tmpLevel;

  tmpLevel.XMLtagName = elementName;
  tmpLevel.XMLtagNumber = elementNumber;

  tXMLpathList.push_back(tmpLevel);
  originalXMLpathList=tXMLpathList;
  ifTriedToPosition=false;
  return *this;
}

//----------------------
simpleXMLparser& simpleXMLparser::reset(){
  originalXMLpathList.clear(); 
  tXMLpathList.clear(); 
  ifTriedToPosition=false;
  return *this;
}


//----------------------
bool simpleXMLparser::Check()
{
  returnTextStr="NULL";
  simpleXMLparser::CherchezLaTag();
  ifTriedToPosition=true;
  currAtrName="NULL";

  tXMLpathList=originalXMLpathList;

  if (returnTextStr=="NULL")
    return false;
  else
    return true;
}


//----------------------
const char* simpleXMLparser::getTextValue()
{
  if ( !(ifTriedToPosition) )
    if (this->Check() == false)              //------begin error message
      {
	errS << "\nError: required path does not exist: [" << this->getPath() << "]\n\n";
	std::cout << errS.str();
	exit(-1);                                //------end error message
      }

  return returnTextStr.c_str();
}


//----------------------
const char* simpleXMLparser::getAttrValue(const char* attrName)
{
  if ( !(ifTriedToPosition) )
    if (this->Check() == false)           //------begin error message
      {
	errS << "\nError: required path does not exist: [" << this->getPath() << "]\n\n";
	std::cout << errS.str();
	exit(-1);                             //------end error message
      }

  currAtrName=attrName;

  std::string tmpReturnStr="NULL";
  std::list<tagAttribute> tmpAttr=returnAttr;

  for (int i=0;(tmpAttr.size()>0) && (tmpReturnStr=="NULL");i++)
    {
      tagAttrTmp=tmpAttr.front();

      if (tagAttrTmp.attrName == attrName)
	tmpReturnStr=tagAttrTmp.attrValueStr;

      tmpAttr.pop_front();
    }

  if (tmpReturnStr=="NULL")                    //------begin error message
    {
      errS << "\nError: required path does not exist: [" << this->getPath() << "]\n\n";
      std::cout << errS.str();
      exit(-1);                                //------end error message
    }

  return tmpReturnStr.c_str();
}

/*
const char* simpleXMLparser::getPath()
{
  std::list<XMLlevel> errCopy=originalXMLpathList;// Copy;
  XMLlevel tmpLevel;
  tmpS=xmlFileName;
  
  
  while ( errCopy.size() > 0 )
    {
      tmpLevel=errCopy.front();
      errCopy.pop_front();
      tmpS += "->";
      tmpS += tmpLevel.XMLtagName;
      if ( tmpLevel.XMLtagNumber > 1 )
	{
	  tmpS += '(';
	  tmpS += tmpLevel.XMLtagNumber;
	  tmpS += ')';
	}
    }
  
    if ( currAtrName != "NULL" )
      {
        tmpS += "(attribute=\"";
        tmpS += currAtrName;
        tmpS += "\")";
      }
  
    return tmpS.c_str();
}
*/

const char* simpleXMLparser::getPath()
{
  std::list<XMLlevel> errCopy=originalXMLpathList;// Copy;
  XMLlevel tmpLevel;
  tmpS=xmlFileName;
  
  while ( errCopy.size() > 0 )
    {
      tmpLevel=errCopy.front();
      errCopy.pop_front();
      tmpS += "->";
      tmpS += tmpLevel.XMLtagName;

      if ( tmpLevel.XMLtagNumber > 1 )
	{
	  tmpS += '(';
	  tmpS += formString(tmpLevel.XMLtagNumber);
	  tmpS += ')';
	}
   }

  if ( currAtrName != "NULL" )
    {
      tmpS += "(attribute=\"";
      tmpS += currAtrName;
      tmpS += "\")";
    }
  
  return tmpS.c_str();
}


//----------------------
void simpleXMLparser::exitOnFormatError(bool ifExit)
{
  if (ifExit)
    {
      std::cout << "\nError: wrong format in [" << this->getPath() <<"]\n\n";
      exit(-1);
    }
}
//----------------------

simpleXMLparser& simpleXMLparser::stepBack(int nSteps)
{
  for (int i=0; i<nSteps; i++)
    originalXMLpathList.pop_back();

  tXMLpathList=originalXMLpathList;
  ifTriedToPosition=false;
  return *this;
}


//----------------------
double simpleXMLparser::getDoubleValue()
{
  double returnDouble;
  std::istringstream inpStr;
  inpStr.str(this->value()); 
  inpStr.clear();
  inpStr >> returnDouble;
  this->exitOnFormatError(inpStr.fail());
  return returnDouble;
}



//----------------------
double simpleXMLparser::getDoubleValue(const char* attrName)
{
  double returnDouble;
  std::istringstream inpStr;
  inpStr.str(this->value(attrName)); 
  inpStr.clear();
  inpStr >> returnDouble;
  this->exitOnFormatError(inpStr.fail());
  return returnDouble;
}

//----------------------
int simpleXMLparser::getIntValue()
{
  int returnInt;
  std::istringstream inpStr;
  inpStr.str(this->value()); 
  inpStr.clear();
  inpStr >> returnInt;
  this->exitOnFormatError(inpStr.fail());
  return returnInt;
}

//----------------------
int simpleXMLparser::getIntValue(const char* attrName)
{
  int returnInt;
  std::istringstream inpStr;
  inpStr.str(this->value(attrName)); 
  inpStr.clear();
  inpStr >> returnInt;
  this->exitOnFormatError(inpStr.fail());
  return returnInt;
}

//----------------------
bool simpleXMLparser::getBoolValue()
{
  bool returnBool;
  std::string tmp_str;
  tmp_str=this->value(); 
  if ((tmp_str=="y") or (tmp_str=="true"))
    returnBool = true;
  else
    {
      if ((tmp_str=="n") or (tmp_str=="false"))
	returnBool = false;
      else
	this->exitOnFormatError(true);
    }
  return returnBool;
}

//----------------------
bool simpleXMLparser::getBoolValue(const char* attrName)
{
  bool returnBool;
  std::string tmp_str;
  tmp_str=this->value(attrName); 
  if ((tmp_str=="y") or (tmp_str=="true"))
    returnBool=true;
  else
    {
      if ((tmp_str=="n") or (tmp_str=="false"))
	returnBool=false;
      else
	this->exitOnFormatError(true);
    }
  return returnBool;
}
//----------------------
void simpleXMLparser::getWordValue(std::string& returnWord)
{
  My_istringstream iStr(this->value()); 
  iStr.getNextWord(returnWord);
  exitOnFormatError(iStr.fail());
}

//----------------------
void simpleXMLparser::getWordValue(std::string& returnWord, const char* attrName)
{
  My_istringstream iStr(this->value(attrName)); 
  iStr.getNextWord(returnWord);
  exitOnFormatError(iStr.fail());
}



//----------------------
////////////////////////////////////
// handlers for eXpat XML parser  //
////////////////////////////////////

//----------------------
// Handler for an opening tag:
void XMLCALL 
simpleXMLparser::openingHndl (void *parserData, const char *currentOpeningTagName, const char **currentAttribute)
{ 
  simpleXMLparser* PData = (simpleXMLparser*)parserData;                // PData: pointer to 'this'

  (PData->Depth)++;                                                     // For any opening tag -- we are one level deeper;
 
  if ((PData->searchTag.XMLtagName == currentOpeningTagName)            // Current OPENING tag is what we are searching for
      && ((PData->searchDepth - PData->Depth) == PData->tXMLpathList.size()-1) // at the correct depth
      && (PData->tXMLpathList.size() != 0)                               // not at the end of the list (no more oppening tags required)
      && !(PData->ifStop))                                              // and it is not too late :)  ;
    {
      PData->currentTagNumber++;                                        // One more correct tag at the correct depth;

      if (PData->currentTagNumber == PData->searchTag.XMLtagNumber)     // If there is enough correct tags at this depth
	{
	  PData->prevTag=PData->tXMLpathList.front();                            // store current found tag as a previous
	  PData->prevDepth=PData->Depth;                                // store curent Depth as depth of previuos tag
	  PData->tXMLpathList.pop_front();                               // and delete found tag from the List;

	  if (PData->tXMLpathList.size() != 0)                           // If not the last element in the List
	    {
	      PData->searchTag=PData->tXMLpathList.front();              // take next tag to search from the List
	      PData->currentTagNumber=0;                                // and reset the counter;
	    }
	}

      PData->currentAttr.clear();                                       // ---- begin [Copy current attributes to the currentAttr list]
      for (int i = 0; currentAttribute[i]; i += 2) 
	{
	  PData->tagAttrTmp.attrName=currentAttribute[i];
	  PData->tagAttrTmp.attrValueStr=currentAttribute[i + 1];
	  PData->currentAttr.push_back(PData->tagAttrTmp);             
	}                                                               // ---- end;
    }
}


//----------------------
//Handler for a closing tag:
void XMLCALL 
simpleXMLparser::closingHndl(void *parserData, const char *currentClosingTagName)

{
  simpleXMLparser* PData = (simpleXMLparser*)parserData;               // PData: pointer to 'this'
  
  if ((PData->searchTag.XMLtagName == currentClosingTagName)           // Current CLOSING tag is what we are looking for
      && ( PData->searchTag.XMLtagNumber == PData->currentTagNumber)   // correct number of oppening tags of this type
      && ( PData->searchDepth == PData->Depth)                         // at the correct depth
      && ( PData->tXMLpathList.size()==0)                               // if all elements from the List are found
      && !(PData->ifStop))                                             // and it is not too late :) ;
    {
      PData->returnTextStr=PData->currentTextStr;                      // Returns current text value

      for ( int i=0; PData->currentAttr.size()>0; i++ )                // ---- begin [return attributes of the current tag];
	{
	  PData->returnAttr.push_back(PData->currentAttr.front());     
	  PData->currentAttr.pop_front(); 
	}                                                              // ---- end;

      PData->ifStop=true;                                              // Stop searching (correct tag found);
    }

  if ((PData->prevTag.XMLtagName == currentClosingTagName)             // Current CLOSING tag is the previous from the List
      && (PData->prevDepth == PData->Depth)                            // at the position of the previous one;
      && !(PData->ifStop))                                             // and not yet stoped;
    {
      PData->ifStop=true;                                              // Stop searching (nothing found);
    }


  (PData->Depth)--;                                                    // For any closing tag -- one leve up;
  
}


//----------------------
//Handler for text in between opening and closing tags:
void XMLCALL
simpleXMLparser::textHndl(void *parserData, const char *textStr, const int length)
{
  simpleXMLparser* PData = (simpleXMLparser*)parserData;               // PData: pointer to 'this'
  
  PData->currentTextStr=textStr;
  PData->currentTextStr=PData->currentTextStr.substr(0,length);        
}



//----------------------
///////////////////////////////////////////////////////
//
// ChercheLaTag() -- better to rewrite this!!
//
///////////////////////////////////////////////////////
void simpleXMLparser::CherchezLaTag()
{
  XML_ParserReset(XMLp, NULL);

  XML_SetElementHandler(XMLp, simpleXMLparser::openingHndl, simpleXMLparser::closingHndl);
  XML_SetCharacterDataHandler(XMLp, simpleXMLparser::textHndl);
  XML_SetUserData(XMLp, this);   

  FILE * inputFile;
  if  ( ! (inputFile = fopen(xmlFileName,"r")) )
    {
      std::cout << "Error: File \"" << xmlFileName << "\" does not exist\n" ;
      exit(-1);
    }
	

  searchTag=tXMLpathList.front();     // tag for search
  currentTagNumber=0;                

  Depth = 0;
  prevDepth = 0;
  searchDepth=  tXMLpathList.size(); 

  currentAttr.clear();
  returnAttr.clear();

  ifStop=false;
 
  for (;;) {
    int done;
    int len;
  
    len = fread(Buff, 1, BUFFSIZE, inputFile);
    if (ferror(inputFile)) {
      std::cout << "\nError: read error\n\n";
      exit(-1);
    }
    done = feof(inputFile);

    if (XML_Parse(XMLp, Buff, len, done) == XML_STATUS_ERROR) {
      std::cout << "\nError: parse error at line " 
		<< XML_GetCurrentLineNumber(XMLp) << ": "
		<< XML_ErrorString(XML_GetErrorCode(XMLp)) << "\n\n";
      exit(-1);
    }

    if (done)
      break;
  }

  fclose(inputFile);
  tXMLpathList.clear();
}

const char* simpleXMLparser::formString(int iNumber)
{
  std::string retStr;
  std::ostringstream retStream;
  
  retStream<<iNumber;
  retStr=retStream.str();
  return retStr.c_str();
}
