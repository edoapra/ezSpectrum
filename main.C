#include "genincludes.h"
#include "genutil.h"
#include "harmonic_pes_main.h"
#include "simple_xml_parser.h"


main(int argc, char* argv[])
{

  if ((argc-1)<1) {
    std::cerr << argv[0] << ": one and only one argument required, <Job.xml>.\n";
    exit(1);
  }
  std::string arg=argv[1];

  std::cout << "Job \"" << argv[0] << ' ' << argv[1] << "\" has been started: " << GetTime() << '\n';

  simpleXMLparser xmlF;
  xmlF.assignFile(arg.c_str()); 


  xmlF.reset().node("input");
  bool if_web_version=false;
  if ( xmlF.CheckSubNode("if_web_version") )
    {
      xmlF.node("if_web_version");
      if_web_version = xmlF.getBoolValue("flag"); 
    }


  if (not(if_web_version))
    {
      std::cout << "A copy of the \"" << argv[1] << "\" input:\n";
      std::cout << "------------------------------------------------------------------------------\n";
      xmlF.printInputFile();
      std::cout << "------------------------------------------------------------------------------\n \n";
    }

  std::string job;
  job=xmlF.reset().node("input").value("job");

  bool done = false;
  if (job == "harmonic_pes")
    done=harmonic_pes_main(arg.c_str());

  if ( !done )
    std::cout << "Method \"" << job <<"\" is unknown, or it has been failed. \n";
  
  std::cout << '\n' << "Job \"" << argv[0] << ' ' << argv[1] << "\" has been finished: " <<  GetTime() << '\n';
  
}


