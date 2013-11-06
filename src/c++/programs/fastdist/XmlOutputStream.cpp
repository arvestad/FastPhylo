#include "XmlOutputStream.hpp"
#include <cstdio>
#include <string>
#include <math.h>


using namespace std;

XmlOutputStream::XmlOutputStream(char * filename = 0 ) : DataOutputStream(filename) 
{
  fprintf(fp,"<?xml version=\"1.0\"?>\n<root>\n <runs>\n");
};

XmlOutputStream::~XmlOutputStream() 
{
  fprintf(fp," </runs>\n</root>\n");
};

void
XmlOutputStream::print( StrDblMatrix & dm ) 
{
	//mehmood changes here... email: malagori@kth.se
	//setXmlFlag( true);
	printPHYLIPfast(dm, fp, true );
}
//mehmood changes here... email: malagori@kth.se
void
XmlOutputStream::printRow( StrFloRow & dm, string name, int row, bool mem_eff_flag)
{

	const size_t numNodes = dm.getColumns();

	char defstr[11];// = "   .      ";
	defstr[0]=' ';
	defstr[3] = '.';
	defstr[10] = 0;
  	//the names PENDING NAME LENGTH

  	int entriesPerRow = numNodes;

	fprintf(fp,"    <row>\n" );
	if (mem_eff_flag == true){
		row=0;
	}
	for ( size_t j = row ; j < entriesPerRow ; j++ ){

		float f = dm.getDistance(j);

	   if ( ! isfinite(f) ){
			USER_WARNING("warning float not finite (use fix factor) " << f );
		   fprintf(fp,"     <entry>-1</entry>\n" );
	      continue;
	   }
	   //warning: this isn't enough to get the correct rounding but it is close
		f += 0.0000005;
		defstr[1]=' ';
		int intpart = (int) f;
		if ( intpart > 99 ){
			if ( f-intpart*1.0 < 0.000001 ){
				fprintf(fp,"     <entry>%d</entry>\n", intpart );
				continue;
			}
		   fprintf(fp,"     <entry>%f</entry>\n", f );
	      continue;
	      }

	      //      printf("F:%10.6f\n",f);
	      float decimalpart = f-1.0*intpart;
	      //warning: this isn't enough to get the correct rounding but it is close
	      //decimalpart += 0.0000005;
	      //write intpart
	      if ( intpart == 0 ) {
	      	defstr[2] = '0';
	      } else {
				defstr[2] = ONEDIGIT[intpart];
	        	intpart = intpart /10;
	        	if ( intpart != 0 ) {
	         	defstr[1] = ONEDIGIT[intpart];
	        	}
	      }

	      //write 6 decimals part
	      int deci = 4;
	      while ( deci <= 9 ){
	        	decimalpart = decimalpart*100.0;
	        	int index = (int) decimalpart;
	        	decimalpart = decimalpart-index;
	        	defstr[deci++] = TENDIGIT[index];
	        	defstr[deci++] = ONEDIGIT[index];
	      }
	      //      cfp << defstr << endl;

			// skip leading spaces
			int i = 0;
			while ( defstr[i] == ' ' ) {
		  		i++;
			}
			//	xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) &defstr[i] );
			fprintf(fp,"     <entry>%s</entry>\n", &defstr[i] );
    }

    fprintf(fp,"    </row>\n");
}

void
XmlOutputStream::printHeader( size_t numNodes ){}

void
XmlOutputStream::printStartRun( std::vector<string> & names, std::string & runId, Extrainfos &extrainfos ) 
{
  fprintf(fp,"  <run id=\"%s\" dim=\"%i\">\n   <identities>\n", runId.c_str(), names.size() );
      for(size_t namei=0 ; namei<names.size() ; namei++ )
	{
	  // This only works if we constrain the input by a schema to not have "<", "&" and such.
          // Otherwise we need to use xmlEncodeSpecialChars(xmlDocPtr doc,  const xmlChar * input)
	  //	  xmlChar * str = xmlEncodeSpecialChars( 0 ,( const xmlChar * ) names[namei].c_str()  );
	  fprintf(fp,"    <identity name=\"%s\"", names[namei].c_str() );

	  if ( extrainfos.size() > namei  && extrainfos[namei].size() > 0 ) {
	    fprintf(fp,">\n     %s\n    </identity>\n", (char *) extrainfos[namei].c_str() );
          } 
          else {
	    fprintf(fp,"/>\n");
         } 

		  // fprintf(fp,"     <entry>%s</entry>\n", ( char * ) str );
		  //  xmlFree(str);
	}
  fprintf(fp,"   </identities>\n   <dms>\n");
  //fprintf(fp,"   <dm>\n"); //mehmood changes here... email: malagori@kth.se
}


void
XmlOutputStream::printEndRun() 
{
//	 fprintf(fp,"   </dm>\n"); //mehmood changes here... email: malagori@kth.se
  fprintf(fp,"   </dms>\n  </run>\n");
}

void XmlOutputStream::printBootstrapSpliter(size_t numNodes){
  fprintf(fp,"   </dm>\n");
  fprintf(fp,"   <dm>\n");
}

