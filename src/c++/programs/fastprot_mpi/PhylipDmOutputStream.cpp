#include "PhylipDmOutputStream.hpp"
#include <cstdio>
#include <libxml/xmlreader.h>

using namespace std;

void PhylipDmOutputStream::print( StrDblMatrix & dm ) {
	printPHYLIPfastSD(dm,fp,false,false);
}

void PhylipDmOutputStream::printSD( StrDblMatrix & dm ) {
	printPHYLIPfastSD(dm,fp,true,false);
}

void
printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out, bool writeXml, bool writeXmlSD ){

  const size_t numNodes = dm.getSize();

  //  xmlNodePtr dmNode;
  if ( writeXml ) {
    //    dmNode = xmlNewNode(0, ( const xmlChar * ) "dm");
    fprintf(out,"   <dm>\n");
  } else if (writeXmlSD) {
    fprintf(out,"   <sdm>\n");
  }
  else {
    fprintf(out,"%5lu\n",numNodes);
  }


  char defstr[12];// = "   .      ";
  defstr[0]=' ';
  defstr[4] = '.';
  defstr[11] = 0;
  //the names PENDING NAME LENGTH

  int entriesPerRow;


  if ( writeXml || writeXmlSD ) {
    entriesPerRow = 0;
  } else
  {
    entriesPerRow = numNodes;
  }
  for ( size_t i = 0 ; i < numNodes ; i++ ){

    //    xmlNodePtr rowNode;
    if (  writeXml || writeXmlSD )   {
      //      rowNode = xmlNewChild(dmNode,0, ( const xmlChar * ) "row",0);
      // xmlSetProp(rowNode, ( const xmlChar * ) "species",( const xmlChar * ) dm.getIdentifier(i).c_str() );

      //      fprintf(out,"    <row species=\"%s\">\n", dm.getIdentifier(i).c_str() );
      fprintf(out,"    <row>\n" );
    }
    else {
      fprintf(out,"%-11s", dm.getIdentifier(i).c_str());
    }

    if ( writeXml || writeXmlSD )  ( entriesPerRow++ );

    for ( size_t j = 0 ; j < entriesPerRow ; j++ ){
      float f = dm.getDistance(i,j);
      if ( ! isfinite(f) ){
        USER_WARNING("warning float not finite (use fix factor) " << f );

        if (  writeXml || writeXmlSD ) {
          //	  xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) "-1" );
          fprintf(out,"     <entry>-1</entry>\n" );
        }
        else {
          fprintf(out,"        -1");
        }
        continue;
      }
      //warning: this isn't enough to get the correct rounding but it is close
      f += 0.0000005;
      defstr[1]=' ';
      defstr[2]=' ';
      int intpart = (int) f;
      if ( intpart > 99 ){
        if ( f-intpart*1.0 <0.000001 ){
          if (  writeXml || writeXmlSD ) {
            // I guess 20 should be more than enough. Please lower this number if you know how it all works. /Erik Sjolund
            // char str[20];
            //	    snprintf(str,20,"%10d",intpart);
            //	    xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) str );
            fprintf(out,"     <entry>%11d</entry>\n", intpart );
          }
          else {
            fprintf(out,"%11d",intpart);

          }
          continue;
        }
        if (  writeXml || writeXmlSD ) {
          // I guess 20 should be more than enough. Please lower this number if you know how it all works. /Erik Sjolund
          //	  char str[20];
          //  snprintf(str,20,"%10f",f);
          // xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) str );
          fprintf(out,"     <entry>%11f</entry>\n", f );
        }
        else {
          fprintf(out,"%11f",f);
        }
        continue;
      }
      //      printf("F:%10.6f\n",f);
      float decimalpart = f-1.0*intpart;
      //warning: this isn't enough to get the correct rounding but it is close
      //decimalpart += 0.0000005;
      //write intpart
      if ( intpart == 0 )
        defstr[3] = '0';
      else {
        defstr[3] = DataOutputStream::ONEDIGIT[intpart];
	//std::cerr<<"inpart"<<intpart<<endl;
        intpart = intpart /10;
        if ( intpart != 0 )
          defstr[2] = DataOutputStream::ONEDIGIT[intpart];
      }

      //write 6 decimals part
      int deci = 5;
      while ( deci <= 10 ){
        decimalpart = decimalpart*100.0;
        int index = (int) decimalpart;
        decimalpart = decimalpart-index;
	//std::cerr<<"decimal part"<<decimalpart<<endl;
        defstr[deci++] = DataOutputStream::TENDIGIT[index];
        defstr[deci++] = DataOutputStream::ONEDIGIT[index];
      }
      //      cout << defstr << endl;

      if (  writeXml || writeXmlSD ) {
        // skip leading spaces
        int i = 0;
        while ( defstr[i] == ' ' ) {
          i++;
        }
        //	xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) &defstr[i] );
        fprintf(out,"     <entry>%s</entry>\n", &defstr[i] );


      }
      else {
	char a[2];
	//a[0]=' ';
        fwrite(defstr,sizeof(char),11,out);
	//std::cerr<<"defstr"<<defstr<<endl;
	//fwrite(a,sizeof(char),1,out);
      }
    }

    if (  writeXml || writeXmlSD )    {
      fprintf(out,"    </row>\n");
    } else {
      fprintf(out,"\n");
    }


  }
  if (  writeXml ) {
    //    xmlElemDump(out, 0, dmNode);
    fprintf(out,"   </dm>\n");
    //  xmlFreeNode(dmNode);
  } else if ( writeXmlSD) {
    fprintf(out,"   </sdm>\n");
  }
}

