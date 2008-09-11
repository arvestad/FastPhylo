#include <float.h>
#include <math.h>
#include <iostream>
#include "log_utils.hpp"
#include "DistanceMatrix.hpp"
#include "config.h"


bool
applyFixFactor(StrDblMatrix &dm, double fixFactor){
  double biggest = 0;
  int size =dm.getSize();
  
  for ( int i = 0 ; i < size ; i++ ){
    for ( int j = 0 ; j < size ; j++ ){
      double d = dm.getDistance(i,j);
      if ( isfinite(d) && d>biggest ){
        biggest = d;
      }
    }
  }

  bool changed = false;
  biggest = biggest*fixFactor;
  for ( int i = 0 ; i < size ; i++ ){
    for ( int j = 0 ; j < size ; j++ ){
      double d = dm.getDistance(i,j);
      if ( !isfinite(d) || d<0 ){
        dm.setDistance(i,j,biggest);
	changed  = true;
      }
    }
  }

  return changed;
}



//-------------------------------------------------------------
// FAST PRINTING OF FLOATS
static const char ONEDIGIT[128]={
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' ,
  '0' , '1' , '2' , '3' , '4' , '5' , '6' , '7'};
static const char  TENDIGIT[128] ={
  '0' , '0' , '0' , '0' , '0' , '0' , '0' , '0' , '0' , '0' ,
  '1' , '1' , '1' , '1' , '1' , '1' , '1' , '1' , '1' , '1' ,
  '2' , '2' , '2' , '2' , '2' , '2' , '2' , '2' , '2' , '2' ,
  '3' , '3' , '3' , '3' , '3' , '3' , '3' , '3' , '3' , '3' ,
  '4' , '4' , '4' , '4' , '4' , '4' , '4' , '4' , '4' , '4' ,
  '5' , '5' , '5' , '5' , '5' , '5' , '5' , '5' , '5' , '5' ,
  '6' , '6' , '6' , '6' , '6' , '6' , '6' , '6' , '6' , '6' ,
  '7' , '7' , '7' , '7' , '7' , '7' , '7' , '7' , '7' , '7' ,
  '8' , '8' , '8' , '8' , '8' , '8' , '8' , '8' , '8' , '8' ,
  '9' , '9' , '9' , '9' , '9' , '9' , '9' , '9' , '9' , '9' ,
  '0' , '0' , '0' , '0' , '0' , '0' , '0' , '0' , '0' , '0' ,
  '1' , '1' , '1' , '1' , '1' , '1' , '1' , '1' , '1' , '1' ,
  '2' , '2' , '2' , '2' , '2' , '2' , '2' , '2' };


void
printPHYLIPfast(const StrDblMatrix &dm, FILE *out, bool writeXml ){

  const size_t numNodes = dm.getSize();

  //  xmlNodePtr dmNode;
  if ( writeXml ) {
    //    dmNode = xmlNewNode(0, ( const xmlChar * ) "dm"); 
    fprintf(out,"   <dm>\n"); 
  } 
  else {
    fprintf(out,"%5lu\n",numNodes); 
  }


  char defstr[11];// = "   .      ";
  defstr[0]=' ';
  defstr[3] = '.';
  defstr[10] = 0;
  //the names PENDING NAME LENGTH

  int entriesPerRow;


  if ( writeXml ) {
    entriesPerRow = 0;
  } else
    {
  entriesPerRow = numNodes;
    }
  for ( size_t i = 0 ; i < numNodes ; i++ ){

    //    xmlNodePtr rowNode;
    if (  writeXml )   { 
      //      rowNode = xmlNewChild(dmNode,0, ( const xmlChar * ) "row",0);
      // xmlSetProp(rowNode, ( const xmlChar * ) "species",( const xmlChar * ) dm.getIdentifier(i).c_str() );

      //      fprintf(out,"    <row species=\"%s\">\n", dm.getIdentifier(i).c_str() ); 
      fprintf(out,"    <row>\n" ); 
    }
    else {
      fprintf(out,"%-10s", dm.getIdentifier(i).c_str());
    }
    
    if ( writeXml )  ( entriesPerRow++ );

    for ( size_t j = 0 ; j < entriesPerRow ; j++ ){
      float f = dm.getDistance(i,j);
      if ( ! isfinite(f) ){
	USER_WARNING("warning float not finite (use fix factor) " << f );

	if (  writeXml ) { 
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
      int intpart = (int) f;
      if ( intpart > 99 ){
        if ( f-intpart*1.0 <0.000001 ){

	  if (  writeXml ) { 
            // I guess 20 should be more than enough. Please lower this number if you know how it all works. /Erik Sjolund
	    // char str[20];   
	    //	    snprintf(str,20,"%10d",intpart);  
	    //	    xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) str ); 
                 fprintf(out,"     <entry>%10d</entry>\n", intpart ); 

	  }
	  else {
	    fprintf(out,"%10d",intpart);
	  }

          continue;
        } 

	if (  writeXml ) { 
          // I guess 20 should be more than enough. Please lower this number if you know how it all works. /Erik Sjolund
	  //	  char str[20];  
	  //  snprintf(str,20,"%10f",f);  
	  // xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) str ); 
                 fprintf(out,"     <entry>%10f</entry>\n", f ); 

	} 
	else {
	  fprintf(out,"%10f",f);
	}

        continue;
      } 
      //      printf("F:%10.6f\n",f);
      float decimalpart = f-1.0*intpart;
      //warning: this isn't enough to get the correct rounding but it is close
      //decimalpart += 0.0000005;
      //write intpart
      if ( intpart == 0 )
        defstr[2] = '0';
      else {
	defstr[2] = ONEDIGIT[intpart];
        intpart = intpart /10;
        if ( intpart != 0 )
          defstr[1] = ONEDIGIT[intpart];
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
      //      cout << defstr << endl;

      if (  writeXml ) { 
        // skip leading spaces 
	int i = 0; 
	while ( defstr[i] == ' ' ) { 
	  i++; 
        }
	//	xmlNodePtr entryNode = xmlNewChild(rowNode,0, ( const xmlChar * ) "entry",  ( const xmlChar * ) &defstr[i] ); 
                 fprintf(out,"     <entry>%s</entry>\n", &defstr[i] ); 


      } 
      else {
	fwrite(defstr,sizeof(char),10,out); 
      }
    }
    
  if (  writeXml )    {
      fprintf(out,"    </row>\n"); 
    } else {
      fprintf(out,"\n"); 
    }
    
   
  }
  if (  writeXml ) { 
    //    xmlElemDump(out, 0, dmNode);
    fprintf(out,"   </dm>\n");
    //  xmlFreeNode(dmNode);
  } 
}














