/*
 * PhylipDmOutputStream.cpp
 *
 *  Created on: Dec 1, 2011
 *      Author: Mehmood Alam Khan
 *      Email: malagori@kth.se
 */
#include "PhylipDmOutputStream.hpp"
//#include <cstdio>
//#include <string>
#include <math.h>

using namespace std;
void
PhylipDmOutputStream::print( StrDblMatrix & dm )
{
	//setXmlFlag( false);
	printPHYLIPfast(dm, fp, false );
}

void
PhylipDmOutputStream::printRow( StrFloRow & dm, string name, int row, bool mem_eff_flag) {

   const size_t numNodes = dm.getColumns();

   char defstr[11];// = "   .      ";
   defstr[0]=' ';
   defstr[3] = '.';
   defstr[10] = 0;
   //the names PENDING NAME LENGTH

   int entriesPerRow = numNodes;
   fprintf(fp,"%-10s", name.c_str());

	float f = 0.0;
	if (mem_eff_flag == false){
		for(size_t i = 0; i < row; ++i){
			fprintf(fp, "%10f",f);
		}
	}else{
		row=0;
	}

   for( size_t j = row ; j < entriesPerRow ; j++ ) {

		float f = dm.getDistance(j);

		if ( ! isfinite(f) ){
			USER_WARNING("warning float not finite (use fix factor) " << f );
			fprintf(fp,"        -1");
			continue;
		}
		//warning: this isn't enough to get the correct rounding but it is close
		f += 0.0000005;
		defstr[1]=' ';
		int intpart = (int) f;
		if ( intpart > 99 ){
			if ( f-intpart*1.0 < 0.000001 ){
				fprintf(fp,"%10d",intpart);
				continue;
			}
			fprintf(fp,"%10f",f);
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

		fwrite(defstr,sizeof(char),10,fp);

   }

	fprintf(fp,"\n");
}

void
PhylipDmOutputStream::printHeader( size_t numNodes ) {
  fprintf(fp,"%5d\n",numNodes);
}

void PhylipDmOutputStream::printBootstrapSpliter(size_t numNodes){
  //fprintf(fp,"%5d\n",numNodes);
}

