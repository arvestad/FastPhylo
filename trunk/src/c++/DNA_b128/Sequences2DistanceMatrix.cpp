//--------------------------------------------------
//                                        
// File: Sequences2DistanceMatrix.cpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Sequences2DistanceMatrix.cpp,v 1.51 2006/12/31 11:17:52 isaac Exp $                                 
//
//--------------------------------------------------

#include "Sequences2DistanceMatrix.hpp"
#include "DNA_b128_String.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "stl_utils.hpp"
#include "log_utils.hpp"
#include "file_utils.hpp"
#include <float.h>

using namespace std;


void
Sequences2DNA_b128(std::vector<Sequence> &seqs, std::vector<DNA_b128_String> &b128){

  b128.resize(seqs.size());
  for(size_t i=0;i<seqs.size();i++){
    size_t cap = seqs[i].seq.length()+1;
    b128[i].reInitiate(cap);
    b128[i].append(seqs[i].seq);
  }
}


void
DNA_b128_StringsFromPHYLIP(ifstream &fin, std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings){

  int numSequences;
  int seqlen;
  const int MAXLINE = 16384;
  char line[MAXLINE];
  do{//skip lines that does not contain two integers
    fin.getline(line,MAXLINE);
  } while( sscanf(line,"%d %d",&numSequences,&seqlen) != 2 );


  names.clear();
  names.reserve(numSequences);
  b128_strings.resize(numSequences);
  for ( int i = 0 ; i < numSequences ; i++ ){
    b128_strings[i].reInitiate(seqlen);
    names.push_back(string());
  }
  //phylip has name lenght 10.
  //read the names and map the sequences onto the tree
  char tmpName[11];
  for ( int i = 0 ; i < numSequences ; i++ ){
    DNA_b128_String &s = b128_strings[i];
  
    fin.getline(tmpName,11);//reads atmost 10 chars
    //skip lines without 10 chars per line
    if( !fin.fail() ){
      if( fin.eof() ) THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
      i--;
      continue;
    }
    appendUntil(names[i],tmpName,fin.gcount(), ' ');
    
    fin.clear();
    fin.getline(line,MAXLINE);

    while( fin.fail() && fin.gcount()==MAXLINE-1 ){//didn't read all the line
      s.append(line);
      fin.clear();
      fin.getline(line,MAXLINE);
    }
    if( !fin.fail()) {//we read it all including the newline char unless it ended with eof
      s.append(line);
    }
    else //fail
      THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof()); 
  }//end for loop

  //The sequences aren't neccesarily on one line but my be spread out interleaving
  //over several lines. Therefore we read until seqlen chars have been read.
  while ( b128_strings[0].getNumChars() < seqlen ){
    for ( int i = 0 ; i < numSequences ; i++ ){
      DNA_b128_String &s = b128_strings[i];
      
      fin.getline(tmpName,11);//reads atmost 10 chars
      //skip lines without 10 chars per line
      if( !fin.fail() ){
	if( fin.eof() ) THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
	i--;
	continue;
      }
      
      fin.clear();
      fin.getline(line,MAXLINE);

      while( fin.fail() && fin.gcount()==MAXLINE-1 ){//didn't read all the line
	s.append(line);
	fin.clear();
	fin.getline(line,MAXLINE);
      }
      if( !fin.fail()) {//we read it all including the newline char unless it ended with eof
	s.append(line);
      }
      else //fail
	THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof()); 
    }//end for loop    
  }

  // CHECK THAT ALL STRINGS HAVE THE SAME LENGTH
  for ( int i = 0 ; i < numSequences ; i++ ){
    if ( b128_strings[i].getNumChars() != seqlen ){
       THROW_EXCEPTION("Sequence not of correct length: " << names[i] 
		       << "    length is " << b128_strings[i].getNumChars()); 
    }
  }
  //-----------------------------

}
//---------------------------------------------------------
void 
bootstrapSequences(const std::vector<Sequence> &seqs, std::vector<DNA_b128_String> &b128_strings){


  //ensure capacity in bootsequences
  b128_strings.resize(seqs.size());

  const size_t seqlen = seqs[0].seq.length();
  for ( size_t i = 0 ; i < seqs.size() ; i++ )
    b128_strings[i].reInitiate(seqlen);
  
  // Do the bootstrapping
  size_t pos=0;
  size_t seq;  
  vector<int> samplePositions(seqlen);
  const size_t BUFFSIZE = (16383>seqlen ? seqlen : 16383); //2^14=16384
  char buff[BUFFSIZE+1];
  buff[BUFFSIZE]='\0';
  size_t const stride = 32;

  if( stride < seqlen){
    for( pos=0; pos<(seqlen-stride); pos+= stride)
      for( size_t i=0; i<stride; i++)
	samplePositions[pos+i] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));
    for (; pos<seqlen; pos++)
      samplePositions[pos] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));
    
    
    for( seq=0; seq<seqs.size(); seq++){
      const string & s = seqs[seq].seq;
      //-----------
      pos = 0;
      for( ; pos<seqlen-BUFFSIZE; pos+=BUFFSIZE){
	for( size_t i=0; i<BUFFSIZE; i++){
	  buff[i] = s[samplePositions[pos+i]]; 
	}
	b128_strings[seq].append(buff);
      }
      size_t i;
      for ( i=0; pos<seqlen; i++,pos++)
	buff[i] = s[samplePositions[pos]]; 
      buff[i] = '\0';
      //----------
      b128_strings[seq].append(buff);
    }
  }
  else{//seqlen<stride
    for (pos=0; pos<seqlen; pos++)
      samplePositions[pos] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));    
    
    for( seq=0; seq<seqs.size(); seq++){
      const string & s = seqs[seq].seq;
      for (pos=0; pos<seqlen; pos++)
	buff[pos] = s[samplePositions[pos]]; 
      buff[pos]='\0';
      b128_strings[seq].append(buff);
    } 
  }
}






//----------------------------------------------------------
        
void
fillMatrix(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
		sequence_translation_model trans_model){

  const size_t numSequences = seqs.size();
  dm.resize(numSequences);
  

  //base frequences
  DNA_b128_String::base_frequences freqs = seqs[0].getBaseFrequences();
  for ( size_t i = 1 ; i < numSequences ; i++ ){
    DNA_b128_String::base_frequences tmpfreqs = seqs[i].getBaseFrequences();
    freqs.num_As_+= tmpfreqs.num_As_;
    freqs.num_Cs_+= tmpfreqs.num_Cs_;
    freqs.num_Gs_+= tmpfreqs.num_Gs_;
    freqs.num_Ts_+= tmpfreqs.num_Ts_;
    freqs.num_unknowns_+= tmpfreqs.num_unknowns_;
    freqs.num_ambiguities_+= tmpfreqs.num_ambiguities_;
  }

  //compute ambiguities probabilities 
  if ( ! trans_model.no_ambiguities ){
    if ( trans_model.use_base_freqs ){
      for ( size_t i = 0 ; i < numSequences ; i++ ){
        seqs[i].calcAmbiguityProbabilities( freqs.num_As_,freqs.num_Cs_,freqs.num_Gs_,freqs.num_Ts_);
      }
    }
    else {
      for ( size_t i = 0 ; i < numSequences ; i++ ){
        seqs[i].calcAmbiguityProbabilitiesUNIFORM();
      }
    }
  }

   //CALL DISTANCE COMPUTATION
   switch( trans_model.model ){
   case HAMMING_DISTANCE: fillMatrix_Hamming(dm,seqs,trans_model);break;
   case JC:  fillMatrix_JC(dm,seqs,trans_model);break;
   case K2P:  fillMatrix_K2P(dm,seqs,trans_model);break;
   case TN93: fillMatrix_TN93(dm,seqs,freqs,trans_model);break;
   default:
     PROG_ERROR("Non handled model: " << trans_model.model);
   }


}




void 
fillMatrix_Hamming(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
		   sequence_translation_model trans_model){

  const size_t numSequences = seqs.size();
  const size_t strlen = seqs[0].getNumChars();

  dm.resize(numSequences);

  //extended information used to compute ambiguities
  static DistanceMatrix<int,simple_string_distance, 
    empty_Data_init<int>,empty_Data_printOn<int>, 
    empty_Data_init<simple_string_distance>, empty_Data_printOn<simple_string_distance> > 
    extendedDistanceInfo(numSequences);
  extendedDistanceInfo.resize(numSequences);

  //
  // The loop will resolve the ambiguities according to the translation model
  // 
  for ( size_t i = 0 ; i < numSequences; i++ ){
    dm.setDistance(i,i,0);
    size_t closestNeig = i;
    float closestDist =  FLT_MAX;

    DNA_b128_String &si = seqs[i];
    //only if it has ambiguities does it need to be resolved
    //start by checking the allready computed distances for si
    if ( si.hasAmbiguities() ){
        for ( size_t k = 0 ; k < i ; k++ ){
          float dist = dm.getDistance(k,i);
          if (  dist < closestDist && dist >= 0 ){
            closestDist = dist;
            closestNeig = k;
          }
        }
    }
    
    // compute the remaining distances for si
    for ( size_t j = i+1 ; j < numSequences ; j++ ){
      //compute distance without using the ambiguities
      simple_string_distance sd = DNA_b128_String::computeDistance(si,seqs[j]);
      float hamdist = compute_Hamming_distance(sd);

      dm.setDistance(i,j,hamdist);
      extendedDistanceInfo.setDistance(i,j,sd);

      if ( hamdist < closestDist && hamdist >= 0 ){//update the closest neighbor
        closestDist = hamdist;
        closestNeig = j;
      }
    }
    //all distances for si have been compted

    //Resolve the ambiguities according to the closest neighbor
    if ( si.hasAmbiguities() && ! trans_model.no_ambig_resolve ){
      si.resolveAmbiguities(seqs[closestNeig]);
    }
  }

  
  //Update the computed distances with the ambiguities
  if ( !trans_model.no_ambiguities ){
    for ( size_t i = 0 ; i < numSequences ; i++ ){
      DNA_b128_String &si = seqs[i];
      for ( size_t j = i+1 ; j < numSequences ; j++ ){
        if ( si.hasAmbiguities() || seqs[j].hasAmbiguities() ){
          simple_string_distance sd = extendedDistanceInfo.getDistance(i,j);
	  ML_string_distance ml_dist = compute_JC(strlen,sd);
          sd = DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd,ml_dist,si,seqs[j]);

          //hamming
          dm.setDistance(i,j,compute_Hamming_distance(sd));
           
        }
      }
    }
  }
   
}


void 
fillMatrix_JC(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
	      sequence_translation_model trans_model){

  const size_t numSequences = seqs.size();
  const size_t strlen = seqs[0].getNumChars();

  dm.resize(numSequences);
  
   //extended information used to compute ambiguities
  static DistanceMatrix<int,pair<simple_string_distance,ML_string_distance>, 
    empty_Data_init<int>,empty_Data_printOn<int>, 
    empty_Data_init<pair<simple_string_distance,ML_string_distance> >,
    empty_Data_printOn<pair<simple_string_distance,ML_string_distance> > > 
    extendedDistanceInfo(numSequences);
  
  extendedDistanceInfo.resize(numSequences);
 
  for ( size_t i = 0 ; i < numSequences ; i++ ){
    dm.setDistance(i,i,0);
    size_t closestNeig = i;
    float closestDist =  FLT_MAX;
    
    DNA_b128_String &si = seqs[i];

    //if has ambig find closest string and resolve.
    if ( si.hasAmbiguities() ){
        for ( size_t k = 0 ; k < i ; k++ ){
          float dist = dm.getDistance(k,i);
          if (  dist < closestDist && dist >= 0 ){
            closestDist = dist;
            closestNeig = k;
          }
        }
    }

    
    ML_string_distance  ml_dist;
    for ( size_t j = i+1 ; j < numSequences ; j++ ){
      simple_string_distance sd = DNA_b128_String::computeDistance(si,seqs[j]);
      //     cout << sd << endl;
      ml_dist = compute_JC(strlen,sd);

      dm.setDistance(i,j,ml_dist.distance);
      extendedDistanceInfo.setDistance(i,j,pair<simple_string_distance,ML_string_distance>(sd,ml_dist));
      
      if ( ml_dist.distance < closestDist && ml_dist.distance >= 0 ){
        closestDist = ml_dist.distance;
        closestNeig = j;
      }
    }

    if ( si.hasAmbiguities() && ! trans_model.no_ambig_resolve ){
      ml_dist = extendedDistanceInfo.getDistance(i,closestNeig).second;
      si.resolveAmbiguitiesUsingTransitionProbabilities(seqs[closestNeig],ml_dist);
    }
  }
  
  //UPDATE USING AMBIGUITIES
  if ( !trans_model.no_ambiguities ){
    for (size_t i = 0 ; i < numSequences ; i++ ){
      DNA_b128_String &si = seqs[i];
      for ( size_t j = i+1 ; j < numSequences ; j++ ){
        if ( si.hasAmbiguities() || seqs[j].hasAmbiguities() ){
          simple_string_distance sd = extendedDistanceInfo.getDistance(i,j).first;
          ML_string_distance ml_dist = extendedDistanceInfo.getDistance(i,j).second;
          sd = DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd,ml_dist,si,seqs[j]);

          ml_dist = compute_JC(strlen,sd);      
          dm.setDistance(i,j,ml_dist.distance);
        }
      }
    }
  }
}



void 
fillMatrix_K2P(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
	       sequence_translation_model trans_model){

  const size_t numSequences = seqs.size();
  const size_t strlen = seqs[0].getNumChars();

  dm.resize(numSequences);
  

  //extended information used to compute ambiguities
  static DistanceMatrix<int,pair<simple_string_distance,ML_string_distance>, 
    empty_Data_init<int>,empty_Data_printOn<int>, 
    empty_Data_init<pair<simple_string_distance,ML_string_distance> >,
    empty_Data_printOn<pair<simple_string_distance,ML_string_distance> > >
    extendedDistanceInfo(numSequences);
  extendedDistanceInfo.resize(numSequences);

  for ( size_t i = 0 ; i < numSequences ; i++ ){
    dm.setDistance(i,i,0);
    int closestNeig = -1;
    float closestDist =  FLT_MAX;
    DNA_b128_String &si = seqs[i];
    
    //if has ambig find closest string and resolve.
    if ( si.hasAmbiguities() ){
      for ( size_t k = 0 ; k < i ; k++ ){
	float dist = dm.getDistance(k,i);
	if (  dist < closestDist && dist >= 0 ){
	  closestDist = dist;
	  closestNeig = k;
	}
      }
    }
    
    
    for ( size_t j = i+1 ; j < numSequences ; j++ ){
      simple_string_distance sd = DNA_b128_String::computeDistance(si,seqs[j]);
      ML_string_distance  ml_dist;
      if ( trans_model.no_tstvratio )
	ml_dist = compute_K2P(strlen,sd);
      else
	ml_dist = compute_K2P_fixratio(strlen,sd,trans_model.tstvratio);
	
      dm.setDistance(i,j,ml_dist.distance);
      extendedDistanceInfo.setDistance(i,j,pair<simple_string_distance,ML_string_distance>(sd,ml_dist));
	
      if ( ml_dist.distance < closestDist && ml_dist.distance >= 0 ){
	closestDist = ml_dist.distance;
	closestNeig = j;
      }
    }
      
    if ( si.hasAmbiguities() && ! trans_model.no_ambig_resolve ){
      ML_string_distance  ml_dist = extendedDistanceInfo.getDistance(i,closestNeig).second;
      si.resolveAmbiguitiesUsingTransitionProbabilities(seqs[closestNeig],ml_dist);
    }
  }
    
  //UPDATE USING AMBIGUITIES
  if( ! trans_model.no_ambiguities ){
    for ( size_t i = 0 ; i < numSequences ; i++ ){
      DNA_b128_String &si = seqs[i];
      for ( size_t j = i+1 ; j < numSequences ; j++ ){
	if ( si.hasAmbiguities() || seqs[j].hasAmbiguities() ){
	  simple_string_distance sd = extendedDistanceInfo.getDistance(i,j).first;
	  ML_string_distance ml_dist = extendedDistanceInfo.getDistance(i,j).second;
	  sd = DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd,ml_dist,si,seqs[j]);
	  if ( trans_model.no_tstvratio )
	    ml_dist = compute_K2P(strlen,sd);
	  else
	    ml_dist = compute_K2P_fixratio(strlen,sd,trans_model.tstvratio);
	  
	  dm.setDistance(i,j,ml_dist.distance);           
	}
      }
    }
  }
  
}



void fillMatrix_TN93(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, 
		     DNA_b128_String::base_frequences freqs,
		     sequence_translation_model trans_model){

  const size_t numSequences = seqs.size();
  const size_t strlen = seqs[0].getNumChars();

  dm.resize(numSequences);
  
  
  //extended information used to compute ambiguities
  static DistanceMatrix<int,pair<TN_string_distance,ML_string_distance>, 
    empty_Data_init<int>,empty_Data_printOn<int>, 
    empty_Data_init<pair<TN_string_distance,ML_string_distance> >,
    empty_Data_printOn<pair<TN_string_distance,ML_string_distance> > >
    extendedDistanceInfo(numSequences);
  extendedDistanceInfo.resize(numSequences);

  for ( size_t i = 0 ; i < numSequences ; i++ ){
    dm.setDistance(i,i,0);
    size_t closestNeig = i;
    float closestDist =  FLT_MAX;
      
    DNA_b128_String &si = seqs[i];
      
    //if has ambig find closest string and resolve.
    if ( si.hasAmbiguities() ){
      for ( size_t k = 0 ; k < i ; k++ ){
	float dist = dm.getDistance(k,i);
	if (  dist < closestDist && dist >= 0 ){
	  closestDist = dist;
	  closestNeig = k;
	}
      }
    }
      
    
    ML_string_distance  ml_dist;
    for ( size_t j = i+1 ; j < numSequences ; j++ ){
      TN_string_distance tn = DNA_b128_String::computeTAMURANEIDistance(si,seqs[j]);
	
      if ( trans_model.no_tstvratio )
	ml_dist = compute_Tamura_Nei(strlen,tn,freqs.num_As_,freqs.num_Cs_,freqs.num_Gs_,freqs.num_Ts_);
      else
	ml_dist = compute_Tamura_Nei_fixratio(strlen,tn,
					      freqs.num_As_,freqs.num_Cs_,freqs.num_Gs_,freqs.num_Ts_,
					      trans_model.tstvratio, trans_model.pyrtvratio);
      dm.setDistance(i,j,ml_dist.distance);
      extendedDistanceInfo.setDistance(i,j,pair<TN_string_distance,ML_string_distance>(tn,ml_dist));
	
      if ( ml_dist.distance < closestDist && ml_dist.distance >= 0 ){
	closestDist = ml_dist.distance;
	closestNeig = j;
      }
    }
      
    if ( si.hasAmbiguities() && ! trans_model.no_ambig_resolve ){
      ml_dist = extendedDistanceInfo.getDistance(i,closestNeig).second;
      si.resolveAmbiguitiesUsingTransitionProbabilities(seqs[closestNeig],ml_dist);
    }
  }
    
  //UPDATE USING AMBIGUITIES
  if( !trans_model.no_ambiguities ){
    for ( size_t i = 0 ; i < numSequences ; i++ ){
      DNA_b128_String &si = seqs[i];
      for ( size_t j = i+1 ; j < numSequences ; j++ ){
        if ( si.hasAmbiguities() || seqs[j].hasAmbiguities() ){
          TN_string_distance tn = extendedDistanceInfo.getDistance(i,j).first;
          ML_string_distance ml_dist = extendedDistanceInfo.getDistance(i,j).second;
          if ( trans_model.no_transition_probs ){
            tn = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(tn,si,seqs[j]);
          }
          else {
            tn = DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(tn,ml_dist,si,seqs[j]);
          }
          if ( trans_model.no_tstvratio )
            ml_dist = compute_Tamura_Nei(strlen,tn,freqs.num_As_,freqs.num_Cs_,freqs.num_Gs_,freqs.num_Ts_);
          else
            ml_dist = compute_Tamura_Nei_fixratio(strlen,tn,
						  freqs.num_As_,freqs.num_Cs_,freqs.num_Gs_,freqs.num_Ts_,
						  trans_model.tstvratio, trans_model.pyrtvratio);
	  
          dm.setDistance(i,j,ml_dist.distance);
	  
        }
      }
    }
  }
}























