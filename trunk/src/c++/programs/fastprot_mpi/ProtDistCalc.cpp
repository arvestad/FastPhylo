#include "ProtDistCalc.hpp"
#include <cctype> // for toupper()
#include <cmath>  // for log() and pow()
#include <string>
#include "MaximumLikelihood.hpp"
#include "ProtSeqUtils.hpp"
#include "Matrix.hpp"
//#include "omp.h"

  /*
   * Method to call to calculate distances  
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, prot_sequence_translation_model t_model){
    switch (t_model.model){
      case id: 
        calculate_id_dists(sv, dm); 
        break;
      case jc: 
        calculate_jc_dists(sv, dm); 
        break;
      case jck: 
        calculate_kimura_dists(sv, dm); 
        break;
      case jcss: 
        calculate_stormsonnhammer_dists(sv, dm); 
        break;
      case wag:
      case day:
      case arve:
      case jtt:
      case mvr:
        if (t_model.ml)
          calculate_ml_dists(sv, dm, t_model.model); 
        else 
          calculate_ed_dists(sv, dm, t_model); 
        break;
    }
  }

  /*
   * Calculates the expected distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param tm Translation model that contains parameters for the calculations
   */
  void calculate_ed_dists(const SeqVec &sv, StrDblMatrix &dm, 
      const prot_sequence_translation_model &tm){

    Matrix Q(get_model_matrix(tm.model)); 
    DblVec eq = get_model_vec(tm.model);
    initialize_ed(tm.tp, tm.step_size, Q, eq); 

//#pragma omp parallel for default(none) shared(dm, sv)   
//#pragma omp schedule(runtime)
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = count_replacements(sv[i], sv[j]);
        double distance = calculate_ed(N);
        dm.setDistance(i, j, distance); 
      }
    }
  }

  /*
   * Calculates the expected distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param tm Translation model that contains parameters for the calculations
   */
  void calculate_ed_dists_mpi(const SeqVec &sv, std::vector<double> &dv, 
      const prot_sequence_translation_model &tm, int rank, int size){

    Matrix Q(get_model_matrix(tm.model)); 
    DblVec eq = get_model_vec(tm.model);
    initialize_ed(tm.tp, tm.step_size, Q, eq); 
    
    ////////////////////////////////////////////////////////////////////
    // Section below is used to find out which distances to calculate
    // Or more specific, what the starting index of 'i' and 'j' should be
    /////////////////////////////////////////////////////////////////// 
    int start_i = 0;
    int start_j = 1;
    int nr_seqs = sv.size();
    int seq_length = sv[0].seq.size();
    // Total number of distances to be calculated by all processes
    unsigned long total_nr_dists = (nr_seqs*nr_seqs-nr_seqs)/2;
    
    // If the number of distances isn't evenly dividible by the
    // the number of processes, some processes should calculate
    // one more distance than the others.
    // The processes to do this is the x first processes, where
    // x is the remainder of total_nr_dists / size
    unsigned long local_nr_dists = total_nr_dists/size;
    if (rank < total_nr_dists % size)
      local_nr_dists++;
    

    // Start with finding out at what point in the distance vector
    // this process should begin calculating, if we think of the
    // distance vector as a vector where ALL distances to be calculated
    // is laid out in a long row
    
    // For example, if I am process 0 of a total of 2 processes
    // and the number of sequences are 4 => total_nr_dists = 6,
    // Then I calculate from distance 0 (to distance 2, which is
    // the first row in the distance matrix). Process 1 on the 
    // other hand calculates from distance 3 to distance 5,
    // which is row two and three.


    unsigned long my_start = rank*total_nr_dists/size;

    if (total_nr_dists % size != 0) {
      if (rank < total_nr_dists % size) {
        // This process is one of the processes calculating one extra
        // element, thus the offset is (total_nr_dists/size+1). +1 for
        // the extra element..
        my_start = (total_nr_dists/size+1)*rank;
      } else {
        // Note that if process i is the last one to calculate one
        // extra distance, then process i+k, k >= 1, needs to take  
        // that into account when calculating where to start. 

        // (total_nr_dists & size) below is the number of processes
        // that are calculating one extra element
        my_start = (total_nr_dists % size) * (total_nr_dists/size+1);

        // Now we need the offset from the last process which did calculate
        // one extra element. that is (rank - (total_nr_dists % size))
        
        // For example, if I am process 2 of a total of 3 processes, and
        // the number of distances to be calculated is 10, then process
        // 0 calculates 4 distances, process 1 calculates 3 and I calculate
        // 3 distances. 
        // my_start above would be (10 % 3) * (10/3+1) = 1 * 4 = 4
        // the line below would be  (2 - (10 % 3)) *(10/3) = (2-1)*3 = 3
        // and my_start would be 4 + 3 = 7 where 4 = distances calculated
        // by process 0, 3 = distances calculated by process 1
        my_start += (rank - (total_nr_dists % size))*(total_nr_dists/size);
      }
    }
    
    // Now we know which distance we should start to calculate
    // But we also need to know the indexes in the sequence vector
    // for that distance. 
    // we want to emulate this loop:
    // for i = 0 to i<nr_sequences
    //    for j = i+1 to j<nr_sequences
    // which creates all pairs of sequences that we are interested in

    // We start at distance 0, and then step through the distance matrix
    // until we reach the number (my_start) we should start at
    for (unsigned long i=0; i<my_start; i++) {
      start_j++;
      if (start_j == nr_seqs) {
        start_i++;
        start_j = start_i+1;
      }
    }

    ////////////////////////////////////////////////////////////////////
    // Calculate actual distances below
    
    // limit is used to calculate the right amount of distances, for
    // every distance calculated it is increased, until we break the
    // outer (and inner) loop when local_nr_dists is reached
    unsigned long limit = 0;
    
    for (int i=start_i; limit<local_nr_dists; i++) {
      for (int j=start_j; j<nr_seqs && limit < local_nr_dists; j++) {
        Matrix N = count_replacements(sv[i], sv[j]);
        dv[limit] = calculate_ed(N);
        limit++;
      }
      // we want j to start at i+1 next round, thus start_j = i+2 
      // (since i++ isn't done yet)
      start_j = i+2;
    }
  }

  /*
   * Calculates the maximum likelihood distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param mt Specifies which model to use for the calculations
   */
  void calculate_ml_dists(const SeqVec &sv, StrDblMatrix &dm, model_type mt){
    Matrix Q = get_model_matrix(mt); 

    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = count_replacements(sv[i], sv[j]);
        double distance = likelihood_calc(N, Q);
        dm.setDistance(i, j, distance);
      }
    }
  }

  /*
   * Computes the identity based distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_id_dists(const SeqVec &sv, StrDblMatrix &dm){
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double distance = 100 * count_id_dist(sv[i], sv[j]);
        dm.setDistance(i, j, distance); 
      }
    }
  }

  /*
   * Calculates the Jukes-Cantor corrected identity distance
   * d = -(1900/20)*log(1-20/19)*(1-count_id_dist(s1, s2))
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_jc_dists(const SeqVec &sv, StrDblMatrix &dm){
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double diff = 1 - count_id_dist(sv[i], sv[j]);
        double distance = -(1900/20.0)*log(1-(20.0/19) * diff);
        dm.setDistance(i, j, distance); 
      }
    }

  }
  /*
   * Calculates the Kimura corrected distance between two sequences
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_kimura_dists(const SeqVec &sv, StrDblMatrix &dm){
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double diff = 1 - count_id_dist(sv[i], sv[j]);
        double adj_distance = diff + 0.2*diff*diff;
        if (adj_distance > 0.854)
          adj_distance = 0.854;
        double distance = -100*log(1-adj_distance);  
        dm.setDistance(i, j, distance); 
      }
    }
  }

  /*
   * Calculates Storm-Sonnhammer distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
void calculate_stormsonnhammer_dists(const SeqVec &sv, StrDblMatrix &dm){
  for (int i=0; i<sv.size(); i++) {
    for (int j=i+1; j<sv.size(); j++){
      double diff = 1 - count_id_dist(sv[i], sv[j]);
      double adj_distance;
      if (diff > 0.916)
        adj_distance = 1000.0;
      else {
        adj_distance = -100 * log(1 - 0.95844*diff
            - 0.69957 * pow(diff, 2)
            + 2.4955 * pow(diff, 3)
            - 4.6353 * pow(diff, 4)
            + 2.8076 * pow(diff, 5));
      }
      dm.setDistance(i, j, adj_distance); 
    }
  }
}
