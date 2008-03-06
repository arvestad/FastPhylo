<?xml version="1.0"?>
<article  
 id="manual">
  <articleinfo>
    <title>fastphylo</title>
  </articleinfo>
  <sect1 id="introduction">
    <title>Introduction</title>
    <para>
fastphylo is software project containing the implementations of the algorithms "Fast Computation of Distance Estimators" and "Fast Neighbor Joining".
 The software is licensed under the MIT license.  </para>
    <para>
        The primary URL for this document is
         <ulink url="http://fastphylo.sourceforge.net">http://fastphylo.sourceforge.net</ulink>.
    </para>
  </sect1>




  <sect1 id="Fast_Computation_of_Distance Estimators">
    <title>Fast Computation of Distance Estimators</title>


  <sect2 id="about_the_article_Fast_Computation_of_Distance Estimators">
    <title>About the published article</title>

<para>
Isaac Elias and Jens Lagergren published the algorithm in the journal <ulink url="http://www.biomedcentral.com/bmcbioinformatics/">BMC Bioinformatics</ulink> in 2007.
      </para>
      <para>

 BibTex
<programlisting>
@Article{EliasLagergren_fastdist
  author =      {Isaac Elias and Jens Lagergren},
  title =	{Fast Computation of Distance Estimators},
  journal =	{BMC Bioinformatics},
  year =        {2007},
  pages =       {89},
  volume =      {8}
}
</programlisting>

      </para>



  </sect2>


  <sect2 id="abstract_of_the_article_Fast_Computation_of_Distance Estimators">
    <title>Abstract of the published article</title>

<para>

<emphasis>Background:</emphasis> Some distance methods are among the most commonly used methods for reconstructing phylogenetic trees from sequence data. The input to a distance method is a distance matrix, containing estimated pairwise distances between all pairs of taxa. Distance methods themselves are often fast, e.g., the famous and popular Neighbor Joining (NJ) algorithm reconstructs a phylogeny of <emphasis>n</emphasis> taxa in time <emphasis>O(n<superscript>3</superscript>)</emphasis>. Unfortunately, the fastest practical algorithms known for computing the distance matrix, from <emphasis>n</emphasis>  sequences of length <emphasis>l</emphasis>, takes time proportional to <emphasis>l·n<superscript>2</superscript></emphasis>. Since the sequence length typically is much larger than the number of taxa, the distance estimation is the bottleneck in phylogeny reconstruction. This bottleneck is especially apparent in reconstruction of large phylogenies or in applications where many trees have to be reconstructed, e.g., bootstrapping and genome wide applications.

      </para>

      <para>
<emphasis>Results:</emphasis> We give an advanced algorithm for computing the number of mutational events between DNA sequences which is significantly faster than both Phylip and Paup. Moreover, we give a new method for estimating pairwise distances between sequences which contain ambiguity symbols. This new method is shown to be more accurate as well as faster than earlier methods.
      </para>

      <para>
<emphasis>Conclusions: </emphasis> Our novel algorithm for computing distance estimators provides a valuable tool in phylogeny reconstruction. Since the running time of our distance estimation algorithm is comparable to that of most distance methods, the previous bottleneck is removed. All distance methods, such as NJ, require a distance matrix as input and, hence, our novel algorithm significantly improves the overall running time of all distance methods. In particular, we show for real world biological applications how the running time of phylogeny reconstruction using NJ is improved from a matter of hours to a matter of seconds.

      </para>



  </sect2>



  <sect2 id="Supplementary_Material">
    <title>Supplementary Material</title>

<para>
Supplementary Material - Fast Computation of Distance Estimators.

Contains additional figures for the tests run on the ambiguity approaches.
(<ulink url="http://www.nada.kth.se/~isaac/publications/fastdist/fastdist_supplementary_material.pdf" >PDF</ulink>)


      </para>


<para>
Simulated Test Data for Ambiguities 
(<ulink url="http://www.nada.kth.se/~isaac/publications/fastdist/fastdist_simulated_data.tar.gz" >Tar archive</ulink>)
      </para>

<para>
Biological Test Data 
(<ulink url="http://www.nada.kth.se/~isaac/publications/fastdist/fastdist_biological_data.tar.gz" >Tar archive</ulink>)
      </para>

<para>
Command file used for running Paup 
(<ulink url="http://www.nada.kth.se/~isaac/publications/fastdist/paup_command.nex" >Nexus file</ulink>)
      </para>



  </sect2>




  </sect1>


  <sect1 id="Fast_Neighbor_Joining">
    <title>Fast Neighbor Joining</title>


  <sect2 id="about_the_article_Fast_Neighbor_Joining">
    <title>About the published article</title>

<para>
Isaac Elias and Jens Lagergren published the algorithm in the book "Proc. of the 32nd International Colloquium on Automata, 
                Languages and Programming ({ICALP}'05)" in 2005.
      </para>
      <para>
 BibTex
<programlisting>@InProceedings{ICALP05:EliasLagergren_FNJ,
  author =      {Isaac Elias and Jens Lagergren},
  title =	{Fast Neighbor Joining},
  booktitle =	{Proc. of the 32nd International Colloquium on Automata, 
                Languages and Programming ({ICALP}'05)},
  pages =	{1263--1274},
  year =	{2005},
  volume =	{3580},
  series =	{Lecture Notes in Computer Science},
  month =	{July},
  publisher =	{Springer-Verlag},
  ISBN =	{3-540-27580-0},
}
</programlisting>

      </para>



  </sect2>


  <sect2 id="abstract_of_the_article_Fast_Neighbor_Joining">
    <title>Abstract of the published article</title>

<para>
Reconstructing the evolutionary history of a set of species is a fundamental problem in biology and methods for solving this problem are gaged based on two characteristics: accuracy and efficiency. Neighbor Joining (NJ) is a so-called distance-based method that, thanks to its good accuracy and speed, has been embraced by the phylogeny community. It takes the distances between <emphasis>n</emphasis> taxa and produces in <emphasis>&#920;(n<superscript>3</superscript>)</emphasis> time a phylogenetic tree, i.e., a tree which aims to describe the evolutionary history of the taxa. In addition to performing well in practice, the NJ algorithm has optimal reconstruction radius.
</para><para>
The contribution of this paper is twofold: (1) we present an algorithm called Fast Neighbor Joining (FNJ) with optimal reconstruction radius and optimal run time complexity <emphasis>O(n<superscript>2</superscript>)</emphasis> and (2) we present a greatly simplified proof for the correctness of NJ. Initial experiments show that FNJ in practice has almost the same accuracy as NJ, indicating that the property of optimal reconstruction radius has great importance to their good performance. Moreover, we show how improved running time can be achieved for computing the so-called correction formulas.
      </para>



  </sect2>



  <sect2 id="Supplementary_Material_Fast_Neighbor_Joining">
    <title>Supplementary Material</title>

<para>

<emphasis>
In Proc. of the 32nd Int. Coll. on Automata, Languages and Programming (ICALP'05)</emphasis>, volume 3580 of <emphasis>Lecture Notes in Computer Science</emphasis>, pages 1263-1274. Springer-Verlag, July 2005. 

(<ulink url="http://www.nada.kth.se/~isaac/publications/FNJ_final_icalp.pdf" >PDF</ulink>,<ulink url="http://www.springerlink.com/openurl.asp?genre=article&amp;issn=0302-9743&amp;volume=3580&amp;spage=1263" >Springer</ulink>)


      </para>


<para>
Slides from presentation at Technion, Israel 2006 

(<ulink url="http://www.nada.kth.se/~isaac/publications/FNJ_technion_presentation.pdf" >PDF</ulink>)
      </para>

<para>

Slides from presentation at ICALP 2005 

(<ulink url="http://www.nada.kth.se/~isaac/publications/FNJ_icalp_presentation.pdf" >PDF</ulink>)
      </para>
      <para>
<ulink url="http://scholar.google.com/scholar?hl=en&amp;lr=&amp;q=Elias+Lagergren+Fast+Neighbor+Joining+" >Google scholar citations</ulink> 

<ulink url="http://citeseer.ist.psu.edu/744489.html" >Go Citeseer</ulink> 
      </para>
  </sect2>




  </sect1>



  <sect1 id="download">
    <title>Download</title>
    <para>      
       Download the software from the <ulink url="http://sourceforge.net/projects/fastphylo">sourceforge</ulink> project page. 

( The software is currently in the subversion repository )

<!-- The latest version of fastphylo is @PACKAGE_VERSION@. -->
    </para>
  </sect1>
  <sect1 id="installation">
    <title>Installation</title>
    <sect2 id="building_from_source">
      <title>Building from source</title>
      <para>To build fastphylo you need to have this installed
        <itemizedlist mark="bullet">
          <listitem>
            <para>
              <ulink url="http://www.cmake.org">cmake</ulink>
            </para>
          </listitem>
<!--
          <listitem>
            <para>
              <ulink url="http://xmlsoft.org/">libxml2</ulink>

            </para>
          </listitem>

-->

        </itemizedlist>
      </para>
      <para>

If you have the fastphylo source code in the directory <filename>/tmp/fastphylo</filename> and you want to install fastphylo into the directory <filename>/tmp/install</filename>, you 

First run <command>cmake</command> then <command>make</command> and then <command>make install</command>
<programlisting><![CDATA[
$ mkdir /tmp/build
$ cd /tmp/build
$ cmake -DCMAKE_INSTALL_PREFIX=/tmp/install /tmp/source && make && make install
-- Configuring done
-- Generating done
-- Build files have been written to: /tmp/build
Scanning dependencies of target fastdist
[  3%] Building CXX object src/c++/CMakeFiles/fastdist.dir/programs/fastdist.o
[  6%] Building CXX object src/c++/CMakeFiles/fastdist.dir/BitVector.o
[  9%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Exception.o
[ 12%] Building CXX object src/c++/CMakeFiles/fastdist.dir/InitAndPrintOn_utils.o
[ 15%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Object.o
[ 18%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Sequence.o
[ 21%] Building CXX object src/c++/CMakeFiles/fastdist.dir/SequenceTree.o
[ 25%] Building CXX object src/c++/CMakeFiles/fastdist.dir/SequenceTree_MostParsimonious.o
[ 28%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Simulator.o
[ 31%] Building CXX object src/c++/CMakeFiles/fastdist.dir/arg_utils_ext.o
[ 34%] Building CXX object src/c++/CMakeFiles/fastdist.dir/file_utils.o
[ 37%] Building CXX object src/c++/CMakeFiles/fastdist.dir/stl_utils.o
[ 40%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/DNA_b128_String.o
[ 43%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/Sequences2DistanceMatrix.o
[ 46%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_LeafLifting.o
[ 50%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_given_edge_probabilities.o
[ 53%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_local_improve.o
[ 56%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_star.o
[ 59%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/Big_AML.o
[ 62%] Building CXX object src/c++/CMakeFiles/fastdist.dir/distance_methods/LeastSquaresFit.o
[ 65%] Building CXX object src/c++/CMakeFiles/fastdist.dir/distance_methods/NeighborJoining.o
[ 68%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/Kimura2parameter.o
[ 71%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/TamuraNei.o
[ 75%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/ambiguity_nucleotide.o
[ 78%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/dna_pairwise_sequence_likelihood.o
[ 81%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/string_compare.o
[ 84%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DistanceMatrix.o
[ 87%] Building C object src/c++/CMakeFiles/fastdist.dir/arg_utils.o
cc1: warning: command line option "-fno-default-inline" is valid for C++/ObjC++ but not for C
[ 90%] Building C object src/c++/CMakeFiles/fastdist.dir/std_c_utils.o
cc1: warning: command line option "-fno-default-inline" is valid for C++/ObjC++ but not for C
[ 93%] Building C object src/c++/CMakeFiles/fastdist.dir/DNA_b128/sse2_wrapper.o
[ 96%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/computeTAMURANEIDistance_DNA_b128_String.o
[100%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/computeDistance_DNA_b128_String.o
Linking CXX executable fastdist
[100%] Built target fastdist
[100%] Built target fastdist
Linking CXX executable CMakeFiles/CMakeRelink.dir/fastdist
Install the project...
-- Install configuration: ""
-- Install configuration: ""
-- Installing /tmp/install/bin/fastdist
-- Install configuration: ""
]]></programlisting>
      </para>
    </sect2>
  </sect1>

</article>
