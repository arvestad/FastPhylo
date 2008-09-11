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

  <sect1 id="algorithms">
    <title>Algorithms</title>


  <sect2 id="Fast_Computation_of_Distance_Estimators">
    <title>Fast Computation of Distance Estimators</title>


  <sect3 id="about_the_article_Fast_Computation_of_Distance_Estimators">
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



  </sect3>


  <sect3 id="abstract_of_the_article_Fast_Computation_of_Distance Estimators">
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



  </sect3>



  <sect3 id="Supplementary_Material">
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



  </sect3>




  </sect2>


  <sect2 id="Fast_Neighbor_Joining">
    <title>Fast Neighbor Joining</title>


  <sect3 id="about_the_article_Fast_Neighbor_Joining">
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



  </sect3>


  <sect3 id="abstract_of_the_article_Fast_Neighbor_Joining">
    <title>Abstract of the published article</title>

<para>
Reconstructing the evolutionary history of a set of species is a fundamental problem in biology and methods for solving this problem are gaged based on two characteristics: accuracy and efficiency. Neighbor Joining (NJ) is a so-called distance-based method that, thanks to its good accuracy and speed, has been embraced by the phylogeny community. It takes the distances between <emphasis>n</emphasis> taxa and produces in <emphasis>&#920;(n<superscript>3</superscript>)</emphasis> time a phylogenetic tree, i.e., a tree which aims to describe the evolutionary history of the taxa. In addition to performing well in practice, the NJ algorithm has optimal reconstruction radius.
</para><para>
The contribution of this paper is twofold: (1) we present an algorithm called Fast Neighbor Joining (FNJ) with optimal reconstruction radius and optimal run time complexity <emphasis>O(n<superscript>2</superscript>)</emphasis> and (2) we present a greatly simplified proof for the correctness of NJ. Initial experiments show that FNJ in practice has almost the same accuracy as NJ, indicating that the property of optimal reconstruction radius has great importance to their good performance. Moreover, we show how improved running time can be achieved for computing the so-called correction formulas.
      </para>



  </sect3>



  <sect3 id="Supplementary_Material_Fast_Neighbor_Joining">
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
  </sect3>




  </sect2>
</sect1>


  <sect1 id="software">
    <title>Software</title>


  <sect2 id="download">
    <title>Download</title>
    <para>      
       Download the software from the <ulink url="http://sourceforge.net/projects/fastphylo">sourceforge</ulink> project page. 

( The software is currently in the subversion repository )

<!-- The latest version of fastphylo is @PACKAGE_VERSION@. -->
    </para>
  </sect2>
  <sect2 id="installation">
    <title>Installation</title>
    <sect3 id="building_from_source">
      <title>Building from source</title>
      <para>To build fastphylo you need to have this installed
        <itemizedlist mark="bullet">
          <listitem>
            <para>
              <ulink url="http://www.cmake.org">cmake</ulink>
            </para>
          </listitem>

          <listitem>
            <para>
              <ulink url="http://xmlsoft.org/">libxml2</ulink> ( You can remove this requirement by passing the option -DWITH_LIBXML=OFF to <application>cmake</application> )

            </para>
          </listitem>



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

If you want to build the html documentation ( i.e. this page ) you need to pass the option -DBUILD_DOCBOOK=ON to <application>cmake</application>.
      </para>
    </sect3>
  </sect2>



  <sect2 id="usage">
    <title>Usage</title>


  <sect3 id="fastdist">
    <title>fastdist</title>

<para><application>fastdist</application> implements the algorithm Fast Computation of Distance Estimators ( see <xref linkend="Fast_Computation_of_Distance_Estimators"/> )
</para>


  <sect4 id="fastdist_command_line_options">
    <title>Command line options</title>


<para>

Type <userinput>fastdist --help</userinput> to see the command line options

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>
</para></sect4>




  <sect4 id="fastdist_input_file_formats">
    <title>fastdist input file formats</title>

<para>

<table frame='all'><title>fastdist input file formats</title>
<tgroup cols='3' align='left' colsep='1' rowsep='1'>
<thead>
<row>
<entry>file format</entry>
<entry>short option</entry>
<entry>description</entry>
</row>
</thead>
<tbody>
<row>
<entry>fasta format</entry>
<entry>-I xml</entry>
<entry><xref linkend="fasta_format"/></entry>
</row>
<row>
<entry>phylip multialignment format</entry>
<entry>-I phylip_multialignment_format</entry>
<entry><xref linkend="phylip_multialignment_format"/></entry>

</row>
<row>
<entry>fastphylo sequence XML format</entry>
<entry>-I fastphylo_sequence_xml_format</entry>
<entry><xref linkend="fastphylo_sequence_xml_format"/></entry>

</row>
</tbody>
</tgroup>
</table>


</para>





  </sect4>

  <sect4 id="fastdist_output_formats">
    <title>fastdist output file formats</title>


<para>

<table frame='all'><title>fastdist output file formats</title>
<tgroup cols='3' align='left' colsep='1' rowsep='1'>
<thead>
<row>
<entry>file format</entry>
<entry>short option</entry>
<entry>description</entry>
</row>
</thead>
<tbody>
<row>
<entry>fastphylo sequence XML format</entry>
<entry>-O xml</entry>
<entry><xref linkend="fastphylo_distance_matrix_xml_format"/></entry>

</row>
<row>
<entry>phylip distance matrix format</entry>
<entry>-O phylip_dm</entry>
<entry><xref linkend="phylip_distance_matrix_format"/></entry>

</row>
</tbody>
</tgroup>
</table>


</para>







  </sect4>

  <sect4 id="fastdist_examples">
    <title>Examples</title>
<para>
<example id="fastdist_phylip_multialignment"><title>fastdist with input in file Phylip multialignment format</title><para>
We use the file described in <xref linkend="seq.phylip_multialignment"/> as input file. 
The file has two datasets so we pass the option <userinput>-r 2</userinput> to <application>fastdist</application>. Per default the output is given in XML format 

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_seq.phylip_multialignment">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_seq.phylip_multialignment">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para>
</example>
<example id="fastdist_seq.fasta"><title>fastdist with input in file fasta format</title><para>
We use the file described in <xref linkend="seq.fasta"/> as input file. 
Per default the output is given in XML format 

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_seq.fasta">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_seq.fasta">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para>
</example>

<example id="example_fastdist_seq_xml"><title>fastdist with input file in XML format</title><para>



We use the file described in <xref linkend="seq.xml"/> as input file. 

<note><para>The -r option can only be used if the input is in phylip_multialignment format. <application>fastdist</application> will for XML files compute all data sets ( runs ). Fasta files can only contain one data set so the -r option does not make any sense there.</para></note>

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_seq.xml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_seq.xml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>


</para>
</example>

<example id="example_fastdist_with_an_xml_stream_on_stdin"><title>fastdist with an XML stream on stdin</title><para>

If you leave out the input filename, the input will be read from stdin. <application>fastdist</application> doesn't wait for the whole xml file to be read before it starts. It starts a computation as soon as an ending <![CDATA[ </run>   ]]>
has been read. The memory consumption will not grow over time so the input can be arbitrarily large. A never ending input stream only works in the fastphylo sequence XML format, because the phylip multialignment input format needs you to specify in advance how many data sets are to be sent to <application>fastdist</application> ( the -r option ).



<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_xml_stdin_xml_output">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_xml_stdin_xml_output">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para>
</example>



<example id="example_reading_the_fastdist_xml_output_stream_with_python"><title>reading the fastdist XML output stream with python</title><para>

If the XML output is very large you might want to use an XML parser that doesn't hold the whole file in memory. This python script is an example of how to do this

<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/fastdist_lxml.py">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>


For each distance matrix the script counts the number of elements with a value below 0.1

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_lxml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_lxml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

Read more about <ulink url="http://codespeak.net/lxml/parsing.html">lxml</ulink> and <ulink url="http://www.zvon.org/xxl/XPathTutorial/General/examples.html">xpath</ulink>.


</para>
</example>


      </para>

    </sect4>

    </sect3>


  <sect3 id="fnj">
    <title>fnj</title>
<para><application>fnj</application> implements the algorithm Fast Neighbor Joining ( see <xref linkend="Fast_Neighbor_Joining"/> )
</para>



  <sect4 id="fnj_command_line_options">
    <title>Command line options</title>


<para>

Type <userinput>fnj --help</userinput> to see the command line options

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fnj_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fnj_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para></sect4>




  <sect4 id="fnj_input_file_formats">
    <title>fnj input file formats</title>





<table frame='all'><title>fnj input file formats</title>
<tgroup cols='3' align='left' colsep='1' rowsep='1'>
<thead>
<row>
<entry>file format</entry>
<entry>short option</entry>
<entry>description</entry>
</row>
</thead>
<tbody>
<row>
<entry>fastphylo sequence XML format</entry>
<entry>-I xml</entry>
<entry><xref linkend="fastphylo_distance_matrix_xml_format"/></entry>

</row>
<row>
<entry>phylip distance matrix format</entry>
<entry>-I phylip_dm</entry>
<entry><xref linkend="phylip_distance_matrix_format"/></entry>

</row>
</tbody>
</tgroup>
</table>



  </sect4>

  <sect4 id="fnj_output_formats">
    <title>fnj output file formats</title>



<table frame='all'><title>fnj output file formats</title>
<tgroup cols='3' align='left' colsep='1' rowsep='1'>
<thead>
<row>
<entry>file format</entry>
<entry>short option</entry>
<entry>description</entry>
</row>
</thead>
<tbody>
<row>
<entry>fastphylo count tree XML format</entry>
<entry>-O xml</entry>
<entry><xref linkend="fastphylo_tree_count_xml_format"/></entry>

</row>
</tbody>
</tgroup>
</table>





  </sect4>




  <sect4 id="fnj_examples">
    <title>Examples</title>
<para>
<example id="fnj_phylip.dm"><title>fnj with input file in Phylip distance matrix format</title><para>
We use the file described in <xref linkend="dm.phylip_dm"/> as input file. The file has two datasets so we pass the option <userinput>-r 2</userinput> to <application>fnj</application>. Per default the output is given in the "fastphylo count tree XML format" ( -O xml ).

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fnj_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fnj_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>
</para>
</example>

<example id="fnj_dm.xml"><title>fnj with input file in XML format</title><para>

We use the file described in <xref linkend="dm.xml"/> as input file. Per default the output is given in the "fastphylo count tree XML format" ( -O xml ).


<note><para>The -r option is not available and also not needed when the input is in XML format. <application>fnj</application> computes all data sets ( runs ).</para></note>

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fnj_dm.xml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fnj_dm.xml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>


</para>
</example>

<example id="example_connecting_fastdist_to_fnj_with_a_pipe"><title>connecting fastdist to fnj with a pipe</title><para>

We use the file described in <xref linkend="seq.phylip_multialignment"/> as input file. 
The file has two data sets. We will bootstrap 3 times.

First we send the data in  phylip_dm format through the pipe:

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/connecting_fnj_to_fastdist_with_a_pipe1">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/connecting_fnj_to_fastdist_with_a_pipe1">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>


We could also send the data in XML format through the pipe:
<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/connecting_fnj_to_fastdist_with_a_pipe2">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/connecting_fnj_to_fastdist_with_a_pipe2">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

As the the XML format is more descriptive, the flags -d and -r are no longer needed by <application>fnj</application>.



</para>
</example>



<example id="example_reading_the_fnj_xml_output_stream_with_python"><title>reading the fnj XML output stream with python</title><para>


If the XML output is very large you might want to use an XML parser that doesn't hold the whole file in memory. This python script is an example of how to do this

<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/fnj_lxml.py">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>


The script prints the maximum count ( just as an example ).

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fnj_lxml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fnj_lxml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

Read more about <ulink url="http://codespeak.net/lxml/parsing.html">lxml</ulink> and <ulink url="http://www.zvon.org/xxl/XPathTutorial/General/examples.html">xpath</ulink>.


</para>
</example>


      </para>

    </sect4>

    </sect3>



    </sect2>

  <sect2 id="file_formats">
    <title>File formats</title>
    <para>      
This software package handles the following file formats
</para>


  <sect3 id="fastphylo_sequence_xml_format">
    <title>Fastphylo sequence XML format</title>
<para>
The Fastphylo sequence XML format is chosen by the option <userinput>-I xml</userinput> to fastdist. 

Type <userinput>fastdist --print-relaxng-input</userinput> to see its relaxng schema

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_print_relaxng_input">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_print_relaxng_input">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>



The Relax NG schema specifies that the extrainfo element is optional and can be inserted as a child to a seq element. The extrainfo element may contain any content and will be passed on to the output XML format.


<example id="seq.xml"><title>seq.xml, an example file in  Fastphylo sequence XML format</title><para>

The example file <ulink url="example_files/seq.xml">seq.xml</ulink> contains 

<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.xml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>
</para></example>

</para></sect3>


  <sect3 id="phylip_multialignment_format">
    <title>Phylip multialignment format</title>
<para>
The Phylip multialignment input format is chosen by the option <userinput>-I phylip_multialignment</userinput> to fastdist. 

This format has an extra option -r ( --number-of-run ) for specifying the number of data sets in the input.





<example id="seq.phylip_multialignment"><title>seq.phylip_multialignment, an example file in Phylip multialignment format</title><para>
The example file <ulink url="example_files/seq.phylip_multialignment">seq.phylip_multialignment</ulink> contains

<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.phylip_multialignment">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>
The file has two data sets.



</para></example>
</para></sect3>


  <sect3 id="fasta_format">
    <title>Fasta format</title>
<para>
The Fasta input format is chosen by the option <userinput>-I fasta</userinput> to fastdist. 
Fasta files can only contain one data set. Read more about the <ulink url="http://en.wikipedia.org/wiki/Fasta_format">Fasta format</ulink> on Wikipedia.
The parser will take the whole header line as the sequence identifier name, i.e. all characters after the greater-than character ( ">" ).




<example id="seq.fasta"><title>seq.fasta, an example file in fasta format</title><para>
The example file <ulink url="example_files/seq.fasta">seq.fasta</ulink> contains


<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.fasta">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>

</para></example>


</para></sect3>

  <sect3 id="fastphylo_distance_matrix_xml_format">
    <title>Fastphylo distance matrix XML format</title>
<para>
The Fastphylo sequence XML format is chosen by the option <userinput>-O xml</userinput> to fastdist and the option <userinput>-I xml</userinput> to fnj. 

Type <userinput>fastdist --print-relaxng-output</userinput> to see its relaxng schema

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_print_relaxng_output">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_print_relaxng_output">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

The Relax NG schema specifies that the extrainfo element is optional and can be inserted as a child to a seq element. The extrainfo element may contain any content.





<example id="dm.xml"><title>dm.xml, an example file in  Fastphylo distance matrix XML format</title><para>
The example file <ulink url="example_files/dm.xml">dm.xml</ulink> contains

<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/dm.xml">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>
</para></example>


</para></sect3>

  <sect3 id="phylip_distance_matrix_format">
    <title>Phylip distance matrix format</title>
<para>
The Phylip distance matrix format is chosen by the option <userinput>-O phylip_dm</userinput> to fastdist or the option  <userinput>-I phylip_dm</userinput> to fnj.


<example id="dm.phylip_dm"><title>dm.phylip_dm, an example file in phylip distance matrix format</title><para>
The example file <ulink url="example_files/dm.phylip_dm">dm.phylip_dm</ulink> contains



<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>

It contains two data sets.
</para>
</example>

</para></sect3>

  <sect3 id="fastphylo_tree_count_xml_format">
    <title>Fastphylo tree count XML format</title>
<para>
The Fastphylo tree count XML format is chosen by the option <userinput>-O xml</userinput> to fnj. 
You can see an example of the format in the example <xref linkend="fnj_dm.xml"/>.

Type <userinput>fnj --print-relaxng-output</userinput> to see the formats relaxng schema. 
<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fnj_print_relaxng_output">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fnj_print_relaxng_output">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>


</para></sect3>



  </sect2>


  </sect1>






</article>




