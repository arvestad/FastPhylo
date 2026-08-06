#include "fastphylo/io/XmlOutputStream.hpp"
#include <array>
#include <cstdio>
#include <cmath>

using namespace std;

XmlOutputStream::XmlOutputStream(char *filename) : DataOutputStream(filename) {
  fprintf(fp, "<?xml version=\"1.0\"?>\n<root>\n <runs>\n");
}

XmlOutputStream::~XmlOutputStream() {
  fprintf(fp, " </runs>\n</root>\n");
}

void XmlOutputStream::print(StrDblMatrix &dm) {
  PhylipDmOutputStream::printPHYLIPfastSD(dm, fp, PhylipDmOutputStream::Format::Xml);
}

void XmlOutputStream::printSD(StrDblMatrix &dm) {
  PhylipDmOutputStream::printPHYLIPfastSD(dm, fp, PhylipDmOutputStream::Format::XmlSD);
}

void XmlOutputStream::printStartRun(std::vector<string> &names, std::string &runId, Extrainfos &extrainfos) {
  fprintf(fp, "  <run id=\"%s\" dim=\"%i\">\n   <identities>\n", runId.c_str(), static_cast<int>(names.size()));
  for (size_t namei = 0; namei < names.size(); namei++) {
    // This only works if we constrain the input by a schema to not have "<", "&" and such.
    // Otherwise we need to use xmlEncodeSpecialChars(xmlDocPtr doc, const xmlChar * input)
    fprintf(fp, "    <identity name=\"%s\"", names[namei].c_str());

    if (extrainfos.size() > namei && !extrainfos[namei].empty()) {
      fprintf(fp, ">\n     %s\n    </identity>\n", const_cast<char *>(extrainfos[namei].c_str()));
    } else {
      fprintf(fp, "/>\n");
    }
  }
  fprintf(fp, "   </identities>\n   <dms>\n");
}

void XmlOutputStream::printEndRun() {
  fprintf(fp, "   </dms>\n  </run>\n");
}

// fastdist's memory-efficient bootstrap-streaming path (see header comment).
void XmlOutputStream::printRow(StrFloRow &dm, string name, size_t row, bool mem_eff_flag) {

  const size_t numNodes = dm.getColumns();

  std::array<char, 11> defstr{};
  defstr[0] = ' ';
  defstr[3] = '.';
  defstr[10] = 0;

  const size_t entriesPerRow = numNodes;

  fprintf(fp, "    <row>\n");
  if (mem_eff_flag) {
    row = 0;
  }
  for (size_t j = row; j < entriesPerRow; j++) {

    float f = dm.getDistance(j);

    if (!isfinite(f)) {
      USER_WARNING("warning float not finite (use fix factor) " << f);
      fprintf(fp, "     <entry>-1</entry>\n");
      continue;
    }
    f += 0.0000005;
    defstr[1] = ' ';
    int intpart = static_cast<int>(f);
    if (intpart > 99) {
      if (f - (intpart * 1.0) < 0.000001) {
        fprintf(fp, "     <entry>%d</entry>\n", intpart);
        continue;
      }
      fprintf(fp, "     <entry>%f</entry>\n", f);
      continue;
    }

    float decimalpart = f - (1.0 * intpart);
    if (intpart == 0) {
      defstr[2] = '0';
    } else {
      defstr[2] = ONEDIGIT[intpart];
      intpart = intpart / 10;
      if (intpart != 0) {
        defstr[1] = ONEDIGIT[intpart];
      }
    }

    int deci = 4;
    while (deci <= 9) {
      decimalpart = decimalpart * 100.0;
      int index = static_cast<int>(decimalpart);
      decimalpart = decimalpart - index;
      defstr[deci++] = TENDIGIT[index];
      defstr[deci++] = ONEDIGIT[index];
    }

    // skip leading spaces
    int i = 0;
    while (defstr[i] == ' ') {
      i++;
    }
    fprintf(fp, "     <entry>%s</entry>\n", &defstr[i]);
  }

  fprintf(fp, "    </row>\n");
}

void XmlOutputStream::printBootstrapSpliter(size_t numNodes) {
  fprintf(fp, "   </dm>\n");
  fprintf(fp, "   <dm>\n");
}
