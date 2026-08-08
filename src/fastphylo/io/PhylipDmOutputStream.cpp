#include "fastphylo/io/PhylipDmOutputStream.hpp"
#include <array>
#include <cstdio>
#include <cmath>

using namespace std;

void PhylipDmOutputStream::print(StrDblMatrix &dm)
{
    printPHYLIPfastSD(dm, fp, Format::Plain, headerWritten);
}

void PhylipDmOutputStream::printSD(StrDblMatrix &dm)
{
    printPHYLIPfastSD(dm, fp, Format::XmlSD, headerWritten);
}

namespace
{

// %f on a finite float needs at most ~39 integer digits (near FLT_MAX) +
// '.' + 6 decimals; 128 is generous headroom over that plus the XML tag
// text. snprintf() is bounds-safe regardless (unlike the fixed-width
// defstr fast path below, this is not the hot path - values this large
// are already unusual, matching the original fprintf()'s "rare" framing).
void appendRareEntry(std::string &row, float f, int intpart, PhylipDmOutputStream::Format format)
{
    std::array<char, 128> rare{};
    int n;
    if (f - (intpart * 1.0) < 0.000001)
    {
        if (format == PhylipDmOutputStream::Format::Plain)
        {
            n = snprintf(rare.data(), rare.size(), "%10d", intpart);
        }
        else
        {
            n = snprintf(rare.data(), rare.size(), "     <entry>%10d</entry>\n", intpart);
        }
    }
    else if (format == PhylipDmOutputStream::Format::Plain)
    {
        n = snprintf(rare.data(), rare.size(), "%10f", f);
    }
    else
    {
        n = snprintf(rare.data(), rare.size(), "     <entry>%10f</entry>\n", f);
    }
    if (n > 0)
    {
        row.append(rare.data(), n < static_cast<int>(rare.size()) ? n : static_cast<int>(rare.size()) - 1);
    }
}

// Fixed-width fast path for |value| <= 99: builds "   d.dddddd" (phylip)
// or the same digits wrapped in <entry> tags (xml), using the same
// digit-lookup tables as printRow()'s streaming equivalent.
void appendFastEntry(std::string &row, float f, int intpart, PhylipDmOutputStream::Format format)
{
    std::array<char, 11> defstr{}; // = "   .      ";
    defstr[0] = ' ';
    defstr[1] = ' ';
    defstr[3] = '.';
    defstr[10] = 0;

    float decimalpart = f - (1.0 * intpart);
    // warning: this isn't enough to get the correct rounding but it is close
    if (intpart == 0)
    {
        defstr[2] = '0';
    }
    else
    {
        defstr[2] = DataOutputStream::ONEDIGIT[intpart];
        intpart = intpart / 10;
        if (intpart != 0)
        {
            defstr[1] = DataOutputStream::ONEDIGIT[intpart];
        }
    }

    // write 6 decimals part
    int deci = 4;
    while (deci <= 9)
    {
        decimalpart = decimalpart * 100.0;
        int index = static_cast<int>(decimalpart);
        decimalpart = decimalpart - index;
        defstr[deci++] = DataOutputStream::TENDIGIT[index];
        defstr[deci++] = DataOutputStream::ONEDIGIT[index];
    }

    if (format == PhylipDmOutputStream::Format::Plain)
    {
        row.append(defstr.data(), 10);
    }
    else
    {
        // skip leading spaces
        int skip = 0;
        while (defstr[skip] == ' ')
        {
            skip++;
        }
        row += "     <entry>";
        row += &defstr[skip];
        row += "</entry>\n";
    }
}

// One matrix entry: the not-finite warning case, then dispatch to the
// fast (|value| <= 99) or rare path. This is what printPHYLIPfastSD()'s
// innermost loop originally did inline, re-checking "writeXml ||
// writeXmlSD" at every step - the main source of that function's
// cognitive complexity (86).
void appendEntry(std::string &row, float f, PhylipDmOutputStream::Format format)
{
    if (!isfinite(f))
    {
        USER_WARNING("warning float not finite (use fix factor) " << f);
        row += format == PhylipDmOutputStream::Format::Plain ? "        -1" : "     <entry>-1</entry>\n";
        return;
    }
    // warning: this isn't enough to get the correct rounding but it is close
    f += 0.0000005;
    int intpart = static_cast<int>(f);
    if (intpart > 99)
    {
        appendRareEntry(row, f, intpart, format);
    }
    else
    {
        appendFastEntry(row, f, intpart, format);
    }
}

} // namespace

// Builds each row into an in-memory buffer and issues one fwrite() per
// row, instead of one fwrite()/fprintf() per matrix entry. Profiling a
// 2500-sequence run (see phase0_audit.md's "Post-Phase-3 finding") found
// this function dominated by stdio's per-call lock/unlock overhead
// (phylip) and fprintf's format-string parsing (XML) - both scale with
// the number of I/O calls, which this reduces from O(numNodes^2) to
// O(numNodes). Every formatted byte written is otherwise identical to
// before; only how it reaches `out` has changed. Originally written for
// fastprot only; now shared with fastdist (Layout Phase C), which was
// still on the pre-optimization per-entry approach.
void PhylipDmOutputStream::printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out, Format format, bool skipLeadingHeader)
{

    const size_t numNodes = dm.getSize();
    const bool xml = format != Format::Plain;

    if (format == Format::Xml)
    {
        fprintf(out, "   <dm>\n");
    }
    else if (format == Format::XmlSD)
    {
        fprintf(out, "   <sdm>\n");
    }
    else if (!skipLeadingHeader)
    {
        fprintf(out, "%5lu\n", numNodes);
    }

    size_t entriesPerRow = xml ? 0 : numNodes;

    string row; // reused per row; one fwrite() flushes it at the end of each row
    row.reserve((xml ? numNodes : 0) * 26 + 32);

    for (size_t i = 0; i < numNodes; i++)
    {
        row.clear();

        if (xml)
        {
            row += "    <row>\n";
            entriesPerRow++;
        }
        else
        {
            // Matches the original fprintf(out,"%-10s", name.c_str()): left-
            // justify, pad to a MINIMUM width of 10 - longer names are not
            // truncated.
            const string &name = dm.getIdentifier(i);
            row += name;
            if (name.size() < 10)
            {
                row.append(10 - name.size(), ' ');
            }
        }

        for (size_t j = 0; j < entriesPerRow; j++)
        {
            appendEntry(row, dm.getDistance(i, j), format);
        }

        row += xml ? "    </row>\n" : "\n";

        fwrite(row.data(), sizeof(char), row.size(), out);
    }
    if (format == Format::Xml)
    {
        fprintf(out, "   </dm>\n");
    }
    else if (format == Format::XmlSD)
    {
        fprintf(out, "   </sdm>\n");
    }
}

// fastdist's bootstrap-streaming path: writes one row at a time as it's
// computed, rather than building the whole matrix first (mem_eff_flag).
// Not exercised by RunExamples.sh, so kept behaviorally identical to the
// original per-entry implementation rather than risking an unverifiable
// rewrite - see the header's comment.
void PhylipDmOutputStream::printRow(StrFloRow &dm, string name, size_t row, bool mem_eff_flag)
{

    const size_t numNodes = dm.getColumns();

    std::array<char, 11> defstr{};
    defstr[0] = ' ';
    defstr[3] = '.';
    defstr[10] = 0;

    const size_t entriesPerRow = numNodes;
    fprintf(fp, "%-10s", name.c_str());

    float f = 0.0;
    if (!mem_eff_flag)
    {
        for (size_t i = 0; i < row; ++i)
        {
            fprintf(fp, "%10f", f);
        }
    }
    else
    {
        row = 0;
    }

    for (size_t j = row; j < entriesPerRow; j++)
    {

        float f = dm.getDistance(j);

        if (!isfinite(f))
        {
            USER_WARNING("warning float not finite (use fix factor) " << f);
            fprintf(fp, "        -1");
            continue;
        }
        f += 0.0000005;
        defstr[1] = ' ';
        int intpart = static_cast<int>(f);
        if (intpart > 99)
        {
            if (f - (intpart * 1.0) < 0.000001)
            {
                fprintf(fp, "%10d", intpart);
                continue;
            }
            fprintf(fp, "%10f", f);
            continue;
        }

        float decimalpart = f - (1.0 * intpart);
        if (intpart == 0)
        {
            defstr[2] = '0';
        }
        else
        {
            defstr[2] = ONEDIGIT[intpart];
            intpart = intpart / 10;
            if (intpart != 0)
            {
                defstr[1] = ONEDIGIT[intpart];
            }
        }

        int deci = 4;
        while (deci <= 9)
        {
            decimalpart = decimalpart * 100.0;
            int index = static_cast<int>(decimalpart);
            decimalpart = decimalpart - index;
            defstr[deci++] = TENDIGIT[index];
            defstr[deci++] = ONEDIGIT[index];
        }

        fwrite(defstr.data(), sizeof(char), 10, fp);
    }

    fprintf(fp, "\n");
}

void PhylipDmOutputStream::printHeader(size_t numNodes)
{
    fprintf(fp, "%5d\n", static_cast<int>(numNodes));
    headerWritten = true;
}
