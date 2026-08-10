#include "fastphylo/io/XmlInputStream.hpp"
#include "fastphylo/core/Exception.hpp"
#include <cstdio>
#include <cstring>
#include <libxml/relaxng.h>

using namespace std;

XmlSequenceReader::~XmlSequenceReader()
{
    if (reader != nullptr)
    {
        xmlFreeTextReader(reader);
    }
    xmlCleanupParser();
}

XmlSequenceReader::XmlSequenceReader(char *filename, const char *relaxngSchemaStr)
{
    //
    //  Re: [xml] Why does "-" read from stdin?
    //  http://mail.gnome.org/archives/xml/2007-February/msg00005.html
    //
    // Comment: The special treatment of the filename "-" in the libxml api, is not
    // a good designed api, but now when it is there let us use it.

    if ((filename != nullptr) && (strncmp(filename, "-", 1) == 0))
    {
        THROW_EXCEPTION("file name \"-\" is not allowed. See \n "
                        "http://mail.gnome.org/archives/xml/2007-February/msg00005.html \n  ");
    }

    if (filename == nullptr)
    {
        filename = const_cast<char *>("-");
    }
    LIBXML_TEST_VERSION

    reader = xmlReaderForFile(filename, nullptr, XML_PARSE_COMPACT | XML_PARSE_NONET);
    if (reader == nullptr)
    {
        THROW_EXCEPTION("Could not open file");
    };

    l.in_root = false;
    l.in_runs = false;
    l.in_run = false;
    l.in_seq = false;

    xmlRelaxNGParserCtxtPtr parserctxt;
    size_t len = strlen(relaxngSchemaStr);
    parserctxt = xmlRelaxNGNewMemParserCtxt(relaxngSchemaStr, len);
    xmlRelaxNGSetParserErrors(parserctxt, reinterpret_cast<xmlRelaxNGValidityErrorFunc>(fprintf),
                              reinterpret_cast<xmlRelaxNGValidityWarningFunc>(fprintf), stderr);
    xmlRelaxNGPtr schema = nullptr;
    schema = xmlRelaxNGParse(parserctxt);
    xmlRelaxNGFreeParserCtxt(parserctxt);

    if (xmlTextReaderRelaxNGSetSchema(reader, schema) != 0)
    {
        THROW_EXCEPTION("failed to set relax ng schema");
    }
}

std::optional<bool> XmlSequenceReader::handleSeqNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                     std::vector<Sequence> &seqs, Extrainfos &extrainfos,
                                                     int &numSequences)
{
    if (!(l.in_root && l.in_runs && l.in_run && depth == 3 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("seq")) != 0))
    {
        return std::nullopt;
    }
    if (type == XML_READER_TYPE_ELEMENT)
    {
        l.in_seq = true;
        numSequences++;
        seqs.resize(numSequences);
        extrainfos.emplace_back();
        Sequence &s = seqs[numSequences - 1];

        xmlChar *name = xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("name"));
        xmlChar *seq = xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("seq"));

        if (name == nullptr)
        {
            THROW_EXCEPTION("failed to read attribute \"name\"");
        }
        if (seq == nullptr)
        {
            THROW_EXCEPTION("failed to read attribute \"seq\"");
        }

        s.name = reinterpret_cast<const char *>(name);
        s.seq = reinterpret_cast<const char *>(seq);
        xmlFree(name);
        xmlFree(seq);
        return std::nullopt;
    }
    if (type == XML_READER_TYPE_END_ELEMENT)
    {
        l.in_seq = false;
        return std::nullopt;
    }
    return std::nullopt;
}

std::optional<bool> XmlSequenceReader::handleExtrainfoNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                           Extrainfos &extrainfos)
{
    if (!(l.in_root && l.in_runs && l.in_run && l.in_seq && depth == 4 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("extrainfo")) != 0))
    {
        return std::nullopt;
    }
    if (type == XML_READER_TYPE_ELEMENT)
    {
        xmlChar *outerStr = xmlTextReaderReadOuterXml(reader);
        extrainfos.back() = reinterpret_cast<char *>(outerStr);
        xmlFree(outerStr);
        return std::nullopt;
    }
    if (type == XML_READER_TYPE_END_ELEMENT)
    {
        return std::nullopt;
    }
    return std::nullopt;
}

std::optional<bool> XmlSequenceReader::handleRootNode(int depth, xmlReaderTypes type, const xmlChar *name)
{
    if (!(depth == 0 && xmlStrEqual(name, reinterpret_cast<const xmlChar *>("root")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT:
        l.in_root = true;
        return std::nullopt;
    case XML_READER_TYPE_END_ELEMENT:
        l.in_root = false;
        return std::nullopt;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<bool> XmlSequenceReader::handleRunsNode(int depth, xmlReaderTypes type, const xmlChar *name)
{
    if (!(l.in_root && depth == 1 && xmlStrEqual(name, reinterpret_cast<const xmlChar *>("runs")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT:
        l.in_runs = true;
        return std::nullopt;
    case XML_READER_TYPE_END_ELEMENT:
        l.in_runs = false;
        return std::nullopt;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<bool> XmlSequenceReader::handleRunNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                     std::string &runId, Extrainfos &extrainfos)
{
    if (!(l.in_root && l.in_runs && depth == 2 && xmlStrEqual(name, reinterpret_cast<const xmlChar *>("run")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT: {
        l.in_run = true;
        extrainfos.clear();
        xmlChar *id = xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("id"));
        if (id == nullptr)
        {
            THROW_EXCEPTION("failed to read attribute \"id\"");
        }
        runId = reinterpret_cast<const char *>(id);
        xmlFree(id);
        return std::nullopt;
    }
    case XML_READER_TYPE_END_ELEMENT:
        l.in_run = false;
        return true;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

// Lint Phase 5 (lint_phase5_refactors.md): this was a single 114-line
// function (cognitive complexity 55, threshold 25) - a chain of five
// "if (right depth/parents/name) {...}" blocks, one per XML element
// this reader recognizes. Split into the five handleXNode() methods
// above; readSequences() itself is now just dispatch.
bool XmlSequenceReader::readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos)
{
    int ret;
    int numSequences = 0;
    while ((ret = xmlTextReaderRead(reader)) == 1)
    {
        if (xmlTextReaderIsValid(reader) != 1)
        {
            THROW_EXCEPTION("xml input does not validate");
        }

        int depth = xmlTextReaderDepth(reader);
        auto type = static_cast<xmlReaderTypes>(xmlTextReaderNodeType(reader));
        const xmlChar *name = xmlTextReaderConstName(reader);

        if (auto result = handleSeqNode(depth, type, name, seqs, extrainfos, numSequences))
        {
            return *result;
        }
        if (auto result = handleExtrainfoNode(depth, type, name, extrainfos))
        {
            return *result;
        }
        if (auto result = handleRootNode(depth, type, name))
        {
            return *result;
        }
        if (auto result = handleRunsNode(depth, type, name))
        {
            return *result;
        }
        if (auto result = handleRunNode(depth, type, name, runId, extrainfos))
        {
            return *result;
        }
    }
    if (ret == 0)
    {
        if (l.in_root)
        {
            THROW_EXCEPTION("failed to parse");
        }
        else
        {
            return false;
        }
    }

    return false;
}
