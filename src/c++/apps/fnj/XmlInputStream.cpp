#include "XmlInputStream.hpp"
#include "fastphylo/core/Exception.hpp"
#include <cstdio>

using namespace std;

XmlInputStream::~XmlInputStream()
{
    if (reader != nullptr)
    {
        xmlFreeTextReader(reader);
    }
    xmlCleanupParser();
}

XmlInputStream::XmlInputStream(char *filename)
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
        THROW_EXCEPTION("Could not open file");
    l.in_root = false;
    l.in_runs = false;
    l.in_run = false;
    l.in_identities = false;
    l.in_identity = false;
    l.row_nr = -1;
    l.entry_nr = -1;
    dmSize = 0;
    xmlRelaxNGParserCtxtPtr parserctxt;
    size_t len = strlen(fastphylo_distance_matrix_xml_relaxngstr);
    parserctxt = xmlRelaxNGNewMemParserCtxt(fastphylo_distance_matrix_xml_relaxngstr, len);
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

std::optional<readstatus> XmlInputStream::handleEntryNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                          StrDblMatrix &dm)
{
    if (!(l.in_root && l.in_runs && l.in_run && l.in_dms && l.in_dm && l.in_row && depth == 6 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("entry")) != 0 && type == XML_READER_TYPE_ELEMENT))
    {
        return std::nullopt;
    }
    l.entry_nr++;
    xmlChar *distanceStr = xmlTextReaderReadString(reader);
    char *endptr = nullptr;
    float distance = strtof(reinterpret_cast<char *>(distanceStr), &endptr);
    if (endptr == reinterpret_cast<char *>(distanceStr))
    {
        THROW_EXCEPTION("expected a floating-point distance, got \"" << distanceStr << "\"");
    }
    dm.setDistance(l.row_nr, l.entry_nr, distance);
    xmlFree(distanceStr);
    return std::nullopt;
}

std::optional<readstatus> XmlInputStream::handleRowNode(int depth, xmlReaderTypes type, const xmlChar *name)
{
    if (!(l.in_root && l.in_runs && l.in_run && l.in_dms && l.in_dm && depth == 5 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("row")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT:
        l.in_row = true;
        l.row_nr++;
        l.entry_nr = -1;
        return std::nullopt;
    case XML_READER_TYPE_END_ELEMENT:
        l.in_row = false;
        return std::nullopt;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<readstatus> XmlInputStream::handleDmNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                       std::vector<std::string> &names, StrDblMatrix &dm)
{
    if (!(l.in_root && l.in_runs && l.in_run && l.in_dms && depth == 4 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("dm")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT:
        dm.resize(dmSize);
        l.in_dm = true;
        l.row_nr = -1;
        return std::nullopt;
    case XML_READER_TYPE_END_ELEMENT:
        l.in_dm = false;
        for (size_t namei = 0; namei < names.size(); namei++)
        {
            dm.setIdentifier(namei, names[namei]);
        }
        return DM_READ;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<readstatus> XmlInputStream::handleIdentityNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                             std::vector<std::string> &names, Extrainfos &extrainfos,
                                                             int &nr_of_ids)
{
    if (!(l.in_root && l.in_runs && l.in_run && l.in_identities && depth == 4 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("identity")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT: {
        nr_of_ids++;
        l.in_identity = true;
        extrainfos.emplace_back();
        xmlChar *nameStr = xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("name"));
        names.emplace_back(reinterpret_cast<char *>(nameStr));
        xmlFree(nameStr);
        return std::nullopt;
    }
    case XML_READER_TYPE_END_ELEMENT:
        l.in_identity = false;
        return std::nullopt;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<readstatus> XmlInputStream::handleExtrainfoNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                              Extrainfos &extrainfos)
{
    if (!(l.in_root && l.in_runs && l.in_run && l.in_identities && l.in_identity && depth == 5 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("extrainfo")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT: {
        xmlChar *outerStr = xmlTextReaderReadOuterXml(reader);
        extrainfos.back() = reinterpret_cast<char *>(outerStr);
        xmlFree(outerStr);
        return std::nullopt;
    }
    case XML_READER_TYPE_END_ELEMENT:
        return std::nullopt;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<readstatus> XmlInputStream::handleDmsNode(int depth, xmlReaderTypes type, const xmlChar *name)
{
    if (!(l.in_root && l.in_runs && l.in_run && depth == 3 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("dms")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT:
        l.in_dms = true;
        return std::nullopt;
    case XML_READER_TYPE_END_ELEMENT:
        l.in_dms = false;
        return std::nullopt;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<readstatus> XmlInputStream::handleRunNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                        std::string &runId)
{
    if (!(l.in_root && l.in_runs && depth == 2 && xmlStrEqual(name, reinterpret_cast<const xmlChar *>("run")) != 0))
    {
        return std::nullopt;
    }
    if (type == XML_READER_TYPE_ELEMENT)
    {
        l.in_run = true;
        xmlChar *dimStr = xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("dim"));
        xmlChar *idStr = xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("id"));
        if (dimStr == nullptr)
        {
            THROW_EXCEPTION("failed to read attribute \"dim\"");
        }
        if (idStr == nullptr)
        {
            THROW_EXCEPTION("failed to read attribute \"id\"");
        }
        {
            char *dimEndptr = nullptr;
            dmSize = strtol(reinterpret_cast<const char *>(dimStr), &dimEndptr, 10);
            if (dimEndptr == reinterpret_cast<const char *>(dimStr))
            {
                THROW_EXCEPTION("expected an integer \"dim\" attribute, got \"" << dimStr << "\"");
            }
        }
        runId = reinterpret_cast<const char *>(idStr);
        xmlFree(idStr);
        xmlFree(dimStr);
        return std::nullopt;
    }
    if (type == XML_READER_TYPE_END_ELEMENT)
    {
        l.in_run = false;
        return END_OF_RUN;
    }
    return std::nullopt;
}

std::optional<readstatus> XmlInputStream::handleIdentitiesNode(int depth, xmlReaderTypes type, const xmlChar *name,
                                                               std::vector<std::string> &names, Extrainfos &extrainfos)
{
    if (!(l.in_root && l.in_runs && l.in_run && depth == 3 &&
          xmlStrEqual(name, reinterpret_cast<const xmlChar *>("identities")) != 0))
    {
        return std::nullopt;
    }
    switch (type)
    {
    case XML_READER_TYPE_ELEMENT:
        l.in_identities = true;
        names.clear();
        extrainfos.clear();
        return std::nullopt;
    case XML_READER_TYPE_END_ELEMENT:
        l.in_identities = false;
        return std::nullopt;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<readstatus> XmlInputStream::handleRunsNode(int depth, xmlReaderTypes type, const xmlChar *name)
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
        return END_OF_RUNS;
    default:
        return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
    }
}

std::optional<readstatus> XmlInputStream::handleRootNode(int depth, xmlReaderTypes type, const xmlChar *name)
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

// Lint Phase 5 (lint_phase5_refactors.md): this was a single 172-line
// function (cognitive complexity 82, threshold 25) - a chain of ten
// "if (right depth/parents/name) { switch(type) {...} }" blocks, one
// per XML element this reader recognizes. Split into the ten
// handleXNode() methods above; readDM() itself is now just dispatch.
readstatus XmlInputStream::readDM(StrDblMatrix &dm, std::vector<std::string> &names, std::string &runId,
                                  Extrainfos &extrainfos)
{
    int nr_of_ids = 0;
    int ret;

    while ((ret = xmlTextReaderRead(reader)) == 1)
    {
        if (xmlTextReaderIsValid(reader) != 1)
        {
            THROW_EXCEPTION("xml input does not validate");
        }
        int depth = xmlTextReaderDepth(reader);
        auto type = static_cast<xmlReaderTypes>(xmlTextReaderNodeType(reader));
        const xmlChar *name = xmlTextReaderConstName(reader);

        if (auto status = handleEntryNode(depth, type, name, dm))
        {
            return *status;
        }
        if (auto status = handleRowNode(depth, type, name))
        {
            return *status;
        }
        if (auto status = handleDmNode(depth, type, name, names, dm))
        {
            return *status;
        }
        if (auto status = handleIdentityNode(depth, type, name, names, extrainfos, nr_of_ids))
        {
            return *status;
        }
        if (auto status = handleExtrainfoNode(depth, type, name, extrainfos))
        {
            return *status;
        }
        if (auto status = handleDmsNode(depth, type, name))
        {
            return *status;
        }
        if (auto status = handleRunNode(depth, type, name, runId))
        {
            return *status;
        }
        if (auto status = handleIdentitiesNode(depth, type, name, names, extrainfos))
        {
            return *status;
        }
        if (auto status = handleRunsNode(depth, type, name))
        {
            return *status;
        }
        if (auto status = handleRootNode(depth, type, name))
        {
            return *status;
        }
    }
    return ERROR;
}
