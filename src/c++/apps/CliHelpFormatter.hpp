#pragma once

#include "CLI11.hpp"

// gengetopt's --help always led with "Usage: <name> [OPTION]...";
// CLI11's default Formatter drops the "Usage: " label whenever a
// program name is available (which it always is for these apps), so
// this override reinstates it. Shared across every app migrated off
// gengetopt (gengetopt_migration_plan.md) rather than copy-pasted per
// app - fnj was first to need it, fastdist is the second.
class FastphyloHelpFormatter : public CLI::Formatter {
  public:
	std::string make_usage(const CLI::App *app, std::string name) const override {
		return CLI::Formatter::make_usage(app, "Usage: " + name);
	}
};
