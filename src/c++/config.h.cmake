#cmakedefine WITH_LIBXML

// github_actions_release_builds_plan.md's Phase E: single source of
// truth for every app's --version output, flowed from CMakeLists.txt's
// PACKAGE_VERSION (tag-derived for a release build, git-describe for
// an ordinary dev/CI build) via the same @VAR@ substitution
// configure_file() already does for this file.
#define PACKAGE_VERSION "@PACKAGE_VERSION@"
