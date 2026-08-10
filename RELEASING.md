# Making a release

This is the short, practical version. For the design rationale behind
any of this (why Ubuntu 22.04 instead of latest, why two workflow
triggers, etc.), see `planning/github_actions_release_builds_plan.md`.

## The short version

```
git tag 1.0.4
git push origin 1.0.4
```

That's it. Pushing a tag matching `X.Y.Z` (no `v` prefix — check
`git tag -l` if unsure, this repo's convention is bare numbers like
`1.0.2`, `1.0.3`) triggers `.github/workflows/release.yml`, which:

1. Builds and tests fastphylo on Linux, macOS, and Windows.
2. Packages each platform's binaries into an archive
   (`fastphylo-<version>-<os>.tar.gz`/`.zip`).
3. Publishes a GitHub Release named after the tag, with all three
   archives attached and auto-generated release notes.

No manual version bump anywhere in the source: the version number in
`--version` output and in `CPACK_PACKAGE_VERSION` comes directly from
the tag you pushed (via `-DPACKAGE_VERSION`). There is nothing else to
edit before tagging.

Watch progress under the repo's **Actions** tab, then check the
**Releases** page once it's done (should be 10-15 minutes for all
three platforms).

### Prereleases (beta/rc)

For a prerelease, append `-` and a label to the version, e.g.
`2.0.0-beta.1`. This still matches the release workflow's tag trigger
(it just needs two dots present somewhere in the tag, which any
`X.Y.Z`-prefixed prerelease tag has), and `release.yml` automatically
marks the GitHub Release as a **prerelease** whenever the tag contains
a `-` — no extra flag to remember. A plain `X.Y.Z` tag (no `-`)
publishes as a normal release, same as always.

## Do a dry run first

Before pushing a real tag, you can run the entire pipeline —
build, test, package, upload — without publishing anything public:

1. Go to **Actions → Release** (left sidebar).
2. Click **Run workflow**.
3. Pick the branch you want to build from (e.g. `modernize-cpp17`).
4. Click **Run workflow**.

This runs everything a real release does, except the final "publish a
GitHub Release" step is skipped — you get downloadable build artifacts
from the run page instead. Safe to run as many times as you like.

## Which branch to tag

Tag `master`. The plan for 2.0.0-beta.1 is to merge `modernize-cpp17`
into `master` first, then tag `master` as usual — unlike the 1.0.x
line, there's no reason to keep shipping from a long-lived side branch
once it's ready, and `master` staying frozen on the old layout only
invites confusion for the (small) existing user base.

If you ever do need to tag a branch other than `master` again (a
dry run of work still in progress, say): `master` carries its own
frozen copy of `release.yml`, but only so the **Run workflow** button
is visible in the Actions UI at all — GitHub requires a
`workflow_dispatch` workflow to exist on the default branch to list it
there. That copy is not what actually runs for a dry run against
another branch (the branch you pick supplies its own copy of the
workflow), and it plays no role in a real tag-push release either — a
tag push always uses whatever `release.yml` exists at the tagged
commit. Keep the `master` copy in sync with the real one if you rely
on this.

## If something goes wrong mid-release

If a release run fails partway (say, one platform's build breaks),
nothing has been published yet — the `release` job that creates the
GitHub Release only runs after all three platform builds succeed. Fix
the problem, delete the tag, and re-tag:

```
git tag -d 1.0.4
git push origin :refs/tags/1.0.4
# fix the issue, commit it
git tag 1.0.4
git push origin 1.0.4
```

If a GitHub Release was already published and needs correcting, delete
it from the **Releases** page in the GitHub UI, then delete and re-push
the tag as above.
