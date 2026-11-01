# GitHub and Docker Hub setup

Actions only a repository owner can do. Nothing here is done by agents, and no credential belongs in Git.

## Branch and tag protection

Current `protect_main` ruleset, read from the public API: no deletion or force-push, pull request with 1 approval, and a required `R-CMD-check` check with strict up-to-date branches. Bypass actors are not visible publicly.

Two rules would stall the merge from `dev`:

- **Required approval.** GitHub does not count an author's approval of their own pull request.
- **Strict up-to-date.** `main` has two commits that `dev` lacks (they add and delete one file, with no net effect), so GitHub would ask to merge `main` into `dev` first.

Recommended, under **Settings > Rules > Rulesets**:

| Ruleset | Target | Rules |
|---|---|---|
| `protect_main` (edit) | `main` | Keep pull request required, no deletion, no force-push. Set required approvals to **0** while there is one maintainer; raise it when a second reviewer exists. Keep required checks, adding `R-CMD-check`, `validate`, `conda`, and `build-and-test`. Keep "require branches to be up to date" |
| `protect_tags` (new) | Tags `v*` | Restrict creation, updates, and deletion. Add yourself as a bypass actor, so that only you can create release tags |
| `protect_dev` (new) | `dev` | No deletion, no force-push. Do not require pull requests |

This keeps every merge going through a pull request with passing checks, without a general bypass that could also skip checks. Setting approvals to 0 removes only a review that one person cannot give anyway.

**Merging `main` into `dev`:** the up-to-date rule needs it once, and project policy requires your authorization. It adds a merge commit with no change to `dev`'s content, since `main`'s two extra commits cancel out. After you authorize it, use the "Update branch" button on the pull request.

## Docker Hub publishing

`publish.yaml` is off until you complete these steps. Do this only when ready to release.

1. **Create a token.** On Docker Hub: Account settings > Personal access tokens > Generate. Description `dndsr-github-publish`; permission **Read & Write**, for this purpose only. Copy it once.
2. **Create an environment.** On GitHub: Settings > Environments > New environment, named exactly `dockerhub-release`.
   - Enable **Required reviewers**, and add yourself.
   - Under **Deployment branches and tags**, choose selected tags, with the rule `v*`.
3. **Add secrets to that environment**, not to the repository: `DOCKERHUB_USERNAME` is `scrumpis`, and `DOCKERHUB_TOKEN` is the token.
4. **Enable it last.** Settings > Secrets and variables > Actions > Variables > New repository variable: `DOCKERHUB_PUBLISH_ENABLED` = `true`.

When a tag such as `v0.3.0` is pushed, the workflow:

1. Refuses unless the tag is on `main`, matches `DESCRIPTION`, the skill, and `NEWS.md`, and does not already exist on Docker Hub.
2. Waits for your approval in the environment.
3. Builds with no cache, labels the image with its version and commit, and runs the full container checks.
4. Pushes only that version tag. It never touches `latest` or an older tag.
5. Records the digest, and uploads a manifest and the versioned skill.

To disable it again, set the variable to anything else, or delete the secrets.

Consider also enabling **tag immutability** for the Docker Hub repository, if your plan offers it.

The first publish tests the push step for real. It cannot be fully tested earlier without publishing.

## Reading CI

GitHub requires a login to read job logs. `container-checks.sh` also reports failures as annotations, which are public, so that agents without credentials can read them.
