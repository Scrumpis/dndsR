#!/usr/bin/env bash
# PreToolUse hook for the audit agents: allows only commands on a read-only list.
# Anything else, including running R, is denied; auditors report such commands
# for the calling session to run.
set -u
command -v jq >/dev/null || { echo "read-only-bash: jq is required" >&2; exit 2; }
input=$(cat)
case $(jq -r '.agent_type // empty' <<<"$input") in
  scientific-auditor | software-auditor) ;;
  *) exit 0 ;;
esac
cmd=$(jq -r '.tool_input.command // empty' <<<"$input")
[ -n "$cmd" ] || exit 0
deny() { printf 'Read-only audit agent: %s. Report the command for the calling session to run.\n' "$1" >&2; exit 2; }

# No output redirection to files, command substitution, or inline code.
stripped=$(printf '%s' "$cmd" | sed -E 's/[0-9]*>&[0-9]//g; s/[0-9]*>[[:space:]]*\/dev\/null//g')
case $stripped in *'>'* | *'$('* | *'`'* | *'<<'*) deny 'redirection and substitution are not allowed' ;; esac

# Check the first word of every segment, without being misled by & in "2>&1".
while IFS= read -r segment; do
  set -f
  # shellcheck disable=SC2086
  set -- $segment
  set +f
  [ $# -gt 0 ] || continue
  case $1 in
    ls | cat | head | tail | grep | rg | wc | sort | uniq | cut | diff | stat | file | du | pwd | echo | which | basename | dirname | realpath | tr | nl | comm | column | sha256sum | true) ;;
    jq) case " $* " in *' -i '* | *' --in-place '*) deny 'in-place editing is not allowed' ;; esac ;;
    find) case " $* " in *' -delete '* | *' -exec '* | *' -execdir '* | *' -ok '* | *' -fprint'*) deny 'find actions are not allowed' ;; esac ;;
    git)
      case ${2:-} in
        status | diff | log | show | ls-files | ls-tree | rev-parse | rev-list | merge-base | blame | grep | shortlog | describe | cat-file | for-each-ref) ;;
        *) deny "git ${2:-} is not on the read-only list" ;;
      esac
      case " $* " in *' --output'* | *' --ext-diff'*) deny 'that git option is not allowed' ;; esac
      ;;
    *) deny "$1 is not on the read-only list" ;;
  esac
done < <(printf '%s\n' "$stripped" | sed -E 's/&&|\|\||[;|&]/\n/g')
exit 0
