#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_PORT="${VO_SIMULATOR_PORT:-8770}"
if ! [[ "$BASE_PORT" =~ ^[0-9]+$ ]] || (( BASE_PORT < 1024 || BASE_PORT > 65530 )); then
  echo "VO_SIMULATOR_PORT must be an integer between 1024 and 65530."
  exit 2
fi

is_python_310_or_newer() {
  "$1" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1
}

is_current_simulator() {
  local page
  page="$(/usr/bin/curl --silent --fail --max-time 2 "http://127.0.0.1:$1/" 2>/dev/null || true)"
  [[ "$page" == *'viewBox="0 0 900 530"'* && "$page" == *'formula concentration (M)'* ]]
}

PORT=""
for offset in 0 1 2 3 4 5; do
  candidate=$((BASE_PORT + offset))
  if is_current_simulator "$candidate"; then
    /usr/bin/open "http://127.0.0.1:${candidate}/"
    exit 0
  fi
  if ! /usr/bin/curl --silent --fail --max-time 2 "http://127.0.0.1:${candidate}/" >/dev/null; then
    PORT="$candidate"
    break
  fi
done
if [[ -z "$PORT" ]]; then
  echo "No free local port found between $BASE_PORT and $((BASE_PORT + 5))."
  exit 1
fi

cd "$PROJECT_DIR"
if [[ -x "$PROJECT_DIR/.venv/bin/python" ]]; then
  PYTHON="$PROJECT_DIR/.venv/bin/python"
  if ! is_python_310_or_newer "$PYTHON"; then
    echo "The existing .venv uses Python older than 3.10. Remove this package's .venv and install Python 3.10+ first."
    exit 1
  fi
else
  PYTHON=""
  for candidate in "${VO_PYTHON:-}" python3.13 python3.12 python3.11 python3.10 python3 /opt/anaconda3/bin/python; do
    [[ -n "$candidate" ]] || continue
    if command -v "$candidate" >/dev/null 2>&1 && is_python_310_or_newer "$(command -v "$candidate")"; then
      PYTHON="$(command -v "$candidate")"
      break
    elif [[ -x "$candidate" ]] && is_python_310_or_newer "$candidate"; then
      PYTHON="$candidate"
      break
    fi
  done
  if [[ -z "$PYTHON" ]]; then
    echo "Python 3.10+ is required. Install it, then run this launcher again."
    echo "For example: conda create -n vo-simulator python=3.11"
    exit 1
  fi
  "$PYTHON" -m venv .venv
  PYTHON="$PROJECT_DIR/.venv/bin/python"
fi
if ! "$PYTHON" -c 'import flask, numpy, scipy' >/dev/null 2>&1; then
  "$PYTHON" -m pip install -e "$PROJECT_DIR"
fi
URL="http://127.0.0.1:${PORT}/"
/usr/bin/nohup "$PYTHON" -m llps_vo_simulator.vo_web --port "$PORT" --no-browser >> "$PROJECT_DIR/vo_simulator.log" 2>&1 </dev/null &
for attempt in $(seq 1 35); do
  if /usr/bin/curl --silent --fail --max-time 1 "$URL" >/dev/null; then
    /usr/bin/open "$URL"
    exit 0
  fi
  sleep 1
done
/usr/bin/tail -n 40 "$PROJECT_DIR/vo_simulator.log"
exit 1
