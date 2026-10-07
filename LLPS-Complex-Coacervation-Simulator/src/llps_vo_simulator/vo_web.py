"""Standalone, local Voorn--Overbeek phase-diagram simulator."""

from __future__ import annotations

import argparse
import threading
import webbrowser
from pathlib import Path

from flask import Flask, jsonify, render_template, request, send_from_directory

from .vo_asymmetric import (
    AsymmetricVOParameters,
    simulate_2d as simulate_asymmetric_2d,
    simulate_high_dim as simulate_asymmetric_high_dim,
)
from .vo_theory import VOParameters, profile_at_mu, simulate


app = Flask(__name__)


@app.get("/")
def index() -> str:
    return render_template("vo_simulator.html")


@app.get("/mathjax/<path:filename>")
def mathjax_asset(filename: str):
    return send_from_directory(Path(__file__).parent / "static" / "vendor" / "mathjax", filename)


@app.post("/api/simulate")
def calculate():
    try:
        params = VOParameters.from_mapping(request.get_json(silent=True) or {})
        return jsonify(simulate(params))
    except (ValueError, TypeError, OverflowError) as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/profile")
def calculate_profile():
    try:
        payload = request.get_json(silent=True) or {}
        params = VOParameters.from_mapping(payload)
        return jsonify(profile_at_mu(params, float(payload["mu_ion_pair"])))
    except (ValueError, TypeError, KeyError, OverflowError) as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/asymmetric/2d")
def calculate_asymmetric_2d():
    try:
        payload = request.get_json(silent=True) or {}
        params = AsymmetricVOParameters.from_mapping(payload)
        return jsonify(simulate_asymmetric_2d(
            params,
            charge_ratio=float(payload.get("charge_ratio", 1.0)),
            polymer_max_m=float(payload.get("polymer_max_m", 0.025)),
            salt_max_m=float(payload.get("salt_max_m", 1.0)),
            polymer_count=int(payload.get("polymer_count", 48)),
            salt_count=int(payload.get("salt_count", 42)),
            feed_polymer_m=float(payload.get("feed_polymer_m", 0.01)),
            feed_salt_m=float(payload.get("feed_salt_m", 0.05)),
            salt_type=str(payload.get("salt_type", "NaCl")),
        ))
    except (ValueError, TypeError, OverflowError) as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/asymmetric/high-dimensional")
def calculate_asymmetric_high_dimensional():
    try:
        payload = request.get_json(silent=True) or {}
        params = AsymmetricVOParameters.from_mapping(payload)
        return jsonify(simulate_asymmetric_high_dim(
            params,
            rq_min=float(payload.get("rq_min", 0.25)),
            rq_max=float(payload.get("rq_max", 4.0)),
            polymer_min_m=float(payload.get("polymer_min_m", 0.0001)),
            polymer_max_m=float(payload.get("polymer_max_m", 0.025)),
            salt_min_m=float(payload.get("salt_min_m", 0.0)),
            salt_max_m=float(payload.get("salt_max_m", 1.0)),
            ratio_count=int(payload.get("ratio_count", 19)),
            polymer_count=int(payload.get("polymer_count", 24)),
            salt_count=int(payload.get("salt_count", 61)),
            salt_type=str(payload.get("salt_type", "NaCl")),
        ))
    except (ValueError, TypeError, OverflowError) as exc:
        return jsonify({"error": str(exc)}), 400


def run_gui(port: int = 8770, open_browser: bool = True) -> None:
    if open_browser:
        threading.Timer(0.8, lambda: webbrowser.open(f"http://127.0.0.1:{port}/")).start()
    app.run(host="127.0.0.1", port=port, debug=False, use_reloader=False, threaded=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Standalone VO coacervation simulator")
    parser.add_argument("--port", type=int, default=8770)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    run_gui(port=args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    main()
