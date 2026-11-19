"""MAGNA++ command-line interface: download, LEDA (.gw) input files, runs and statistics."""

import subprocess
import zipfile

import pandas as pd
import requests

from cmn.paths import DATA, SRC
from cmn.string_db import STRING_DIR, edges_file

URL = "https://academicweb.nd.edu/~cone/MAGNA++/downloads/MAGNA++_CLI_Linux64_01_30_2015.zip"
TOOL_DIR = SRC.parent / "tools" / "magna"
LINUX_BINARY = TOOL_DIR / "cli_linux" / "magnapp_cli_linux64"
IMAGE = "ubuntu:22.04"
RUNS_DIR = DATA / "magna"


def download_magna_linux():
    """Download and unpack the Linux MAGNA++ CLI (skipped if already present)."""
    if not LINUX_BINARY.exists():
        TOOL_DIR.mkdir(parents=True, exist_ok=True)
        archive = TOOL_DIR / "cli_linux.zip"
        archive.write_bytes(requests.get(URL, timeout=120).content)
        with zipfile.ZipFile(archive) as z:
            z.extractall(TOOL_DIR / "cli_linux")
        LINUX_BINARY.chmod(0o755)
    return LINUX_BINARY


def gw_file(name):
    return STRING_DIR / f"{name}.gw"


def write_gw(name):
    """Convert an edge-list file into LEDA .gw format (nodes numbered from 1 in the order listed)."""
    edges = pd.read_csv(edges_file(name), sep=" ", header=None, names=["u", "v"])
    nodes = pd.unique(edges[["u", "v"]].values.ravel())
    index = {n: i for i, n in enumerate(nodes, 1)}
    lines = ["LEDA.GRAPH", "void", "void", "-2", str(len(nodes))]
    lines += [f"|{{{n}}}|" for n in nodes]
    lines.append(str(len(edges)))
    lines += [f"{index[u]} {index[v]} 0 |{{}}|" for u, v in edges.itertuples(index=False)]
    path = gw_file(name)
    path.write_text("\n".join(lines) + "\n")
    return path, len(nodes), len(edges)


def run(small, large, population, generations, measure="S3", threads=None, tag=""):
    """Align network `small` (fewer nodes) to `large` with the Linux CLI in a sandboxed Docker container.

    The container has no network, sees the binary and the .gw files read-only, and writes only to RUNS_DIR.
    No node-similarity file is given, so alpha keeps its default of 1 (edge conservation only).
    Returns the parsed final statistics.
    """
    download_magna_linux()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    prefix = f"{small}_{large}{tag}"
    cmd = ["docker", "run", "--rm", "--platform", "linux/amd64", "--network", "none",
           "-v", f"{LINUX_BINARY.parent}:/magna:ro", "-v", f"{STRING_DIR}:/nets:ro", "-v", f"{RUNS_DIR}:/out",
           "-w", "/out", IMAGE, f"/magna/{LINUX_BINARY.name}",
           "-G", f"/nets/{gw_file(small).name}", "-H", f"/nets/{gw_file(large).name}", "-o", prefix,
           "-m", measure, "-p", str(population), "-n", str(generations)]
    if threads:
        cmd += ["-t", str(threads)]
    with open(RUNS_DIR / f"{prefix}.log", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    return read_stats(prefix)


def read_stats(prefix):
    """Key/value pairs of a MAGNA++ <prefix>_final_stats.txt file."""
    lines = (RUNS_DIR / f"{prefix}_final_stats.txt").read_text().split("\n")
    return {k: v for k, v in (line.split()[:2] for line in lines if line.strip())}
