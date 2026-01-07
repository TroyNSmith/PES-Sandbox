import os
import subprocess
import time
from math import ceil
from pathlib import Path

from hyperqueue import Client, Job
from hyperqueue.ffi.protocol import ResourceRequest
from hyperqueue.task.function import PythonEnv

from .ref import CustomTypes as CT

HQ_SERVER_DIR = Path(os.environ["HQ_SERVER_DIR"])
HQ_LOG = HQ_SERVER_DIR / "server.log"
DEFAULT_RESOURCES = ResourceRequest(cpus=1, resources={"mem": 500})

def pixi_activation_hook() -> str:
    """Returns the activation script for pixi env."""
    return subprocess.check_output(["pixi", "shell-hook"], text=True)

def start_server():
    """Starts HQ server in background if not already running."""
    if HQ_LOG.exists():
        HQ_LOG.unlink()  # clear old log

    # Clear stale access token or lock
    for fname in ["access-token", "lock", "server.pid"]:
        stale = HQ_SERVER_DIR / fname
        if stale.exists():
            stale.unlink()

    HQ_SERVER_DIR.mkdir(parents=True, exist_ok=True)
    log_file = HQ_LOG.open("w")

    subprocess.Popen(
        ["hq", "server", "start", "--server-dir", str(HQ_SERVER_DIR)],
        stdout=log_file,
        stderr=subprocess.STDOUT,
    )


def submit_tasks_orca(task_graph: CT.NetworkXGraph):
    """Submits ORCA calculations from a graph whose nodes contain submit paths and edges describe job dependencies."""

    def _bash(filename: str):
        cmd = f"bash {filename}.sh"
        def fn():
            subprocess.run(cmd, shell=True, executable="/bin/bash")

        return fn

    def _slurm_alloc(pars):
        mem_mib = ceil(pars.max_memory / 1.049)
        return [
            "hq",
            "alloc",
            "add",
            "slurm",
            "--time-limit",
            pars.time_limit,
            f"--cpus={pars.processors}",
            f"--resource=mem=sum({mem_mib})",
            "--",
            "--partition=batch",
            f"--ntasks={pars.processors}",
            f"--mem-per-cpu={pars.max_memory}",
            f"--gres=lscratch:{pars.lscratch_size}",
        ]

    start_server()
    wait_for_server()

    all_tasks = {}
    allocation_requests = set()

    job = Job()
    # py_env = PythonEnv(prologue=pixi_activation_hook)
    # client = Client(HQ_SERVER_DIR, python_env=py_env)
    client = Client(HQ_SERVER_DIR)

    for n, d in task_graph.nodes(data=True):
        pars = d["pars"]

        dependent_tasks = [all_tasks[dep] for dep in list(task_graph.predecessors(n))]

        allocation = _slurm_alloc(pars)
        if "".join(allocation) not in allocation_requests:
            subprocess.run(allocation, check=True)
            allocation_requests.add("".join(allocation))

        mem_mib = ceil(pars.max_memory / 1.049)

        task = job.function(
            fn=_bash(pars.name_out),
            cwd=n,
            deps=dependent_tasks,
            resources=ResourceRequest(cpus=pars.processors, resources={"mem": mem_mib}),
            stderr=n / "stderr.log",
            stdout=n / "stdout.log",
        )

        all_tasks[n] = task

    submitted = client.submit(job)
    client.wait_for_jobs([submitted])


def wait_for_server():
    """Waits until HQ server is responsive."""
    for _ in range(10):
        try:
            subprocess.run(
                ["hq", "alloc", "list"], check=True, stdout=subprocess.DEVNULL
            )
            return
        except subprocess.CalledProcessError:
            time.sleep(2)

    raise RuntimeError("HQ server not responding after multiple attempts.")
