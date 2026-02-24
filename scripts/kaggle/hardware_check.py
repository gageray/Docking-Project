"""Hardware isolation and GPU management utilities for Kaggle notebooks.

Test GPU visibility, manage CUDA_VISIBLE_DEVICES, and verify parallel execution.
"""

import subprocess
import os
import queue
import concurrent.futures


def check_gpu_info():
    """Query GPU hardware information.

    Returns:
        dict: GPU count, names, and memory info

    Example:
        >>> info = check_gpu_info()
        >>> print(f"GPUs available: {info['count']}")
    """
    result = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
        capture_output=True, text=True, check=False
    )

    if result.returncode != 0:
        return {"count": 0, "available": False, "error": result.stderr}

    gpus = result.stdout.strip().split('\n')
    return {
        "count": len(gpus),
        "available": True,
        "gpus": [{"info": gpu} for gpu in gpus]
    }


def test_gpu_isolation(gpu_id):
    """Test if a specific GPU can be isolated via CUDA_VISIBLE_DEVICES.

    Args:
        gpu_id (int): GPU index to test

    Returns:
        dict: Isolation test results

    Example:
        >>> result = test_gpu_isolation(0)
        >>> print(result['visible_gpus'])  # Should be 1
    """
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": str(gpu_id)}

    # Check CUDA visibility via PyTorch
    py_cmd = "import torch; print(torch.cuda.device_count())"

    gpu_result = subprocess.run(
        ["python", "-c", py_cmd],
        env=env, capture_output=True, text=True
    )

    visible_gpus = int(gpu_result.stdout.strip()) if gpu_result.returncode == 0 else 0

    # Get CPU count
    cpu_result = subprocess.run(
        ["nproc"], capture_output=True, text=True
    )
    cpu_count = int(cpu_result.stdout.strip()) if cpu_result.returncode == 0 else 0

    return {
        "gpu_id": gpu_id,
        "visible_gpus": visible_gpus,
        "cpu_count": cpu_count,
        "isolated": visible_gpus == 1
    }


def test_parallel_gpu_isolation(num_workers=2, gpu_ids=None):
    """Test hardware isolation across parallel workers.

    Args:
        num_workers (int): Number of parallel workers. Default 2.
        gpu_ids (list, optional): GPU IDs to test. Defaults to [0, 1].

    Returns:
        list: Test results for each worker

    Example:
        >>> results = test_parallel_gpu_isolation(num_workers=2)
        >>> for r in results:
        ...     print(f"Worker GPU {r['gpu_id']}: {r['visible_gpus']} visible")
    """
    if gpu_ids is None:
        gpu_ids = list(range(num_workers))

    print(f"=== Testing hardware isolation with {num_workers} workers ===")

    # Thread-safe GPU queue
    gpu_queue = queue.Queue()
    for gpu_id in gpu_ids:
        gpu_queue.put(gpu_id)

    def worker(_):
        gpu_id = gpu_queue.get()
        try:
            result = test_gpu_isolation(gpu_id)
            return result
        finally:
            gpu_queue.put(gpu_id)

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(worker, range(num_workers)))

    # Print results
    for r in results:
        isolated_str = "✓ ISOLATED" if r["isolated"] else "✗ NOT ISOLATED"
        print(f"  Worker GPU {r['gpu_id']}: {r['visible_gpus']} visible, "
              f"{r['cpu_count']} CPUs | {isolated_str}")

    return results


def create_gpu_queue(gpu_ids=None):
    """Create a thread-safe GPU queue for parallel docking.

    Args:
        gpu_ids (list, optional): GPU IDs to add to queue. Defaults to [0, 1].

    Returns:
        queue.Queue: GPU queue for worker assignment

    Example:
        >>> gpu_q = create_gpu_queue([0, 1])
        >>> gpu_id = gpu_q.get()  # Worker gets GPU
        >>> # ... do work ...
        >>> gpu_q.put(gpu_id)  # Return GPU to queue
    """
    if gpu_ids is None:
        gpu_ids = [0, 1]

    gpu_queue = queue.Queue()
    for gpu_id in gpu_ids:
        gpu_queue.put(gpu_id)

    return gpu_queue


def get_optimal_workers():
    """Determine optimal number of parallel workers based on available GPUs.

    Returns:
        int: Recommended number of workers

    Example:
        >>> workers = get_optimal_workers()
        >>> print(f"Use {workers} parallel workers")
    """
    gpu_info = check_gpu_info()

    if not gpu_info["available"]:
        return 1

    # Use all available GPUs
    return gpu_info["count"]


def verify_hardware_setup(test_parallel=True):
    """Complete hardware verification: GPU info and isolation tests.

    Args:
        test_parallel (bool): Run parallel isolation tests. Default True.

    Returns:
        dict: Complete hardware status

    Example:
        >>> status = verify_hardware_setup()
        >>> if status['gpus']['available']:
        ...     print(f"Ready with {status['gpus']['count']} GPUs")
    """
    print("=== Hardware Verification ===\n")

    status = {
        "gpus": check_gpu_info(),
        "optimal_workers": get_optimal_workers()
    }

    print(f"GPUs detected: {status['gpus']['count']}")
    print(f"Optimal workers: {status['optimal_workers']}")

    if test_parallel and status["gpus"]["available"] and status["gpus"]["count"] > 1:
        print()
        isolation_results = test_parallel_gpu_isolation(
            num_workers=status["optimal_workers"]
        )
        status["isolation_test"] = isolation_results
        status["all_isolated"] = all(r["isolated"] for r in isolation_results)
    else:
        status["isolation_test"] = None
        status["all_isolated"] = None

    return status
