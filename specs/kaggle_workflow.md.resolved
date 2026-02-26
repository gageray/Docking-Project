# Kaggle CLI GNINA Workflow

Here is the condensed guide on leveraging the official Kaggle CLI to run your GNINA docking jobs. Since your `gnina_worker.py` handles its own Google Drive uploads, the Kaggle side just needs to receive the input queue and execute the worker.

### 1. Initial Setup
1. **API Token**: Go to Kaggle Settings -> Create New Token. Save the `kaggle.json` to `~/.kaggle/kaggle.json`.
2. **Install CLI**: `pip install kaggle`

### 2. Uploading Inputs (Datasets)
Kaggle Kernels read input data from Kaggle Datasets. Your `local_dispatcher.py` should output the prepared receptors, ligands, and the `work_queue.json` into a single directory (e.g., `data/kaggle_input`).

1. **Initialize Metadata**:
   ```bash
   kaggle datasets init -p data/kaggle_input
   ```
   *Edit the generated `dataset-metadata.json` to set your desired `title` and `id` (e.g., `your-username/gnina-docking-queue`).*
2. **Create/Update the Dataset**:
   ```bash
   # First time creation:
   kaggle datasets create -p data/kaggle_input -r zip
   # Subsequent updates (new docking batches):
   kaggle datasets version -p data/kaggle_input -m "New batch of jobs" -r zip
   ```

### 3. Setting Up the Worker Kernel
Your worker script (`gnina_worker.py`) will run as a Kaggle Kernel (Script type) attached to the dataset you just created.

1. **Initialize Kernel Metadata** in the directory containing `gnina_worker.py`:
   ```bash
   kaggle kernels init -p scripts/kaggle_worker
   ```
2. **Configure `kernel-metadata.json`**:
   Edit the file to look something like this:
   ```json
   {
     "id": "your-username/gnina-gpu-worker",
     "title": "GNINA GPU Worker",
     "code_file": "gnina_worker.py",
     "language": "python",
     "kernel_type": "script",
     "enable_gpu": "true",
     "enable_internet": "true",
     "dataset_sources": ["your-username/gnina-docking-queue"]
   }
   ```
   *Note: `enable_internet: true` is crucial for your worker to upload results back to Google Drive.*

### 4. Running the Jobs
Whenever a new `work_queue.json` is ready:

1. Update the dataset (see step 2).
2. Push the kernel to trigger execution. You can explicitly request T4 GPUs:
   ```bash
   kaggle kernels push -p scripts/kaggle_worker --accelerator NvidiaTeslaT4
   ```

### 5. Monitoring
You can monitor the status of your worker kernel directly from the CLI:
```bash
# Check if it's still running or finished
kaggle kernels status your-username/gnina-gpu-worker

# Tail the logs / output files if needed
kaggle kernels output your-username/gnina-gpu-worker
```
