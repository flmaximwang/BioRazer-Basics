import os, re

def dump_AF2_jobs(af2_jobs: list[set], target_file: str):
    with open(target_file, 'w') as f:
        for job_name, full_seq in af2_jobs:
            final_seq = "\n".join([full_seq[i:i+80] for i in range(0, len(full_seq), 80)])
            f.write(f">{job_name}\n")
            f.write(f"{final_seq}\n")

def remove_predicted_AF2_jobs(af2_jobs: list[set], res_dir: str):
    jobs_to_remove = []
    for i, job_request in enumerate(af2_jobs):
        job_name = job_request[0]
        job_dir = os.path.join(res_dir, f"{job_name}")
        if os.path.exists(job_dir):
            jobs_to_remove.append(i)
    for i in reversed(jobs_to_remove):
        af2_jobs.pop(i)
