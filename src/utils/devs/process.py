
def process_snakemake_standard_output(snakemake_cmd, outfilelog):
    time_stamp_pattern = r'\[\w{3} \w{3} \d{2} \d{2}:\d{2}:\d{2} \d{4}\]'
    rule_name_pattern = r'rule (\w+):$'
    jobid_pattern = r'jobid: (\d+)'

    running_jobs = {}
    finished_jobs = {}

    proc = subprocess.Popen(f"{snakemake_cmd} ", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    with open(outfilelog, "w") as outfile:

        def process_next_line():
            # Read next line
            line = proc.stdout.readline().decode('utf-8').rstrip()
            # writes line to log file
            outfile.write(line + "\n")
            # Check if line matches timestamp pattern
            timestamp_match = re.search(time_stamp_pattern, line)
            if timestamp_match:
                # A new timestamp means a new job is starting or has finished
                rule_name = None
                job_number = None

                # Check next line for rule name
                next_line = proc.stdout.readline().decode('utf-8').rstrip()
                # writes line to log file
                outfile.write(next_line + "\n")
                rule_name_match = re.search(rule_name_pattern, next_line)

                if rule_name_match:
                    rule_name = rule_name_match.group(1)

                    # Read lines until we find the job id or reach the next timestamp
                    while True:
                        following_line = proc.stdout.readline().decode('utf-8').rstrip()
                        outfile.write(following_line + "\n")
                        # Check if line matches timestamp pattern
                        if re.search(time_stamp_pattern, following_line):
                            # We reached the next timestamp, stop reading lines
                            process_next_line()
                            return

                        job_id_match = re.search(jobid_pattern, following_line)

                        if job_id_match:
                            job_number = job_id_match.group(1)
                            print(f"Job {job_number} of rule {rule_name} is currently running.")
                            running_jobs[job_number] = rule_name
                else:
                    # Check if the line indicates a job has finished
                    finished_match = re.search(r'Finished job (\d+)\.', next_line)
                    if finished_match:
                        finished_job_number = finished_match.group(1)
                        print(f"Finished job Number: {finished_job_number}")
                        finished_jobs[finished_job_number] = rule_name

                        while True:
                            following_following_line = proc.stdout.readline().decode('utf-8').rstrip()
                            outfile.write(following_following_line + "\n")
                            # Check if line matches timestamp pattern
                            if re.search(time_stamp_pattern, following_following_line):
                                # We reached the next timestamp, stop reading lines
                                process_next_line()
                                return

                            progress_match = re.search(r'\d+ of \d+ steps \((\d+)%\) done', following_following_line)

                            if progress_match:
                                progress = progress_match.group(1)
                                print(f"Job {finished_job_number} of rule {rule_name} has finished with {progress}% progress.")
            else:
                # Line didn't match timestamp pattern, check next line
                process_next_line()

        # Start processing output
        process_next_line()

        return running_jobs, finished_jobs
    

def parse_snakemake_output(snakemake_cmd, outlogfile):
    time_stamp_pattern = r'\[\w{3} \w{3} \d{2} \d{2}:\d{2}:\d{2} \d{4}\]'
    rule_name_pattern = r'rule (\w+):$'
    jobid_pattern = r'jobid: (\d+)'
    running_jobs = {}
    finished_jobs = {}

    # Define the recursive function
    def process_next_line(line):
        nonlocal running_jobs, finished_jobs
        if re.search(time_stamp_pattern, line):
            # A new timestamp means a new job is starting or has finished
            rule_name = None
            job_number = None
            while True:
                next_line = proc.stdout.readline().decode('utf-8').rstrip()
                if re.search(time_stamp_pattern, next_line):
                    # We reached the next timestamp, stop reading lines and return
                    return next_line
                outfile.write(next_line + "\n")
                rule_name_match = re.search(rule_name_pattern, next_line)
                if rule_name_match:
                    # We found the name of the rule, which means a new job is starting
                    rule_name = rule_name_match.group(1)
                    job_id_match = re.search(jobid_pattern, next_line)
                    if job_id_match:
                        # We also found the ID of the job
                        job_number = job_id_match.group(1)
                        print(f"Job {job_number} of rule {rule_name} is currently running.")
                        running_jobs[job_number] = rule_name
                else:
                    # Check if the line indicates a job has finished
                    finished_match = re.search(r'Finished job (\d+)\.', next_line)
                    if finished_match:
                        finished_job_number = finished_match.group(1)
                        print(f"Finished job Number: {finished_job_number}")
                        finished_jobs[rule_name] = finished_job_number
                        while True:
                            following_line = process_next_line(proc.stdout.readline().decode('utf-8').rstrip())
                            if not following_line:
                                break
                            progress_match = re.search(r'\d+ of \d+ steps \((\d+)%\) done', following_line)
                            if progress_match:
                                progress = progress_match.group(1)
                                print(f"Job {finished_job_number} of rule {rule_name} has finished with {progress}% progress.")
            return next_line
        else:
            # The line doesn't match the timestamp pattern, keep reading lines
            outfile.write(line + "\n")
            return process_next_line(proc.stdout.readline().decode('utf-8').rstrip())

    # Open output file
    with open("output.log", "w") as outfile:
        # Start the Snakemake process
        proc = subprocess.Popen(f"{snakemake_cmd} ", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # Process the output lines recursively
        while True:
            line = process_next_line(proc.stdout.readline().decode('utf-8').rstrip())
            if not line:
                break

    return running_jobs, finished_jobs